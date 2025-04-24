#include <atrip/Reader.hpp>
#include <cstring>
#include <atrip/Complex.hpp>
#include <thread>
#include <chrono>

namespace atrip {

template <typename F>
static void permute_copy(std::vector<size_t> ranges,
                         const F* reorder_buffer,
                         F* source_buffer) {
    size_t dims = ranges.size();
    std::vector<size_t> permutation(dims);
    std::iota(permutation.begin(), permutation.end(), 0);
    std::reverse(permutation.begin(), permutation.end());

    std::vector<size_t> strides(dims, 1), permuted_strides(dims, 1);

    // Compute strides (column-major style, last index varies fastest)
    for (int i = dims - 2; i >= 0; --i) {
        strides[i] = strides[i + 1] * ranges[i + 1];
        permuted_strides[i] = permuted_strides[i + 1] * ranges[permutation[i + 1]];
    }

    // Multi-loop iteration (recursive lambda)
    std::vector<size_t> indices(dims, 0);
    std::function<void(size_t)> loop = [&](size_t depth) {
        if (depth == dims) {
            // Compute source & destination indices using strides
            size_t src_idx = 0, dest_idx = 0;
            for (size_t d = 0; d < dims; ++d) {
                src_idx  += indices[d] * strides[d];           // Normal order
                dest_idx += indices[permutation[d]] * permuted_strides[d];  // Permuted order
            }
            //std::cout << "source[ " << dest_idx << " ] <-- original[ "
            //          << src_idx << "] :" << reorder_buffer[src_idx] << "\n";
            source_buffer[dest_idx] = reorder_buffer[src_idx]; // Fast write!
            return;
        }

        // Recursively loop over dimensions
        for (indices[depth] = 0; indices[depth] < ranges[depth]; ++indices[depth]) {
            loop(depth + 1);
        }
    };

    loop(0); // Start recursion

}




template <typename F>
Sources<F> diskReader(const std::string &file_path,
                      std::vector<size_t> tensor_dimension,
                      std::vector<size_t> slice_mapping,
                      bool rowMajor) {

  size_t Nv = [&]{
      size_t v = 0;
      for(size_t i = 0; i < tensor_dimension.size(); ++i)
          if(slice_mapping[i] == 1)
              v = std::max(v, tensor_dimension[i]);
      return v;
  }();
  // check if the slicing is valid
  int transitions = 0;
  for (size_t i = 1; i < slice_mapping.size(); i++) {
    if (slice_mapping[i] != slice_mapping[i-1]) ++transitions;
  }
  assert(slice_mapping.size() && transitions == 1); // {0,1,0,1} is not allowed!
  // if we have {1,0,0,0} we have to work with column major storage!
  assert(static_cast<bool>(slice_mapping.front()) == rowMajor);


  std::vector<size_t> slice_dimension, source_dimension;

  for (auto i(0UL); i < slice_mapping.size(); i++) {
    slice_mapping[i] ? slice_dimension.push_back(tensor_dimension[i])
                     : source_dimension.push_back(tensor_dimension[i]);
  }

  RankMap<F> rank_map(slice_dimension);

  size_t s_sources = std::accumulate(source_dimension.begin(),
                                     source_dimension.end(),
                                     1UL,
                                     std::multiplies<size_t>());

  size_t number_slices = std::accumulate(slice_dimension.begin(),
                                         slice_dimension.end(),
                                         1UL,
                                         std::multiplies<size_t>());

  auto np(Atrip::np);
  size_t n_sources(0);
  if (np > 0) {
    n_sources = number_slices / np;
    if (number_slices % np > 0 && Atrip::rank < (number_slices % np)) n_sources++;
  }


  MPI_File handle;
  MPI_File_open(MPI_COMM_WORLD,
                file_path.c_str(),
                MPI_MODE_RDONLY,
                MPI_INFO_NULL,
                &handle);

  std::vector<std::vector<F>> sources(n_sources, std::vector<F>(s_sources));

  // in case we have to reorder the buffer we create a new array
  // otherwise we write to the result array directly
  F *buffer;
  if (rowMajor) buffer = new F[s_sources];
  size_t a,b;
  for (auto s(0UL); s < n_sources; s++) {
    // if we have a 1D map the result is smaller Nv, so b will always b=0
    std::tie(a,b) = orbitalMap(rank_map.find_element({Atrip::rank, s}), Nv);
    size_t off = (slice_dimension.size() == 1 || !rowMajor) ? a + Nv*b : a*Nv + b;
    MPI_Offset offset = s_sources * off * sizeof(F);
    if (!rowMajor) buffer = sources[s].data();
    if (MPI_SUCCESS != MPI_File_read_at(handle,
                                        offset,
                                        buffer,
                                        s_sources,
                                        traits::mpi::datatype_of<F>(),
                                        MPI_STATUS_IGNORE)) {
      throw "error reading!";
    }

    if (rowMajor) {
      permute_copy(source_dimension, buffer, sources[s].data());
    }
  }

  MPI_File_close(&handle);
  return {n_sources, s_sources, sources};
}



static std::vector<int> largest_factors(int N) {
    for (int i = N / 2; i >= 1; --i) if (N % i == 0) return {i, N / i};
    return {1, N};  // For N = 1 case
}
// This can only be a incomplete documentation of the algorithm as both,
// CTF and Atrip have complicated data-layouts.
// The general concept is as follows:
// 1.) bring the CTF tensor from any given distribution in an atrip-favorable one
//     (this means that the tensor is distributed along the virtual indices (Nv)
//     which are distributed in atrip (V_ABPH means distributed among first two indices)
// 2.) now we use a simple sendrecv pattern to send the local batches of the CTF
//     tensor to the atrip 'sources':
//     atrip: el = a + Nv * b for 2d slices; or el = a for 1d slices
//            elements el are distributed in node-aware round robin fashion
//            el=0: rank0 node 0; el=1: rank0 node1; ....
//     ctf:   processor mesh p=
//
//
// Notes: 1.) within CTF, the local tensor every ranks holds is stored in column major.
//            this implies that we have to rearrange the tensor data such that the
//            indices to slice are the slowest (V_PHAB)
//        2.)

#if defined(HAVE_CTF)
template <typename F>
Sources<F> ctfReader(CTF::Tensor<F>& tensor,
                     std::vector<size_t> tensor_dimension,
                     std::vector<size_t> slice_mapping) {

  size_t Nv = [&]{
      size_t v = 0;
      for(size_t i = 0; i < tensor_dimension.size(); ++i)
          if(slice_mapping[i] == 1)
              v = std::max(v, tensor_dimension[i]);
      return v;
  }();

  // We have to work with the communicator of the world's tensor.
  // the atrip::communicator is a subworld (or identical) of it

  MPI_Comm world = tensor.wrld->comm;
  int ctf_rank, ctf_np, atrip_np(Atrip::np);
  MPI_Comm_size(world, &ctf_np);
  MPI_Comm_rank(world, &ctf_rank);
  if (ctf_rank == 0)  assert(atrip_np);
  MPI_Bcast(&atrip_np, 1, MPI_INT, 0, world);

  // TODO: right now the logic works only for more than one core
  assert(ctf_np > 1);
  auto order(tensor_dimension.size());

  std::vector<size_t> slice_dimension, source_dimension, ptensor_dimension(tensor_dimension);
  assert(order == slice_mapping.size() && tensor.order == (int) order);

  bool is_reversed;
  auto reorder(slice_mapping);
  std::sort(reorder.begin(), reorder.end());
  if (reorder == slice_mapping) {
    is_reversed = false;
    std::iota(reorder.begin(), reorder.end(), 0);
    LOG(0, "Atrip") << "CitfReader: " << tensor.name << " | tensor is in correct order." << std::endl;
  } else {
    LOG(0, "Atrip") << "CitfReader: " << tensor.name << " | tensor will be reverted." << std::endl;
    is_reversed = true;
    std::iota(reorder.begin(), reorder.end(), 0);
    std::sort(reorder.begin(), reorder.end(), [&slice_mapping](size_t a, size_t b) {
      return slice_mapping[a] < slice_mapping[b];
    });
    for (size_t i(0); i < order; i++) ptensor_dimension[i] = tensor_dimension[reorder[i]];
    std::sort(slice_mapping.begin(), slice_mapping.end());
  }

  for (auto i(0UL); i < slice_mapping.size(); i++) {
    slice_mapping[i] ? slice_dimension.push_back(ptensor_dimension[i])
                     : source_dimension.push_back(ptensor_dimension[i]);
  }


  RankMap<F> rank_map(slice_dimension);

  size_t s_sources = std::accumulate(source_dimension.begin(),
                                     source_dimension.end(),
                                     1UL,
                                     std::multiplies<size_t>());

  size_t number_slices = std::accumulate(slice_dimension.begin(),
                                         slice_dimension.end(),
                                         1UL,
                                         std::multiplies<size_t>());

  size_t n_sources(0);
  // This means: you are not MPI_UNDEFINED in the atrip communicator!
  if (Atrip::np > 0) {
    n_sources = number_slices / Atrip::np;
    if (number_slices % Atrip::np > 0
        && Atrip::rank < (number_slices % Atrip::np)) n_sources++;
  }
  std::vector<std::vector<F>> sources(n_sources, std::vector<F>(s_sources));

  assert(tensor.order == tensor_dimension.size());
  for (auto i(0); i < tensor.order; i++) assert(tensor.lens[i] == tensor_dimension[i]);

  // identify the slice-mapping in order to create the CTF::Partition
  // one dimensional case is easy - we have to find the non-zero entry
  std::vector<int> plens(tensor.order, 1);
  if (slice_dimension.size() == 1) {
    for (auto i(0UL); i < slice_mapping.size(); i++) {
      if (slice_mapping[i] > 0) plens[i] = atrip_np;
    }
  } else if (slice_dimension.size() == 2) {
    auto facs(largest_factors(atrip_np));
    size_t u(0);
    for (auto i(0UL); i < slice_mapping.size(); i++) {
      if (slice_mapping[i] > 0) plens[i] = facs[u++];
    }
  } else {
    assert(0);
  }
  CTF::Partition part(tensor.order, plens.data());
  std::string albet{"abcdefghijklmnopqrstuvwxyz"};
  auto s1 = albet.substr(0, tensor.order);
  auto s2 = albet.substr(0, tensor.order);

  std::vector<int64_t> lens(tensor.lens, tensor.lens + order);
  for (size_t i(0); i < order; i++) lens[i] = tensor.lens[reorder[i]];
  for (size_t i(0); i < order; i++) s2[i] = s1[reorder[i]];

  //TODO: this is hardcoded for the Vpphh
  // we have to transpose the processor grid
  auto s3(s1);
  if (is_reversed && slice_dimension.size() == 2) {
    std::swap(s3[2], s3[3]);
    std::swap(plens[2], plens[3]);
  }

  CTF::Tensor<F> ptensor(order, lens.data(), tensor.sym, *tensor.wrld, s3.c_str(), part[s1.c_str()]);
  ptensor[s2.c_str()] = tensor[s1.c_str()];

  ptensor.set_name("pVpphh");

  // we run into problems if the local buffer gets too large!
  assert(s_sources*sizeof(F) < std::numeric_limits<int>::max());

  std::vector<int> list_target,list_origin;

  // Atrip source distribution
  for (int i(0); i < n_sources; i++) {
    list_target.push_back(rank_map.find_element({Atrip::rank, i}));
  }
  // Ctf distribution
  std::vector<size_t> grid;
  for (auto i(0); i < slice_mapping.size(); i++) {
    if (slice_mapping[i] > 0) grid.push_back(plens[i]);
  }
  // if we have reversed the tensor the grid is also revered
  // the grid is also reversed?! result is wrong if we do not reverse once more
  if (is_reversed) std::reverse(grid.begin(), grid.end());


  if (slice_dimension.size() == 1) {
    assert(grid.size() == 1);
    for (auto _v(0); _v < Nv; _v++) {
      if (ctf_rank == _v % grid[0]) list_origin.push_back(_v);
    }
  } else {
    for (auto _b(0); _b < Nv; _b++)
    for (auto _a(0); _a < Nv; _a++) {
      auto proc =  _a % grid[0] + (_b % grid[1]) * grid[0];
      if (ctf_rank == proc) {
        auto el = (is_reversed) ?  _b + _a*Nv : _a + _b*Nv;
        // this is the tag for the element atrip is waiting for.
        list_origin.push_back(el);
      }
    }
  }
  auto n_list = (slice_dimension.size() == 1 ) ? Nv : Nv*Nv;
  // this is a poor man's database which tells to origin and target of the data
  std::vector<int> send_recv_list(n_list*2, -1), glb_send_recv_list(n_list*2, -1);

  for (auto _v(0); _v < n_list; _v++) {
    if (std::find(list_origin.begin(), list_origin.end(), _v) != list_origin.end())
      send_recv_list[_v*2] = static_cast<int>(ctf_rank);
    if (std::find(list_target.begin(), list_target.end(), _v) != list_target.end())
      send_recv_list[_v*2+1] = static_cast<int>(ctf_rank);
  }
  //distribute the result on all ranks
  MPI_Allreduce(send_recv_list.data(),
                glb_send_recv_list.data(),
                2*n_list,
                MPI_INT,
                MPI_MAX,
                world);

  auto get_matching_indices = [&glb_send_recv_list](int rank, bool count_even) {
  std::vector<int> matching_indices;
    for (size_t index = 0; index < glb_send_recv_list.size(); ++index) {
        bool is_target_index = count_even ? (index % 2 == 0) : (index % 2 != 0);
        if (is_target_index && glb_send_recv_list[index] == rank) {
            matching_indices.push_back(index/2);
        }
    }
    return matching_indices;
  };

  auto origin_indices = get_matching_indices(ctf_rank, 0);
  auto target_indices = get_matching_indices(ctf_rank, 1);

  size_t n_origins = origin_indices.size();
  size_t n_targets = target_indices.size();

  //++++++++++++++++++
  // +++ MPI PHASE +++
  //++++++++++++++++++
  std::vector<MPI_Request> send_requests(n_targets);
  std::vector<MPI_Request> recv_requests(n_origins);

  std::vector<int64_t> loc_lens(ptensor.pad_edge_len, ptensor.pad_edge_len + order);
  loc_lens.erase(loc_lens.begin(), loc_lens.begin() + order - grid.size());
  // we know that the ctf tensor is distribtued by the n most right dimensions
  if (is_reversed) std::reverse(loc_lens.begin(), loc_lens.end());
  // this is wrong for general tensors...but every 2d slice is currently reversed
  if (is_reversed && slice_dimension.size() > 1) {
    for (size_t i = grid.size(); i-- > 0;) {
      loc_lens[i] /= grid[i];
    }
  }
  // Recv phase
  for (size_t i(0); i < n_origins; i++) {
    auto &t = origin_indices[i];
    int recv_rank = glb_send_recv_list[t*2];
    char *dest = reinterpret_cast<char *>(sources[i].data());
    MPI_Irecv(dest,
              static_cast<int>(s_sources*sizeof(F)),
              MPI_CHAR,
              recv_rank,
              t,
              world,
              &recv_requests[i]);
  }
  // Send phase
  for (size_t i(0); i < n_targets; i++) {
    auto &t = target_indices[i];
    int j(i);
    if (is_reversed && slice_dimension.size() > 1) {
      int aa(t%Nv), bb(t/Nv);
      j = ((aa / grid[1])% loc_lens[1]) + ((bb/grid[0]) % loc_lens[0]) * loc_lens[1];
    }
    int send_rank = glb_send_recv_list[2*t+1];
    char *ctf_data_pointer = ptensor.data + j * s_sources * sizeof(F);
    MPI_Isend(ctf_data_pointer,
              static_cast<int>(s_sources*sizeof(F)),
              MPI_CHAR,
              send_rank,
              t,
              world,
              &send_requests[i]);
  }

  MPI_Waitall(n_origins, recv_requests.data(), MPI_STATUSES_IGNORE);
  MPI_Waitall(n_targets, send_requests.data(), MPI_STATUSES_IGNORE);

  if (0) {
    for (auto &s: sources) {
      F* buffer = new F[s.size()];
      std::copy(s.begin(), s.end(), buffer);
      permute_copy(source_dimension, buffer, s.data());
    }
  }


  return {n_sources, s_sources, sources};
}


#endif /* defined(HAVE_CTF) */

template <typename F>
std::vector<F> read_all(std::vector<size_t> lengths,
                        std::string const &ctf_file_path,
                        MPI_Comm comm) {
  MPI_File handle;
  MPI_Offset offset = 0;
  const size_t count = std::accumulate(lengths.begin(),
                                       lengths.end(),
                                       1UL,
                                       std::multiplies<size_t>());
  std::vector<F> buffer(count);

  MPI_File_open(comm,
                ctf_file_path.c_str(),
                MPI_MODE_RDONLY,
                MPI_INFO_NULL,
                &handle);

  char *dest = reinterpret_cast<char *>(buffer.data());
  LOG(0, "Atrip") << "Reading file " << ctf_file_path << " from disk" << std::endl;
  if (MPI_SUCCESS
      != MPI_File_read_at(handle,
                          offset,
                          dest,
                          count * sizeof(F),
                          MPI_CHAR,
                          MPI_STATUS_IGNORE)) {
    throw "error reading!";
  }

  MPI_File_close(&handle);
  return buffer;
}






#define INSTANTIATE_DISK_READER(T) \
template Sources<T> diskReader<T>(const std::string&, \
                                  std::vector<size_t>, \
                                  std::vector<size_t>, \
                                  bool);


INSTANTIATE_DISK_READER(Complex)
INSTANTIATE_DISK_READER(double)
INSTANTIATE_DISK_READER(float)


#if defined(HAVE_CTF)

#define INSTANTIATE_CTF_READER(T) \
template Sources<T> ctfReader<T>(CTF::Tensor<T>&, \
                                 std::vector<size_t>, \
                                 std::vector<size_t>);


INSTANTIATE_CTF_READER(double)
INSTANTIATE_CTF_READER(float)
INSTANTIATE_CTF_READER(Complex)
#endif  /* HAVE CTF */

template std::vector<float> read_all<float>(std::vector<size_t> lengths,
                                            std::string const &ctf_file_path,
                                            MPI_Comm comm);
template std::vector<double> read_all<double>(std::vector<size_t> lengths,
                                              std::string const &ctf_file_path,
                                              MPI_Comm comm);
template std::vector<Complex>
read_all<Complex>(std::vector<size_t> lengths,
                  std::string const &ctf_file_path,
                  MPI_Comm comm);



template static void permute_copy(std::vector<size_t> ranges,
                                  const float* reorder_buffer,
                                  float* source_buffer);
template static void permute_copy(std::vector<size_t> ranges,
                                  const double* reorder_buffer,
                                  double* source_buffer);
template static void permute_copy(std::vector<size_t> ranges,
                                  const Complex* reorder_buffer,
                                  Complex* source_buffer);





template struct Sources<Complex>;
template struct Sources<double>;
template struct Sources<float>;
} // namespace atrip
