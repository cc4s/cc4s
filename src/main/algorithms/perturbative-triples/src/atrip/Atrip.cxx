// Copyright 2022 Alejandro Gallo
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// [[file:~/cuda/atrip/atrip.org::*Main][Main:1]]
#include <iomanip>

#include <atrip/Atrip.hpp>
#include <atrip/Utils.hpp>
#include <atrip/Equations.hpp>
#include <atrip/SliceUnion.hpp>
#include <atrip/Unions.hpp>
#include <atrip/Checkpoint.hpp>
#include <atrip/DatabaseCommunicator.hpp>
#include <atrip/Malloc.hpp>

#if defined(HAVE_CUDA)
#  include <nvToolsExt.h>
#  include <atrip/CUDA.hpp>
#endif

using namespace atrip;

template <typename F>
bool RankMap<F>::RANK_ROUND_ROBIN;
template bool RankMap<double>::RANK_ROUND_ROBIN;
template bool RankMap<Complex>::RANK_ROUND_ROBIN;
size_t Atrip::rank;
size_t Atrip::np;
ClusterInfo *Atrip::cluster_info;
#if defined(HAVE_CUDA)
typename Atrip::CudaContext Atrip::cuda;
typename Atrip::KernelDimensions Atrip::kernel_dimensions;
#endif
MPI_Comm Atrip::communicator;
Timings Atrip::chrono;
size_t Atrip::ppn;
size_t Atrip::network_send;
size_t Atrip::local_send;
double Atrip::bytes_sent;

// user printing block
IterationDescriptor IterationDescription::descriptor;
void atrip::register_iteration_descriptor(IterationDescriptor d) {
  IterationDescription::descriptor = d;
}

void Atrip::init(MPI_Comm world) {
  Atrip::communicator = world;
  MPI_Comm_rank(world, (int *)&Atrip::rank);
  MPI_Comm_size(world, (int *)&Atrip::np);
  Atrip::cluster_info = new ClusterInfo(get_cluster_info(world));
  Atrip::network_send = 0UL;
  Atrip::local_send = 0UL;
  Atrip::bytes_sent = 0.0;
  Atrip::ppn = Atrip::cluster_info->ranks_per_node;
}

template <typename F>
Atrip::Output Atrip::run(Atrip::Input<F> const &in) {

  const size_t np = Atrip::np;
  const size_t rank = Atrip::rank;
  MPI_Comm universe = Atrip::communicator;

  const size_t No = in.ei->lens[0];
  const size_t Nv = in.ea->lens[0];
  LOG(0, "Atrip") << "No: " << No << "\n";
  LOG(0, "Atrip") << "Nv: " << Nv << "\n";
  LOG(0, "Atrip") << "np: " << np << "\n";

#if defined(HAVE_CUDA)
  int ngcards;
  _CHECK_CUDA_SUCCESS("initializing cuda", cuInit(0));
  _CHECK_CUDA_SUCCESS("getting device count", cuDeviceGetCount(&ngcards));
  const auto cluster_info = *Atrip::cluster_info;
  LOG(0, "Atrip") << "ngcards: " << ngcards << "\n";
  if (cluster_info.ranks_per_node > ngcards) {
    const auto msg = _FORMAT(
        "ATRIP: You are running on more ranks per node than the number of "
        "graphic cards\n"
        "You have %d cards at your disposal\n",
        ngcards);
    std::cerr << msg;
    throw msg;
  } else if (cluster_info.ranks_per_node < ngcards) {
    const auto msg = _FORMAT(
        "You have %d cards at your disposal.\n"
        "You will be only using %d, i.e, the number of ranks\n",
        ngcards,
        cluster_info.ranks_per_node);
    std::cerr << msg;
  }

  for (size_t _rank = 0; _rank < np; _rank++) {
    if (rank == _rank) {
      CUcontext ctx;
      CUdevice dev;
      CUdevprop prop;
      struct {
        struct {
          size_t free, total;
        } avail;
        size_t total;
      } memory;
      char *name = (char *)malloc(256);

      // - TODO :: we should check that the Zuweisung of graphic cards
      //           to nodes works as expected, i.e., node k should get from 0
      //           to ngcards with the formula =rank % ngcards=.

      // set current device
      _CHECK_CUDA_SUCCESS("getting device for index <rank>",
                          cuDeviceGet(&dev, rank % ngcards));
      _CHECK_CUDA_SUCCESS("creating a cuda context", cuCtxCreate(&ctx, 0, dev));
      _CHECK_CUDA_SUCCESS("setting the context", cuCtxSetCurrent(ctx));

      // get information of the device
      _CHECK_CUDA_SUCCESS("getting  properties of current device",
                          cuDeviceGetProperties(&prop, dev));
      _CHECK_CUDA_SUCCESS(
          "getting memory information",
          cuMemGetInfo(&memory.avail.free, &memory.avail.total));
      _CHECK_CUDA_SUCCESS("getting name", cuDeviceGetName(name, 256, dev));
      _CHECK_CUDA_SUCCESS("getting total memory",
                          cuDeviceTotalMem(&memory.total, dev));

      printf(
          "\n"
          "CUDA CARD RANK %d\n"
          "=================\n"
          "\tnumber: %1$ld\n"
          "\tname: %s\n"
          "\tMem. clock rate (KHz): %ld\n"
          "\tShared Mem Per Block (KB): %f\n"
          "\tAvail. Free/Total mem (GB): %f/%f\n"
          "\tFree memory (GB): %f\n"
          "\n",
          Atrip::rank,
          name,
          prop.clockRate,
          prop.sharedMemPerBlock / 1024.0,
          memory.avail.free / 1024.0 / 1024.0 / 1024.0,
          memory.avail.total / 1024.0 / 1024.0 / 1024.0,
          memory.total / 1024.0 / 1024.0 / 1024.0);
      std::free((void *)name);

      _CHECK_CUBLAS_SUCCESS("creating a cublas handle",
                            cublasCreate(&Atrip::cuda.handle));
    }
    MPI_Barrier(universe);
  }

  if (in.ooo_threads > 0) {
    Atrip::kernel_dimensions.ooo.threads = in.ooo_threads;
  }
  if (in.ooo_blocks > 0) {
    Atrip::kernel_dimensions.ooo.blocks = in.ooo_blocks;
  }

  if (Atrip::kernel_dimensions.ooo.threads <= 0
      || Atrip::kernel_dimensions.ooo.blocks <= 0) {
    Atrip::kernel_dimensions.ooo.blocks = No / 32 + No % 32;
    Atrip::kernel_dimensions.ooo.threads = 32;
  }

  LOG(0, "Atrip") << "ooo blocks: " << Atrip::kernel_dimensions.ooo.blocks
                  << "\n";
  LOG(0, "Atrip") << "ooo threads per block: "
                  << Atrip::kernel_dimensions.ooo.threads << "\n";
#endif

  // allocate the three scratches, see piecuch
  // we need local copies of the following tensors on every
  // rank
  std::vector<F> _epsi(No), _epsa(Nv), _Tai(No * Nv);

  // copy the data from the tensors into the vectors
  in.ei->read_all(_epsi.data());
  in.ea->read_all(_epsa.data());
  in.Tph->read_all(_Tai.data());

  if (in.ijkabc) {
    std::transform(
      _Tai.begin(), _Tai.end(), _Tai.begin(), [](F& c){return -c;}
    );
  }
  // TODO: free memory pointers in the end of the algorithm
  DataPtr<F> Tijk, Zijk;

#if defined(HAVE_CUDA)
  DataPtr<F> Tai, epsi, epsa;

  // TODO: free memory pointers in the end of the algorithm

  _CHECK_CUDA_SUCCESS("Tai", cuMemAlloc(&Tai, sizeof(F) * _Tai.size()));
  _CHECK_CUDA_SUCCESS("epsi", cuMemAlloc(&epsi, sizeof(F) * _epsi.size()));
  _CHECK_CUDA_SUCCESS("epsa", cuMemAlloc(&epsa, sizeof(F) * _epsa.size()));

  _CHECK_CUDA_SUCCESS(
      "memcpy Tai",
      cuMemcpyHtoD(Tai, (void *)_Tai.data(), sizeof(F) * _Tai.size()));
  _CHECK_CUDA_SUCCESS(
      "memcpy epsi",
      cuMemcpyHtoD(epsi, (void *)_epsi.data(), sizeof(F) * _epsi.size()));
  _CHECK_CUDA_SUCCESS(
      "memcpy epsa",
      cuMemcpyHtoD(epsa, (void *)_epsa.data(), sizeof(F) * _epsa.size()));

#else
  DataPtr<F> Tai = _Tai.data(), epsi = _epsi.data();
#endif

  MALLOC_DATA_PTR("Zijk", &Zijk, sizeof(DataFieldType<F>) * No * No * No);
  MALLOC_DATA_PTR("Tijk", &Tijk, sizeof(DataFieldType<F>) * No * No * No);

  RankMap<F>::RANK_ROUND_ROBIN = in.rank_round_robin;
  if (RankMap<F>::RANK_ROUND_ROBIN) {
    LOG(0, "Atrip") << "Doing rank round robin slices distribution\n";
  } else {
    LOG(0, "Atrip")
        << "Doing node > local rank round robin slices distribution\n";
  }

  // COMMUNICATOR CONSTRUCTION ========================================={{{1
  //
  // Construct a new communicator living only on a single rank
  int child_size = 1, child_rank;
  const int color = rank / child_size, crank = rank % child_size;
  MPI_Comm child_comm;
  if (np == 1) {
    child_comm = universe;
  } else {
    MPI_Comm_split(universe, color, crank, &child_comm);
    MPI_Comm_rank(child_comm, &child_rank);
    MPI_Comm_size(child_comm, &child_size);
  }

  // a, b, c, d, e, f and P => Nv
  // H                      => No
  // total_source_sizes contains a list of the number of elements
  // in all sources of every tensor union, therefore n_slices * slice_size

  // TODO: remove the alignment of sources in one big block or
  //        add it as an option

  // const std::vector<size_t> total_source_sizes = {
  //     // ABPH
  //     SliceUnion<F>::get_size({Nv, No}, {Nv, Nv}, (size_t)np, universe),
  //     // ABHH
  //     SliceUnion<F>::get_size({No, No}, {Nv, Nv}, (size_t)np, universe),
  //     // TABHH
  //     SliceUnion<F>::get_size({No, No}, {Nv, Nv}, (size_t)np, universe),
  //     // TAPHH
  //     SliceUnion<F>::get_size({Nv, No, No}, {Nv}, (size_t)np, universe),
  //     // HHHA
  //     SliceUnion<F>::get_size({No, No, No}, {Nv}, (size_t)np, universe),
  // };

  // const size_t total_source_size = sizeof(DataFieldType<F>)
  //                                * std::accumulate(total_source_sizes.begin(),
  //                                                  total_source_sizes.end(),
  //                                                  0UL);

  // #if defined(HAVE_CUDA)
  //   DataPtr<F> all_sources_pointer;
  //   cuMemAlloc(&all_sources_pointer, total_source_size);
  // #else
  //   DataPtr<F> all_sources_pointer = (DataPtr<F>)malloc(total_source_size);
  // #endif
  //   size_t _source_pointer_idx = 0;

  // BUILD SLICES PARAMETRIZED BY NV x NV =============================={{{1
  WITH_CHRONO(
      "nv-nv-slices", LOG(0, "Atrip") << "building NV x NV slices\n";
      // TODO
      // DataPtr<F> offseted_pointer = all_sources_pointer
      //                             * total_source_sizes[_source_pointer_idx++];
      ABPH<F> abph(*in.Vppph,
                   (size_t)No,
                   (size_t)Nv,
                   (size_t)np,
                   child_comm,
                   universe);

      // TODO
      // DataPtr<F> offseted_pointer = all_sources_pointer
      //                             * total_source_sizes[_source_pointer_idx++];
      ABHH<F> abhh(*in.Vpphh,
                   (size_t)No,
                   (size_t)Nv,
                   (size_t)np,
                   child_comm,
                   universe);

      // TODO
      // DataPtr<F> offseted_pointer = all_sources_pointer
      //                             * total_source_sizes[_source_pointer_idx++];
      TABHH<F> tabhh(*in.Tpphh,
                     (size_t)No,
                     (size_t)Nv,
                     (size_t)np,
                     child_comm,
                     universe);)

  // delete the Vppph so that we don't have a HWM situation for the NV slices
  if (in.delete_Vppph) { delete in.Vppph; }

  // BUILD SLICES PARAMETRIZED BY NV ==================================={{{1
  WITH_CHRONO(
      "nv-slices", LOG(0, "Atrip") << "building NV slices\n";
      // TODO
      // DataPtr<F> offseted_pointer = all_sources_pointer
      //                             * total_source_sizes[_source_pointer_idx++];
      TAPHH<F> taphh(*in.Tpphh,
                     (size_t)No,
                     (size_t)Nv,
                     (size_t)np,
                     child_comm,
                     universe);
      // TODO
      // DataPtr<F> offseted_pointer = all_sources_pointer
      //                             * total_source_sizes[_source_pointer_idx++];
      HHHA<F> hhha(*in.Vhhhp,
                   (size_t)No,
                   (size_t)Nv,
                   (size_t)np,
                   child_comm,
                   universe);)

  // all tensors
  std::vector<SliceUnion<F> *> unions = {&taphh, &hhha, &abph, &abhh, &tabhh};

  // IF (cT) IS USED: HANDLE TWO FURTHER SLICES==========================={{{1
  HHHA<F> *jhhha = nullptr;
  ABPH<F> *jabph = nullptr;
  if (in.Jhhhp) {
    WITH_CHRONO("Jhhha-slice", LOG(0, "Atrip") << "slicing Jijka" << std::endl;
                jhhha = new HHHA<F>(*in.Jhhhp,
                                    (size_t)No,
                                    (size_t)Nv,
                                    (size_t)np,
                                    child_comm,
                                    universe);)
    jhhha->name = Slice<F>::Name::JIJKA;
    unions.push_back(jhhha);
  }

  if (in.Jppph) {
    WITH_CHRONO("Jabph-slice", LOG(0, "Atrip") << "slicing Jabci" << std::endl;
                jabph = new ABPH<F>(*in.Jppph,
                                    (size_t)No,
                                    (size_t)Nv,
                                    (size_t)np,
                                    child_comm,
                                    universe);)
    jabph->name = Slice<F>::Name::JABCI;
    unions.push_back(jabph);
  }

  DataFieldType<F> *_t_buffer, *_vhhh;
  MALLOC_DATA_PTR("Tijk buffer for doubles equations",
                  &_t_buffer,
                  sizeof(DataFieldType<F>) * No * No * No);
  MALLOC_DATA_PTR("Vijk buffer for doubles equations",
                  &_vhhh,
                  sizeof(DataFieldType<F>) * No * No * No);

#ifdef HAVE_CUDA
  WITH_CHRONO(
      "double:cuda:alloc",
      _CHECK_CUDA_SUCCESS("Allocating _t_buffer",
                          cuMemAlloc((CUdeviceptr *)&_t_buffer,
                                     No * No * No * sizeof(DataFieldType<F>)));
      _CHECK_CUDA_SUCCESS("Allocating _vhhh",
                          cuMemAlloc((CUdeviceptr *)&_vhhh,
                                     No * No * No * sizeof(DataFieldType<F>)));)
#endif

  // get tuples for the current rank
  TuplesDistribution *distribution;

  if (in.tuples_distribution == Atrip::Input<F>::TuplesDistribution::NAIVE) {
    LOG(0, "Atrip") << "Using the naive distribution\n";
    distribution = new NaiveDistribution();
  } else {
    LOG(0, "Atrip") << "Using the group-and-sort distribution\n";
    distribution = new group_and_sort::Distribution();
  }

  LOG(0, "Atrip") << "BUILDING TUPLE LIST\n";
  WITH_CHRONO("tuples:build",
              auto const tuples_list = distribution->get_tuples(Nv, universe);)
  const size_t n_iterations = tuples_list.size();
  {
    LOG(0, "Atrip") << "#iterations: " << n_iterations << "/"
                    << n_iterations * np << "\n";
  }

  const size_t iteration_mod = (in.percentage_mod > 0)
                                 ? n_iterations * in.percentage_mod / 100.0
                                 : in.iteration_mod,
               iteration1Percent = n_iterations * 0.01;

  auto const is_fake_tuple = [&tuples_list, distribution](size_t const i) {
    return distribution->tuple_is_fake(tuples_list[i]);
  };

  double db_last_iteration_time = 0.0;

  using Database = typename Slice<F>::Database;
  auto communicate_database =
      [&unions, &in, Nv, np](ABCTuple const &abc,
                             MPI_Comm const &c,
                             size_t iteration) -> Database {
    if (in.tuples_distribution == Atrip::Input<F>::TuplesDistribution::NAIVE) {

      WITH_CHRONO("db:comm:naive",
                  auto const &db =
                      naive_database<F>(unions, Nv, np, iteration);)
      return db;

    } else {
      WITH_CHRONO("db:comm:type:do",
                  auto MPI_LDB_ELEMENT =
                      Slice<F>::mpi::local_database_element();)

      WITH_CHRONO(
          "db:comm:ldb", typename Slice<F>::LocalDatabase ldb;
          for (auto const &tensor
               : unions) {
            auto const &tensor_db = tensor->build_local_database(abc);
            ldb.insert(ldb.end(), tensor_db.begin(), tensor_db.end());
          })

      Database db(np * ldb.size(), ldb[0]);

      WITH_CHRONO(
          "oneshot-db:comm:allgather",
          WITH_CHRONO("db:comm:allgather",
                      MPI_Allgather(ldb.data(),
                                    /* ldb.size() * sizeof(typename
                                       Slice<F>::LocalDatabaseElement) */
                                    ldb.size(),
                                    MPI_LDB_ELEMENT,
                                    db.data(),
                                    /* ldb.size() * sizeof(typename
                                       Slice<F>::LocalDatabaseElement), */
                                    ldb.size(),
                                    MPI_LDB_ELEMENT,
                                    c);))

      WITH_CHRONO("db:comm:type:free", MPI_Type_free(&MPI_LDB_ELEMENT);)

      return db;
    }
  };

  auto do_io_phase =
      [&unions, &rank, &np, No, Nv, &universe, &db_last_iteration_time](
          Database const &db,
          ABCTuple const abc,
          size_t iteration) {
        IGNORABLE(iteration); // iteration used to print database
        const size_t localDBLength = db.size() / np;

        size_t send_tag = 0, recv_tag = rank * localDBLength;

        // RECIEVE PHASE ======================================================
        {
          // At this point, we have already send to everyone that fits
          auto const &begin = &db[rank * localDBLength],
                     end = begin + localDBLength;
          for (auto it = begin; it != end; ++it) {
            recv_tag++;
            auto const &el = *it;
            auto &u = union_by_name(unions, el.name);

            WITH_DBG std::cout << rank << ":r"
                               << "♯" << recv_tag << " =>"
                               << " «n" << el.name << ", t" << el.info.type
                               << ", s" << el.info.state << "»"
                               << " ⊙ {" << rank << "⇐" << el.info.from.rank
                               << ", " << el.info.from.source << "}"
                               << " ∴ {" << el.info.tuple[0] << ", "
                               << el.info.tuple[1] << "}"
                               << "\n";

            WITH_CHRONO("db:io:recv", u.receive(el.info, recv_tag);)

#if defined(ATRIP_PRINT_DB)
            if (!Atrip::rank)
              std::cout << _FORMAT(
                  "%4s %d %d %d %5d %5s %2s %14s (%ld,%ld) %ld %ld %f\n",
                  "RECV",
                  iteration,
                  el.info.from.rank,
                  rank,
                  recv_tag,
                  name_to_string<double>(el.name).c_str(),
                  type_to_string<double>(el.info.type).c_str(),
                  state_to_string<double>(el.info.state).c_str(),
                  el.info.tuple[0],
                  el.info.tuple[1],
                  name_to_size<double>(el.name, No, Nv),
                  u.free_pointers.size(),
#  if defined(ATRIP_MPI_STAGING_BUFFERS)
                  u.mpi_staging_buffers.size(),
#  else
                  0,
#  endif /* defined(ATRIP_MPI_STAGING_BUFFERS) */
                  db_last_iteration_time);
#endif /* defined(ATRIP_PRINT_DB) */

          } // recv
        }

        MPI_Barrier(universe);

        // SEND PHASE =========================================================
        for (size_t other_rank = 0; other_rank < np; other_rank++) {
          auto const &begin = &db[other_rank * localDBLength],
                     end = begin + localDBLength;
          for (auto it = begin; it != end; ++it) {
            send_tag++;
            typename Slice<F>::LocalDatabaseElement const &el = *it;

            if (el.info.from.rank != rank) continue;

            auto &u = union_by_name(unions, el.name);
            WITH_DBG std::cout
                << rank << ":s"
                << "♯" << send_tag << " =>"
                << " «n" << el.name << ", t" << el.info.type << ", s"
                << el.info.state << "»"
                << " ⊙ {" << el.info.from.rank << "⇒" << other_rank << ", "
                << el.info.from.source << "}"
                << " ∴ {" << el.info.tuple[0] << ", " << el.info.tuple[1] << "}"
                << "\n";

#if defined(ATRIP_PRINT_DB)
            if (!Atrip::rank)
              std::cout << _FORMAT(
                  "%4s %d %d %d %5d %5s %2s %14s (%ld,%ld) %ld %ld %f\n",
                  "SEND",
                  iteration,
                  el.info.from.rank,
                  other_rank,
                  send_tag,
                  name_to_string<double>(el.name).c_str(),
                  type_to_string<double>(el.info.type).c_str(),
                  state_to_string<double>(el.info.state).c_str(),
                  el.info.tuple[0],
                  el.info.tuple[1],
                  name_to_size<double>(el.name, No, Nv),
                  u.free_pointers.size(),
#  if defined(ATRIP_MPI_STAGING_BUFFERS)
                  u.mpi_staging_buffers.size(),
#  else
                  0,
#  endif /* defined(ATRIP_MPI_STAGING_BUFFERS) */
                  db_last_iteration_time);
#endif /* defined(ATRIP_PRINT_DB) */

            WITH_CHRONO("db:io:send", u.send(other_rank, el, send_tag, abc);)

          } // send phase

        } // other_rank
      };

#if defined(HAVE_OCD) || defined(ATRIP_PRINT_TUPLES)
  std::map<ABCTuple, double> tuple_energies;
#endif

  const double doubles_flops =
      double(No) * double(No) * double(No) * (double(No) + double(Nv)) * 2.0
      * (traits::is_complex<F>() ? 4.0 : 1.0) * 6.0 / 1e9;

  // START MAIN LOOP ======================================================{{{1

  MPI_Barrier(universe);
  Output local_output = {0, 0}, global_output = {0, 0};
  size_t first_iteration = 0;
  Checkpoint c;
  const size_t checkpoint_mod =
      in.checkpoint_at_every_iteration != 0
          ? in.checkpoint_at_every_iteration
          : n_iterations * in.checkpoint_at_percentage / 100;
  if (in.read_checkpoint_if_exists) {
    std::ifstream fin(in.checkpoint_path);
    if (fin.is_open()) {
      LOG(0, "Atrip") << "Reading checkpoint from " << in.checkpoint_path
                      << "\n";
      c = read_checkpoint(fin);
      first_iteration = (size_t)c.iteration;
      if (first_iteration > n_iterations) {
        // TODO: throw an error here
        // first_iteration is bigger than n_iterations,
        // you probably started the program with a different number
        // of cores
      }
      if (No != c.no) { /* TODO: write warning */
      }
      if (Nv != c.nv) { /* TODO: write warning */
      }
      // TODO write warnings for nrank and so on
      if (Atrip::rank == 0) {
        // take the negative of the energy to correct for the
        // negativity of the equations, the energy in the checkpoint
        // should always be the correct physical one.
        local_output.energy = -(double)c.energy;
      }
      LOG(0, "Atrip") << "energy from checkpoint " << local_output.energy
                      << "\n";
      LOG(0, "Atrip") << "iteration from checkpoint " << first_iteration
                      << "\n";
    }
  }

  const auto compute_local_energy =
      [Tijk, Zijk, No, epsi, &_epsa](ABCTuple const &abc,
                                     bool const fake_tuple_p) {
#if defined(ATRIP_ONLY_DGEMM)
        if (false)
#endif /* defined(ATRIP_ONLY_DGEMM) */
          if (!fake_tuple_p) {
#if defined(HAVE_CUDA)
            double *tuple_energy;
            cuMemAlloc((DataPtr<double> *)&tuple_energy, sizeof(double));
#else
          double _tuple_energy(0.);
          double *tuple_energy = &_tuple_energy;
#endif /* defined(HAVE_CUDA) */

            int distinct(0);
            if (abc[0] == abc[1]) distinct++;
            if (abc[1] == abc[2]) distinct--;
            const double epsabc =
                std::real(_epsa[abc[0]] + _epsa[abc[1]] + _epsa[abc[2]]);

            DataFieldType<F> _epsabc{epsabc};

            WITH_CHRONO(
                "energy",
                if (distinct == 0) {
                  ACC_FUNCALL(get_energy_distinct<DataFieldType<F>>,
                              1,
                              1, // for cuda
                              _epsabc,
                              No,
                              (DataFieldType<F> *)epsi,
                              (DataFieldType<F> *)Tijk,
                              (DataFieldType<F> *)Zijk,
                              tuple_energy);
                } else {
                  ACC_FUNCALL(get_energy_same<DataFieldType<F>>,
                              1,
                              1, // for cuda
                              _epsabc,
                              No,
                              (DataFieldType<F> *)epsi,
                              (DataFieldType<F> *)Tijk,
                              (DataFieldType<F> *)Zijk,
                              tuple_energy);
                })

#if defined(HAVE_CUDA)
            double host_tuple_energy;
            cuMemcpyDtoH((void *)&host_tuple_energy,
                         (DataPtr<double>)tuple_energy,
                         sizeof(double));
#else
          double host_tuple_energy = *tuple_energy;
#endif /* defined(HAVE_CUDA) */

            return host_tuple_energy;
          }
        return 0.0;
      };

  for (size_t i = first_iteration, iteration = first_iteration + 1;
       i < tuples_list.size();
       i++, iteration++) {
    Atrip::chrono["iterations"].start();
    Atrip::chrono["db-last-iteration"].start();

#if defined(HAVE_CUDA)
    char nvtx_name[60];
    sprintf(nvtx_name, "iteration: %d", i);
    nvtxRangePushA(nvtx_name);
#endif /* defined(HAVE_CUDA) */

    // check overhead from chrono over all iterations
    WITH_CHRONO("start:stop", {})

    // check overhead of doing a barrier at the beginning
    WITH_CHRONO(
        "oneshot-mpi:barrier",
        WITH_CHRONO("mpi:barrier", if (in.barrier) MPI_Barrier(universe);))

    // write checkpoints
    // TODO: ENABLE THIS
    if (iteration % checkpoint_mod == 0 && false) {
      double global_energy = 0;
      MPI_Reduce(&local_output.energy,
                 &global_energy,
                 1,
                 MPI_DOUBLE,
                 MPI_SUM,
                 0,
                 universe);
      Checkpoint out = {No,
                        Nv,
                        Atrip::cluster_info->ranks_per_node,
                        Atrip::cluster_info->n_nodes,
                        -global_energy,
                        iteration - 1,
                        in.rank_round_robin};
      LOG(0, "Atrip") << "Writing checkpoint\n";
      if (Atrip::rank == 0) write_checkpoint(out, in.checkpoint_path);
    }

    // write reporting
    if (iteration % iteration_mod == 0 || iteration == iteration1Percent) {

      if (IterationDescription::descriptor) {
        IterationDescription::descriptor(
            {iteration, n_iterations, Atrip::chrono["iterations"].count()});
      }

      const double _doubles_time = Atrip::chrono["doubles"].count(),
                   _its_time = Atrip::chrono["iterations"].count();

      size_t network_send(0);
      MPI_Reduce(&Atrip::network_send,
                 &network_send,
                 1,
                 MPI_UINT64_T,
                 MPI_SUM,
                 0,
                 universe);

      size_t local_send(0);
      MPI_Reduce(&Atrip::local_send,
                 &local_send,
                 1,
                 MPI_UINT64_T,
                 MPI_SUM,
                 0,
                 universe);

      double bytes_sent(0.0);
      MPI_Reduce(&Atrip::bytes_sent,
                 &bytes_sent,
                 1,
                 MPI_DOUBLE,
                 MPI_SUM,
                 0,
                 universe);

      const size_t total_send = network_send + local_send;

      LOG(0, "Atrip")
          << "iteration " << iteration << " [" << 100 * iteration / n_iterations
          << "%]"
          << " ("
          << (_doubles_time > 0.0 ? doubles_flops * iteration / _doubles_time
                                  : -1)
          << "GF)"
          << " ("
          << (_its_time > 0.0 ? doubles_flops * iteration / _its_time : -1)
          << "GF) :: GB sent per rank: " << bytes_sent / 1073741824.0
          << " :: "
          << (total_send > 0UL ? (double) network_send / total_send : 0UL)
          << " % network communication" << std::endl;

      // PRINT TIMINGS
      if (in.chrono)
        for (auto const &pair : Atrip::chrono)
          LOG(1, " ") << pair.first << " :: " << pair.second.count()
                      << std::endl;
    }

    const ABCTuple abc = is_fake_tuple(i) ? tuples_list[tuples_list.size() - 1]
                                          : tuples_list[i],
                   *abc_next = i == (tuples_list.size() - 1)
                                 ? nullptr
                                 : &tuples_list[i + 1];

    WITH_CHRONO("with_rank",
                WITH_RANK << " :it " << iteration << " :abc "
                          << pretty_print(abc) << " :abcN "
                          << (abc_next ? pretty_print(*abc_next) : "None")
                          << "\n";)

    // COMM FIRST DATABASE
    // ================================================{{{1
    if (i == first_iteration) {
      WITH_RANK << "__first__:first database ............ \n";
      const auto db = communicate_database(abc, universe, i);
      WITH_RANK << "__first__:first database communicated \n";
      WITH_RANK << "__first__:first database io phase \n";
      do_io_phase(db, abc, i);
      WITH_RANK << "__first__:first database io phase DONE\n";
      WITH_RANK << "__first__::::Unwrapping all slices for first database\n";
      for (auto &u : unions) u->unwrap_all(abc);
      WITH_RANK << "__first__::::Unwrapping slices for first database DONE\n";
      MPI_Barrier(universe);
    }

    // COMM NEXT DATABASE
    // ================================================={{{1
    if (abc_next) {
      WITH_RANK << "__comm__:" << iteration << "th communicating database\n";
      WITH_CHRONO("db:comm",
                  const auto db = communicate_database(*abc_next, universe, i);)
      WITH_CHRONO("db:io", do_io_phase(db, abc, i + 1);)
      WITH_RANK << "__comm__:" << iteration << "th database io phase DONE\n";
    }

    // COMPUTE DOUBLES
    // ===================================================={{{1
    OCD_Barrier(universe);
    if (!is_fake_tuple(i)) {
      WITH_RANK << iteration << "-th doubles\n";
      WITH_CHRONO(
          "oneshot-unwrap",
          WITH_CHRONO(
              "unwrap",
              WITH_CHRONO(
                  "unwrap:doubles",
                  for (auto &u
                       : decltype(unions){&abph, &hhha, &taphh, &tabhh}) {
                    u->unwrap_all(abc);
                  })))
      if (in.blocking && abc_next) {
        WITH_CHRONO(
            "blockingCommunication",
            for (auto &u
                 : decltype(unions){&abph, &hhha, &taphh, &tabhh}) {
              u->unwrap_all(*abc_next);
            })
        WITH_CHRONO("blocking-barrier", MPI_Barrier(universe);)
      }

      WITH_CHRONO("oneshot-doubles",
                  WITH_CHRONO("doubles",
                              doubles_contribution<F>(
                                  (size_t)No,
                                  (size_t)Nv,
                                  // -- VABCI
                                  abph.unwrap_slice(Slice<F>::AB, abc),
                                  abph.unwrap_slice(Slice<F>::AC, abc),
                                  abph.unwrap_slice(Slice<F>::BC, abc),
                                  abph.unwrap_slice(Slice<F>::BA, abc),
                                  abph.unwrap_slice(Slice<F>::CA, abc),
                                  abph.unwrap_slice(Slice<F>::CB, abc),
                                  // -- VHHHA,
                                  hhha.unwrap_slice(Slice<F>::A, abc),
                                  hhha.unwrap_slice(Slice<F>::B, abc),
                                  hhha.unwrap_slice(Slice<F>::C, abc),
                                  // -- TA,
                                  taphh.unwrap_slice(Slice<F>::A, abc),
                                  taphh.unwrap_slice(Slice<F>::B, abc),
                                  taphh.unwrap_slice(Slice<F>::C, abc),
                                  // -- TABIJ
                                  tabhh.unwrap_slice(Slice<F>::AB, abc),
                                  tabhh.unwrap_slice(Slice<F>::AC, abc),
                                  tabhh.unwrap_slice(Slice<F>::BC, abc),
                                  // -- TIJK
                                  (DataFieldType<F> *)Tijk,
                                  // -- tmp buffers
                                  (DataFieldType<F> *)_t_buffer,
                                  (DataFieldType<F> *)_vhhh);

                              WITH_RANK << iteration << "-th doubles done\n";))
    }

    // COMPUTE SINGLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // {{{1
    OCD_Barrier(universe);
#if defined(ATRIP_ONLY_DGEMM)
    if (false)
#endif
      if (!is_fake_tuple(i)) {
        WITH_CHRONO(
            "oneshot-unwrap",
            WITH_CHRONO("unwrap",
                        WITH_CHRONO("unwrap:singles", abhh.unwrap_all(abc);)))
        WITH_CHRONO("reorder", /**/
                    int ooo = No * No * No,
                    stride = 1;
                    atrip::xcopy<F>(&ooo,
                                    (DataFieldType<F> *)Tijk,
                                    &stride,
                                    (DataFieldType<F> *)Zijk,
                                    &stride);)
        WITH_CHRONO(
            "singles",
            ACC_FUNCALL(
                singles_contribution<F>,
                1, // gpu
                1, // gpu
                No,
                Nv,
                abc[0],
                abc[1],
                abc[2],
                (DataFieldType<F> *)Tai,
                (DataFieldType<F> *)abhh.unwrap_slice(Slice<F>::AB, abc),
                (DataFieldType<F> *)abhh.unwrap_slice(Slice<F>::AC, abc),
                (DataFieldType<F> *)abhh.unwrap_slice(Slice<F>::BC, abc),
                (DataFieldType<F> *)Zijk);)
      }

    local_output.energy += compute_local_energy(abc, is_fake_tuple(i));

    // COMPUTE (cT) DOUBLES WITH THE J-INTERMEDIATE%%%%%%%%%%%%%%%%%%%%%
    if (!is_fake_tuple(i) && jhhha && jabph) {

      WITH_CHRONO("oneshot-doubles-J",
                  WITH_CHRONO("doubles-J",
                              doubles_contribution<F>(
                                  (size_t)No,
                                  (size_t)Nv,
                                  // -- VABCI
                                  jabph->unwrap_slice(Slice<F>::AB, abc),
                                  jabph->unwrap_slice(Slice<F>::AC, abc),
                                  jabph->unwrap_slice(Slice<F>::BC, abc),
                                  jabph->unwrap_slice(Slice<F>::BA, abc),
                                  jabph->unwrap_slice(Slice<F>::CA, abc),
                                  jabph->unwrap_slice(Slice<F>::CB, abc),
                                  // -- VHHHA,
                                  jhhha->unwrap_slice(Slice<F>::A, abc),
                                  jhhha->unwrap_slice(Slice<F>::B, abc),
                                  jhhha->unwrap_slice(Slice<F>::C, abc),
                                  // -- TA,
                                  taphh.unwrap_slice(Slice<F>::A, abc),
                                  taphh.unwrap_slice(Slice<F>::B, abc),
                                  taphh.unwrap_slice(Slice<F>::C, abc),
                                  // -- TABIJ
                                  tabhh.unwrap_slice(Slice<F>::AB, abc),
                                  tabhh.unwrap_slice(Slice<F>::AC, abc),
                                  tabhh.unwrap_slice(Slice<F>::BC, abc),
                                  // -- TIJK
                                  (DataFieldType<F> *)Tijk,
                                  // -- tmp buffers
                                  (DataFieldType<F> *)_t_buffer,
                                  (DataFieldType<F> *)_vhhh);

                              WITH_RANK << iteration << "-th doubles done\n";))
    }

    local_output.ct_energy += compute_local_energy(abc, is_fake_tuple(i));

    // TODO: remove this
    if (is_fake_tuple(i)) {
      // fake iterations should also unwrap whatever they got
      WITH_RANK << iteration << "th unwrapping because of fake in " << i
                << "\n";
      for (auto &u : unions) u->unwrap_all(abc);
    }

#ifdef HAVE_OCD
    for (auto const &u : unions) {
      WITH_RANK << "__dups__:" << iteration << "-th n" << u->name
                << " checking duplicates\n";
      u->check_for_duplicates();
    }
#endif

    // CLEANUP UNIONS
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%{{{1
    OCD_Barrier(universe);
    if (abc_next) {
      WITH_RANK << "__gc__:" << iteration << "-th cleaning up.......\n";
      for (auto &u : unions) {

        u->unwrap_all(abc);
        WITH_RANK << "__gc__:n" << u->name << " :it " << iteration << " :abc "
                  << pretty_print(abc) << " :abcN " << pretty_print(*abc_next)
                  << "\n";
        // for (auto const& slice: u->slices)
        //   WITH_RANK << "__gc__:guts:" << slice.info << "\n";
        u->clear_unused_slices_for_next_tuple(*abc_next);

        WITH_RANK << "__gc__: checking validity\n";

#ifdef HAVE_OCD
        // check for validity of the slices
        for (auto type : u->slice_types) {
          auto tuple = Slice<F>::subtuple_by_slice(abc, type);
          for (auto &slice : u->slices) {
            if (slice.info.type == type && slice.info.tuple == tuple
                && slice.is_directly_fetchable()) {
              if (slice.info.state == Slice<F>::Dispatched)
                throw std::domain_error(
                    "This slice should not be undispatched! "
                    + pretty_print(slice.info));
            }
          }
        }
#endif
      }
    }

#if defined(ATRIP_MPI_STAGING_BUFFERS)
    // Cleanup mpi staging buffers
    for (auto &u : unions) {
      std::vector<typename SliceUnion<F>::StagingBufferInfo> to_erase;
      for (auto &i : u->mpi_staging_buffers) {
        int completed = i.abc == abc;
        if (completed) { to_erase.push_back(i); }
      }
      for (auto &i : to_erase) {
        u->free_pointers.insert(i.data);
        u->mpi_staging_buffers.erase(i);
      }
    }
#endif /* defined(ATRIP_MPI_STAGING_BUFFERS) */

    WITH_RANK << iteration << "-th cleaning up....... DONE\n";

    Atrip::chrono["iterations"].stop();
    // ITERATION END
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%{{{1

    // AMB: debugging only
    WITH_CHRONO("mpi:barrier", MPI_Barrier(universe););
#if defined(HAVE_CUDA)
    cudaDeviceSynchronize();
    nvtxRangePop();
#endif

    Atrip::chrono["db-last-iteration"].stop();
    db_last_iteration_time = Atrip::chrono["db-last-iteration"].count();
    Atrip::chrono["db-last-iteration"].clear();
    if (in.max_iterations != 0 && i >= in.max_iterations) {
      if (abc_next)
        for (auto &u : unions) u->unwrap_all(*abc_next);
      break;
    }
  }
  // END OF MAIN LOOP

#if defined(HAVE_CUDA)
  _CUDA_FREE(Tai);
  _CUDA_FREE(epsi);
  _CUDA_FREE(epsa);
#endif

  FREE_DATA_PTR("Zijk", Zijk);
  FREE_DATA_PTR("Tijk", Tijk);
  FREE_DATA_PTR("Doubles Tijk buffer", _t_buffer);
  FREE_DATA_PTR("Doubles Vijk buffer", _vhhh);

  if (jhhha) delete jhhha;
  if (jabph) delete jabph;
  MPI_Barrier(universe);

  // PRINT TUPLES
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%{{{1
#if defined(HAVE_OCD) || defined(ATRIP_PRINT_TUPLES)
  LOG(0, "Atrip") << "tuple energies"
                  << "\n";
  for (size_t i = 0; i < np; i++) {
    MPI_Barrier(universe);
    for (auto const &pair : tuple_energies) {
      if (i == rank)
        std::cout << pair.first[0] << " " << pair.first[1] << " "
                  << pair.first[2] << std::setprecision(15) << std::setw(23)
                  << " tuple_energy: " << pair.second << "\n";
    }
  }
#endif

  // COMMUNICATE THE ENERGIES
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%{{{1
  LOG(0, "Atrip") << "COMMUNICATING ENERGIES \n";
  MPI_Reduce(&local_output.energy,
             &global_output.energy,
             1,
             MPI_DOUBLE,
             MPI_SUM,
             0,
             universe);
  MPI_Reduce(&local_output.ct_energy,
             &global_output.ct_energy,
             1,
             MPI_DOUBLE,
             MPI_SUM,
             0,
             universe);
  if (!in.ijkabc) {
    global_output.energy = -global_output.energy;
    global_output.ct_energy = -global_output.ct_energy;
  }
  WITH_RANK << "local energy " << local_output.energy << "\n";
  LOG(0, "Atrip") << "Energy: " << std::setprecision(15) << std::setw(23)
                  << global_output.energy << std::endl;
  LOG(0, "Atrip") << "Energy (cT): " << std::setprecision(15) << std::setw(23)
                  << global_output.ct_energy << std::endl;

  // PRINT TIMINGS {{{1
  if (in.chrono)
    for (auto const &pair : Atrip::chrono)
      LOG(0, "atrip:chrono")
          << pair.first << " " << pair.second.count() << std::endl;

  LOG(0, "atrip:flops(doubles)")
      << n_iterations * doubles_flops / Atrip::chrono["doubles"].count()
      << "\n";
  LOG(0, "atrip:flops(iterations)")
      << n_iterations * doubles_flops / Atrip::chrono["iterations"].count()
      << "\n";

  // TODO: change the sign in  the getEnergy routines
  return global_output;
}
// instantiate
template Atrip::Output Atrip::run(Atrip::Input<double> const &in);
template Atrip::Output Atrip::run(Atrip::Input<Complex> const &in);
// Main:1 ends here
