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
#include <atrip/Checkpoint.hpp>
#include <atrip/DatabaseCommunicator.hpp>
#include <atrip/Malloc.hpp>
#include <atrip/Acc.hpp>

using namespace atrip;

int Atrip::rank = MPI_UNDEFINED;
int Atrip::np = 0;
int Atrip::omp_threads;
int Atrip::omp_granularity;
bool Atrip::rank_round_robin;

#if defined(HAVE_ACC)
typename Atrip::CudaContext Atrip::cuda;
typename Atrip::KernelDimensions Atrip::kernel_dimensions;
#endif
MPI_Comm Atrip::communicator;
Timings Atrip::chrono;
int Atrip::ranks_per_node;
int Atrip::n_nodes;
int Atrip::node_id;
int Atrip::local_rank;
std::vector<int> Atrip::node_ids;
size_t Atrip::network_send;
size_t Atrip::local_send;
double Atrip::bytes_sent;

// user printing block
IterationDescriptor IterationDescription::descriptor;
void atrip::register_iteration_descriptor(IterationDescriptor d) {
  IterationDescription::descriptor = d;
}

template <typename F>
struct LocalOutput {
  F energy;
  F ct_energy;
};

void Atrip::init(MPI_Comm atrip_world, int omp_granularity_) {
  Atrip::communicator = atrip_world;
  MPI_Comm_rank(Atrip::communicator, &Atrip::rank);
  MPI_Comm_size(Atrip::communicator, &Atrip::np);
  Atrip::network_send = 0UL;
  Atrip::local_send = 0UL;
  Atrip::bytes_sent = 0.0;
  Atrip::omp_granularity = omp_granularity_;

  int num_threads = omp_get_max_threads();
  Atrip::omp_threads = num_threads * omp_granularity_;
  // not sure if this is still working correctly
  std::tie(Atrip::local_rank,
           Atrip::node_id,
           Atrip::n_nodes,
           Atrip::ranks_per_node) = get_mpi_info(atrip_world);

  //WATCH OUT: this is the node_id in the default rank distribution!
  Atrip::node_ids.resize(Atrip::np);
  MPI_Allgather(&Atrip::node_id, 1, MPI_INT, Atrip::node_ids.data(), 1, MPI_INT, Atrip::communicator);
}

template <typename F>
Atrip::Output Atrip::run(Atrip::Input<F> const &in) {


  const size_t No = in.epsilon_i->size();
  const size_t Nv = in.epsilon_a->size();
  LOG(0, "Atrip") << "No: " << No << "\n";
  LOG(0, "Atrip") << "Nv: " << Nv << "\n";
  LOG(0, "Atrip") << "np: " << Atrip::np << "\n";

#if defined(HAVE_ACC)
  int ngcards;
  ACC_CHECK_SUCCESS("initializing accelerator", ACC_INIT(0));
  ACC_CHECK_SUCCESS("getting device count", ACC_DEVICE_GET_COUNT(&ngcards));
  LOG(0, "Atrip") << "ngcards: " << ngcards << "\n";
  if (Atrip::ranks_per_node > ngcards) {
    const auto msg = _FORMAT(
        "ATRIP: You are running on more ranks per node than the number of "
        "graphic cards\n"
        "You have %d cards at your disposal\n",
        ngcards);
    std::cerr << msg;
    throw msg;
  } else if (Atrip::ranks_per_node < ngcards) {
    const auto msg = _FORMAT(
        "You have %d cards at your disposal.\n"
        "You will be only using %d, i.e, the number of ranks\n",
        ngcards,
        Atrip::ranks_per_node);
    std::cerr << msg;
  }

  for (size_t _rank = 0; _rank < Atrip::np; _rank++) {
    if (Atrip::rank == _rank) {
      ACC_CONTEXT ctx;
      ACC_DEVICE dev;
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
      ACC_CHECK_SUCCESS("getting device for index <rank>",
                        ACC_DEVICE_GET(&dev, Atrip::rank % ngcards));
      ACC_CHECK_SUCCESS("creating a accelerator context",
                        ACC_CONTEXT_CREATE(&ctx, 0, dev));
      ACC_CHECK_SUCCESS("setting the context", ACC_CONTEXT_SET_CURRENT(ctx));

      // get information of the device
      ACC_CHECK_SUCCESS(
          "getting memory information",
          ACC_MEM_GET_INFO(&memory.avail.free, &memory.avail.total));
      ACC_CHECK_SUCCESS("getting name", ACC_GET_DEVICE_NAME(name, 256, dev));
      ACC_CHECK_SUCCESS("getting total memory",
                        ACC_DEVICE_TOTAL_MEM(&memory.total, dev));

      printf(
          "\n"
          "CARD RANK %d\n"
          "============\n"
          "\tnumber: %1$ld\n"
          "\tname: %s\n"
          "\tAvail. Free/Total mem (GB): %f/%f\n"
          "\tFree memory (GB): %f\n"
          "\n",
          Atrip::rank,
          name,
          memory.avail.free / 1024.0 / 1024.0 / 1024.0,
          memory.avail.total / 1024.0 / 1024.0 / 1024.0,
          memory.total / 1024.0 / 1024.0 / 1024.0);
      std::free((void *)name);

      ACC_CHECK_BLAS("creating a cublas handle",
                     ACC_BLAS_CREATE(&Atrip::cuda.handle));
    }
    MPI_Barrier(Atrip::communicator);
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
  std::vector<F> &_epsi = *in.epsilon_i, &_epsa = *in.epsilon_a,
                 &_Tai = *in.Tph;

  if (in.ijkabc) {
    std::transform(_Tai.begin(), _Tai.end(), _Tai.begin(), [](F &c) {
      return -c;
    });
  }
  // TODO: free memory pointers in the end of the algorithm
  DataPtr<F> Tijk, Zijk;

#if defined(HAVE_ACC)
  DataPtr<F> Tai, epsi, epsa;

  // TODO: free memory pointers in the end of the algorithm

  MALLOC_DATA_PTR("Tai", &Tai, sizeof(F) * _Tai.size());
  MALLOC_DATA_PTR("epsi", &epsi, sizeof(F) * _epsi.size());
  MALLOC_DATA_PTR("epsa", &epsa, sizeof(F) * _epsa.size());

  ACC_CHECK_SUCCESS("memcpy Tai",
                    ACC_MEMCPY_HOST_TO_DEV(Tai,
                                           (void *)_Tai.data(),
                                           sizeof(F) * _Tai.size()));
  ACC_CHECK_SUCCESS("memcpy epsi",
                    ACC_MEMCPY_HOST_TO_DEV(epsi,
                                           (void *)_epsi.data(),
                                           sizeof(F) * _epsi.size()));
  ACC_CHECK_SUCCESS("memcpy epsa",
                    ACC_MEMCPY_HOST_TO_DEV(epsa,
                                           (void *)_epsa.data(),
                                           sizeof(F) * _epsa.size()));

#else
  DataPtr<F> Tai = _Tai.data(), epsi = _epsi.data();
#endif

  MALLOC_DATA_PTR("Zijk", &Zijk, sizeof(DataFieldType<F>) * No * No * No);
  MALLOC_DATA_PTR("Tijk", &Tijk, sizeof(DataFieldType<F>) * No * No * No);

  if (Atrip::rank_round_robin) {
    LOG(0, "Atrip") << "Doing rank round robin slices distribution\n";
  } else {
    LOG(0, "Atrip")
        << "Doing node > local rank round robin slices distribution\n";
  }

  // a, b, c, d, e, f and P => Nv
  // H                      => No
  // total_source_sizes contains a list of the number of elements
  // in all sources of every tensor union, therefore n_slices * slice_size

  // TODO: remove the alignment of sources in one big block or
  //        add it as an option

  // const std::vector<size_t> total_source_sizes = {
  //     // ABPH
  //     SliceUnion<F>::get_size({Nv, No}, {Nv, Nv}, (size_t)np, Atrip::communicator),
  //     // ABHH
  //     SliceUnion<F>::get_size({No, No}, {Nv, Nv}, (size_t)np, Atrip::communicator),
  //     // TABHH
  //     SliceUnion<F>::get_size({No, No}, {Nv, Nv}, (size_t)np, Atrip::communicator),
  //     // TAPHH
  //     SliceUnion<F>::get_size({Nv, No, No}, {Nv}, (size_t)np, Atrip::communicator),
  //     // HHHA
  //     SliceUnion<F>::get_size({No, No, No}, {Nv}, (size_t)np, Atrip::communicator),
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

    SliceUnion<F> abph({Slice<F>::AB, Slice<F>::BC, Slice<F>::AC, Slice<F>::BA, Slice<F>::CB, Slice<F>::CA},
                       Slice<F>::Name::VABCI,
                       *in.sVabph,
                       {Nv, Nv},
                       12
                      );

    SliceUnion<F> abhh({Slice<F>::AB, Slice<F>::BC, Slice<F>::AC},
                        Slice<F>::Name::VABIJ,
                        *in.sVabhh,
                        {Nv, Nv},
                        6);

    SliceUnion<F> tabhh({Slice<F>::AB, Slice<F>::BC, Slice<F>::AC},
                        Slice<F>::Name::TABIJ,
                        *in.sTabhh,
                        {Nv, Nv},
                        6);

  )

  // BUILD SLICES PARAMETRIZED BY NV ==================================={{{1
  WITH_CHRONO(
      "nv-slices", LOG(0, "Atrip") << "building NV slices\n";
      // TODO
      // DataPtr<F> offseted_pointer = all_sources_pointer
      //                             * total_source_sizes[_source_pointer_idx++];
    SliceUnion<F> taphh({Slice<F>::Type::A, Slice<F>::Type::B, Slice<F>::Type::C},
                        Slice<F>::Name::TA,
                        *in.sTaphh,
                        {Nv},
                        6);

    SliceUnion<F> hhha({Slice<F>::Type::A, Slice<F>::Type::B, Slice<F>::Type::C},
                       Slice<F>::Name::VIJKA,
                       *in.sVhhha,
                       {Nv},
                       6);

  )

  // all tensors
  std::vector<SliceUnion<F> *> unions = {&taphh, &hhha, &abph, &abhh, &tabhh};

  // IF (cT) IS USED: HANDLE TWO FURTHER SLICES==========================={{{1
  ///HHHA<F> *jhhha = nullptr;
  SliceUnion<F> *jhhha = nullptr;
  SliceUnion<F> *jabph = nullptr;
  //ABPH<F> *jabph = nullptr;
  if (in.cT) {
    WITH_CHRONO("Jhhha-slice",
                /**/ LOG(0, "Atrip") << "slicing Jijka" << std::endl;

      jhhha = new SliceUnion<F>({Slice<F>::Type::A, Slice<F>::Type::B, Slice<F>::Type::C},
                                Slice<F>::Name::JIJKA,
                                *in.sJhhha,
                                {Nv},
                                6);

    )
    unions.push_back(jhhha);
  }

  if (in.cT) {
    WITH_CHRONO("Jabph-slice",
                /**/ LOG(0, "Atrip") << "slicing Jabci" << std::endl;

      jabph = new SliceUnion<F>({Slice<F>::AB, Slice<F>::BC, Slice<F>::AC, Slice<F>::BA, Slice<F>::CB, Slice<F>::CA},
                                Slice<F>::Name::JABCI,
                                *in.sJabph,
                                {Nv, Nv},
                                12);
    )


    unions.push_back(jabph);
  }

  DataFieldType<F> *_t_buffer, *_vhhh;
  MALLOC_DATA_PTR("Tijk buffer for doubles equations",
                  &_t_buffer,
                  sizeof(DataFieldType<F>) * No * No * No);
  MALLOC_DATA_PTR("Vijk buffer for doubles equations",
                  &_vhhh,
                  sizeof(DataFieldType<F>) * No * No * No);

  // #if defined(HAVE_ACC)
  //   WITH_CHRONO("double:acc:alloc",
  //               MALLOC_DATA_PTR("t_buffer",
  //                               (ACC_DEVICE_PTR *)&_t_buffer,
  //                               No * No * No * sizeof(DataFieldType<F>));
  //               MALLOC_DATA_PTR("vhhh",
  //                               (ACC_DEVICE_PTR *)&_vhhh,
  //                               No * No * No * sizeof(DataFieldType<F>));)
  // #endif

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
              auto const tuples_list = distribution->get_tuples(Nv, Atrip::communicator);)
  const size_t n_iterations = tuples_list.size();
  {
    LOG(0, "Atrip") << "#iterations: " << n_iterations << "/"
                    << n_iterations * Atrip::np << "\n";
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
      [&unions, &in, Nv](ABCTuple const &abc,
                                    MPI_Comm const &c,
                                    size_t iteration) -> Database {
    if (in.tuples_distribution == Atrip::Input<F>::TuplesDistribution::NAIVE) {

      WITH_CHRONO("db:comm:naive",
                  auto const &db =
                      naive_database<F>(unions, Nv, Atrip::np, iteration);)
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

      Database db(Atrip::np * ldb.size(), ldb[0]);

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

      MPI_Barrier(c); // we need a barrier - otherwise mpi oject gets free'd while in use!
      WITH_CHRONO("db:comm:type:free", MPI_Type_free(&MPI_LDB_ELEMENT);)

//      WITH_CHRONO("db:comm:type:free", Slice<F>::mpi::free_local_database_element();)

      return db;
    }
  };

  auto do_io_phase =
      [&unions, No, Nv, &db_last_iteration_time](
          Database const &db,
          ABCTuple const abc,
          size_t iteration) {
        IGNORABLE(iteration); // iteration used to print database
        const size_t localDBLength = db.size() / Atrip::np;

        size_t send_tag = 0, recv_tag = Atrip::rank * localDBLength;

        // RECIEVE PHASE ======================================================
        {
          // At this point, we have already send to everyone that fits
          auto const &begin = &db[Atrip::rank * localDBLength],
                     end = begin + localDBLength;
          for (auto it = begin; it != end; ++it) {
            recv_tag++;
            auto const &el = *it;
            auto &u = union_by_name(unions, el.name);

            WITH_DBG std::cout << Atrip::rank << ":r"
                               << "♯" << recv_tag << " =>"
                               << " «n" << el.name << ", t" << el.info.type
                               << ", s" << el.info.state << "»"
                               << " ⊙ {" << Atrip::rank << "⇐" << el.info.from.rank
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
                  Atrip::rank,
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

        MPI_Barrier(Atrip::communicator);

        // SEND PHASE =========================================================
        for (size_t other_rank = 0; other_rank < Atrip::np; other_rank++) {
          auto const &begin = &db[other_rank * localDBLength],
                     end = begin + localDBLength;
          for (auto it = begin; it != end; ++it) {
            send_tag++;
            typename Slice<F>::LocalDatabaseElement const &el = *it;

            if (el.info.from.rank != Atrip::rank) continue;

            auto &u = union_by_name(unions, el.name);
            WITH_DBG std::cout
                << Atrip::rank << ":s"
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

  MPI_Barrier(Atrip::communicator);
  LocalOutput<EnergyType<F>> local_output = {0, 0};
  Output global_output = {0, 0};
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
        local_output.energy = -EnergyType<F>(c.energy);
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
#if defined(HAVE_ACC)
            EnergyType<F> *tuple_energy;
            MALLOC_DATA_PTR("tuple energy",
                            (DataPtr<EnergyType<F>> *)&tuple_energy,
                            sizeof(EnergyType<F>));
#else
          EnergyType<F> _tuple_energy(0.);
          EnergyType<F> *tuple_energy = &_tuple_energy;
#endif /* defined(HAVE_ACC) */

            int distinct(0);
            if (abc[0] == abc[1]) distinct++;
            if (abc[1] == abc[2]) distinct--;
            const PrecisionType<F> epsabc =
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

#if defined(HAVE_ACC)
            EnergyType<F> host_tuple_energy;
            ACC_MEMCPY_DEV_TO_HOST((void *)&host_tuple_energy,
                                   (DataPtr<EnergyType<F>>)tuple_energy,
                                   sizeof(EnergyType<F>));
#else
          EnergyType<F> host_tuple_energy = *tuple_energy;
#endif /* defined(HAVE_ACC) */

            return host_tuple_energy;
          }
        return EnergyType<F>(0.0);
      };

  // Swich number of omp threads
  omp_set_num_threads(Atrip::omp_threads);
  LOG(0, "atrip") << "number of omp threads " << Atrip::omp_threads << std::endl;
  for (size_t i = first_iteration, iteration = first_iteration + 1;
       i < tuples_list.size();
       i++, iteration++) {
    Atrip::chrono["iterations"].start();
    Atrip::chrono["db-last-iteration"].start();
#if defined(HAVE_ACC)
    char nvtx_name[60];
    sprintf(nvtx_name, "iteration: %ld", i);
#  if defined(HAVE_CUDA)
    nvtxRangePushA(nvtx_name);
#  endif /* defined(HAVE_CUDA) */
#  if defined(HAVE_HIP)
    omnitrace_user_push_region(nvtx_name);
#  endif /* defined(HAVE_CUDA) */
#endif   /* defined(HAVE_ACC) */

    // check overhead from chrono over all iterations
    WITH_CHRONO("start:stop", {})

    // check overhead of doing a barrier at the beginning
    WITH_CHRONO(
        "oneshot-mpi:barrier",
        WITH_CHRONO("mpi:barrier", if (in.barrier) MPI_Barrier(Atrip::communicator);))

    // write checkpoints
    // TODO: ENABLE THIS
    if (iteration % checkpoint_mod == 0 && false) {
      EnergyType<F> global_energy = 0;
      MPI_Reduce(&local_output.energy,
                 &global_energy,
                 1,
                 traits::mpi::datatype_of<EnergyType<F>>(),
                 MPI_SUM,
                 0,
                 Atrip::communicator);
      Checkpoint out = {No,
                        Nv,
                        Atrip::ranks_per_node,
                        Atrip::n_nodes,
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
                 Atrip::communicator);

      size_t local_send(0);
      MPI_Reduce(&Atrip::local_send,
                 &local_send,
                 1,
                 MPI_UINT64_T,
                 MPI_SUM,
                 0,
                 Atrip::communicator);

      double bytes_sent(0.0);
      MPI_Reduce(&Atrip::bytes_sent,
                 &bytes_sent,
                 1,
                 MPI_DOUBLE,
                 MPI_SUM,
                 0,
                 Atrip::communicator);

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
          << "GF) :: GB sent per rank: " << bytes_sent / 1073741824.0 << " :: "
          << (total_send > 0UL ? (double)network_send / total_send : 0UL)
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

    WITH_CHRONO(
        "with_rank",
        WITH_RANK << " :it " << iteration << " :abc " << abc << " :abcN "
                  << (abc_next ? array_to_string(*abc_next) : "None") << "\n";)

    // COMM FIRST DATABASE
    // ================================================{{{1
    if (i == first_iteration) {
      WITH_RANK << "__first__:first database ............ \n";
      const auto db = communicate_database(abc, Atrip::communicator, i);
      WITH_RANK << "__first__:first database communicated \n";
      WITH_RANK << "__first__:first database io phase \n";
      do_io_phase(db, abc, i);
      WITH_RANK << "__first__:first database io phase DONE\n";
      WITH_RANK << "__first__::::Unwrapping all slices for first database\n";
      for (auto &u : unions) u->unwrap_all(abc);
      WITH_RANK << "__first__::::Unwrapping slices for first database DONE\n";
      MPI_Barrier(Atrip::communicator);
    }

    // COMM NEXT DATABASE
    // ================================================={{{1
    if (abc_next) {
      WITH_RANK << "__comm__:" << iteration << "th communicating database\n";
      WITH_CHRONO("db:comm",
                  const auto db = communicate_database(*abc_next, Atrip::communicator, i);)
      WITH_CHRONO("db:io", do_io_phase(db, abc, i + 1);)
      WITH_RANK << "__comm__:" << iteration << "th database io phase DONE\n";
    }

    // COMPUTE DOUBLES
    // ===================================================={{{1
    OCD_Barrier(Atrip::communicator);
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
        WITH_CHRONO("blocking-barrier", MPI_Barrier(Atrip::communicator);)
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
    OCD_Barrier(Atrip::communicator);
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
    if (!is_fake_tuple(i) && in.cT) {

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
    OCD_Barrier(Atrip::communicator);
    if (abc_next) {
      WITH_RANK << "__gc__:" << iteration << "-th cleaning up.......\n";
      for (auto &u : unions) {

        u->unwrap_all(abc);
        WITH_RANK << "__gc__:n" << u->name << " :it " << iteration << " :abc "
                  << abc << " :abcN " << *abc_next << "\n";
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
                    _FORMAT("This slice should not be undispatched! %s",
                            slice.info));
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
    // WITH_CHRONO("mpi:barrier", MPI_Barrier(Atrip::communicator););
    ACC_DEVICE_SYNCHRONIZE();
#if defined(HAVE_ACC)
#  if defined(HAVE_CUDA)
    nvtxRangePop();
#  endif /* defined(HAVE_CUDA) */
#  if defined(HAVE_HIP)
    omnitrace_user_pop_region(nvtx_name);
#  endif /* defined(HAVE_HIP) */
#endif   /* defined(HAVE_ACC) */

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

#if defined(HAVE_ACC)
  ACC_FREE(Tai);
  ACC_FREE(epsi);
  ACC_FREE(epsa);
#endif

  FREE_DATA_PTR("Zijk", Zijk);
  FREE_DATA_PTR("Tijk", Tijk);
  FREE_DATA_PTR("Doubles Tijk buffer", _t_buffer);
  FREE_DATA_PTR("Doubles Vijk buffer", _vhhh);

  if (jhhha) delete jhhha;
  if (jabph) delete jabph;
  MPI_Barrier(Atrip::communicator);

// PRINT TUPLES
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%{{{1
#if defined(HAVE_OCD) || defined(ATRIP_PRINT_TUPLES)
  LOG(0, "Atrip") << "tuple energies"
                  << "\n";
  for (size_t i = 0; i < Atrip::np; i++) {
    MPI_Barrier(Atrip::communicator);
    for (auto const &pair : tuple_energies) {
      if (i == Atrip::rank)
        std::cout << pair.first[0] << " " << pair.first[1] << " "
                  << pair.first[2] << std::setprecision(15) << std::setw(23)
                  << " tuple_energy: " << pair.second << "\n";
    }
  }
#endif

  // COMMUNICATE THE ENERGIES
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%{{{1
  LOG(0, "Atrip") << "COMMUNICATING ENERGIES \n";
  EnergyType<F> _pt_energy = 0, _ct_energy = 0;
  MPI_Reduce(&local_output.energy,
             &_pt_energy,
             1,
             MPI_DOUBLE,
             MPI_SUM,
             0,
             Atrip::communicator);
  global_output.energy = _pt_energy;

  MPI_Reduce(&local_output.ct_energy,
             &_ct_energy,
             1,
             MPI_DOUBLE,
             MPI_SUM,
             0,
             Atrip::communicator);
  global_output.ct_energy = _ct_energy;
  if (!in.ijkabc) {
    global_output.energy = -global_output.energy;
    global_output.ct_energy = -global_output.ct_energy;
  }
  WITH_RANK << "local energy " << local_output.energy << "\n";
  LOG(0, "Atrip") << "Energy (T): " << std::setprecision(15) << std::setw(23)
                  << global_output.energy << std::endl;
  LOG(0, "Atrip") << "Energy (cT): " << std::setprecision(15) << std::setw(23)
                  << global_output.ct_energy << std::endl;

  // PRINT TIMINGS {{{1
  if (in.chrono)
    for (auto const &pair : Atrip::chrono)
      LOG(0, "atrip:chrono")
          << pair.first << " " << pair.second.count() << std::endl;
  size_t _it = (in.max_iterations) ? in.max_iterations : n_iterations;
  LOG(0, "atrip:flops(doubles)")
      << _it * doubles_flops / Atrip::chrono["doubles"].count()
      << "\n";
  LOG(0, "atrip:flops(iterations)")
      << _it * doubles_flops / Atrip::chrono["iterations"].count()
      << "\n";

  // TODO: change the sign in  the getEnergy routines
  return global_output;
}

// instantiate
template Atrip::Output Atrip::run(Atrip::Input<float> const &in);
template Atrip::Output Atrip::run(Atrip::Input<double> const &in);
template Atrip::Output Atrip::run(Atrip::Input<Complex> const &in);
// Main:1 ends here
