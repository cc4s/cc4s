/* Copyright 2021 cc4s.org
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <Cc4s.hpp>
#include <Keywords.hpp>
#include <Log.hpp>
#include <MathFunctions.hpp>
#include <TensorSet.hpp>
#include <algorithms/PerturbativeTriples.hpp>
#include <algorithms/PerturbativeTriplesReference.hpp>
#include <algorithms/PerturbativeTriplesStar.hpp>
#include <atrip.hpp>
#include <atrip/Debug.hpp>
#include <atrip/Reader.hpp>
#include <atrip/Sources.hpp>
#include <fstream>
#include <iomanip>
#include <tcc/Tcc.hpp>


using namespace cc4s;
ALGORITHM_REGISTRAR_DEFINITION(PerturbativeTriples)

template <typename F, typename TE>
Ptr<MapNode> runAtrip(const Ptr<MapNode> &arguments) {
  auto amplitudes = arguments->getPtr<TensorSet<F, TE>>("amplitudes");
  if (!amplitudes)
    return nullptr;

  auto result(New<MapNode>(SOURCE_LOCATION));

  const bool single_precision =
      arguments->getValue<int>("singlePrecision", 0) == 1;

  const size_t omp_granularity = arguments->getValue<int>("ompGranularity", 1);

  const bool delete_tnsr = arguments->getValue<int>("deleteTensors", 1) == 1;

  if (delete_tnsr) {
    OUT() << "Warning: All tensors will be deleted. Ensure this is"
          << " the final algorithm accessing the tensor data." << std::endl;
  }
  /*
    This should be added to the manual for CCSD(T).

    The checkpoint can be written every $x$ iterations
    or every $x$ percentage.

    The default in cc4s is to write a checkpoint every 1.0%.
    The checkpoint variables are controlled by the input parameters

    | Name                         | Default                 |
    |------------------------------+-------------------------|
    | =checkpointAtEveryIteration= | 0                       |
    | =checkpointAtPercentage=     | 1.0                     |
    | =checkpointPath=             | "atrip-checkpoint.yaml" |
    | =useCheckpoint=              | 1                       |

    There is however right now a pitfall, you should
    be careful to perform two calculation in a row, since then
    CCSD(T) will use the automatically generated checkpoint file
    from the previous calculation.

   */

  const auto print_at_percentage =
                 arguments->getValue<int>("printAtPercentage", 10),
             checkpoint_at_every_iteration =
                 arguments->getValue<int>("checkpointAtEveryIteration", 500);
  const auto checkpoint_at_percentage =
      arguments->getValue<double>("checkpointAtPercentage", 1.0);
  const auto checkpoint_path = arguments->getValue<std::string>(
      "checkpointPath", "PerturbativeTriples.checkpoint");

  const bool rankRoundRobin =
      arguments->getValue<std::string>("tuplesRoundRobin", "node") == "node"
          ? false
          : true;
  bool use_checkpoint = arguments->getValue<int>("useCheckpoint", 1) == 1;

  size_t first_iteration = arguments->getValue<int>("firstIteration", 0);
  if (first_iteration)
    OUT() << "Starting atrip from iteration: " << first_iteration << std::endl;

  atrip::Atrip::init(MPI_COMM_WORLD);
  atrip::Atrip::Input<F> in;

#define __V__(_idx)                                                            \
  ([&arguments]() {                                                            \
    COMPILE(                                                                   \
        arguments->getPtr<TensorSet<F, TE>>("coulombIntegrals")->get(#_idx))   \
        ->execute();                                                           \
    return &(arguments->getPtr<TensorSet<F, TE>>("coulombIntegrals")           \
                 ->get(#_idx)                                                  \
                 ->evaluate()                                                  \
                 ->getMachineTensor()                                          \
                 ->tensor);                                                    \
  })()
#define __tccV__(_idx)                                                         \
  ([&arguments]() {                                                            \
    COMPILE(                                                                   \
        arguments->getPtr<TensorSet<F, TE>>("coulombIntegrals")->get(#_idx))   \
        ->execute();                                                           \
    return (arguments->getPtr<TensorSet<F, TE>>("coulombIntegrals")            \
                ->get(#_idx)                                                   \
                ->evaluate());                                                 \
  })()
#define __delete_V__(_idx)                                                     \
  ([&arguments]() {                                                            \
    COMPILE(                                                                   \
        arguments->getPtr<TensorSet<F, TE>>("coulombIntegrals")->get(#_idx))   \
        ->execute();                                                           \
    return arguments->getPtr<TensorSet<F, TE>>("coulombIntegrals")             \
        ->get(#_idx)                                                           \
        ->evaluate()                                                           \
        ->deleteMachineTensor();                                               \
  })()
#define __tccT__(_idx)                                                         \
  (arguments->getPtr<TensorSet<F, TE>>("amplitudes")->get(_idx)->evaluate())
#define __T__(_idx)                                                            \
  &(arguments->getPtr<TensorSet<F, TE>>("amplitudes")                          \
        ->get(_idx)                                                            \
        ->evaluate()                                                           \
        ->getMachineTensor()                                                   \
        ->tensor)
#define __delete_T__(_idx)                                                     \
  do {                                                                         \
    arguments->getPtr<TensorSet<F, TE>>("amplitudes")                          \
        ->get(_idx)                                                            \
        ->evaluate()                                                           \
        ->deleteMachineTensor();                                               \
  } while (0)
#define __eps__(_idx)                                                          \
  &(arguments->getPtr<TensorSet<Real<>, TE>>("slicedEigenEnergies")            \
        ->get(#_idx)                                                           \
        ->evaluate()                                                           \
        ->getMachineTensor()                                                   \
        ->tensor)
#define __CT__(_idx)                                                           \
  &(arguments->getPtr<TensorSet<F, TE>>(kwd::CT_INTERMEDIATE)                  \
        ->get(_idx)                                                            \
        ->evaluate()                                                           \
        ->getMachineTensor()                                                   \
        ->tensor)
#define __delete_CT__(_idx)                                                    \
  do {                                                                         \
    arguments->getPtr<TensorSet<F, TE>>(kwd::CT_INTERMEDIATE)                  \
        ->get(#_idx)                                                           \
        ->evaluate()                                                           \
        ->deleteMachineTensor();                                               \
  } while (0)

  CTF::Tensor<F> epsi(1, (__eps__(h))->lens, (__eps__(h))->sym),
      epsa(1, (__eps__(p))->lens, (__eps__(p))->sym);
  const auto toComplex =
      CTF::Transform<double, F>([](double d, F &f) { f = d; });
  toComplex((*__eps__(h))["i"], epsi["i"]);
  toComplex((*__eps__(p))["a"], epsa["a"]);

  const auto _no = *(epsi.lens), _nv = *(epsa.lens);
  std::vector<F> epsi_v(_no), epsa_v(_nv), Tai_v(_no * _nv);
  epsi.read_all(epsi_v.data());
  epsa.read_all(epsa_v.data());
  (__T__("ph"))->read_all(Tai_v.data());

  {
    // Register printing function
    double _No = *(__eps__(h))->lens, _Nv = *(__eps__(p))->lens,
           doublesFlops = _No * _No * _No * (_No + _Nv) * 2.0 * 6.0 / 1.0e9;
    if (arguments->isGiven(kwd::CT_INTERMEDIATE))
      doublesFlops *= 2.0;
    // this is not working if we have single prec complex!!
    if (typeid(F) == typeid(Complex<>))
      doublesFlops *= 4.0;

    double lastElapsedTime = 0;
    auto const rank = Cc4s::world->getRank();
    bool firstHeaderPrinted = false;
    atrip::register_iteration_descriptor(
        [doublesFlops, &firstHeaderPrinted, rank,
         &lastElapsedTime](atrip::IterationDescription const &d) {
          const char *fmt_header = "%-13s%-10s%-13s",
                     *fmt_nums = "%-13.0f%-10.0f%-13.3f";
          char out[256];
          if (!firstHeaderPrinted) {
            sprintf(out, fmt_header, "Progress(%)", "time(s)", "GFLOP/s");
            firstHeaderPrinted = true;
            OUT() << out << std::endl;
          }
          sprintf(out, fmt_nums,
                  double(d.current_iteration) / double(d.total_iterations) *
                      100,
                  (d.current_elapsed_time - lastElapsedTime),
                  d.current_iteration * doublesFlops / d.current_elapsed_time);
          lastElapsedTime = d.current_elapsed_time;
          OUT() << out << std::endl;
        });
  }
  // atrip's logger which writes into the cc4s logstream
  atrip::getLogger() = [](const std::string &s) {
    Log::getLogStream() << Log::getLogHeaderFunction()(SourceLocation())
                        << "[ATRIP]:" << s << std::endl;
  };

  // currently we will log 500 times. should be in the minutes to couple of
  // minutes range!
  size_t log_iterations =
      _nv * (_nv + 1) * (_nv + 2) / 6 / Cc4s::getProcessesCount() / 500;
  LOG_LOCATION(SOURCE_LOCATION)
      << "Atrip will write a log every " << log_iterations << " iterations\n";

  atrip::Atrip::Output out;

  MPI_Comm atrip_comm = atrip::create_comm(Cc4s::world->getComm(),
                                           omp_granularity, rankRoundRobin);
  // atrip::print_comm_mapping(comm, atrip_comm);
  if (atrip_comm != MPI_COMM_NULL) {
    atrip::Atrip::init(atrip_comm, omp_granularity);
  }

  //////////////////////////////////////////////////////////////////////////////
  // CHECKPOINT CHECKER
  //////////////////////////////////////////////////////////////////////////////
  auto checkpoint_file_exists = [&checkpoint_path]() {
    ifstream checkpoint(checkpoint_path.c_str());
    return checkpoint.good();
  };

  // print warning if checkpoint file exists
  if (checkpoint_file_exists() && use_checkpoint) {
    WARNING_LOCATION(arguments->sourceLocation)
        << ("Checkpoint file found. Calculation will be\n"
            "restarted from the given checkpoint. Make sure that the number\n"
            "of cores & nodes are the same as in the previous run.")
        << endl;

    {
      auto Tpphh = __tccT__("pphh");
      auto Tph = __tccT__("ph");
      Real<> degeneracy(2);
      auto direct(Tcc<TE>::template tensor<F>("D"));
      auto exchange(Tcc<TE>::template tensor<F>("X"));
      auto directSingles(Tcc<TE>::template tensor<F>("D"));
      auto exchangeSingles(Tcc<TE>::template tensor<F>("X"));

      COMPILE((*direct)[""] <<= 0.5 * degeneracy * degeneracy *
                                (*Tpphh)["abij"] * (*(__tccV__(hhpp)))["ijab"],
              (*directSingles)[""] <<= 0.5 * degeneracy * degeneracy *
                                       (*Tph)["ai"] * (*Tph)["bj"] *
                                       (*(__tccV__(hhpp)))["ijab"],
              (*exchange)[""] <<= -0.5 * degeneracy * (*Tpphh)["abij"] *
                                  (*(__tccV__(hhpp)))["ijba"],
              (*exchangeSingles)[""] <<= -0.5 * degeneracy * (*Tph)["ai"] *
                                         (*Tph)["bj"] *
                                         (*(__tccV__(hhpp)))["ijba"])
          ->execute();
      F D(direct->read()), X(exchange->read()), Ds(directSingles->read()),
          Xs(exchangeSingles->read()), e = D + X + Ds + Xs;

      OUT() << "CCSD correlation energy: " << std::real(e) << std::endl;
    }
  }
  //////////////////////////////////////////////////////////////////////////////

  if (!single_precision) {

    const auto distribution =
        arguments->getValue<std::string>("tuplesDistribution", "group") ==
                "group"
            ? atrip::Atrip::Input<F>::TuplesDistribution::GROUP_AND_SORT
            : atrip::Atrip::Input<F>::TuplesDistribution::NAIVE;

    auto &Tpphh = get_ctf_tensor<F, TE>(amplitudes, "pphh");
    size_t nv(_nv), no(_no);
    auto sVabph = atrip::ctfReader<F>(*(__V__(ppph)), {nv, nv, nv, no},
                                      {1, 1, 0, 0}, delete_tnsr);
    auto sVabhh =
        atrip::ctfReader<F>(*(__V__(pphh)), {nv, nv, no, no}, {1, 1, 0, 0});
    auto sTabhh = atrip::ctfReader<F>(Tpphh, {nv, nv, no, no}, {1, 1, 0, 0});
    auto sTaphh =
        atrip::ctfReader<F>(Tpphh, {nv, nv, no, no}, {1, 0, 0, 0}, delete_tnsr);
    auto sVhhha =
        atrip::ctfReader<F>(*(__V__(hhhp)), {no, no, no, nv}, {0, 0, 0, 1});

    atrip::Atrip::Input<F> in;
    in
        // setup tensors
        .with_epsilon_i(&epsi_v)
        .with_epsilon_a(&epsa_v)
        .with_Tph(&Tai_v)
        .with_sVabph(&sVabph)
        .with_sVabhh(&sVabhh)
        .with_sTabhh(&sTabhh)
        .with_sTaphh(&sTaphh)
        .with_sVhhha(&sVhhha)
        // some options
        .with_barrier(false)
        .with_percentage_mod(print_at_percentage)
        .with_tuples_distribution(distribution)
        .with_rank_round_robin(rankRoundRobin)
        // checkpoint options
        .with_checkpoint_at_every_iteration(checkpoint_at_every_iteration)
        .with_checkpoint_at_percentage(checkpoint_at_percentage)
        .with_checkpoint_path(checkpoint_path)
        .with_read_checkpoint_if_exists(use_checkpoint)
        .with_cT(arguments->isGiven(kwd::CT_INTERMEDIATE))
        .with_first_iteration(first_iteration)
        .with_iteration_log(log_iterations);

    atrip::Sources<F> *sJabph = nullptr, *sJhhha = nullptr;
    if (arguments->isGiven(kwd::CT_INTERMEDIATE)) {
      OUT() << "Using ct" << std::endl;
      std::vector<size_t> vvvo{nv, nv, nv, no}, ooov{no, no, no, nv};
      sJabph = new atrip::Sources<F>(atrip::ctfReader<F>(
          *(__CT__("ppph")), vvvo, {1, 1, 0, 0}, delete_tnsr));
      sJhhha = new atrip::Sources<F>(atrip::ctfReader<F>(
          *(__CT__("hhhp")), ooov, {0, 0, 0, 1}, delete_tnsr));
      in.with_sJabph(sJabph);
      in.with_sJhhha(sJhhha);
    }
    OUT() << "All tensors loaded; starting main loop." << std::endl;
    out = atrip::Atrip::run<F>(in);
  } else {

    const auto distribution =
        arguments->getValue<std::string>("tuplesDistribution", "group") ==
                "group"
            ? atrip::Atrip::Input<float>::TuplesDistribution::GROUP_AND_SORT
            : atrip::Atrip::Input<float>::TuplesDistribution::NAIVE;
    atrip::Atrip::Input<float> in;

    std::vector<float> epsi_f(epsi_v.size()), epsa_f(epsa_v.size()),
        Tai_f(Tai_v.size());

    std::transform(epsi_v.begin(), epsi_v.end(), epsi_f.begin(),
                   [](const F &from) { return std::real(from); });
    std::transform(epsa_v.begin(), epsa_v.end(), epsa_f.begin(),
                   [](const F &from) { return std::real(from); });
    std::transform(Tai_v.begin(), Tai_v.end(), Tai_f.begin(),
                   [](const F &from) { return std::real(from); });

    const auto to_float = CTF::Transform<F, float>(
        [](F from, float &to) { to = std::real(from); });

    // source tensors
    CTF::Tensor<F> &F_T = *(__T__("pphh")), &F_Vpphh = *(__V__(pphh)),
                   &F_Vhhhp = *(__V__(hhhp)), &F_Vppph = *(__V__(ppph));

    // deal with Tabij
    CTF::Tensor<float> T(F_T.order, F_T.lens, F_T.sym, *F_T.wrld);
    to_float(F_T["abij"], T["abij"]);
    Cc4s::world->barrier();
    if (delete_tnsr)
      __delete_T__("pphh");

    // deal with Vpphh
    CTF::Tensor<float> Vpphh(F_Vpphh.order, F_Vpphh.lens, F_Vpphh.sym,
                             *F_Vpphh.wrld);
    to_float(F_Vpphh["abij"], Vpphh["abij"]);
    Cc4s::world->barrier();
    if (delete_tnsr)
      __delete_V__(pphh);

    // deal with Vhhhp
    CTF::Tensor<float> Vhhhp(F_Vhhhp.order, F_Vhhhp.lens, F_Vhhhp.sym,
                             *F_Vhhhp.wrld);
    to_float(F_Vhhhp["abij"], Vhhhp["abij"]);
    Cc4s::world->barrier();
    if (delete_tnsr)
      __delete_V__(hhhp);

    // deal with Vppph
    CTF::Tensor<float> Vppph(F_Vppph.order, F_Vppph.lens, F_Vppph.sym,
                             *F_Vppph.wrld);
    to_float(F_Vppph["abij"], Vppph["abij"]);
    Cc4s::world->barrier();
    if (delete_tnsr)
      __delete_V__(ppph);

    size_t nv(_nv), no(_no);
    auto sVabph = atrip::ctfReader<float>(Vppph, {nv, nv, nv, no}, {1, 1, 0, 0},
                                          delete_tnsr);
    auto sVabhh =
        atrip::ctfReader<float>(Vpphh, {nv, nv, no, no}, {1, 1, 0, 0});
    auto sTabhh = atrip::ctfReader<float>(T, {nv, nv, no, no}, {1, 1, 0, 0});
    auto sTaphh =
        atrip::ctfReader<float>(T, {nv, nv, no, no}, {1, 0, 0, 0}, delete_tnsr);
    auto sVhhha = atrip::ctfReader<float>(Vhhhp, {no, no, no, nv}, {0, 0, 0, 1},
                                          delete_tnsr);

    in
        // setup tensors
        .with_epsilon_i(&epsi_f)
        .with_epsilon_a(&epsa_f)
        .with_Tph(&Tai_f)
        .with_sVabph(&sVabph)
        .with_sVabhh(&sVabhh)
        .with_sTabhh(&sTabhh)
        .with_sTaphh(&sTaphh)
        .with_sVhhha(&sVhhha)
        // some options
        .with_barrier(false)
        .with_percentage_mod(print_at_percentage)
        .with_tuples_distribution(distribution)
        .with_rank_round_robin(rankRoundRobin)
        // checkpoint options
        .with_checkpoint_at_every_iteration(checkpoint_at_every_iteration)
        .with_checkpoint_at_percentage(checkpoint_at_percentage)
        .with_checkpoint_path(checkpoint_path)
        .with_read_checkpoint_if_exists(use_checkpoint)
        .with_cT(arguments->isGiven(kwd::CT_INTERMEDIATE))
        .with_first_iteration(first_iteration)
        .with_iteration_log(log_iterations);

    OUT() << "using single precision" << std::endl;

    std::shared_ptr<CTF::Tensor<float>> Jhhhp, Jppph;
    atrip::Sources<float> *sJabph = nullptr, *sJhhha = nullptr;

    if (arguments->isGiven(kwd::CT_INTERMEDIATE)) {
      OUT() << "Using ct" << std::endl;

      // Deal with Jhhhp
      CTF::Tensor<F> &F_Jhhhp = *(__CT__("hhhp")), F_Jppph = *(__CT__("ppph"));
      Jhhhp = std::make_shared<CTF::Tensor<float>>(F_Jhhhp.order, F_Jhhhp.lens,
                                                   F_Jhhhp.sym, *F_Jhhhp.wrld);
      to_float(F_Jhhhp["abij"], (*Jhhhp)["abij"]);
      Cc4s::world->barrier();
      if (delete_tnsr)
        __delete_CT__(hhhp);

      // Deal with Jppph
      Jppph = std::make_shared<CTF::Tensor<float>>(F_Jppph.order, F_Jppph.lens,
                                                   F_Jppph.sym, *F_Jppph.wrld);
      to_float(F_Jppph["abij"], (*Jppph)["abij"]);
      Cc4s::world->barrier();
      if (delete_tnsr)
        __delete_CT__(ppph);

      std::vector<size_t> vvvo{nv, nv, nv, no}, ooov{no, no, no, nv};
      sJabph = new atrip::Sources<float>(
          atrip::ctfReader<float>(*Jppph, vvvo, {1, 1, 0, 0}, delete_tnsr));
      sJhhha = new atrip::Sources<float>(
          atrip::ctfReader<float>(*Jhhhp, ooov, {0, 0, 0, 1}, delete_tnsr));
      in.with_sJabph(sJabph);
      in.with_sJhhha(sJhhha);
    }

    OUT() << "All tensors loaded; starting main loop." << std::endl;
    out = atrip::Atrip::run<float>(in);
  }

  if (first_iteration)
    OUT()
        << "WARNING: the provided energy is NOT the energy of ALL iterations\n";

  OUT() << "(T) correlation energy: " << std::setprecision(15) << std::setw(23)
        << out.energy << std::endl;
  if (arguments->isGiven(kwd::CT_INTERMEDIATE)) {
    OUT() << "(cT) correlation energy: " << std::setprecision(15)
          << std::setw(23) << out.ct_energy << std::endl;
  }
  // if atrip checkpoint file exists.
  // Delete it after a successfully finished calculation
  // if (checkpoint_file_exists()) std::remove(checkpoint_path.c_str());
  /* DUPLICATION WARNING:
   * --------------------
   *
   * If you change something in the output of this algorithm
   * be sure to change accordingly the output of the
   * PerturbativeTriplesReference, as they should be
   * symmetric.
   *
   */

  auto energy(New<MapNode>(SOURCE_LOCATION));
  energy->setValue("correlation", std::real(out.energy));
  energy->setValue("t", std::real(out.energy));

  if (arguments->isGiven(kwd::CT_INTERMEDIATE)) {
    energy->setValue("ct", std::real(out.ct_energy));
  }
  if (arguments->isGiven("mp2PairEnergies")) {
    const F triples_star = computeCssdPtStar<F, TE>(arguments, out.energy);
    OUT() << "(T*)-Bsie energy correction: " << std::setprecision(15)
          << std::setw(23) << std::real(triples_star) << std::endl;
    energy->setValue("starCorrection", std::real(triples_star));
  }

  using TSr = TensorSet<Real<>, TE>;
  energy->setValue("unit", arguments->getPtr<TSr>("slicedEigenEnergies")
                               ->get("h")
                               ->inspect()
                               ->getUnit());

  result->get("energy") = energy;

  return result;
}

template <typename F, typename TE>
Ptr<MapNode> atripDryRun(const Ptr<MapNode> &arguments) {
  auto result(New<MapNode>(SOURCE_LOCATION));
  auto eps(arguments->getPtr<TensorSet<Real<>, TE>>("slicedEigenEnergies"));
  auto epsh(eps->get("h"));
  auto epsp(eps->get("p"));
  const Natural<128> no = epsh->inspect()->getLen(0),
                     nv = epsp->inspect()->getLen(0),
                     np = Cc4s::getProcessesCount(), f = sizeof(F),
                     n_tuples = nv * (nv + 1) * (nv + 2) / 6 - nv;

  // Note: what is detected by the tcc memory stuff is: Tpphh and Jppph+Jhhhp
  //       in normal operation Vppph is still a recipe and will first be created
  DryMemory::allocate(nv * nv * no * no, f, "tpphh2", SOURCE_LOCATION);
  DryMemory::allocate(nv * nv * nv * no, f, "vppph", SOURCE_LOCATION);
  // this will be done before the main algorithm of atrip starts....
  DryMemory::allocate(nv * nv * nv * no, f, "vppph - hwm", SOURCE_LOCATION);
  DryMemory::free(nv * nv * nv * no, f, "vppph - hwm");

  DryMemory::allocate(n_tuples, sizeof(size_t), "tuples", SOURCE_LOCATION);
  DryMemory::allocate(n_tuples * 2, sizeof(size_t), "createTuples",
                      SOURCE_LOCATION);
  DryMemory::free(n_tuples * 2, sizeof(size_t), "createTuples");
  // these buffers are allocated for the main loop of atrip
  DryMemory::allocate(np * no * no * no, f, "tijk", SOURCE_LOCATION);
  DryMemory::allocate(np * no * no * no, f, "zijk", SOURCE_LOCATION);
  DryMemory::allocate(np * (no + nv), f, "epsp", SOURCE_LOCATION);
  DryMemory::allocate(np * no * nv, f, "tai", SOURCE_LOCATION);

  DryMemory::allocate(6 * np * nv * no * no, f, "taphh", SOURCE_LOCATION);
  DryMemory::allocate(6 * np * no * no * no, f, "vhhha", SOURCE_LOCATION);
  DryMemory::allocate(12 * np * nv * no, f, "vabph", SOURCE_LOCATION);
  DryMemory::allocate(6 * np * no * no, f, "vabhh", SOURCE_LOCATION);
  DryMemory::allocate(6 * np * no * no, f, "tabhh", SOURCE_LOCATION);

  if (arguments->isGiven(kwd::CT_INTERMEDIATE)) {
    DryMemory::allocate(6 * np * no * no * no, f, "jhhha", SOURCE_LOCATION);
    DryMemory::allocate(12 * np * nv * no, f, "jabph", SOURCE_LOCATION);
  }

  // Deallocate
  DryMemory::free(n_tuples, sizeof(size_t), "tuples");

  DryMemory::free(nv * nv * no * no, f, "tpphh2");
  DryMemory::free(nv * nv * nv * no, f, "vppph");

  DryMemory::free(np * no * no * no, f, "tijk");
  DryMemory::free(np * no * no * no, f, "zijk");
  DryMemory::free(np * (no + nv), f, "epsp");
  DryMemory::free(np * no * nv, f, "tai");

  DryMemory::free(6 * np * nv * no * no, f, "taphh");
  DryMemory::free(6 * np * no * no * no, f, "vhhha");
  DryMemory::free(12 * np * nv * no, f, "vabph");
  DryMemory::free(6 * np * no * no, f, "vabhh");
  DryMemory::free(6 * np * no * no, f, "tabhh");

  if (arguments->isGiven(kwd::CT_INTERMEDIATE)) {
    DryMemory::free(6 * np * no * no * no, f, "jhhha");
    DryMemory::free(12 * np * nv * no, f, "jabph");
  }

  //  DryMemory::allocate(atrip_memory, f, "Atrip-memory", SOURCE_LOCATION);
  //  DryMemory::free(atrip_memory, f, "Atrip-memory");
  double flopsPer = (typeid(F) == typeid(Complex<>)) ? 8.0 : 2.0;
  Operation<TE>::addFloatingPointOperations(no * no * no * nv * nv * nv *
                                            (no + nv) * flopsPer);
  if (arguments->isGiven(kwd::CT_INTERMEDIATE)) {
    Operation<TE>::addFloatingPointOperations(no * no * no * nv * nv * nv *
                                              (no + nv) * flopsPer);
  }

  return result;
}

Ptr<MapNode> PerturbativeTriples::run(const Ptr<MapNode> &arguments) {

  Ptr<MapNode> result;
  if (Cc4s::dryRun) {
    using TE = DefaultDryTensorEngine;
    if (arguments->getPtr<TensorSet<Real<>, TE>>("amplitudes") != nullptr) {
      result = atripDryRun<Real<>, TE>(arguments);
    } else {
      result = atripDryRun<Complex<>, TE>(arguments);
    }
  } else {
    using TE = DefaultTensorEngine;
    (result = runAtrip<Real<>, TE>(arguments)) ||
        (result = runAtrip<Complex<>, TE>(arguments));
  }

  return result;
}
