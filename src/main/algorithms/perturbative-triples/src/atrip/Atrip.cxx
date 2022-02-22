// [[file:~/cc4s/src/atrip/bbbfb30/atrip.org::*Main][Main:1]]
#include <iomanip>

#include <atrip/Atrip.hpp>
#include <atrip/Utils.hpp>
#include <atrip/Equations.hpp>
#include <atrip/SliceUnion.hpp>
#include <atrip/Unions.hpp>

using namespace atrip;

bool RankMap<Complex>::RANK_ROUND_ROBIN;
bool RankMap<double>::RANK_ROUND_ROBIN;
int Atrip::rank;
int Atrip::np;
Timings Atrip::chrono;

// user printing block
IterationDescriptor IterationDescription::descriptor;
void atrip::registerIterationDescriptor(IterationDescriptor d) {
  IterationDescription::descriptor = d;
}

void Atrip::init()  {
  MPI_Comm_rank(MPI_COMM_WORLD, &Atrip::rank);
  MPI_Comm_size(MPI_COMM_WORLD, &Atrip::np);
}

template <typename F>
Atrip::Output Atrip::run(Atrip::Input<F> const& in) {

  const int np = Atrip::np;
  const int rank = Atrip::rank;
  MPI_Comm universe = in.ei->wrld->comm;

  const size_t No = in.ei->lens[0];
  const size_t Nv = in.ea->lens[0];
  LOG(0,"Atrip") << "No: " << No << "\n";
  LOG(0,"Atrip") << "Nv: " << Nv << "\n";
  LOG(0,"Atrip") << "np: " << np << "\n";

  // allocate the three scratches, see piecuch
  std::vector<F> Tijk(No*No*No) // doubles only (see piecuch)
               , Zijk(No*No*No) // singles + doubles (see piecuch)
               // we need local copies of the following tensors on every
               // rank
               , epsi(No)
               , epsa(Nv)
               , Tai(No * Nv)
               ;

  in.ei->read_all(epsi.data());
  in.ea->read_all(epsa.data());
  in.Tph->read_all(Tai.data());

  RankMap<F>::RANK_ROUND_ROBIN = in.rankRoundRobin;
  if (RankMap<F>::RANK_ROUND_ROBIN) {
    LOG(0,"Atrip") << "Doing rank round robin slices distribution" << "\n";
  } else {
    LOG(0,"Atrip")
      << "Doing node > local rank round robin slices distribution" << "\n";
  }


  // COMMUNICATOR CONSTRUCTION ========================================={{{1
  //
  // Construct a new communicator living only on a single rank
  int child_size = 1
    , child_rank
    ;
  const
  int color = rank / child_size
    , crank = rank % child_size
    ;
  MPI_Comm child_comm;
  if (np == 1) {
    child_comm = universe;
  } else {
    MPI_Comm_split(universe, color, crank, &child_comm);
    MPI_Comm_rank(child_comm, &child_rank);
    MPI_Comm_size(child_comm, &child_size);
  }


  // BUILD SLICES PARAMETRIZED BY NV ==================================={{{1
  WITH_CHRONO("nv-slices",
    LOG(0,"Atrip") << "BUILD NV-SLICES\n";
    TAPHH<F> taphh(*in.Tpphh, (size_t)No, (size_t)Nv, (size_t)np, child_comm, universe);
    HHHA<F>  hhha(*in.Vhhhp, (size_t)No, (size_t)Nv, (size_t)np, child_comm, universe);
  )

  // BUILD SLICES PARAMETRIZED BY NV x NV =============================={{{1
  WITH_CHRONO("nv-nv-slices",
    LOG(0,"Atrip") << "BUILD NV x NV-SLICES\n";
    ABPH<F> abph(*in.Vppph, (size_t)No, (size_t)Nv, (size_t)np, child_comm, universe);
    ABHH<F> abhh(*in.Vpphh, (size_t)No, (size_t)Nv, (size_t)np, child_comm, universe);
    TABHH<F> tabhh(*in.Tpphh, (size_t)No, (size_t)Nv, (size_t)np, child_comm, universe);
  )

  // all tensors
  std::vector< SliceUnion<F>* > unions = {&taphh, &hhha, &abph, &abhh, &tabhh};

  // get tuples for the current rank
  TuplesDistribution *distribution;

  if (in.tuplesDistribution == Atrip::Input<F>::TuplesDistribution::NAIVE) {
    LOG(0,"Atrip") << "Using the naive distribution\n";
    distribution = new NaiveDistribution();
  } else {
    LOG(0,"Atrip") << "Using the group-and-sort distribution\n";
    distribution = new group_and_sort::Distribution();
  }

  LOG(0,"Atrip") << "BUILDING TUPLE LIST\n";
  WITH_CHRONO("tuples:build",
    auto const tuplesList = distribution->getTuples(Nv, universe);
    )
  const size_t nIterations = tuplesList.size();

  {
    const size_t _all_tuples = Nv * (Nv + 1) * (Nv + 2) / 6 - Nv;
    LOG(0,"Atrip") << "#iterations: "
                  << nIterations
                  << "/"
                  << nIterations * np
                  << "\n";
  }

  const size_t
      iterationMod = (in.percentageMod > 0)
                  ? nIterations * in.percentageMod / 100
                  : in.iterationMod

    , iteration1Percent = nIterations * 0.01
    ;



  auto const isFakeTuple
    = [&tuplesList, distribution](size_t const i) {
      return distribution->tupleIsFake(tuplesList[i]);
    };


  using Database = typename Slice<F>::Database;
  using LocalDatabase = typename Slice<F>::LocalDatabase;
  auto communicateDatabase
    = [ &unions
      , np
      ] (ABCTuple const& abc, MPI_Comm const& c) -> Database {

        WITH_CHRONO("db:comm:type:do",
          auto MPI_LDB_ELEMENT = Slice<F>::mpi::localDatabaseElement();
        )

        WITH_CHRONO("db:comm:ldb",
          typename Slice<F>::LocalDatabase ldb;
          for (auto const& tensor: unions) {
            auto const& tensorDb = tensor->buildLocalDatabase(abc);
            ldb.insert(ldb.end(), tensorDb.begin(), tensorDb.end());
          }
        )

        Database db(np * ldb.size(), ldb[0]);

        WITH_CHRONO("oneshot-db:comm:allgather",
        WITH_CHRONO("db:comm:allgather",
          MPI_Allgather( ldb.data()
                       , ldb.size()
                       , MPI_LDB_ELEMENT
                       , db.data()
                       , ldb.size()
                       , MPI_LDB_ELEMENT
                       , c);
        ))

        WITH_CHRONO("db:comm:type:free",
          MPI_Type_free(&MPI_LDB_ELEMENT);
        )

        return db;
      };

  auto doIOPhase
    = [&unions, &rank, &np, &universe] (Database const& db) {

    const size_t localDBLength = db.size() / np;

    size_t sendTag = 0
         , recvTag = rank * localDBLength
         ;

    // RECIEVE PHASE ======================================================
    {
      // At this point, we have already send to everyone that fits
      auto const& begin = &db[rank * localDBLength]
                , end   = begin + localDBLength
                ;
      for (auto it = begin; it != end; ++it) {
        recvTag++;
        auto const& el = *it;
        auto& u = unionByName(unions, el.name);

        WITH_DBG std::cout
          << rank << ":r"
          << "♯" << recvTag << " =>"
          << " «n" << el.name
          << ", t" << el.info.type
          << ", s" << el.info.state
          << "»"
          << " ⊙ {" << rank << "⇐" << el.info.from.rank
                    << ", "
                    << el.info.from.source << "}"
          << " ∴ {" << el.info.tuple[0]
                    << ", "
                    << el.info.tuple[1]
                    << "}"
          << "\n"
          ;

        WITH_CHRONO("db:io:recv",
          u.receive(el.info, recvTag);
        )

      } // recv
    }

    // SEND PHASE =========================================================
    for (size_t otherRank = 0; otherRank<np; otherRank++) {
      auto const& begin = &db[otherRank * localDBLength]
                , end = begin + localDBLength
                ;
      for (auto it = begin; it != end; ++it) {
        sendTag++;
        typename Slice<F>::LocalDatabaseElement const& el = *it;

        if (el.info.from.rank != rank) continue;

        auto& u = unionByName(unions, el.name);
        WITH_DBG std::cout
          << rank << ":s"
          << "♯" << sendTag << " =>"
          << " «n" << el.name
          << ", t" << el.info.type
          << ", s" << el.info.state
          << "»"
          << " ⊙ {" << el.info.from.rank << "⇒" << otherRank
                    << ", "
                    << el.info.from.source << "}"
          << " ∴ {" << el.info.tuple[0]
                    << ", "
                    << el.info.tuple[1]
                    << "}"
          << "\n"
          ;

        WITH_CHRONO("db:io:send",
          u.send(otherRank, el, sendTag);
        )

      } // send phase

    } // otherRank


  };

#if defined(HAVE_OCD) || defined(ATRIP_PRINT_TUPLES)
  std::map<ABCTuple, double> tupleEnergies;
#endif

  const double doublesFlops
    = double(No)
    * double(No)
    * double(No)
    * (double(No) + double(Nv))
    * 2.0
    * (traits::isComplex<F>() ? 2.0 : 1.0)
    * 6.0
    / 1e9
    ;

  // START MAIN LOOP ======================================================{{{1

  double energy(0.);

  for ( size_t i = 0, iteration = 1
      ; i < tuplesList.size()
      ; i++, iteration++
      ) {
    Atrip::chrono["iterations"].start();

    // check overhead from chrono over all iterations
    WITH_CHRONO("start:stop", {})

    // check overhead of doing a barrier at the beginning
    WITH_CHRONO("oneshot-mpi:barrier",
    WITH_CHRONO("mpi:barrier",
      if (in.barrier) MPI_Barrier(universe);
    ))

    if (iteration % iterationMod == 0 || iteration == iteration1Percent) {

      if (IterationDescription::descriptor) {
        IterationDescription::descriptor({
          iteration,
          nIterations,
          Atrip::chrono["iterations"].count()
        });
      }

      LOG(0,"Atrip")
        << "iteration " << iteration
        << " [" << 100 * iteration / nIterations << "%]"
        << " (" << doublesFlops * iteration / Atrip::chrono["doubles"].count()
        << "GF)"
        << " (" << doublesFlops * iteration / Atrip::chrono["iterations"].count()
        << "GF)"
        << " ===========================\n";

      // PRINT TIMINGS
      if (in.chrono)
      for (auto const& pair: Atrip::chrono)
        LOG(1, " ") << pair.first << " :: "
                    << pair.second.count()
                    << std::endl;

    }

    const ABCTuple abc = isFakeTuple(i)
                       ? tuplesList[tuplesList.size() - 1]
                       : tuplesList[i]
                 , *abcNext = i == (tuplesList.size() - 1)
                            ? nullptr
                            : &tuplesList[i + 1]
                 ;

    WITH_CHRONO("with_rank",
      WITH_RANK << " :it " << iteration
                << " :abc " << pretty_print(abc)
                << " :abcN "
                << (abcNext ? pretty_print(*abcNext) : "None")
                << "\n";
    )


    // COMM FIRST DATABASE ================================================{{{1
    if (i == 0) {
      WITH_RANK << "__first__:first database ............ \n";
      const auto db = communicateDatabase(abc, universe);
      WITH_RANK << "__first__:first database communicated \n";
      WITH_RANK << "__first__:first database io phase \n";
      doIOPhase(db);
      WITH_RANK << "__first__:first database io phase DONE\n";
      WITH_RANK << "__first__::::Unwrapping all slices for first database\n";
      for (auto& u: unions) u->unwrapAll(abc);
      WITH_RANK << "__first__::::Unwrapping slices for first database DONE\n";
      MPI_Barrier(universe);
    }

    // COMM NEXT DATABASE ================================================={{{1
    if (abcNext) {
      WITH_RANK << "__comm__:" << iteration << "th communicating database\n";
      WITH_CHRONO("db:comm",
        const auto db = communicateDatabase(*abcNext, universe);
      )
      WITH_CHRONO("db:io",
        doIOPhase(db);
      )
      WITH_RANK << "__comm__:" <<  iteration << "th database io phase DONE\n";
    }

    // COMPUTE DOUBLES ===================================================={{{1
    OCD_Barrier(universe);
    if (!isFakeTuple(i)) {
      WITH_RANK << iteration << "-th doubles\n";
      WITH_CHRONO("oneshot-unwrap",
      WITH_CHRONO("unwrap",
      WITH_CHRONO("unwrap:doubles",
        for (auto& u: decltype(unions){&abph, &hhha, &taphh, &tabhh}) {
          u->unwrapAll(abc);
        }
      )))
      WITH_CHRONO("oneshot-doubles",
      WITH_CHRONO("doubles",
        doublesContribution<F>( abc, (size_t)No, (size_t)Nv
                              // -- VABCI
                              , abph.unwrapSlice(Slice<F>::AB, abc)
                              , abph.unwrapSlice(Slice<F>::AC, abc)
                              , abph.unwrapSlice(Slice<F>::BC, abc)
                              , abph.unwrapSlice(Slice<F>::BA, abc)
                              , abph.unwrapSlice(Slice<F>::CA, abc)
                              , abph.unwrapSlice(Slice<F>::CB, abc)
                              // -- VHHHA
                              , hhha.unwrapSlice(Slice<F>::A, abc)
                              , hhha.unwrapSlice(Slice<F>::B, abc)
                              , hhha.unwrapSlice(Slice<F>::C, abc)
                              // -- TA
                              , taphh.unwrapSlice(Slice<F>::A, abc)
                              , taphh.unwrapSlice(Slice<F>::B, abc)
                              , taphh.unwrapSlice(Slice<F>::C, abc)
                              // -- TABIJ
                              , tabhh.unwrapSlice(Slice<F>::AB, abc)
                              , tabhh.unwrapSlice(Slice<F>::AC, abc)
                              , tabhh.unwrapSlice(Slice<F>::BC, abc)
                              // -- TIJK
                              , Tijk.data()
                              );
        WITH_RANK << iteration << "-th doubles done\n";
      ))
    }

    // COMPUTE SINGLES =================================================== {{{1
    OCD_Barrier(universe);
    if (!isFakeTuple(i)) {
      WITH_CHRONO("oneshot-unwrap",
      WITH_CHRONO("unwrap",
      WITH_CHRONO("unwrap:singles",
        abhh.unwrapAll(abc);
      )))
      WITH_CHRONO("reorder",
        for (size_t I(0); I < Zijk.size(); I++) Zijk[I] = Tijk[I];
      )
      WITH_CHRONO("singles",
      singlesContribution<F>( No, Nv, abc
                            , Tai.data()
                            , abhh.unwrapSlice(Slice<F>::AB, abc)
                            , abhh.unwrapSlice(Slice<F>::AC, abc)
                            , abhh.unwrapSlice(Slice<F>::BC, abc)
                            , Zijk.data());
      )
    }


    // COMPUTE ENERGY ==================================================== {{{1
    if (!isFakeTuple(i)) {
      double tupleEnergy(0.);

      int distinct(0);
      if (abc[0] == abc[1]) distinct++;
      if (abc[1] == abc[2]) distinct--;
      const F epsabc(epsa[abc[0]] + epsa[abc[1]] + epsa[abc[2]]);

      WITH_CHRONO("energy",
        if ( distinct == 0)
          tupleEnergy = getEnergyDistinct<F>(epsabc, epsi, Tijk, Zijk);
        else
          tupleEnergy = getEnergySame<F>(epsabc, epsi, Tijk, Zijk);
      )

#if defined(HAVE_OCD) || defined(ATRIP_PRINT_TUPLES)
      tupleEnergies[abc] = tupleEnergy;
#endif

      energy += tupleEnergy;

    }

    // TODO: remove this
    if (isFakeTuple(i)) {
      // fake iterations should also unwrap whatever they got
      WITH_RANK << iteration
                << "th unwrapping because of fake in "
                << i << "\n";
      for (auto& u: unions) u->unwrapAll(abc);
    }

#ifdef HAVE_OCD
    for (auto const& u: unions) {
      WITH_RANK << "__dups__:"
                << iteration
                << "-th n" << u->name << " checking duplicates\n";
      u->checkForDuplicates();
    }
#endif


    // CLEANUP UNIONS ===================================================={{{1
    OCD_Barrier(universe);
    if (abcNext) {
      WITH_RANK << "__gc__:" << iteration << "-th cleaning up.......\n";
      for (auto& u: unions) {

        u->unwrapAll(abc);
        WITH_RANK << "__gc__:n" << u->name  << " :it " << iteration
                  << " :abc " << pretty_print(abc)
                  << " :abcN " << pretty_print(*abcNext)
                  << "\n";
        // for (auto const& slice: u->slices)
        //   WITH_RANK << "__gc__:guts:" << slice.info << "\n";
        u->clearUnusedSlicesForNext(*abcNext);

        WITH_RANK << "__gc__: checking validity\n";

#ifdef HAVE_OCD
        // check for validity of the slices
        for (auto type: u->sliceTypes) {
          auto tuple = Slice<F>::subtupleBySlice(abc, type);
        for (auto& slice: u->slices) {
          if ( slice.info.type == type
             && slice.info.tuple == tuple
             && slice.isDirectlyFetchable()
             ) {
            if (slice.info.state == Slice<F>::Dispatched)
              throw std::domain_error( "This slice should not be undispatched! "
                                     + pretty_print(slice.info));
          }
        }
        }
#endif


      }
    }

      WITH_RANK << iteration << "-th cleaning up....... DONE\n";

    Atrip::chrono["iterations"].stop();
    // ITERATION END ====================================================={{{1

  }
    // END OF MAIN LOOP

  MPI_Barrier(universe);

  // PRINT TUPLES ========================================================={{{1
#if defined(HAVE_OCD) || defined(ATRIP_PRINT_TUPLES)
  LOG(0,"Atrip") << "tuple energies" << "\n";
  for (size_t i = 0; i < np; i++) {
    MPI_Barrier(universe);
    for (auto const& pair: tupleEnergies) {
      if (i == rank)
        std::cout << pair.first[0]
                  << " " << pair.first[1]
                  << " " << pair.first[2]
                  << std::setprecision(15) << std::setw(23)
                  << " tupleEnergy: " << pair.second
                  << "\n"
                  ;
    }
  }
#endif

  // COMMUNICATE THE ENERGIES ============================================={{{1
  LOG(0,"Atrip") << "COMMUNICATING ENERGIES \n";
  double globalEnergy = 0;
  MPI_Reduce(&energy, &globalEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, universe);

  WITH_RANK << "local energy " << energy << "\n";
  LOG(0, "Atrip") << "Energy: "
    << std::setprecision(15) << std::setw(23)
    << (- globalEnergy) << std::endl;

  // PRINT TIMINGS {{{1
  if (in.chrono)
  for (auto const& pair: Atrip::chrono)
    LOG(0,"atrip:chrono") << pair.first << " "
                          << pair.second.count() << std::endl;


  LOG(0, "atrip:flops(doubles)")
    << nIterations * doublesFlops / Atrip::chrono["doubles"].count() << "\n";
  LOG(0, "atrip:flops(iterations)")
    << nIterations * doublesFlops / Atrip::chrono["iterations"].count() << "\n";

  // TODO: change the sign in  the getEnergy routines
  return { - globalEnergy };

}
// instantiate
template Atrip::Output Atrip::run(Atrip::Input<double> const& in);
template Atrip::Output Atrip::run(Atrip::Input<Complex> const& in);
// Main:1 ends here
