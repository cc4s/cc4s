/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <Cc4s.hpp>
#include <Parser.hpp>
#include <algorithms/Algorithm.hpp>
#include <util/Timer.hpp>
#include <tcc/DryTensor.hpp>
#include <util/FlopsCounter.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <fstream>

// TODO: to be removed from the main class
#include <math/MathFunctions.hpp>

using namespace cc4s;
using namespace CTF;

Cc4s::Cc4s() {
}

Cc4s::~Cc4s() {
}

void Cc4s::run() {
  printBanner();
  Parser parser(options->file);
  std::vector<Algorithm *> algorithms(parser.parse());
  LOG(0, "root") <<
    "execution plan read, steps=" << algorithms.size() << std::endl;

  int64_t rootFlops, totalFlops;
  Time totalTime;
  {
    FlopsCounter rootCounter(&rootFlops);
    FlopsCounter totalCounter(&totalFlops, world->comm);
    Timer totalTimer(&totalTime);

    for (unsigned int i(0); i < algorithms.size(); ++i) {
      LOG(0, "root") << "step=" << (i+1) << ", " << algorithms[i]->getName() << std::endl;

      int64_t flops;
      Time time;
      {
        FlopsCounter flopsCounter(&flops);
        Timer timer(&time);
        algorithms[i]->run();
        delete algorithms[i];
      }

      LOG(1, "root") << "step=" << (i+1) << ", realtime=" << time << " s"
        << ", operations=" << flops / 1e9 << " GFLOPS/core"
        << ", speed=" << flops / 1e9 / time.getFractionalSeconds() << " GFLOPS/s/core" << std::endl;
    }
  }

  printStatistics(rootFlops, totalFlops, totalTime);
}

void Cc4s::dryRun() {
  printBanner();
  LOG(0, "root") <<
    "DRY RUN - nothing will be calculated" << std::endl;
  OUT() << std::endl;
  Parser parser(options->file);
  std::vector<Algorithm *> algorithms(parser.parse());
  LOG(0, "root") <<
    "execution plan read, steps=" << algorithms.size() << std::endl;

  for (unsigned int i(0); i < algorithms.size(); ++i) {
    LOG(0, "root") << "step=" << (i+1) << ", " << algorithms[i]->getName() << std::endl;
    algorithms[i]->dryRun();
  }

  LOG(0, "root")
    << "estimated memory=" << DryMemory::maxTotalSize / (1024.0*1024.0*1024.0)
    << " GB" << std::endl;
}


void Cc4s::printBanner() {
  OUT() << "                __ __      " << std::endl
        << "     __________/ // / _____" << std::endl
        << "    / ___/ ___/ // /_/ ___/" << std::endl
        << "   / /__/ /__/__  __(__  ) " << std::endl
        << "   \\___/\\___/  /_/ /____/  " << std::endl
        << "  Coupled Cluster for Solids" << std::endl << std::endl;
  LOG(0, "root") << "version=" << CC4S_VERSION <<
    ", date=" << CC4S_DATE << std::endl;
  LOG(0, "root") << "build date=" << __DATE__ << " " << __TIME__ << std::endl;
  LOG(0, "root") << "compiler=" << COMPILER_VERSION << std::endl;
  LOG(0, "root") << "number of processes=" << Cc4s::world->np << std::endl;
  OUT() << std::endl;
}

void Cc4s::printStatistics(
  int64_t rootFlops, int64_t totalFlops, Time const &totalTime
) {
  std::string pid, comm, state, ppid, pgrp, session, ttyNr,
    tpgid, flags, minflt, cminflt, majflt, cmajflt,
    utime, stime, cutime, cstime, priority, nice,
    O, itrealvalue, starttime;
  int64_t vsize, rss;
  // assuming LINUX 
  std::ifstream statStream("/proc/self/stat", std::ios_base::in);
  statStream >> pid >> comm >> state >> ppid >> pgrp >> session >> ttyNr
    >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
    >> utime >> stime >> cutime >> cstime >> priority >> nice
    >> O >> itrealvalue >> starttime >> vsize >> rss;
  statStream.close();
  // in case x86-64 is configured to use 2MB pages
  int64_t pageSize = sysconf(_SC_PAGE_SIZE);
  LOG(0, "root") << "total realtime=" << totalTime << " s" << std::endl;
  LOG(0, "root") << "total operations=" << rootFlops / 1e9 << " GFLOPS/core"
    << " speed=" << rootFlops/1e9 / totalTime.getFractionalSeconds() << " GFLOPS/s/core" << std::endl;
  LOG(0, "root") << "physical memory=" << rss * pageSize / 1e9 << " GB/core"
    << ", virtual memory: " << vsize / 1e9 << " GB/core" << std::endl;

  int64_t globalVSize, globalRss;
  MPI_Reduce(&vsize, &globalVSize, 1, MPI_LONG_LONG, MPI_SUM, 0, world->comm);
  MPI_Reduce(&rss, &globalRss, 1, MPI_LONG_LONG, MPI_SUM, 0, world->comm);
  LOG(0, "root") << "overall operations=" << totalFlops / 1.e9 << " GFLOPS"
    << std::endl;
  LOG(0, "root") << "overall physical memory="
    << globalRss * pageSize / 1e9 << " GB"
    << ", overall virtual memory=" << globalVSize / 1e9 << " GB" << std::endl;
}


World *Cc4s::world;
Options *Cc4s::options;


int main(int argumentCount, char **arguments) {
  MPI_Init(&argumentCount, &arguments);

  Cc4s::world = new World(argumentCount, arguments);
  Cc4s::options = new Options(argumentCount, arguments);
  Log::setRank(Cc4s::world->rank);
  LogStream logStream(
    Cc4s::options->logFile, Cc4s::options->logLevel
  );
  Log::setLogStream(&logStream);

//  try {
    Cc4s cc4s;
    if (Cc4s::options->dryRun) cc4s.dryRun();
    else cc4s.run();
//  } catch (DetailedException *cause) {
//    LOG(0) << std::endl << cause->getMessage() << std::endl;
//  }
  MPI_Finalize();
  return 0;
}

