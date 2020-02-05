/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <Cc4s.hpp>
#include <Parser.hpp>
#include <algorithms/Algorithm.hpp>
#include <util/Timer.hpp>
#include <tcc/DryTensor.hpp>
#include <util/FlopsCounter.hpp>
#include <util/MpiCommunicator.hpp>
#include <util/Log.hpp>
#include <util/Emitter.hpp>
#include <util/Exception.hpp>
#include <fstream>
#include <string>
#include <sstream>

// TODO: to be removed from the main class
#include <math/MathFunctions.hpp>

using namespace cc4s;
using namespace CTF;

Cc4s::Cc4s() {
}

Cc4s::~Cc4s() {
}

void Cc4s::run() {
  EMIT() << YAML::BeginMap;
  printBanner();
  Parser parser(options->file);
  std::vector<Algorithm *> algorithms(parser.parse());
  LOG(0, "root") <<
    "execution plan read, steps=" << algorithms.size() << std::endl;
  EMIT() <<
    YAML::Key << "execution-plan-size" << YAML::Value << algorithms.size();

  EMIT() << YAML::Key << "steps" << YAML::Value << YAML::BeginSeq;

  int64_t rootFlops, totalFlops;
  Time totalTime;
  {
    FlopsCounter rootCounter(&rootFlops);
    FlopsCounter totalCounter(&totalFlops, world->comm);
    Timer totalTimer(&totalTime);

    for (unsigned int i(0); i < algorithms.size(); ++i) {
      EMIT() << YAML::BeginMap;
      LOG(0, "root") << "step=" << (i+1) << ", " << algorithms[i]->getName() << std::endl;
      EMIT() << YAML::Key << "step" << YAML::Value << (i+1)
        << YAML::Key << "name" << YAML::Value << algorithms[i]->getName();

      int64_t flops;
      Time time;
      {
        FlopsCounter flopsCounter(&flops);
        Timer timer(&time);
        algorithms[i]->run();
        delete algorithms[i];
      }

      std::stringstream realtime;
      realtime << time;
      LOG(1, "root") << "step=" << (i+1) << ", realtime=" << realtime.str() << " s"
        << ", operations=" << flops / 1e9 << " GFLOPS/core"
        << ", speed=" << flops / 1e9 / time.getFractionalSeconds() << " GFLOPS/s/core" << std::endl;
      EMIT() << YAML::Key << "realtime" << YAML::Value << realtime.str()
        << YAML::Comment(" seconds")
        << YAML::Key << "floating-point-operations" << YAML::Value << flops
        << YAML::Comment("on root process")
        << YAML::Key << "flops" << YAML::Value << flops / time.getFractionalSeconds();
      printStatistics();
      EMIT() << YAML::EndMap;
    }
  }

  EMIT() << YAML::EndSeq;

  OUT() << std::endl;
  std::stringstream totalRealtime;
  totalRealtime << totalTime;
  LOG(0, "root") << "total realtime=" << totalRealtime.str() << " s" << std::endl;
  LOG(0, "root") << "total operations=" << rootFlops / 1e9 << " GFLOPS/core"
    << " speed=" << rootFlops/1e9 / totalTime.getFractionalSeconds() << " GFLOPS/s/core" << std::endl;
  LOG(0, "root") << "overall operations=" << totalFlops / 1.e9 << " GFLOPS"
    << std::endl;
  EMIT() << YAML::Key << "realtime" << YAML::Value << totalRealtime.str()
    << YAML::Key << "floating-point-operations" << YAML::Value << rootFlops
    << YAML::Comment("on root process")
    << YAML::Key << "flops" << YAML::Value << rootFlops / totalTime.getFractionalSeconds()
    << YAML::Key << "total-floating-point-operations" << totalFlops
    << YAML::Comment("of all processes");

  EMIT() << YAML::EndMap;
}

void Cc4s::dryRun() {
  EMIT() << YAML::BeginMap;
  printBanner();
  LOG(0, "root") <<
    "DRY RUN - nothing will be calculated" << std::endl;
  OUT() << std::endl;
  Parser parser(options->file);
  std::vector<Algorithm *> algorithms(parser.parse());
  LOG(0, "root") <<
    "execution plan read, steps=" << algorithms.size() << std::endl;
  EMIT() <<
    YAML::Key << "execution-plan-size" << YAML::Value << algorithms.size();

  EMIT() << YAML::Key << "steps" << YAML::Value << YAML::BeginSeq;

  for (unsigned int i(0); i < algorithms.size(); ++i) {
    EMIT() << YAML::BeginMap;
    LOG(0, "root") << "step=" << (i+1) << ", " << algorithms[i]->getName() << std::endl;
    EMIT() << YAML::Key << "step" << YAML::Value << (i+1)
      << YAML::Key << "name" << YAML::Value << algorithms[i]->getName();
    algorithms[i]->dryRun();
    LOG(0, "root")
      << "estimated memory=" << DryMemory::maxTotalSize / (1024.0*1024.0*1024.0)
      << " GB" << std::endl;
    EMIT()
      << YAML::Key << "estimated-total-memory"
      << YAML::Value << DryMemory::maxTotalSize / (1024.0*1024.0*1024.0)
      << YAML::Comment("GB");
    EMIT() << YAML::EndMap;
  }
  EMIT() << YAML::EndSeq;
  EMIT() << YAML::EndMap;
}


void Cc4s::printBanner() {
  std::stringstream buildDate;
  buildDate << __DATE__ << " " << __TIME__;

  OUT() << "                __ __      " << std::endl
        << "     __________/ // / _____" << std::endl
        << "    / ___/ ___/ // /_/ ___/" << std::endl
        << "   / /__/ /__/__  __(__  ) " << std::endl
        << "   \\___/\\___/  /_/ /____/  " << std::endl
        << "  Coupled Cluster for Solids" << std::endl << std::endl;
  LOG(0, "root") << "version=" << CC4S_VERSION <<
    ", date=" << CC4S_DATE << std::endl;
  LOG(0, "root") << "build date=" << buildDate.str() << std::endl;
  LOG(0, "root") << "compiler=" << COMPILER_VERSION << std::endl;
  LOG(0, "root") << "number of processes=" << Cc4s::world->np << std::endl;
  OUT() << std::endl;

  EMIT()
    << YAML::Key << "version" << YAML::Value << CC4S_VERSION
    << YAML::Key << "build-date" << YAML::Value << buildDate.str()
    << YAML::Key << "compiler" << YAML::Value << COMPILER_VERSION
    << YAML::Key << "number-of-processes" << Cc4s::world->np;
}

void Cc4s::printStatistics() {
  std::string fieldName;
  int64_t peakVirtualSize, peakPhysicalSize;
  // assuming LINUX
  std::ifstream statusStream("/proc/self/status", std::ios_base::in);
  std::string line;
  while (std::getline(statusStream, line)) {
    std::istringstream lineStream(line);
    lineStream >> fieldName;
    if (fieldName == "VmPeak:") {
      lineStream >> peakVirtualSize;
    } else if (fieldName == "VmHWM:") {
      lineStream >> peakPhysicalSize;
    }
    // TODO: check memory unit, currently assumed to be kB
  }
  statusStream.close();
  real unitsPerGB(1024.0*1024.0);
  LOG(0, "root") << "peak physical memory=" << peakPhysicalSize / unitsPerGB << " GB/core"
    << ", peak virtual memory: " << peakVirtualSize / unitsPerGB << " GB/core" << std::endl;
  EMIT()
    << YAML::Key << "peak-physical-memory" << YAML::Value << peakPhysicalSize / unitsPerGB
    << YAML::Comment("GB/core")
    << YAML::Key << "peak-virtual-memory" << YAML::Value <<  peakVirtualSize / unitsPerGB
    << YAML::Comment("GB/core");

  int64_t globalPeakVirtualSize, globalPeakPhysicalSize;
  MpiCommunicator communicator(world->rank, world->np, world->comm);
  communicator.reduce(peakPhysicalSize, globalPeakPhysicalSize);
  communicator.reduce(peakVirtualSize, globalPeakVirtualSize);
  LOG(0, "root") << "overall peak physical memory="
    << globalPeakPhysicalSize / unitsPerGB << " GB"
    << ", overall virtual memory=" << globalPeakVirtualSize / unitsPerGB << " GB" << std::endl;
  EMIT()
    << YAML::Key << "total-peak-physical-memory" << YAML::Value << globalPeakPhysicalSize / unitsPerGB
    << YAML::Comment("GB")
    << YAML::Key << "peak-virtual-memory" << YAML::Value <<  globalPeakVirtualSize / unitsPerGB
    << YAML::Comment("GB");
}

bool Cc4s::isDebugged() {
  // assuming LINUX
  std::ifstream statusStream("/proc/self/status", std::ios_base::in);
  std::string line;
  while (std::getline(statusStream, line)) {
    std::string pidField("TracerPid:");
    size_t position(line.find(pidField));
    if (position != std::string::npos) {
      std::stringstream pidStream(line.substr(position + pidField.length()));
      size_t pid; pidStream >> pid;
      if (pid > 0) {
        LOG(0, "root") << "Debugger present" << std::endl;
      }
      return pid > 0;
    }
  }
  return false;
}


World *Cc4s::world;
Options *Cc4s::options;


int main(int argumentCount, char **arguments) {
  MPI_Init(&argumentCount, &arguments);

  Cc4s::world = new World(argumentCount, arguments);
  Cc4s::options = new Options(argumentCount, arguments);
  Log::setRank(Cc4s::world->rank);
  Log::setFileName(Cc4s::options->logFile);
  Log::setLogLevel(Cc4s::options->logLevel);
  Emitter::setRank(Cc4s::world->rank);

  Cc4s cc4s;
  if (Cc4s::isDebugged()) {
    // run without try-catch in debugger to allow tracing throwing code
    if (Cc4s::options->dryRun) cc4s.dryRun();
    else cc4s.run();
  } else {
    // without debugger: catch and write exception cause
    try {
      if (Cc4s::options->dryRun) cc4s.dryRun();
      else cc4s.run();
    } catch (DetailedException *cause) {
      LOG(0) << std::endl << cause->getMessage() << std::endl;
    }
  }
  MPI_Finalize();
  return 0;
}

