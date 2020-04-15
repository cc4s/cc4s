/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <Cc4s.hpp>
#include <Parser.hpp>
#include <algorithms/Algorithm.hpp>
#include <util/Timer.hpp>
#include <util/MpiCommunicator.hpp>
#include <util/Log.hpp>
#include <util/Emitter.hpp>
#include <util/Exception.hpp>
#include <fstream>
#include <string>
#include <sstream>

using namespace cc4s;

void Cc4s::run() {
  EMIT() << YAML::BeginMap;
  printBanner();
  listHosts();

  // parse input
  Parser parser(options->inFile);
  auto root(parser.parse()->map());
  Assert(root, "expecting map as input root");
  auto algorithms(root->getMap("algorithms"));
  LOG(0, "root") <<
    "execution plan read, steps=" << algorithms->size() << std::endl;
  EMIT() <<
    YAML::Key << "execution-plan-size" << YAML::Value << algorithms->size();

  EMIT() << YAML::Key << "steps" << YAML::Value << YAML::BeginSeq;

  size_t rootFlops, totalFlops;
  Time totalTime;
  {
    // TODO: flops counter
    Timer totalTimer(&totalTime);

    for (unsigned int i(0); i < algorithms->size(); ++i) {
      EMIT() << YAML::BeginMap;
      auto algorithmNode(algorithms->getMap(i));
      auto algorithmName(algorithmNode->getSymbol("name"));
      LOG(0, "root") << "step=" << (i+1) << ", " << algorithmName << std::endl;
      EMIT() << YAML::Key << "step" << YAML::Value << (i+1)
        << YAML::Key << "name" << YAML::Value << algorithmName;

      // create algorithm
      auto algorithm(AlgorithmFactory::create(algorithmName));
      Assert(algorithm, "unknown algorithm: " + algorithmName);

      // get input arguments
      auto inputArguments(algorithmNode->getMap("in"));
      // FIXME: move input arguments from storage

      size_t flops;
      Time time;
      {
        // TODO: flops counter
        Timer timer(&time);
        auto output(
          options->dryRun ?
            algorithm->dryRun(inputArguments) : algorithm->run(inputArguments)
        );
      }

      // get output arguments
      if (algorithmNode->get("out")) {
        auto outputArguments(algorithmNode->getMap("out"));
        // FIXME: move output arguments to storage
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
      // resources which maybe held by the algorithm automatically released
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
  // TODO: flops counter
  EMIT() << YAML::Key << "realtime" << YAML::Value << totalRealtime.str()
    << YAML::Key << "floating-point-operations" << YAML::Value << rootFlops
    << YAML::Comment("on root process")
    << YAML::Key << "flops" << YAML::Value << rootFlops / totalTime.getFractionalSeconds()
    << YAML::Key << "total-floating-point-operations" << totalFlops
    << YAML::Comment("of all processes");

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
  LOG(0, "root") << "total processes=" << world->getProcesses() << std::endl;
  OUT() << std::endl;

  EMIT()
    << YAML::Key << "version" << YAML::Value << CC4S_VERSION
    << YAML::Key << "build-date" << YAML::Value << buildDate.str()
    << YAML::Key << "compiler" << YAML::Value << COMPILER_VERSION
    << YAML::Key << "total-processes" << world->getProcesses();

  if (options->dryRun) {
    LOG(0, "root") <<
      "DRY RUN - nothing will be calculated" << std::endl;
  }
}

void Cc4s::printStatistics() {
  if (options->dryRun) {
    LOG(0, "root")
      << "estimated memory=" << DryMemory::maxTotalSize / (1024.0*1024.0*1024.0)
      << " GB" << std::endl;
    EMIT()
      << YAML::Key << "estimated-total-memory"
      << YAML::Value << DryMemory::maxTotalSize / (1024.0*1024.0*1024.0)
      << YAML::Comment("GB");
  } else {
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
    Real<> unitsPerGB(1024.0*1024.0);
    LOG(0, "root") << "peak physical memory=" << peakPhysicalSize / unitsPerGB << " GB/core"
      << ", peak virtual memory: " << peakVirtualSize / unitsPerGB << " GB/core" << std::endl;
    EMIT()
      << YAML::Key << "peak-physical-memory" << YAML::Value << peakPhysicalSize / unitsPerGB
      << YAML::Comment("GB/core")
      << YAML::Key << "peak-virtual-memory" << YAML::Value <<  peakVirtualSize / unitsPerGB
      << YAML::Comment("GB/core");

    int64_t globalPeakVirtualSize, globalPeakPhysicalSize;
    world->reduce(peakPhysicalSize, globalPeakPhysicalSize);
    world->reduce(peakVirtualSize, globalPeakVirtualSize);
    LOG(0, "root") << "overall peak physical memory="
      << globalPeakPhysicalSize / unitsPerGB << " GB"
      << ", overall virtual memory=" << globalPeakVirtualSize / unitsPerGB << " GB" << std::endl;
    EMIT()
      << YAML::Key << "total-peak-physical-memory" << YAML::Value << globalPeakPhysicalSize / unitsPerGB
      << YAML::Comment("GB")
      << YAML::Key << "peak-virtual-memory" << YAML::Value <<  globalPeakVirtualSize / unitsPerGB
      << YAML::Comment("GB");
  }
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

void Cc4s::listHosts() {
  // TODO: list planned execution environment from options
  if (!options->dryRun) {
    char ownName[MPI_MAX_PROCESSOR_NAME];
    int nameLength;
    MPI_Get_processor_name(ownName, &nameLength);
    ownName[nameLength] = 0;

    if (world->getRank() == 0) {
      std::map<std::string,std::vector<int>> ranksOfHosts;
      // enter hostname/rank
      ranksOfHosts[ownName].push_back(0);
      // receive all names but own name from remote ranks
      for (int remoteRank(1); remoteRank < world->getProcesses(); ++remoteRank) {
        char remoteName[MPI_MAX_PROCESSOR_NAME];
        MPI_Recv(
          remoteName, MPI_MAX_PROCESSOR_NAME, MPI_BYTE,
          remoteRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE
        );
        // and enter remote name/rank into map
        ranksOfHosts[remoteName].push_back(remoteRank);
      }
      EMIT() << YAML::Key << "hosts" << YAML::Value;
      EMIT() << YAML::BeginSeq;
      for (auto &ranksOfHost: ranksOfHosts) {
        EMIT() << YAML::BeginMap
          << YAML::Key << "host" << YAML::Value << ranksOfHost.first
          << YAML::Key << "ranks" << YAML::Value;
        EMIT() << YAML::Flow << YAML::BeginSeq;
        for (auto &rank: ranksOfHost.second) {
          EMIT() << rank;
        }
        EMIT() << YAML::EndSeq;
        EMIT() << YAML::EndMap;
      }
      EMIT() << YAML::EndSeq;
    } else {
      // send own name 
      MPI_Send(
        ownName, MPI_MAX_PROCESSOR_NAME, MPI_BYTE,
        0, 0, MPI_COMM_WORLD
      );
    }
  }
}


Ptr<MpiCommunicator> Cc4s::world;
Ptr<Options> Cc4s::options;


int main(int argumentCount, char **arguments) {
  MPI_Init(&argumentCount, &arguments);

  Cc4s::world = New<MpiCommunicator>();
  Cc4s::options = New<Options>(argumentCount, arguments);
  Log::setRank(Cc4s::world->getRank());
  Log::setFileName(Cc4s::options->logFile);
  Log::setLogLevel(Cc4s::options->logLevel);
  Emitter::setFileName(Cc4s::options->outFile);
  Emitter::setRank(Cc4s::world->getRank());

  Cc4s cc4s;
  if (Cc4s::isDebugged()) {
    // run without try-catch in debugger to allow tracing throwing code
    cc4s.run();
  } else {
    // without debugger: catch and write exception cause
    try {  
      cc4s.run();
    } catch (DetailedException *cause) {
      LOG(0) << std::endl << cause->getMessage() << std::endl;
    }
  }

  EMIT_FLUSH();
  MPI_Finalize();
  return 0;
}

