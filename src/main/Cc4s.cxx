/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <Cc4s.hpp>
#include <Parser.hpp>
#include <Emitter.hpp>
#include <algorithms/Algorithm.hpp>
#include <util/Timer.hpp>
#include <util/MpiCommunicator.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <fstream>
#include <string>
#include <sstream>

using namespace cc4s;

void Cc4s::run() {
  // TODO: truncate output.yaml or write it also in case of failure
  auto job(New<MapNode>());
  printBanner(job);
  listHosts(job);

  // start with empty storage
  storage = New<MapNode>();

  // parse input
  Parser parser(options->inFile);
  auto steps(parser.parse()->map());
  Assert(steps, "expecting map as input");
  job->get("steps") = steps;
  LOG(0, "Cc4s") <<
    "execution plan read, steps=" << steps->size() << std::endl;

  size_t rootFlops, totalFlops;
  Time totalTime;
  {
    // TODO: flops counter
    Timer totalTimer(&totalTime);

    for (unsigned int i(0); i < steps->size(); ++i) {
      auto step(steps->getMap(i));
      auto algorithmName(step->getSymbol("name"));
      LOG(0, "Cc4s") << "step=" << (i+1) << ", " << algorithmName << std::endl;

      // create algorithm
      auto algorithm(AlgorithmFactory::create(algorithmName));
      Assert(algorithm, "unknown algorithm: " + algorithmName);

      // get input arguments
      auto inputArguments(step->getMap("in"));
      fetchSymbols(inputArguments);

      size_t flops;
      Ptr<MapNode> output;
      Time time;
      {
        // TODO: flops counter
        Timer timer(&time);
        output = algorithm->run(inputArguments);
      }

      // get output variables, if given
      if (step->get("out")) {
        auto outputVariables(step->getMap("out"));
        storeSymbols(output, outputVariables);
      }
      // store output in job data
      step->get("out") = output;

      std::stringstream realtime;
      realtime << time;
      LOG(1, "Cc4s") << "step=" << (i+1) << ", realtime=" << realtime.str() << " s"
        << ", operations=" << flops / 1e9 << " GFLOPS/core"
        << ", speed=" << flops / 1e9 / time.getFractionalSeconds() << " GFLOPS/s/core" << std::endl;
      // TODO: enter data in node tree
/*
      EMIT() << YAML::Key << "realtime" << YAML::Value << realtime.str()
        << YAML::Comment(" seconds")
        << YAML::Key << "floating-point-operations" << YAML::Value << flops
        << YAML::Comment("on root process")
        << YAML::Key << "flops" << YAML::Value << flops / time.getFractionalSeconds();
      EMIT() << YAML::EndMap;
*/
      printStatistics(job);
      // resources which maybe held by the algorithm automatically released
    }
  }

  OUT() << std::endl;
  std::stringstream totalRealtime;
  totalRealtime << totalTime;
  LOG(0, "Cc4s") << "total realtime=" << totalRealtime.str() << " s" << std::endl;
  LOG(0, "Cc4s") << "total operations=" << rootFlops / 1e9 << " GFLOPS/core"
    << " speed=" << rootFlops/1e9 / totalTime.getFractionalSeconds() << " GFLOPS/s/core" << std::endl;
  LOG(0, "Cc4s") << "overall operations=" << totalFlops / 1.e9 << " GFLOPS"
    << std::endl;
  // TODO: flops counter
/*
  EMIT() << YAML::Key << "realtime" << YAML::Value << totalRealtime.str()
    << YAML::Key << "floating-point-operations" << YAML::Value << rootFlops
    << YAML::Comment("on root process")
    << YAML::Key << "flops" << YAML::Value << rootFlops / totalTime.getFractionalSeconds()
    << YAML::Key << "total-floating-point-operations" << totalFlops
    << YAML::Comment("of all processes");

  EMIT() << YAML::EndMap;
*/

  // emit job output
  Emitter emitter(options->name + ".yaml");
  emitter.emit(job);
}

void Cc4s::fetchSymbols(const Ptr<MapNode> &arguments) {
  for (auto key: arguments->getKeys()) {
//    auto mapNode(arguments->get(key)->map());
//    if (mapNode) {
//      fetchSymbols(mapNode);
//      break;
//    }
    auto symbolNode(arguments->get(key)->symbol());
    if (symbolNode) {
      // search symbol in storage
      auto storedNode(storage->get(symbolNode->value));
      // if found, replace symbol with stored node
      if (storedNode) {
        arguments->get(key) = storedNode;
      }
    }
  }
}

void Cc4s::storeSymbols(const Ptr<MapNode> &result, const Ptr<MapNode> &variables) {
  for (auto key: variables->getKeys()) {
    // search key in result
    if (result->get(key)) {
      // store key's value node in storage under given symbol name
      auto symbolName(variables->getSymbol(key));
      storage->get(symbolName) = result->get(key);
    } else {
      LOG(0,"Storage") << "Symbol '" << key
        << "' should be stored but is not present in output" << std::endl;
    }
  }
}


void Cc4s::printBanner(const Ptr<MapNode> &job) {
  std::stringstream buildDate;
  buildDate << __DATE__ << " " << __TIME__;

  OUT() << "                __ __      " << std::endl
        << "     __________/ // / _____" << std::endl
        << "    / ___/ ___/ // /_/ ___/" << std::endl
        << "   / /__/ /__/__  __(__  ) " << std::endl
        << "   \\___/\\___/  /_/ /____/  " << std::endl
        << "  Coupled Cluster for Solids" << std::endl << std::endl;
  LOG(0, "Cc4s") << "version=" << CC4S_VERSION <<
    ", date=" << CC4S_DATE << std::endl;
  LOG(0, "Cc4s") << "build date=" << buildDate.str() << std::endl;
  LOG(0, "Cc4s") << "compiler=" << COMPILER_VERSION << std::endl;
  LOG(0, "Cc4s") << "total processes=" << world->getProcesses() << std::endl;
  OUT() << std::endl;

  job->setValue<std::string>("version", CC4S_VERSION);
  job->setValue("buildDate", buildDate.str());
  job->setValue<std::string>("compiler", COMPILER_VERSION);
  job->setValue("dryRun", options->dryRun);
/*
  EMIT()
    << YAML::Key << "version" << YAML::Value << CC4S_VERSION
    << YAML::Key << "build-date" << YAML::Value << buildDate.str()
    << YAML::Key << "compiler" << YAML::Value << COMPILER_VERSION
    << YAML::Key << "total-processes" << world->getProcesses();
*/
  if (options->dryRun) {
    LOG(0, "Cc4s") <<
      "DRY RUN - nothing will be calculated" << std::endl;
  }
}

void Cc4s::printStatistics(const Ptr<MapNode> &job) {
  if (options->dryRun) {
    LOG(0, "Cc4s")
      << "estimated memory=" << DryMemory::maxTotalSize / (1024.0*1024.0*1024.0)
      << " GB" << std::endl;
/*
    EMIT()
      << YAML::Key << "estimated-total-memory"
      << YAML::Value << DryMemory::maxTotalSize / (1024.0*1024.0*1024.0)
      << YAML::Comment("GB");
*/
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
    LOG(0, "Cc4s") << "peak physical memory=" << peakPhysicalSize / unitsPerGB << " GB/core"
      << ", peak virtual memory: " << peakVirtualSize / unitsPerGB << " GB/core" << std::endl;
/*
    EMIT()
      << YAML::Key << "peak-physical-memory" << YAML::Value << peakPhysicalSize / unitsPerGB
      << YAML::Comment("GB/core")
      << YAML::Key << "peak-virtual-memory" << YAML::Value <<  peakVirtualSize / unitsPerGB
      << YAML::Comment("GB/core");
*/
    int64_t globalPeakVirtualSize, globalPeakPhysicalSize;
    world->reduce(peakPhysicalSize, globalPeakPhysicalSize);
    world->reduce(peakVirtualSize, globalPeakVirtualSize);
    LOG(0, "Cc4s") << "overall peak physical memory="
      << globalPeakPhysicalSize / unitsPerGB << " GB"
      << ", overall virtual memory=" << globalPeakVirtualSize / unitsPerGB << " GB" << std::endl;
/*
    EMIT()
      << YAML::Key << "total-peak-physical-memory" << YAML::Value << globalPeakPhysicalSize / unitsPerGB
      << YAML::Comment("GB")
      << YAML::Key << "peak-virtual-memory" << YAML::Value <<  globalPeakVirtualSize / unitsPerGB
      << YAML::Comment("GB");
*/
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
        LOG(0, "Cc4s") << "Debugger present" << std::endl;
      }
      return pid > 0;
    }
  }
  return false;
}

void Cc4s::listHosts(const Ptr<MapNode> &job) {
  // TODO: list planned execution environment from options
  if (!options->dryRun) {
    char ownName[MPI_MAX_PROCESSOR_NAME];
    int nameLength;
    MPI_Get_processor_name(ownName, &nameLength);
    ownName[nameLength] = 0;

    std::map<std::string,std::vector<int>> ranksOfHosts;
    if (world->getRank() == 0) {
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
      auto hosts(New<MapNode>());
      for (auto &ranksOfHost: ranksOfHosts) {
        auto host(New<MapNode>());
        host->setValue("host", ranksOfHost.first);
        auto ranks(New<MapNode>());
        for (auto &rank: ranksOfHost.second) {
          ranks->push_back(New<AtomicNode<int64_t>>(rank));
        }
        host->get("ranks") = ranks;
        hosts->push_back(host);
      }
      job->get("hosts") = hosts;
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
  Log::setFileName(Cc4s::options->name + ".log");
  Log::setLogLevel(Cc4s::options->logLevel);

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

  MPI_Finalize();
  return 0;
}

