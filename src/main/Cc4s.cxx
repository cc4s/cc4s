/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <Cc4s.hpp>
#include <Parser.hpp>
#include <Emitter.hpp>
#include <algorithms/Algorithm.hpp>
#include <util/Timer.hpp>
#include <util/MpiCommunicator.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <tcc/Tcc.hpp>
#include <fstream>
#include <string>
#include <sstream>

using namespace cc4s;

size_t Cc4s::getFlops() {
  return options->dryRun ?
    Operation<DryTensorEngine>::flops : Operation<DefaultTensorEngine>::flops;
}

void Cc4s::run(const Ptr<MapNode> &report) {
  printBanner(report);
  listHosts(report);

  // start with empty storage
  storage = New<MapNode>(SOURCE_LOCATION);

  // parse input
  Parser parser(options->inFile);
  auto input(parser.parse());
  auto steps(input->toMap());
  ASSERT_LOCATION(steps, "expecting map as input", input->sourceLocation);
  report->get("steps") = steps;
  OUT() << "execution plan read, steps=" << steps->size() << std::endl;

  size_t totalFlops(-getFlops());
  Time totalTime;
  {
    // TODO: flops counter
    Timer totalTimer(&totalTime);

    for (unsigned int i(0); i < steps->size(); ++i) {
      auto step(steps->getMap(i));
      auto algorithmName(step->getSymbol("name"));
      OUT() << "step=" << (i+1) << ", " << algorithmName << std::endl;

      // create algorithm
      auto algorithm(AlgorithmFactory::create(algorithmName));
      ASSERT_LOCATION(
        algorithm, "unknown algorithm: " + algorithmName,
        step->get("name")->sourceLocation
      );

      // get input arguments
      auto inputArguments(step->getMap("in"));
      fetchSymbols(inputArguments);

      Ptr<MapNode> output;
      // wait until all processes finished previous work
      Cc4s::world->barrier();
      Time time;
      size_t flops(-getFlops());
      {
        // TODO: flops counter
        Timer timer(&time);
        output = algorithm->run(inputArguments);
        flops += getFlops();
      }

      // get output variables, if given
      if (step->get("out")) {
        auto outputVariables(step->getMap("out"));
        storeSymbols(output, outputVariables);
      }
      // store output in report
      step->get("out") = output;

      std::stringstream realtime;
      realtime << time;
      LOG() << "step=" << (i+1) << ", realtime=" << realtime.str() << " s"
        << ", operations=" << flops / 1e9 << " GFLOPS"
        << ", speed=" << flops / 1e9 / time.getFractionalSeconds() << " GFLOPS/s" << std::endl;
      step->setValue<std::string>("realtime", realtime.str());
      step->setValue<size_t>("floatingPointOperations", flops);
      step->setValue<Real<>>("flops", flops / time.getFractionalSeconds());
      // resources which maybe held by the algorithm automatically released
    }
  }
  totalFlops += getFlops();

  OUT() << std::endl;
  std::stringstream totalRealtime;
  totalRealtime << totalTime;
  OUT() << "total realtime=" << totalRealtime.str() << " s" << std::endl;
  OUT() << "total operations=" << totalFlops / 1e9 << " GFLOPS"
    << " speed=" << totalFlops/1e9 / totalTime.getFractionalSeconds() << " GFLOPS/s" << std::endl;
  LOG() << "total realtime=" << totalRealtime.str() << " s" << std::endl;
  LOG() << "total operations=" << totalFlops / 1e9 << " GFLOPS"
    << " speed=" << totalFlops/1e9 / totalTime.getFractionalSeconds() << " GFLOPS/s" << std::endl;
  report->setValue<std::string>("realtime", totalRealtime.str());
  report->setValue<size_t>("floatingPointOperations", totalFlops);
  report->setValue<Real<>>("flops", totalFlops/totalTime.getFractionalSeconds());
}

void Cc4s::fetchSymbols(const Ptr<MapNode> &arguments) {
  for (auto key: arguments->getKeys()) {
//    auto mapNode(arguments->get(key)->toMap());
//    if (mapNode) {
//      fetchSymbols(mapNode);
//      break;
//    }
    auto symbolNode(arguments->get(key)->toSymbol());
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
      WARNING_LOCATION(variables->sourceLocation) << "Symbol '" << key
        << "' should be stored but is not present in output" << std::endl;
    }
  }
}


void Cc4s::printBanner(const Ptr<MapNode> &report) {
  std::stringstream buildDate;
  buildDate << __DATE__ << " " << __TIME__;
  time_t rawtime;
  time (&rawtime);

  OUT() << std::endl
        << "                __ __      " << std::endl
        << "     __________/ // / _____" << std::endl
        << "    / ___/ ___/ // /_/ ___/" << std::endl
        << "   / /__/ /__/__  __(__  ) " << std::endl
        << "   \\___/\\___/  /_/ /____/  " << std::endl
        << "  Coupled Cluster for Solids" << std::endl << std::endl;
  OUT() << "version= " << CC4S_VERSION <<
    ", date=" << CC4S_DATE << std::endl;
  OUT() << "build date= " << buildDate.str() << std::endl;
  OUT() << "compiler= " << COMPILER_VERSION << std::endl;
  OUT() << "total processes= " << world->getProcesses() << std::endl;
  OUT() << "calculation started on: " << ctime (&rawtime) << std::endl << std::endl;
  report->setValue<std::string>("version", CC4S_VERSION);
  report->setValue("buildDate", buildDate.str());
  report->setValue<std::string>("compiler", COMPILER_VERSION);
  report->setValue("dryRun", options->dryRun);
/*
  EMIT()
    << YAML::Key << "version" << YAML::Value << CC4S_VERSION
    << YAML::Key << "build-date" << YAML::Value << buildDate.str()
    << YAML::Key << "compiler" << YAML::Value << COMPILER_VERSION
    << YAML::Key << "total-processes" << world->getProcesses();
*/
  if (options->dryRun) {
    OUT() << "DRY RUN - nothing will be calculated" << std::endl;
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
        LOG() << "Debugger present" << std::endl;
      }
      return pid > 0;
    }
  }
  return false;
}

void Cc4s::listHosts(const Ptr<MapNode> &report) {
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
      auto hosts(New<MapNode>(SOURCE_LOCATION));
      for (auto &ranksOfHost: ranksOfHosts) {
        auto host(New<MapNode>(SOURCE_LOCATION));
        host->setValue("host", ranksOfHost.first);
        auto ranks(New<MapNode>(SOURCE_LOCATION));
        for (auto &rank: ranksOfHost.second) {
          ranks->push_back(New<AtomicNode<int64_t>>(rank, SOURCE_LOCATION));
        }
        host->get("ranks") = ranks;
        hosts->push_back(host);
      }
      report->get("hosts") = hosts;
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
  Log::setFileName(Cc4s::options->name + ".log");
  Log::setRank(Cc4s::world->getRank());
  Time startTime(Time::getCurrentRealTime());
  Log::setLogHeaderFunction(
    [startTime](const SourceLocation &location) {
      std::stringstream header;
      header << (Time::getCurrentRealTime() - startTime) << ":";
      if (location.isValid()) {
        header << location;
      }
      header << ":";
      return header.str();
    }
  );
  auto report(New<MapNode>(SOURCE_LOCATION));
  bool errorHappened(false);

  Cc4s cc4s;
  if (Cc4s::isDebugged()) {
    // run without try-catch in debugger to allow tracing throwing code
    cc4s.run(report);
  } else {
    // without debugger: catch and write list of causes
    try {
      cc4s.run(report);
    } catch (Ptr<Exception> cause) {
      errorHappened = true;
      auto sourceLocation(cause->getSourceLocation());
      if (!sourceLocation.isValid()) sourceLocation = SOURCE_LOCATION;
      ERROR_LOCATION(sourceLocation) << cause->what() << std::endl;
      cause = cause->getCause();
      while (cause) {
        if (cause->getSourceLocation().isValid()) {
          sourceLocation = cause->getSourceLocation();
        }
        ERROR_LOCATION(sourceLocation) <<
          "Caused by: " << cause->what() << std::endl;
        cause = cause->getCause();
      }
    } catch (std::exception &cause) {
      errorHappened = true;
      OUT() << "unhandled exception encountered (std::exception):" << std::endl;
      OUT() << cause.what() << std::endl;
    } catch (const char *message) {
      errorHappened = true;
      OUT() << "unhandled exception encountered (const char *):" << std::endl;
      OUT() << message << std::endl;
    } catch (...) {
      errorHappened = true;
      OUT() << "unhandled exception encountered (...)." << std::endl;
    }
  }

  // emit report, also in case of error
  Emitter emitter(Cc4s::options->name + ".out");
  emitter.emit(report);

  MPI_Finalize();
  return errorHappened ? 1 : 0;
}

