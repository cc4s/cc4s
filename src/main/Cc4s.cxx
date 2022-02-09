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
#include <algorithms/Algorithm.hpp>
#include <tcc/Tcc.hpp>
#include <Parser.hpp>
#include <Emitter.hpp>
#include <Timer.hpp>
#include <MpiCommunicator.hpp>
#include <Log.hpp>
#include <Exception.hpp>

#include <fstream>
#include <string>
#include <sstream>

using namespace cc4s;

Natural<128> Cc4s::getFloatingPointOperations() {
  return dryRun ?
    Operation<DefaultDryTensorEngine>::getFloatingPointOperations() :
    Operation<DefaultTensorEngine>::getFloatingPointOperations();
}

void Cc4s::addFloatingPointOperations(const Natural<128> ops) {
  return dryRun ?
    Operation<DefaultDryTensorEngine>::addFloatingPointOperations(ops) :
    Operation<DefaultTensorEngine>::addFloatingPointOperations(ops);
}

void Cc4s::runStep(Natural<> i, const Ptr<MapNode> &step) {
  auto algorithmName(step->getSymbol("name"));
  OUT() << "step: " << (i+1) << ", " << algorithmName << std::endl;

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
  Natural<128> operations;
  Time time;
  {
    OperationsCounter operationsCounter(&operations);
    Timer timer(&time);
    output = algorithm->run(inputArguments);
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
  OUT() << "realtime " << realtime.str() << " s" << std::endl;
  OUT() << "--" << std::endl;
  LOG() << "step: " << (i+1) << ", realtime: " << realtime.str() << " s"
    << ", operations: " << operations / 1e9 << " GFLOP"
    << ", speed: " << operations / 1e9 / time.getFractionalSeconds() << " GFLOP/s" << std::endl;
  step->setValue("realtime", realtime.str());
  step->setValue("floatingPointOperations", operations);
  step->setValue("flops", operations / time.getFractionalSeconds());
  // resources which maybe held by the algorithm automatically released
}

void Cc4s::run(const Ptr<MapNode> &report) {
  printBanner(report);
  listHosts(report);

  // start with empty storage
  storage = New<MapNode>(SOURCE_LOCATION);

  // parse input
  Parser parser(options->inFile);
  auto input(parser.parse());
  auto steps(input->toPtr<MapNode>());
  ASSERT_LOCATION(steps, "expecting map as input", input->sourceLocation);
  report->get("steps") = steps;
  OUT() << "execution plan read, steps: " << steps->getSize()
    << std::endl << std::endl;

  if (options->dryRunOnly) dryRun = true;
  Natural<128> totalOperations;
  Time totalTime;
  {
    OperationsCounter totalOperationsCounter(&totalOperations);
    Timer totalTimer(&totalTime);

    for (Natural<> i(0); i < steps->getSize(); ++i) {
      runStep(i, steps->getMap(i));
      // emit report, overwrite from previous step
      Emitter emitter(options->yamlOutFile);
      emitter.emit(report);
    }
  }
  dryRun = false;

  if (dryRun)
    OUT() << "\nMemory Estimate: " <<  DryMemory::maxTotalSize / (1024.0*1024.0*1024.0) << " GB\n";

  std::stringstream totalRealtime;
  totalRealtime << totalTime;
  if (dryRun){
    OUT() << "total operations: " << totalOperations / 1e9 << " GFLOP" << std::endl;
  }
  else {
    OUT() << "total realtime: " << totalRealtime.str() << " s" << std::endl;
    OUT() << "total operations: "
      << totalOperations / 1e9 << " GFLOPS, "
      << "speed: "
      << totalOperations/1e9 / totalTime.getFractionalSeconds()
        / world->getProcesses()
      << " GFLOPS/s/core" << std::endl;
  }
  LOG() << "total realtime: " << totalRealtime.str() << " s" << std::endl;
  LOG() << "total operations: " << totalOperations / 1e9 << " GFLOPS, "
    << "speed: "
    << totalOperations/1e9 / totalTime.getFractionalSeconds()
    << " GFLOP/s" << std::endl;
  report->setValue("realtime", totalRealtime.str());
  report->setValue("floatingPointOperations", totalOperations);
  report->setValue("flops", totalOperations/totalTime.getFractionalSeconds());
}

void Cc4s::fetchSymbols(const Ptr<MapNode> &arguments) {
  for (auto key: arguments->getKeys()) {
//    auto mapNode(arguments->get(key)->toPtr<MapNode>());
//    if (mapNode) {
//      fetchSymbols(mapNode);
//      break;
//    }
    auto symbolNode(arguments->get(key)->toPtr<SymbolNode>());
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
  OUT() << "version: " << CC4S_VERSION <<
    ", date: " << CC4S_DATE << std::endl;
  OUT() << "build date: " << buildDate.str() << std::endl;
  OUT() << "compiler: " << COMPILER_VERSION << std::endl;
  OUT() << "total processes: " << world->getProcesses() << std::endl;
  OUT() << "calculation started on: " << ctime (&rawtime) << std::endl << std::endl;
  report->setValue("version", std::string(CC4S_VERSION));
  report->setValue("buildDate", buildDate.str());
  report->setValue("compiler", std::string(COMPILER_VERSION));
  report->setValue("dry run only", options->dryRunOnly);
  if (options->dryRunOnly) {
    OUT() << "DRY RUN ONLY - nothing will be calculated" << std::endl;
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
  if (!options->dryRunOnly) {
    char ownName[MPI_MAX_PROCESSOR_NAME];
    int nameLength;
    MPI_Get_processor_name(ownName, &nameLength);
    ownName[nameLength] = 0;

    std::map<std::string,std::vector<int>> ranksOfHosts;
    if (world->getRank() == 0) {
      // enter hostname/rank
      ranksOfHosts[ownName].push_back(0);
      // receive all names but own name from remote ranks
      for (
        Natural<> remoteRank(1);
        remoteRank < world->getProcesses();
        ++remoteRank
      ) {
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
          ranks->push_back(New<AtomicNode<Natural<>>>(rank, SOURCE_LOCATION));
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
bool Cc4s::dryRun = false;

int main(int argumentCount, char **arguments) {
  MPI_Init(&argumentCount, &arguments);

  Cc4s::world = New<MpiCommunicator>();
  Cc4s::options = New<Options>(argumentCount, arguments);
  if (int errcode = Cc4s::options->parse()) std::exit(errcode);

  Log::setFileName(Cc4s::options->logFile);
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

  MPI_Finalize();
  return errorHappened ? 1 : 0;
}

