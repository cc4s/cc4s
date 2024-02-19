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
#include <GC.hpp>
#include <algorithms/Algorithm.hpp>
#include <tcc/Tcc.hpp>
#include <Parser.hpp>
#include <Emitter.hpp>
#include <Timer.hpp>
#include <MpiCommunicator.hpp>
#include <Log.hpp>
#include <Exception.hpp>

#include <TensorSet.hpp>

#include <fstream>
#include <string>
#include <sstream>

using namespace cc4s;

void Cc4s::run() {
  printBanner();
  LOG() << "Start dry Run" << std::endl;
  runSteps(true);
  LOG() << "Finished dry Run" << std::endl;
  if (options->dryRanks > 0) return;
  //clear log file if demanded
  if (options->persistLog < 2) {
    Log::closeLogFile();
    Log::setRank(Cc4s::world->getRank());
  }
  LOG() << "Start calculation" << std::endl;
  runSteps(false);
  LOG() << "Finished calculation" << std::endl;
  printEndBanner();
}

static std::map<std::string, std::string> check_tensor_type(Ptr<Node> &value) {
#define _RETURN_IF_IS_(_T, _F, _TE)                                            \
  do {                                                                         \
    Ptr<PointerNode<_T<_F, _TE>>> ptr =                                        \
        value->toPtr<PointerNode<_T<_F, _TE>>>();                              \
    if (ptr) {                                                                 \
      auto t = std::dynamic_pointer_cast<_T<_F, _TE>>(ptr->value);             \
      std::stringstream ptr_name;                                              \
      ptr_name << t.get();                                                     \
      return {                                                                 \
          {"field", #_F},                                                      \
          {"tensor-engine", #_TE},                                             \
          {"type", #_T},                                                       \
          {"pointer", ptr_name.str()},                                         \
      };                                                                       \
    }                                                                          \
  } while (0)
  _RETURN_IF_IS_(TensorSet, Real<>, DefaultTensorEngine);
  _RETURN_IF_IS_(TensorSet, Complex<>, DefaultTensorEngine);
  _RETURN_IF_IS_(Tensor, Real<>, DefaultTensorEngine);
  _RETURN_IF_IS_(Tensor, Complex<>, DefaultTensorEngine);
  return {{"type", "none"}};
#undef _RETURN_IF_IS_
}

static Ptr<MapNode> make_step_for_output(Ptr<MapNode> const &step) {
  Ptr<MapNode> result = New<MapNode>(SOURCE_LOCATION);
  result->get("name") = step->get("name");
  for (const auto &in_or_out_key : {"in", "out"}) {
    auto result_in_or_out = New<MapNode>(SOURCE_LOCATION);
    result->get(in_or_out_key) = result_in_or_out;
    Ptr<MapNode> in_or_out = step->getMap(in_or_out_key);
    for (const auto &key : in_or_out->getKeys()) {
      Ptr<Node> value = in_or_out->get(key);
      // If it is not an atomic node, it means it might be energy and so on,
      // or a mixer node information, so let us leave it like it is
      if (!value->isAtomic()) {
        result_in_or_out->get(key) = value;
        continue;
      } else {
        // We first check if we have a tensor node, which will
        // be done using the function check_tensor_type,
        // if we do have this, then we will just have some information
        // and will write it in the output
        auto tensor_info = check_tensor_type(value);
        if (tensor_info["type"] != "none") {
          auto _result = New<MapNode>(SOURCE_LOCATION);
          for (auto const &__pair : tensor_info) {
            auto __key = __pair.first, __val = __pair.second;
            _result->get(__key) =
                New<AtomicNode<std::string>>(__val, SOURCE_LOCATION);
          }
          result_in_or_out->get(key) = _result;
        } else {
          // if it is not a tensor, just leave it as it was before, this
          // will be strings, integers, real numbers and so on
          result_in_or_out->get(key) = value;
        }
      }
    }
  }
  return result;
}

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

void Cc4s::runSteps(const bool dry) {
  auto output(New<MapNode>(SOURCE_LOCATION));
  output->get("executionEnvironment") = executionEnvironment;

  Cc4s::dryRun = dry;
  // parse input
  Parser parser(options->inFile);
  auto input(parser.parse());
  auto steps(input->toPtr<MapNode>());
  ASSERT_LOCATION(steps, "expecting map as input", input->sourceLocation);

  // start with empty storage
  storage = New<MapNode>(SOURCE_LOCATION);

  auto executedSteps(New<MapNode>(SOURCE_LOCATION));
  output->get("steps") = executedSteps;

  Natural<128> totalOperations;
  Time totalTime;
  {
    OperationsCounter totalOperationsCounter(&totalOperations);
    Timer totalTimer(&totalTime);

    for (Natural<> i(0); i < steps->getSize(); ++i) {
      auto step(steps->getMap(i));
      runStep(i, step);
      if (Cc4s::dryRun) executedSteps->get(i) = step;
      else              executedSteps->get(i) = make_step_for_output(step);
      // executedSteps->get(i) = step;
      if (!options->no_gc)
        collect_garbage(i, steps, storage);
      // emit output, overwrite from previous step
      if ( options->outfiles == 2 || ( options->outfiles == 1 && !dry )) {
        std::string stagePrefix(dry ? "dry-" : "");
        Emitter emitter(stagePrefix + options->yamlOutFile);
        emitter.emit(output);
      }
    }
  }
  Cc4s::dryRun = false;

  auto statistics(New<MapNode>(SOURCE_LOCATION));
  std::stringstream totalRealtime;
  totalRealtime << totalTime;
  LOG() << "total realtime: " << totalRealtime.str() << " s" << std::endl;
  LOG() << "total operations: " << totalOperations / 1e9 << " GFLOPS, "
    << "speed: "
    << totalOperations/1e9 / totalTime.getFractionalSeconds() / getProcessesCount()
    << " GFLOP/core/s" << std::endl;
  if (dry) {
    auto GB(1024.0*1024.0*1024.0);
    auto assumedGflops(10);
    OUT() << "Dry run finished. Estimates provided for "
      << getProcessesCount() << " ranks.\n";
    OUT() << "Memory estimate (per Rank/Total): ";
    OUT() <<  DryMemory::maxTotalSize / GB / getProcessesCount() << " / "
          <<  DryMemory::maxTotalSize / GB << " GB\n";
    OUT() << "Operations estimate (per Rank/Total): ";
    OUT() << totalOperations / 1e9 / getProcessesCount() << " / "
          << totalOperations / 1e9 << " GFLOPS" << std::endl;
    OUT() << "Time estimate with assumed performance of "
      << assumedGflops << " GFLOPS/core/s: ";
    OUT() << totalOperations / 1e9 / getProcessesCount() / assumedGflops
          << " s "
          << "(" << totalOperations / 1e9 / getProcessesCount() / assumedGflops / 3600 << " h)\n";
    OUT() << "--" << std::endl;
    LOG() << "memory estimate: " << DryMemory::maxTotalSize / GB << " GB"
      << std::endl;
  }
  statistics->setValue("realtime", totalRealtime.str());
  statistics->setValue("floatingPointOperations", totalOperations);
  statistics->setValue("flops", totalOperations/totalTime.getFractionalSeconds());
  output->get("statistics") = statistics;

  // emit final stage output
  if ( options->outfiles == 2 || ( options->outfiles == 1 && !dry )) {
    std::string stagePrefix( dry ? "dry-" : "" );
    Emitter emitter(stagePrefix + options->yamlOutFile);
    emitter.emit(output);
  }
}

void Cc4s::runStep(Natural<> i, const Ptr<MapNode> &step) {
  auto algorithmName(step->getSymbol("name"));
  OUT() << "step: " << (i+1) << ", " << algorithmName << std::endl;
  LOG() << "begin step: " << (i+1)
        << ", name: " << algorithmName << std::endl;

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

  auto statistics(New<MapNode>(SOURCE_LOCATION));
  std::stringstream realtime;
  realtime << time;
  OUT() << "realtime " << realtime.str() << " s" << std::endl;
  OUT() << "--" << std::endl;
  LOG() << "end step: " << (i+1)
        << ", name: " << algorithmName
        << ", realtime: " << realtime.str() << " s"
        << ", operations: " << operations / 1e9 << " GFLOP"
        << ", speed: "
        << operations / 1e9 / time.getFractionalSeconds() / getProcessesCount()
        << " GF/rank/s" << std::endl;
  statistics->setValue("realtime", realtime.str());
  statistics->setValue("floatingPointOperations", operations);
  statistics->setValue("flops", operations / time.getFractionalSeconds());
  step->get("statistics") = statistics;
  // resources held by the algorithm are released when it goes out of scope
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


void Cc4s::printBanner() {
  OUT() << std::endl
        << "                __ __      " << std::endl
        << "     __________/ // / _____" << std::endl
        << "    / ___/ ___/ // /_/ ___/" << std::endl
        << "   / /__/ /__/__  __(__  ) " << std::endl
        << "   \\___/\\___/  /_/ /____/  " << std::endl
        << "  Coupled Cluster for Solids" << std::endl << std::endl;

  executionEnvironment = New<MapNode>(SOURCE_LOCATION);
  std::stringstream buildDate;
  buildDate << __DATE__ << " " << __TIME__;
  time_t rawtime;
  time (&rawtime);
  OUT() << "version: " << CC4S_VERSION <<
    ", date: " << CC4S_DATE << std::endl;
  OUT() << "build date: " << buildDate.str() << std::endl;
  OUT() << "compiler: " << COMPILER_VERSION << std::endl;
  OUT() << "mpi: " << CC4S_MPI_VERSION << std::endl;
  OUT() << "total processes: " << world->getProcesses() << std::endl;
  OUT() << "calculation started on: " << ctime (&rawtime) << std::endl;
  executionEnvironment->setValue("version", std::string(CC4S_VERSION));
  executionEnvironment->setValue("buildDate", buildDate.str());
  executionEnvironment->setValue("compiler", std::string(COMPILER_VERSION));
  executionEnvironment->setValue("totalProcesses", world->getProcesses());
  executionEnvironment->setValue("startTime", std::string(ctime (&rawtime)));
  executionEnvironment->setValue("dryRanks", options->dryRanks);
  if (options->dryRanks > 0) {
    OUT() << "DRY RUN ONLY - nothing will be calculated" << std::endl;
  }
  executionEnvironment->get("hosts") = getHostList();
}

void Cc4s::printEndBanner() {
  time_t rawtime;
  time (&rawtime);
  OUT() << "Calculation finished successfully.\n";
  OUT() << "Calculation finished at: " << ctime (&rawtime) << std::endl;
  //check if we have to remove the log files
  if (options->persistLog == 0) {
    Log::closeLogFile();
    remove(Log::getFileName().c_str());
  }
}

Natural<> Cc4s::getProcessesCount() {
  return options->dryRanks > 0 ? options->dryRanks : Cc4s::world->getProcesses();
}


Ptr<MapNode> Cc4s::getHostList() {
  auto hosts(New<MapNode>(SOURCE_LOCATION));
  if (options->dryRanks == 0) {
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
    } else {
      // send own name
      MPI_Send(
        ownName, MPI_MAX_PROCESSOR_NAME, MPI_BYTE,
        0, 0, MPI_COMM_WORLD
      );
    }
  } else {
    // TODO: list planned execution environment from options
  }
  return hosts;
}

bool isDebugged() {
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
  bool isSuccessful(true);

  Cc4s cc4s;
  if (isDebugged()) {
    // run without try-catch in debugger to allow tracing throwing code
    cc4s.run();
  } else {
    // without debugger: catch and write list of causes
    try {
      cc4s.run();
    } catch (Ptr<Exception> cause) {
      isSuccessful = false;
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
      isSuccessful = false;
      OUT() << "unhandled exception encountered (std::exception):" << std::endl;
      OUT() << cause.what() << std::endl;
    } catch (const char *message) {
      isSuccessful = false;
      OUT() << "unhandled exception encountered (const char *):" << std::endl;
      OUT() << message << std::endl;
    } catch (...) {
      isSuccessful = false;
      OUT() << "unhandled exception encountered (...)." << std::endl;
    }
  }

  MPI_Finalize();
  return isSuccessful ? 0 : 1;
}

