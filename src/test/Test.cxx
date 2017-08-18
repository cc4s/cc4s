#define CATCH_CONFIG_RUNNER
#include <test/Test.hpp>
#include <Cc4s.hpp>
#include <Parser.hpp>
#include <util/MpiCommunicator.hpp>
#include <util/Log.hpp>

using namespace cc4s;
using namespace CTF;

World *Cc4s::world;
Options *Cc4s::options;

void printBanner() {
  OUT() << "       UNIT TESTING__      " << std::endl
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


int main(int argumentCount, char **arguments) {
  MPI_Init(&argumentCount, &arguments);

  Cc4s::world = new World(argumentCount, arguments);
  Cc4s::options = new Options(argumentCount, arguments);
  Log::setRank(Cc4s::world->rank);
  LogStream logStream(
    Cc4s::options->logFile, Cc4s::options->logLevel
  );
  Log::setLogStream(&logStream);

  printBanner();

  const int result = Catch::Session().run(argumentCount, arguments);

  MPI_Finalize();

  return result;
}
