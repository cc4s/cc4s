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

#ifndef OPTIONS_DEFINED
#define OPTIONS_DEFINED

#include <string>
#include <extern/CLI11.hpp>

namespace cc4s {

  struct Options {

    std::string inFile, logFile, yamlOutFile;
    int dryRanks;
    int persistLog;
    int outfiles;
    CLI::App app;
    int argc;
    char** argv;

    Options(int _argc, char **_argv)
      : inFile("cc4s.in")
      , logFile("cc4s.log")
      , yamlOutFile("cc4s.out.yaml")
      , dryRanks(0)
      , persistLog(0)
      , outfiles(1)
      , app{"CC4S: Coupled Cluster For Solids"}
      , argc(_argc)
      , argv(_argv)
    {
      app.add_option("-i,--in", inFile, "Input file path")
         ->default_val(inFile);
      app.add_option("-o,--out", yamlOutFile, "Output yaml file")
         ->default_val(yamlOutFile);
      app.add_option("-l,--log", logFile, "Output log file")
         ->default_val(logFile);
      app.add_option("-d,--dry-ranks",
                     dryRanks,
                    "Number of processes for dry run.\n"
                    "If zero, run normally after dry run.\n"
                    "If non-zero, specifies number ranks for resource\n"
                    "estimation and do not run after dry run")
         ->default_val(dryRanks);
      app.add_option("-p, --persist-log",
                      persistLog,
                     "Persistence of the log-file:\n"
                     " 0: Logfile is removed after sucessful run.\n"
                     " 1: Logfile's dryRun-part is deleted after successful dryRun.\n"
                     " 2: Logfile is left untouched.\n"
                    )->default_val(persistLog);
      app.add_option("-O, --outfiles"
                    , outfiles
                    , "Specify the different output files:\n"
                      " 0: No *out.yaml is written.\n"
                      " 1: Write output.out.yaml but suppress dry-output.out.yaml.\n"
                      " 2: Write ouptut.out.yaml and dry-output.out.yaml.\n"
                    )->default_val(outfiles);
    }

    int parse() {
      CLI11_PARSE(app, argc, argv);
      return 0;
    }

  };


}

#endif

