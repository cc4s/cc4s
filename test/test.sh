#! /usr/bin/env bash

function header()   { echo -e "\n\033[1m$@\033[0m"; }
function success()  { echo -e " \033[1;32m==>\033[0m  $@"; }
function error()    { echo -e " \033[1;31mX\033[0m  $@"; }
function arrow()    { echo -e " \033[1;34m==>\033[0m  $@"; }

# READONLY PARAMETERS
declare -r TEST_OUT_FILE=test.out
declare -r TEST_NAME=doTest.sh
declare -r MAIN_TEST_FOLDER=$(dirname $(readlink -f $0))
declare -r GLOBALS_FILE=${MAIN_TEST_FOLDER}/globals.conf
declare -r CLASS_MAGIC_WORD="@CLASS"
declare -r ALL_CLASS=all

# GLOBAL VARIABLES
TEST_RESULT=1
FAILED_TEST_COUNT=0
FAILED_TEST_LIST=()
ALL_SCRIPTS=()

# DEFAULT FLAG VALUES
RUN_COMMAND=
TEST_CLASS=essential
CONFIG=icc
CC4S_PATH=
TEST_DEBUG=

# SCRIPT PARAMETERS
declare -r __SCRIPT_VERSION="0.5"
declare -r __SCRIPT_NAME=$( basename $0 )
declare -r __DESCRIPTION="Test suite for cc4s"
declare -r __OPTIONS=":hvt:r:c:x:ld"

function check_class() {
  local testScript=$1
  for className in $(get_classes ${testScript} | tr "," "\n"); do
    if [[ ${className} = ${TEST_CLASS} ]]; then
      return 0
    fi
  done
  return 1
}

function get_classes() {
  local testScript=$1
  echo ${ALL_CLASS},$(grep "${CLASS_MAGIC_WORD}" ${testScript} | sed "s/.*${CLASS_MAGIC_WORD}=//")
}

function get_description() {
  local testScript=$1
  grep "TEST_DESCRIPTION=" ${testScript} | sed "s/.*TEST_DESCRIPTION=//" | tr -d "\""
}

function print_long_description() {
  local testScript=$1
  if grep LONG_TEST_DESCRIPTION ${testScript} 1>&2 > /dev/null; then
    echo ""
    sed -n "/cat.*LONG_TEST_DESCRIPTION/,/^LONG_TEST_DESCRIPTION/ { /LONG_TEST_DESCRIPTION/!p } " ${testScript}
    echo ""
  fi
}

function filter_scripts() {
  local filteredScripts=()
  local testScript
  for testScript in ${ALL_SCRIPTS[@]} ; do
    if check_class ${testScript}; then
      filteredScripts=(${filteredScripts[@]} ${testScript})
    fi
  done
  ALL_SCRIPTS=(${filteredScripts[@]})
}

function list_tests() {
  local testScript
  local testClasses
  local testDescription
  for testScript in ${ALL_SCRIPTS[@]} ; do
    header ${testScript#$MAIN_TEST_FOLDER/}
    testClasses=$(get_classes ${testScript})
    testDescription=$(get_description ${testScript})
    [[ -n ${testDescription} ]] && arrow "Description: ${testDescription}" || error "Description: No description available..."
    print_long_description ${testScript}
    arrow "Classes:  ${testClasses}"
  done
}

function run_testScript() {
  local testScript
  local testFolder
  local testDescription
  testScript=$1
  testFolder=$(dirname ${testScript})
  testDescription=$(get_description ${testScript})
  TEST_RESULT=1
  header "Testing ${testScript#${MAIN_TEST_FOLDER}/} ... "
  arrow "${testDescription}"
  print_long_description ${testScript}
  cd ${testFolder}
  source ${TEST_NAME} > ${TEST_OUT_FILE}
  if [[ ${TEST_RESULT} = 0 ]]; then
    success "Sucess"
  else
    error "Test FAILED"
    FAILED_TEST_LIST[${FAILED_TEST_COUNT}]=${testScript}
    let FAILED_TEST_COUNT+=1
  fi
  cd ${MAIN_TEST_FOLDER}
}


function usage_head() { echo "Usage :  $__SCRIPT_NAME [-h|-help] [-v|-version] [-c CLASS_NAME] [-r RUN_COMMAND] ... <test files>"; }
function usage ()
{
cat <<EOF
$(usage_head)

    $__DESCRIPTION

    Options:
      -h|help       Display this message
      -v|version    Display script version
      -l            List available test scripts for a given test class.
      -d            Enable debug messages
      -r            Run command
                    (e.g. "mpirun", "mpiexec.hydra -bootstrap ll" ...)
                    Default: ${RUN_COMMAND}
      -t            Test class (e.g. all essential extended)
                    Default: ${TEST_CLASS}
      -c            Compiled configuration for the cc4s binary
                    Default: ${CONFIG}
      -x            Path to cc4s executable, this overrides -c flag

    Examples:

      List all tests:
        ./${__SCRIPT_NAME} -l -t all

      List default tests (${TEST_CLASS}):
        ./${__SCRIPT_NAME} -l

      Run tests with class "silicon" using gxx configuration
        ./${__SCRIPT_NAME} -t silicon -c gxx

      Run essential tests with some custom cc4s path
        ./${__SCRIPT_NAME} -t essential -x /path/to/cc4s

      Run test files given by path using icc configuration
        ./${__SCRIPT_NAME} -c icc path/to/custom/testScript.sh path/to/other/doTest.sh

EOF
}    # ----------  end of function usage  ----------

while getopts $__OPTIONS opt
do
  case $opt in

  h|help     )  usage; exit 0   ;;

  v|version  )  echo "$__SCRIPT_NAME -- Version $__SCRIPT_VERSION"; exit 0   ;;

  l  ) LIST_TESTS=TRUE  ;;

  r  ) RUN_COMMAND=${OPTARG} ;;

  d  ) TEST_DEBUG=TRUE ;;

  c  ) CONFIG=${OPTARG} ;;

  x  ) CC4S_PATH=${OPTARG} ;;

  t  ) TEST_CLASS=${OPTARG} ;;

  * )  echo -e "\n  Option does not exist : $OPTARG\n"
      usage_head; exit 1   ;;

  esac    # --- end of case ---
done
shift $(($OPTIND-1))

# Check if some input was given after flag parsing
# if so, consider them to be the scripts to be run
if [[ -n $@ ]]; then
  ALL_SCRIPTS=($@)
else
  ALL_SCRIPTS=($(find ${MAIN_TEST_FOLDER} -name ${TEST_NAME}))
  # Filter scripts by classes
  filter_scripts
fi

# list scripts
if [[ ${LIST_TESTS} = TRUE ]]; then
  header "${#ALL_SCRIPTS[@]} tests with class ${TEST_CLASS}"
  list_tests
  exit 0
fi

# Check if cc4s path was overriden
if [[ -z ${CC4S_PATH} ]] ; then
  CC4S_PATH=$(readlink -f "${MAIN_TEST_FOLDER}/../build/${CONFIG}/bin/Cc4s")
fi

#Title
header "${__DESCRIPTION} version ${__SCRIPT_VERSION}"
# Print out all relevant configuration
arrow "CONFIG      = ${CONFIG}"
arrow "CC4S_PATH   = ${CC4S_PATH}"
arrow "TEST_CLASS  = ${TEST_CLASS}"
arrow "RUN_COMMAND = ${RUN_COMMAND}"


# Check for cc4s executable
if [[ ! -x ${CC4S_PATH} ]] ; then
  error "CC4S not found"
  exit 1
fi

arrow "Sourcing ${GLOBALS_FILE} file"
source ${GLOBALS_FILE}

# run the tests
for TEST_SCRIPT in ${ALL_SCRIPTS[@]}; do
  run_testScript ${TEST_SCRIPT}
done

header "${#ALL_SCRIPTS[@]} tests DONE for class '${TEST_CLASS}'"
header "${FAILED_TEST_COUNT} tests FAILED for class '${TEST_CLASS}'"

#Print out the tests that failed, if any
if [[ ! ${#FAILED_TEST_LIST[@]} = 0 ]]; then
  error "Tests failed:"
  for TEST_SCRIPT in ${FAILED_TEST_LIST[@]} ; do
    echo -e "\t${TEST_SCRIPT}"
  done
else
  success "All tests passed!"
fi

if [[ ${FAILED_TEST_COUNT} != 0 && ${TEST_CLASS} == "essential" ]]; then
  cat <<EOF

############################################
#                                          #
#  DO NOT COMMIT TO MASTER BRANCH, PLEASE  #
#                                          #
############################################

EOF
fi


