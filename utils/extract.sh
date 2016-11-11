#! /usr/bin/env bash

__SCRIPT_VERSION="0.1"
__SCRIPT_NAME=$( basename $0 )
__DESCRIPTION="Extract information from files in an ordered way"
__OPTIONS=":hvRp:m:E:b:d:"

header()   { echo -e "\n\033[1m$@\033[0m"; }
success()  { echo -e " \033[1;32m==>\033[0m  $@"; }
error()    { echo -e " \033[1;31mX\033[0m  $@"; }
arrow()    { echo -e " \033[1;34m==>\033[0m  $@"; }
warning()  { echo -e " \033[0;93m==>\033[0m  $@"; }

# For this script to work you only need awk, find and getopts
declare -r AWK=awk
declare -r FIND=find

# parameters
mainMagicWord=wiki
mainMagicRegex='[ ]+([a-zA-Z0-9._-]+):([0-9]+)'
magicWordStart=">${mainMagicWord}"
magicWordEnd="<${mainMagicWord}"
magicRegexStart="${magicWordStart}${mainMagicRegex}"
magicRegexEnd="${magicWordEnd}${mainMagicRegex}"
buildDir=$(dirname $0)/${mainMagicWord}_build
distDir=${mainMagicWord}
searchPath=$(pwd)
recursiveSearch=

getFiles() {
  local depth
  local singleFiles=("$@")
  if [[ -z ${recursiveSearch} ]]; then
    depth="-maxdepth 1"
  fi
  ${FIND} ${searchPath} ${depth} -type f
  for f in ${singleFiles[@]} ; do
    echo "${f}"
  done
}

usage_head() { echo "Usage :  $__SCRIPT_NAME [-h|-help] [-v|-version]"; }
usage ()
{
cat <<EOF
$(usage_head)

    $__DESCRIPTION

    Options:
      -h|help       Display this message
      -v|version    Display script version
      -p            Path to search for files
      -R            Recursively search the tree
      -m            Main magic word (Default: ${mainMagicWord})
      -E            Regex for parsing (Default: ${mainMagicRegex})
      -b            Build directory (Default: ${buildDir})
      -d            Dist directory (Default: ${distDir})


    Example:

      - Write a simple html page

      (file main.c)
      =============

        /* ${magicWordStart} index.html:0

        <h1> Main program </h1>

        ${magicWordEnd} index.html:0 */

         int main(int argc, char *argv[])
         {
           printf("Hello World\n");
           return 0;
         }


        /* ${magicWordStart} index.html:1

          <b> Here is a <it> MathJax </it> formula

            \$\$ \int_{\mathbb{R}} f(x) = 0 \$\$

        ${magicWordEnd} index.html:1 */

      - Then you parse the contents by doing

        ${__SCRIPT_NAME} -d wiki main.c

        This will create a wiki directory
        with the file index.html:

      (file index.html)
      =================

        <h1> Main program </h1>

        <b> Here is a <it> MathJax </it> formula

          \$\$ \int_{\mathbb{R}} f(x) = 0 \$\$


    This program is maintained by Alejandro Gallo.
EOF
}    # ----------  end of function usage  ----------

while getopts $__OPTIONS opt; do
  case $opt in

  h|help     )  usage; exit 0   ;;

  v|version  )  echo "$__SCRIPT_NAME -- Version $__SCRIPT_VERSION"; exit 0   ;;

  R  )  recursiveSearch=1;;

  p  )  searchPath=${OPTARG};;

  m  )  magicWord=${OPTARG};;

  E  )  magicRegex=${OPTARG};;

  b  )  buildDir=${OPTARG};;

  d  )  distDir=${OPTARG};;

  * )  echo -e "\n  Option does not exist : $OPTARG\n"
      usage_head; exit 1   ;;

  esac    # --- end of case ---
done
shift $(($OPTIND-1))

header "Parsing files"

[[ -d ${buildDir} ]] && rm -r ${buildDir}
mkdir -p ${buildDir}

for fileName in $(getFiles $@ | sort | uniq); do
  arrow $fileName;
  ${AWK} '
    /'"${magicRegexStart}"'/,/'"${magicRegexEnd}"'/ {
      # Save matching elements in array a 
      matched = match($0, "'"${magicRegexStart}"'", a);
      if (matched){
        newFileName = a[2]"-"a[1];
        printf "\tInformation found for <%s>\n",newFileName;
      }
      if (\
          $0 !~ "'"${magicRegexStart}"'"\
          &&\
          $0 !~ "'"${magicRegexEnd}"'"\
        ){
        print $0 > "'"${buildDir}"'/"newFileName
      }
    }
    ' ${fileName};
done

header "Building..."

[[ -d ${distDir} ]] && rm -r ${distDir}
mkdir -p ${distDir}

for fileName in $(ls -v ${buildDir}/*); do
  generalFileName=$(basename ${fileName} | sed "s/[0-9]\+-//")
  cat ${fileName} >> ${distDir}/${generalFileName}
done




#vim-run: bash % -h
#vim-run: clear; bash %
#vim-run: bash % test.wiki
#vim-run: bash % -R -p src/
#vim-run: bash %
