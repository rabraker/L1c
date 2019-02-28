#!/bin/bash

cur_dir=$(pwd)

function cleanup(){
    cd $cur_dir
}
trap cleanup EXIT

if test -z $verbose; then
    verbose=0
fi



LIB_DIR="${ABS_TOP_BUILDDIR}/interfaces"
TEST_SRCDIR="${ABS_TOP_SRCDIR}/test"
ml_script_dir="${TEST_SRCDIR}"
data_path="${ABS_TOP_SRCDIR}/test/test_data/example_img_data.json"


RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color
err_count=0
function parse_status() {
    err_count+=$1
    if test $1 -eq 0 ; then
        parse_return="${GREEN}success${NC}"
    else
        parse_return="${RED}FAILED with status $1${NC}"
    fi
}

#
function ml_cmd() {
    # expect: $1=path we need to add (where mex and .so/.dll are located)
    #         $2=path to test src dir, containing matlab scripts.
    #         $2=function name
    #         $3...$n=function arguments
    #
    if test -z $1; then
        echo "Not enough input arguments"
        exit 1
    fi

    fcn=$1
    shift # on to arguments

    arg_str=
    for arg in $@; do
        arg_str="${arg_str} ${arg},"
    done
    arg_str=$(echo $arg_str|sed -e 's/,$//')  #remove trailing comma

cmd="try; \
          addpath('${LIB_DIR}');
          addpath('${ml_script_dir}');
          addpath('${ABS_TOP_SRCDIR}/interfaces/');
          $fcn($arg_str); \
          exit(0); \
     catch ME;               \
          fprintf('%s\n', ME.message); \
          exit(1);
     end"
echo $cmd
}


# Check the matlab dlopen interface
cmd=$(ml_cmd  "test_l1qc_dct_dlopen" "'${data_path}'" "${verbose}")
matlab -nojvm -r "${cmd}"
parse_status $?
matlab_dlopen_success=$parse_return

# Check the mex interface
cmd=$(ml_cmd "test_l1qc_dct_mex" "'${data_path}'" "${verbose}")
matlab -nojvm -r "${cmd}"
parse_status $?
matlab_mex_success=$parse_return

# check the python dlopen interface
LIB_DIR="${LIB_DIR}" python3 "${ABS_TOP_SRCDIR}/interfaces/dct_example.py" "${data_path}" 0
parse_status $?
python_success=$parse_return



"${LIB_DIR}/l1qc_dct_c" "${data_path}"> /dev/null
parse_status $?
c_success=$parse_return

echo -e "${NC}-------- Checks ------------"
echo -e "Matlab dlopen: $matlab_dlopen_success"
echo -e "Matlab mex interface: $matlab_mex_success"
echo -e "Python dlopen: $python_success"
echo -e "C-code implementation: $c_success"

exit $err_count
