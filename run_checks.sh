#!/bin/bash

cur_dir=$(pwd)

function cleanup(){
    cd $cur_dir
}
trap cleanup EXIT

cd interfaces
# cp libl1qc_dct.so.0 libl1qc_dct.so

verbose=0
data_path="example_img_data.json"

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

function parse_status() {
    if test $1 -eq 0 ; then
        parse_return="${GREEN}success${NC}"
    else
        parse_return="${RED}FAILED with status $1${NC}"
    fi
}

#
function ml_cmd() {
    if test -z $1; then
        echo "Not enough input arguments"
        exit 1
    fi
    fcn=$1
    shift # on to arguements
    arg_str=
    for arg in $@; do
        arg_str="${arg_str} ${arg},"
    done
    arg_str=$(echo $arg_str|sed -e 's/,$//')  #remove trailing comma

cmd="try; \
          $fcn($arg_str); \
          exit(0); \
     catch ME;               \
          fprintf('%s\n', ME.message); \
          exit(1);
     end"
echo $cmd
}




cmd=$(ml_cmd "test_l1qc_dct_dlopen" "'${data_path}'" "${verbose}")
matlab -nojvm -r "${cmd}"
parse_status $?
matlab_dlopen_success=$parse_return

cmd=$(ml_cmd "test_l1qc_dct_mex" "'${data_path}'" "${verbose}")
matlab -nojvm -r "${cmd}"
parse_status $?
matlab_mex_success=$parse_return


python3 dct_mkl_example.py "${data_path}" 0
parse_status $?
python_success=$parse_return



./l1qc_dct_c > /dev/null
parse_status $?
c_success=$parse_return

echo -e "${NC}-------- Checks ------------"
echo -e "Matlab dlopen: $matlab_dlopen_success"
echo -e "Matlab mex interface: $matlab_mex_success"
echo -e "Python dlopen: $python_success"
echo -e "C-code implementation: $c_success"

