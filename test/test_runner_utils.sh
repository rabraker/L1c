#!/bin/bash



function get_os(){
    case $(uname) in
        MINGW*) host_os=windows;;
        *) host_os=unix;;
    esac

    echo $host_os
}


function ml_cmd() {
    # expect: $1=path we need to add (where mex and .so/.dll are located)
    #         $2=path to test src dir, containing matlab scripts.
    #         $2=function name
    #         $3...$n=function arguments
    #
    if test "$#" -le 3; then
        echo "Not enough input arguments"
        exit 1
    fi
    lib_dir=$1
    shift
    matlab_script_dir=$1
    shift
    fcn=$1
    shift # on to arguments

    arg_str=
    for arg in $@; do
        arg_str="${arg_str} ${arg},"
    done
    arg_str=$(echo $arg_str|sed -e 's/,$//')  #remove trailing comma

    cmd="try; \
          addpath('${lib_dir}');
          addpath('${matlab_script_dir}');
          addpath('${ABS_TOP_SRCDIR}/interfaces/');
          $fcn($arg_str); \
          exit(0);\
     catch ME;\
          fprintf('%s\n', ME.message); \
          exit(1);
     end"
    echo $cmd
}
