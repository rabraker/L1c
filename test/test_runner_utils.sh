#!/bin/bash



function get_os(){
    case $(uname) in
        MINGW*) host_os=windows;;
        *) host_os=unix;;
    esac

    echo $host_os
}


function ml_cmd() {
    # expect: $1=path where mex src is, which contains, e.g., l1qc_dct_opts.m
    #         $2=path we need to add (where mex and .so/.dll are located)
    #         $3=path to test src dir, containing matlab scripts.
    #         $4=function name
    #         $5...$n=function arguments
    #
    if test "$#" -le 5; then
        echo "Not enough input arguments"
        exit 1
    fi
    mex_src_dir_=$1
    shift
    lib_dir_=$1
    shift
    matlab_script_dir_=$1
    shift
    fcn=$1
    shift # on to arguments

    arg_str=
    for arg in $@; do
        arg_str="${arg_str} ${arg},"
    done
    arg_str=$(echo $arg_str|sed -e 's/,$//')  #remove trailing comma
    # It is important that the interfaces get added first, and lib_dir last,
    # so that matlab will find the mex function when it searches lib_dir first,
    # not the help file in interfaces (only relevant for VPATH builds).
    cmd="try; \
          addpath('${mex_src_dir_}');
          addpath('${matlab_script_dir_}');
          addpath('${lib_dir_}');
          $fcn($arg_str); \
          exit(0);\
     catch ME;\
          fprintf('%s\n', ME.message); \
          exit(1);
     end"
    echo $cmd
}
