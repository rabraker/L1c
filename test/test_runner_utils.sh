#!/bin/bash



function get_os(){
    case $(uname) in
        MINGW*) host_os=windows;;
        *) host_os=unix;;
    esac

    echo $host_os
}


function ml_cmd() {
    #         $1=function name
    #         $2...$n=function arguments
    #
    if test "$#" -le 2; then
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
    # It is important that the interfaces get added first, and lib_dir last,
    # so that matlab will find the mex function when it searches lib_dir first,
    # not the help file in interfaces (only relevant for VPATH builds).
    cmd="try; \
          init_matlab_paths; \
          $fcn($arg_str); \
          exit(0);\
     catch ME;\
          fprintf('%s\n', ME.message); \
          exit(1);
     end"
    echo $cmd
}
