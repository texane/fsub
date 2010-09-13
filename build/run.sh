#!/usr/bin/env sh

export LD_LIBRARY_PATH=$HOME/install/xkaapi_master/lib:$LD_LIBRARY_PATH
#export KAAPI_CPUSET=0,1,2,3,4,5,6,7
export KAAPI_CPUSET=0,1,2,3
#export KAAPI_CPUSET=0,1
#export KAAPI_CPUSET=0
./a.out