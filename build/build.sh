#!/usr/bin/env sh
gcc -std=gnu99 -Wall -O3 -march=native \
-I. -I$HOME/install/xkaapi_master/include \
../src/main.c \
../src/fsub_gsl.c \
../src/fsub_seq.c \
../src/fsub_pthread.c \
../src/fsub_xkaapi.c \
-L$HOME/install/xkaapi_master/lib \
-lpthread -lgsl -lblas -lxkaapi

