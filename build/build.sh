#!/usr/bin/env sh
gcc -Wall -O3 -march=native -I. ../src/main.c -lpthread -lgsl -lblas
