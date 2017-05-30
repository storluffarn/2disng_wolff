#!/bin/bash

rm *.pbm
rm *.png

g++ -std=c++17 2dising.cpp -o 2dising -O3

./2dising

./postprocess.sh


