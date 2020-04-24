#!/bin/bash

g++ -std=c++11 -Wall $(gsl-config --cflags) $(globes-config --cflags) $(gsl-config --libs) $(globes-config --libs) -o kamlandstd kamlandstd.cpp -lglobes -lgsl -lgslcblas -lm -fopenmp
./kamlandstd
rm kamlandstd
mv chisq.dat  datfiles/chisq.dat 
mv bfp.dat  datfiles/bfp.dat 

cd datfiles
sh plot.sh
cd ../
notify-send "Acabei"





