#!/bin/bash

path=$(pwd)

mkdir -p Hubs
mkdir -p Superimposed
mkdir -p Edgelists
mkdir -p HCTs
mkdir -p PLaw
mkdir -p LGraphs

for i in $(ls)
do cd $i; 
cp $i"_hubs.txt" $path/Hubs;
cp $i"_si.svg" $path/Superimposed;
cp $i"_edgelist.txt" $path/Edgelists;
cp $i"_weighted_edgelist.csv" $path/HCTs;
cp $i".svg" $path/PLaw;
cp $i"_lg.svg" $path/LGraphs;
cd $path; done
