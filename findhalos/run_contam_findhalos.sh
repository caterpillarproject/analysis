#!/bin/bash

scriptpath=~/analysis/findhalos/caterpillar_findhalos.py
halobase=/bigbang/data/AnnaGroup/caterpillar/halos
flags="-w" #"-f

python $scriptpath $flags -c 1 > $halobase/low_mass_halos/findhalos.o
python $scriptpath $flags -c 2 > $halobase/middle_mass_halos/findhalos.o
python $scriptpath $flags -c 3 > $halobase/high_mass_halos/findhalos.o
