#!/bin/bash

scriptpath=/spacebase/data/alexji/analysis/findhalos/caterpillar_findhalos.py
halobase=/bigbang/data/AnnaGroup/caterpillar/halos

#python $scriptpath -w -c 1 > $halobase/low_mass_halos/findhalos.o
#python $scriptpath -w -c 2 > $halobase/middle_mass_halos/findhalos.o
python $scriptpath -w -c 3 > $halobase/high_mass_halos/findhalos.o
