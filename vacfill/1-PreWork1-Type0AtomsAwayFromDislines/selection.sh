#!/bin/bash

echo 'ATOMS SELECTION:'
# go to
cd AtomsAwayOrNearFromDislines_f/
# inputs
cp ../datafiles-input/disline.input ./disline.input
cp ../datafiles-input/input.data ./input.data
# RUN
sh word_count.sh
gfortran AtomsAwayOrNearFromDislines.f90
./a.out
# move OUTPUTs and remove INPUTs
mv output_*.data ../datafiles-output/
rm disline.input input.data
# go out 
cd ../

echo 'Selecting atoms away from dislines, DONE!'

