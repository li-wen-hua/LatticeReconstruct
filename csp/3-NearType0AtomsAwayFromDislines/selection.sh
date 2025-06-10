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
mv output_away.data ../datafiles-output/
rm disline.input input.data output_near.data
# go out 
cd ../
# go to
cd ./datafiles-input 
rm disline.input input.data # rm inputs
cd ../

echo 'Selecting atoms away from dislines, DONE!'

