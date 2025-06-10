#!/bin/bash



echo '----------------------------------------------------------------------------------------'
echo '------STEP-1:Type0 atoms AWAY from disline, and note that this type0 is from CSP--------'
echo '----------------------------------------------------------------------------------------'
# go to
cd 1-Type0AtomsAwayFromDislines/
# go to input file
cd datafiles-input/
# inputs
cp ../../datafile-input/type0.data ./input.data
cp ../../datafile-input/disline.input ./disline.input
# go out of input-file
cd ../
sh selection.sh # inputs and outputs are processed in this shell
# go out
cd ../
echo '----------------------------------------------------------------------------------------'
echo '------------------------------STEP-2:Subtracting----------------------------------------'
echo '----------------------------------------------------------------------------------------'
# go to
cd 2-SubtractType0FromDump0/
# inputs: a.input - b.input = c.output
cp ../datafile-input/dump.0 ./a.input
mv ../1-Type0AtomsAwayFromDislines/datafiles-output/output_away.data ./b.input
gfortran SubtractType0FromDump0.f90
./a.out
# save outputs, rm inputs
mv c.d* ../datafile-output/
rm a.input b.input
# go out
cd ../

echo '==>> SUBTRACTION IS DONE!'

