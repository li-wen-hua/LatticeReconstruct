#!/bin/bash


echo '----------------------------------------------------------------------------------------'
echo '---------------------------STEP-1:Wrap type-0 atoms a Shell-----------------------------'
echo '----------------------------------------------------------------------------------------'
# go to
cd 1-WrapType0AShell/
# input: input.data
cp ../datafile-input/type0.data ./input.data
gfortran WrapAShell.f90
./a.out
# remove input and dump-format output
rm input.data # redundant input
rm output.dump # redundant output
# go out
cd ../
echo '----------------------------------------------------------------------------------------'
echo '------------------------STEP-2:Selecting atoms NEAR type-0 atoms------------------------'
echo '----------------------------------------------------------------------------------------'
# go to
cd 2-AtomsAroundType0/
# move input, so we don't need to delete
mv ../1-WrapType0AShell/output.data ./wrapped.type0.data
cp ../datafile-input/dump.0 ./dump.0
gfortran -fopenmp SelectingAtomsAroundType0.f90
./a.out
# remove inputs, and away-atoms, and near-atom dump file
rm wrapped.type0.data # rm input
rm dump.0 # rm input
rm away.d* near.dump # rm redundant output
# go out
cd ../
echo '----------------------------------------------------------------------------------------'
echo '-------------------STEP-3:Selecting NEAR-TYPR0-atoms AWAY from dislines-----------------'
echo '----------------------------------------------------------------------------------------'
# go to 
cd 3-NearType0AtomsAwayFromDislines/
cd ./datafiles-input/ 
# inputs, in inputfile
mv ../../2-AtomsAroundType0/near.data ./input.data
cp ../../datafile-input/disline.input ./disline.input
# go out inputfile
cd ../
sh selection.sh
cd ../ # go out
echo '----------------------------------------------------------------------------------------'
echo '---------------------------STEP-4:Wrap Dump0-atoms a Shell------------------------------'
echo '----------------------------------------------------------------------------------------'
# go to
cd 4-WrapDump0AShell/
# input
cp ../datafile-input/dump.0 ./input.dump
gfortran WrapAShell.f90
./a.out
# input, dump-format output
rm input.dump output.dump # remove input and dump-format output
# go out
cd ../
echo '----------------------------------------------------------------------------------------'
echo '---------STEP-5:CSP, working only on type0 atoms which are AWAY from dislines, so-------'
echo '---------here the presupposition is, less type0 atoms in the latter configuration-------'
echo '----------------------------------------------------------------------------------------'
# go to
cd 5-SimplifiedCSP/
# inputs, move
mv ../4-WrapDump0AShell/output.data ./wrapped.data.output
mv ../3-NearType0AtomsAwayFromDislines/datafiles-output/output_away.data ./near.data
gfortran -fopenmp SimplifiedCSP.f90
./a.out
# save outputs, remove inputs
mv good.d* bad.d* ../datafile-output/ # move outputs
rm near.data wrapped.data.output # rm inputs
# go out
cd ../

echo '==>> CSP IS DONE!'









