#!/bin/bash

echo '----------------------------------------------------------------------------------------'
echo '--------------------STEP-1:Finding type-0 atoms away from dislines----------------------'
echo '----------------------------------------------------------------------------------------'
# go to
cd ./1-PreWork1-Type0AtomsAwayFromDislines/
# go to
cd ./datafiles-input/
# inputs
cp ../../datafiles-input/type0.data ./input.data
cp ../../datafiles-input/disline.input ./disline.input
# go out
cd ../
# run
sh selection.sh
# outputs and remove
mv datafiles-output/output_away.data ../datafiles-output/
rm datafiles-output/output_near.data
rm datafiles-input/*
# go out
cd ../

echo '----------------------------------------------------------------------------------------'
echo '-------------------------------STEP-2:Wrap dump.0 a shell-------------------------------'
echo '----------------------------------------------------------------------------------------'
# go to
cd ./2-PreWork2-WrapDump0AShell/
# inputs
cp ../datafiles-input/dump.0 ./input.dump
# run
gfortran WrapAShell.f90
./a.out
# outputs and remove
mv output.* ../datafiles-output/
rm input.dump
# go out
cd ../

echo '----------------------------------------------------------------------------------------'
echo '---------------------------------STEP-3:Mirror atoms------------------------------------'
echo '----------------------------------------------------------------------------------------'
# go to
cd 3-FindingMirrorAtoms/
# inputs
cp ../datafiles-output/output_away.data ./input_atomsself.data
cp ../datafiles-output/output.dump ./input_mirror.dump
# run
gfortran -fopenmp FindingMirrorAtoms.f90
./a.out
# outputs and remove
mv output_mirroratoms.d* ../datafiles-output/
rm input_*
# go out
cd ../

echo '----------------------------------------------------------------------------------------'
echo '----------------------------STEP-4:Wrap mirror atoms a shell----------------------------'
echo '----------------------------------------------------------------------------------------'
# go to
cd 4-WrapMirrorAtomsAShell/
# inputs
cp ../datafiles-output/output_mirroratoms.dump ./input.dump
# run
gfortran WrapAShell.f90
./a.out
# outputs and remove
mv output.d* ../datafiles-output/
rm input.dump
# go out
cd ../

echo '----------------------------------------------------------------------------------------'
echo '--------------------------------STEP-5:Cluster Analysis---------------------------------'
echo '----------------------------------------------------------------------------------------'
# go to
cd 5-ClusterAnalysis/
# inputs
cp ../datafiles-output/output_mirroratoms.data ./input.data
cp ../datafiles-output/output.data ./input_wrapped.data
# run
gfortran ClusterAnalysis.f90
./a.out
# outputs and remove
mv output_cluster.dump ../datafiles-output/
rm input*
# go out
cd ../

echo '----------------------------------------------------------------------------------------'
echo '--------------------------------STEP-6:Position Average---------------------------------'
echo '----------------------------------------------------------------------------------------'
# go to
cd 6-AverageMirrorAtomPositions/
# inputs
cp ../datafiles-output/output_cluster.dump ./input_cluster.dump
cp ../datafiles-input/dump.0 ./input_added.dump
# run
gfortran PositionAverage.f90
./a.out
# outputs and remove
mv filledbox.* ../datafiles-output/
rm input_*
# go out
cd ../

echo '==>>VAC FILL IS DONE!'

################# Remove inputs and redundant outputs ###############
# go to
cd datafiles-output/
# remove
rm output*
# go out
cd ../

## go to
#cd datafiles-input/
## remove
#rm disline.input dump.0 type0.data
## go out
#cd ../

echo '==>>Remove inputs and redundant outputs, DONE!'


