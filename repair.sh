#!/bin/bash


echo '@@@@ Repair Start'

echo '######################################################################################## 0 prepare works'
echo '                                 Backups and Uniform Parameters                   '



echo '---------Backup'
# INPUTS
cp ./datafiles-theINI/* ./datafiles-input/

echo '---------Uniforming parameters'
# uniform-parameters
cp repair.params vacfill/repair.params 
cp repair.params csp/repair.params 
cp repair.params subtract/repair.params



echo '######################################################################################## 1 vacfill'
echo '                          Works on single vac and small vac-clusters                    '
# go to
cd vacfill/
# go to input-directory
cd datafiles-input/
# inputs
cp ../../datafiles-input/* ./
# go out of input-directory
cd ../
# run
sh vac_fill.sh
# copy dump-format filled box as new dump.0 !!!
cp datafiles-output/filledbox.dump ../datafiles-input/dump.0
# move all outputs to repair-output-dir
mv datafiles-output/filledbox.* ../datafiles-output/
# go out
cd ../

echo '######################################################################################## 2 csp'
echo '                     Selecting new type0.data file of filledbox for repair '
# go to
cd csp/
# go to input-dir
cd datafile-input/
# inputs, new dump.0 from vacfill
cp ../../datafiles-input/* ./
# go out
cd ../
# run
sh csp.sh
# NEW type0.data file !!!
mv datafile-output/bad.data ../datafiles-input/type0.data
# remove other outputs
rm datafile-output/bad.dump
rm datafile-output/good.d*
# go out
cd ../

echo '######################################################################################## 3 subtract'
echo '   Subtract type0 atoms from dump0 and create a new dump0 with namy vacancies and viods '
# go to
cd subtract/
# inputs
cp ../datafiles-input/* ./datafile-input/
# run
sh subtract.sh
# NEW dump.0 with many vancancies
mv datafile-output/c.dump.output ../datafiles-input/dump.0
# rm data format output
rm datafile-output/c.data.output
# go out
cd ../

echo '######################################################################################## 4 csp'
echo '                     Selecting new type0.data file of filledbox for repair '
# go to
cd csp/
# go to input-dir
cd datafile-input/
# inputs, new dump.0 from vacfill
cp ../../datafiles-input/* ./
# go out
cd ../
# run
sh csp.sh
# NEW type0.data file !!!
mv datafile-output/bad.data ../datafiles-input/type0.data
# remove other outputs
rm datafile-output/bad.dump
rm datafile-output/good.d*
# go out
cd ../

echo '######################################################################################## 5 vacfill'
echo '                          Works on single vac and small vac-clusters                    '
# go to
cd vacfill/
# go to input-directory
cd datafiles-input/
# inputs
cp ../../datafiles-input/* ./
# go out of input-directory
cd ../
# run
sh vac_fill.sh
# copy dump-format filled box as new dump.0 !!!
cp datafiles-output/filledbox.dump ../datafiles-input/dump.0
# move all outputs to repair-output-dir
mv datafiles-output/filledbox.* ../datafiles-output/
# go out
cd ../

echo '######################################################################################## 6 csp'
echo '                     Selecting new type0.data file of filledbox for repair '
# go to
cd csp/
# go to input-dir
cd datafile-input/
# inputs, new dump.0 from vacfill
cp ../../datafiles-input/* ./
# go out
cd ../
# run
sh csp.sh
# NEW type0.data file !!!
mv datafile-output/bad.data ../datafiles-input/type0.data
# remove other outputs
rm datafile-output/bad.dump
rm datafile-output/good.d*
# go out
cd ../

echo '######################################################################################## 7 vacfill'
echo '                          Works on single vac and small vac-clusters                    '
# go to
cd vacfill/
# go to input-directory
cd datafiles-input/
# inputs
cp ../../datafiles-input/* ./
# go out of input-directory
cd ../
# run
sh vac_fill.sh
# copy dump-format filled box as new dump.0 !!!
cp datafiles-output/filledbox.dump ../datafiles-input/dump.0
# move all outputs to repair-output-dir
mv datafiles-output/filledbox.* ../datafiles-output/
# go out
cd ../

echo '######################################################################################## 8 csp'
echo '                     Selecting new type0.data file of filledbox for repair '
# go to
cd csp/
# go to input-dir
cd datafile-input/
# inputs, new dump.0 from vacfill
cp ../../datafiles-input/* ./
# go out
cd ../
# run
sh csp.sh
# NEW type0.data file !!!
mv datafile-output/bad.data ../datafiles-input/type0.data
# remove other outputs
rm datafile-output/bad.dump
rm datafile-output/good.d*
# go out
cd ../

echo '######################################################################################## 9 vacfill'
echo '                          Works on single vac and small vac-clusters                    '
# go to
cd vacfill/
# go to input-directory
cd datafiles-input/
# inputs
cp ../../datafiles-input/* ./
# go out of input-directory
cd ../
# run
sh vac_fill.sh
# copy dump-format filled box as new dump.0 !!!
cp datafiles-output/filledbox.dump ../datafiles-input/dump.0
# move all outputs to repair-output-dir
mv datafiles-output/filledbox.* ../datafiles-output/
# go out
cd ../

#echo '######################################################################################## 10 csp'
#echo '                     Selecting new type0.data file of filledbox for repair '
## go to
#cd csp/
## go to input-dir
#cd datafile-input/
## inputs, new dump.0 from vacfill
#cp ../../datafiles-input/* ./
## go out
#cd ../
## run
#sh csp.sh
## NEW type0.data file !!!
#mv datafile-output/bad.data ../datafiles-input/type0.data
## remove other outputs
#rm datafile-output/bad.dump
#rm datafile-output/good.d*
## go out
#cd ../
#
#echo '######################################################################################## 11 vacfill'
#echo '                          Works on single vac and small vac-clusters                    '
## go to
#cd vacfill/
## go to input-directory
#cd datafiles-input/
## inputs
#cp ../../datafiles-input/* ./
## go out of input-directory
#cd ../
## run
#sh vac_fill.sh
## copy dump-format filled box as new dump.0 !!!
#cp datafiles-output/filledbox.dump ../datafiles-input/dump.0
## move all outputs to repair-output-dir
#mv datafiles-output/filledbox.* ../datafiles-output/
## go out
#cd ../


echo '######################################################################################## Repair, DONE!'
echo '@@@@ Repair End'



