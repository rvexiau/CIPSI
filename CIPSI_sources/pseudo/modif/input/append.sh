#!/bin/bash

[ $# -ge 1 -a -f "$1" ] && input="$1" || exit

cd ..
./deloc<input/$input
cd ..

cp $CIPSI_ROOT/pseudo/MWBPOT_back ./MWBPOT
cp $CIPSI_ROOT/pseudo/PSNL_back ./PSNL
./read_PSNL
./write_PSNL

rm -f MWBPOT_FORM
./read_MWBPOT
cat "modif/MWBPOT_new" >> "MWBPOT_FORM"
./write_MWBPOT



mv MWBPOT_NEW modif/output/MWBPOT
mv MWBPOT_FORM modif/output/MWBPOT_FORM
mv PSNL_NEW modif/output/PSNL
mv PSNL_FORM modif/output/PSNL_FORM
rm MWBPOT
rm PSNL
rm modif/MWBPOT_new
rm modif/fort.*
rm modif/psnl_gaussian*