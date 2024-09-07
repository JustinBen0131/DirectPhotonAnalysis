#!/bin/bash
mkdir -p FileLists/
mkdir -p FileLists/list_allFEM_clock/

workplace=$(pwd)

cd $workplace
echo " "
echo "1. Generating Golden Calorimeter run list"
echo " "
echo "Processing GoldenCaloRunListGenerator.py..."
python3 GoldenCaloRunListGenerator.py

cd $workplace
echo " "
echo "2. Generating GL1 event number list"
echo " "
echo "Processing get_eventnumberlist.sh..."
sh GenerateEventNumberList.sh

cd $workplace
echo " "
echo "3. Checking GL1 event number list"
root -l -q -b GoldenGL1RunListGenerator.C

cd $workplace
echo " "
echo "4. Generating FEM clock list"
echo " "
echo "Processing get_clocklist.sh..."
sh GenerateClockList.sh

cd $workplace
echo " "
echo "5. Checking FEM clock list"
root -l -q -b GoldenFEMrunListGenerator.C

cd $workplace
echo " "
echo "6. Generating Calo + GL1 + FEM golden run number list"
root -l -q -b FinalizeGoldenRunList_calo_gl1_fem.C
cd $workplace
cp FileLists/FinalGoldenRunList_calo_gl1_femChecked.txt ../GoldenRunNumbers_afterRun46619.txt


cd $workplace
echo " "
echo "7. Generating DST list"
echo " "
sh generateEveryPossibleDSTlist.sh


echo " "
echo "Done. A runnumber.txt is generated with the golden runs. DST files are also generated in dst_list/."



