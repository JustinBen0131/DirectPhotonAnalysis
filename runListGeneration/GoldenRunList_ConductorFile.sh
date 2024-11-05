#!/bin/bash

# Create necessary directories
mkdir -p FileLists/
mkdir -p FileLists/list_allFEM_clock/
mkdir -p dst_list/

# Set the working directory
workplace=$(pwd)

# Step 1: Generate Golden Calorimeter run list
cd "$workplace"
echo " "
echo "1. Generating Golden Calorimeter run list"
echo " "
echo "Processing GoldenCaloRunListGenerator.py..."
python3 GoldenCaloRunListGenerator.py

# Step 2: Generate GL1 event number list
cd "$workplace"
echo " "
echo "2. Generating GL1 event number list"
echo " "
echo "Processing GenerateEventNumberList.sh..."
sh GenerateEventNumberList.sh

# Step 3: Check GL1 event number list
cd "$workplace"
echo " "
echo "3. Checking GL1 event number list"
root -l -q -b GoldenGL1RunListGenerator.C

# Step 4: Generate FEM clock list
cd "$workplace"
echo " "
echo "4. Generating FEM clock list"
echo " "
echo "Processing GenerateClockList.sh..."
sh GenerateClockList.sh

# Step 5: Check FEM clock list
cd "$workplace"
echo " "
echo "5. Checking FEM clock list"
root -l -q -b GoldenFEMrunListGenerator.C

# Step 6: Generate the final Calo + GL1 + FEM golden run number list
cd "$workplace"
echo " "
echo "6. Generating Calo + GL1 + FEM golden run number list"
root -l -q -b FinalizeGoldenRunList_calo_gl1_fem.C

# Copy the initial final golden run list to a separate file
echo "Copying the initial golden run list to Full_ppGoldenRunList.txt..."
cp FileLists/FinalGoldenRunList_calo_gl1_femChecked.txt ../Full_ppGoldenRunList.txt

# Step 7: Filter the final golden run list based on bad tower maps
echo " "
echo "7. Filtering final golden run list based on bad tower maps"
echo " "

# Find run numbers that have a bad tower map
echo "Finding run numbers that have a bad tower map..."
bad_tower_runs=$(find /cvmfs/sphenix.sdcc.bnl.gov/calibrations/sphnxpro/cdb/CEMC_BadTowerMap -name "*p0*" | cut -d '-' -f2 | cut -d c -f1 | sort | uniq)

# Output the runs found to a temporary file
echo "Listing bad tower runs in bad_tower_runs.txt..."
echo "$bad_tower_runs" > bad_tower_runs.txt

# Check that the initial golden run list exists
golden_run_list="../Full_ppGoldenRunList.txt"
if [[ ! -f "$golden_run_list" ]]; then
    echo "[ERROR] The file $golden_run_list does not exist. Please check the path and try again."
    exit 1
fi

# Filter the final golden run list to keep only runs that have a bad tower map
filtered_golden_run_list="/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/Full_ppGoldenRunList_FinalList_withBadTowerMaps.txt"
echo "Filtering the golden run list to keep only runs with a bad tower map..."
grep -Ff bad_tower_runs.txt "$golden_run_list" > "$filtered_golden_run_list"

# Find the removed run numbers by comparing the initial and filtered lists
echo "Finding removed runs..."
removed_runs=$(comm -23 <(sort "$golden_run_list") <(sort "$filtered_golden_run_list"))

# Print the removed run numbers and count
echo "The following run numbers were removed because they do not have a bad tower map:"
echo "$removed_runs"
echo "Number of runs removed: $(echo "$removed_runs" | wc -l)"

# Clean up temporary files
rm bad_tower_runs.txt


echo " "
echo "Final golden run list has been filtered to include only runs with a bad tower map and saved as Full_ppGoldenRunList_FinalList_withBadTowerMaps.txt."
echo " "

# Step 8: Create a .list file for CreateDstList.pl command
echo "Creating .list file for CreateDstList.pl..."
cp "$filtered_golden_run_list" ../Full_ppGoldenRunList_FinalList_withBadTowerMaps.list

# Step 9: Change to the dst_list directory and generate DST list using CreateDstList.pl
echo "Changing to dst_list directory to generate DST lists..."
cd ../dst_list

echo "Running CreateDstList.pl to generate the DST list..."
CreateDstList.pl --build ana437 --cdb 2024p007 DST_CALO_run2pp --list ../Full_ppGoldenRunList_FinalList_withBadTowerMaps.list

echo " "
echo "DST list generated and saved to ../dst_list for the runs present in Full_ppGoldenRunList_FinalList_withBadTowerMaps.txt."

