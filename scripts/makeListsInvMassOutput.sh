#!/usr/bin/env bash

# Define base directories
outputBaseDir="/sphenix/tg/tg01/bulk/jbennett/DirectPhotons/output"
listOutputDir="/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/invMassHist_lists"

# Ensure the output directory exists
mkdir -p "$listOutputDir"

# Loop through each runNumber folder in the output base directory
for runDir in "$outputBaseDir"/*/; do
  # Extract the runNumber from the folder name
  runNumber=$(basename "$runDir")

  # Define the output list file for this runNumber
  listFile="$listOutputDir/invMass_SegmentsHistList-000${runNumber}.list"

  # Ensure the list file is empty before writing
  > "$listFile"

  # Loop through each file in the current runNumber folder
  for segmentFile in "$runDir"CaloTreeGen_*.root; do
    # Write the full path of the segment file to the list file
    echo "$segmentFile" >> "$listFile"
  done

  # Print message confirming list file creation
  echo "List file created: $listFile"
done
