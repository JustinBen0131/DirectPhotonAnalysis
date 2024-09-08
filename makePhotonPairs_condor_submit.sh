#!/usr/bin/env bash

# Check if running locally
if [ "$1" == "local" ]; then
  # Run locally using the first segment file from the first run number
  echo "Running locally on the first file of the first run number..."

  # Read the run numbers from the file
  runNumbers=($(cat GoldenRunNumbers_afterRun46619_part1.txt))

  # Grab the first run number
  runNumber=${runNumbers[0]}

  # Get the first filename from the dst list for that run number
  filename=$(head -n 1 dst_list/dst_calo_run2pp-000${runNumber}.list)

  # Ensure necessary directories exist
  mkdir -p /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/stdout
  mkdir -p /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/error

  # Set the output directory
  outputDir="/sphenix/tg/tg01/bulk/jbennett/DirectPhotons/outputInvMassHistograms/${runNumber}"
  mkdir -p ${outputDir}  # Ensure the directory exists

  # Get just the base name of the file (strip the directory part)
  fileBaseName=$(basename "$filename")

  # Construct the output filename with the base name
  treeOutName="${outputDir}/InvMassDists_${fileBaseName%.*}.root"

  # Run the ROOT macro locally
  root -b -l -q "macro/makePhotonPairs.C(0, \"$filename\", \"$treeOutName\")"

else
  # Read the run numbers from the file and store them in the array
  runNumbers=($(cat GoldenRunNumbers_afterRun46619_part1.txt))

  # Ensure necessary directories exist
  mkdir -p /tmp/patsfan753_condor_logs
  mkdir -p /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/stdout
  mkdir -p /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/error

  # Loop over each run number and create a custom submission file
  for runNumber in "${runNumbers[@]}"; do
    # Create a custom submission file for the current run number
    cat > triggerCondor_${runNumber}.sub <<EOL
universe                = vanilla
executable              = CaloTreeGen_Condor.sh
arguments               = ${runNumber} \$(filename) \$(Cluster)
log                     = /tmp/patsfan753_condor_logs/job.\$(Cluster).\$(Process).log
output                  = /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/stdout/job.\$(Cluster).\$(Process).out
error                   = /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/error/job.\$(Cluster).\$(Process).err
request_memory          = 7GB
queue filename from invMassHist_lists/invMass_SegmentsHistList-000${runNumber}.list
EOL

    # Submit the job
    condor_submit triggerCondor_${runNumber}.sub
  done
fi
