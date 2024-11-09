#!/usr/bin/env bash

# Read run numbers from the input file
inputFile="/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/runListTriggerLUTv1.txt"

# Read run numbers into an array
mapfile -t runNumbers < ${inputFile}

# Ensure directories are set up
mkdir -p /tmp/patsfan753_condor_logs
mkdir -p /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/stdout
mkdir -p /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/error

# Check if 'local' is passed as an argument
if [ "$1" == "local" ]; then
    echo "Running locally for testing purposes."

    testRunNumber=${runNumbers[0]}

    # Execute the processing script locally
    ./mergeSegmentFilesForRuns.sh ${testRunNumber}
else
    # Loop over each run number and create custom submission files
    for runNumber in "${runNumbers[@]}"; do
        cat > condor_${runNumber}.submit <<EOL
universe                = vanilla
executable              = mergeSegmentFilesForRuns.sh
arguments               = ${runNumber}
log                     = /tmp/patsfan753_condor_logs/job.\$(Cluster).\$(Process).log
output                  = /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/stdout/job.\$(Cluster).\$(Process).out
error                   = /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/error/job.\$(Cluster).\$(Process).err
request_memory          = 1.5GB
PeriodicHold            = (NumJobStarts>=1 && JobStatus == 1)
concurrency_limits      = CONCURRENCY_LIMIT_DEFAULT:100
queue
EOL

        # Submit the job
        condor_submit condor_${runNumber}.submit
    done
fi

