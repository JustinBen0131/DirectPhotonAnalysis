#!/usr/bin/env bash

# Ensure directories are set up
mkdir -p /tmp/${USER}_condor_logs
mkdir -p /sphenix/user/${USER}/tutorials/tutorials/CaloDataAnaRun24pp/stdout
mkdir -p /sphenix/user/${USER}/tutorials/tutorials/CaloDataAnaRun24pp/error

# Read run numbers from the input file
inputFile="/sphenix/user/${USER}/tutorials/tutorials/CaloDataAnaRun24pp/Final_RunNumbers_After_All_Cuts.txt"
mapfile -t runNumbers < ${inputFile}

# Check if 'local' is passed as an argument
if [ "$1" == "local" ]; then
    if [ -z "$2" ]; then
        echo "Please provide a run number or a file containing run numbers for local processing."
        exit 1
    fi

    if [ -f "$2" ]; then
        # The second argument is a file containing run numbers
        mapfile -t localRunNumbers < "$2"
        echo "Running locally on run numbers from file ${2}."
        for runNumber in "${localRunNumbers[@]}"; do
            echo "Processing run number ${runNumber} locally."
            ./mergeSegmentFilesForRuns.sh ${runNumber}
        done
    else
        # The second argument is a single run number
        testRunNumber=$2
        echo "Running locally for testing purposes on run number ${testRunNumber}."
        # Execute the processing script locally
        ./mergeSegmentFilesForRuns.sh ${testRunNumber}
    fi
else
    # Existing functionality: Submit jobs to Condor
    for runNumber in "${runNumbers[@]}"; do
        cat > condor_${runNumber}.submit <<EOL
universe                = vanilla
executable              = mergeSegmentFilesForRuns.sh
arguments               = ${runNumber}
log                     = /tmp/${USER}_condor_logs/job.\$(Cluster).\$(Process).log
output                  = /sphenix/user/${USER}/tutorials/tutorials/CaloDataAnaRun24pp/stdout/job.\$(Cluster).\$(Process).out
error                   = /sphenix/user/${USER}/tutorials/tutorials/CaloDataAnaRun24pp/error/job.\$(Cluster).\$(Process).err
request_memory          = 1.5GB
PeriodicHold            = (NumJobStarts>=1 && JobStatus == 1)
concurrency_limits      = CONCURRENCY_LIMIT_DEFAULT:100
queue
EOL

        # Submit the job
        condor_submit condor_${runNumber}.submit
    done
fi
