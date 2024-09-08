#!/usr/bin/env bash
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}
export MYINSTALL="/sphenix/user/patsfan753/install"

source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

# Capture command line arguments
runNumber=$1
filename=$2
clusterID=$3
events=0  # Default number of events

# Create a directory based on the run number
outputDir="/sphenix/tg/tg01/bulk/jbennett/DirectPhotons/outputInvMassHistograms/${runNumber}"
mkdir -p ${outputDir}  # Ensure the directory exists

# Get just the base name of the file (strip the directory part)
fileBaseName=$(basename "$filename")

# Construct the output filename with the base name
treeOutName="${outputDir}/InvMassDists_${fileBaseName%.*}.root"

# Check if filename is provided
if [ -z "$filename" ]; then
  echo "Error: Filename is not provided."
  exit 1
fi

# Run the ROOT macro with dynamically generated output paths and log the output
root -b -l -q "macro/makePhotonPairs.C(0, \"$filename\", \"$treeOutName\")"
