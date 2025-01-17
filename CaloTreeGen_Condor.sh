#!/usr/bin/env bash
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}
export MYINSTALL="/sphenix/user/patsfan753/install"

source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

# Capture command line arguments
runNumber=$1
fileList=$2
clusterID=$3
events=0  # Default number of events



if [ "$runNumber" == "SIMMODE" ]; then
  echo "[INFO] Detected simulation mode."
  # e.g., no runNumber needed
  # we can define a special output directory for sim
  outputDir="/sphenix/tg/tg01/bulk/jbennett/DirectPhotons/output/SimOut"
  mkdir -p ${outputDir}

  # The file list for sim is presumably a single file with 1 line
  mapfile -t rootFiles < "$fileList"
  fileBaseName=$(basename "${rootFiles[0]}")
  treeOutName="${outputDir}/CaloTreeGenSim_${fileBaseName%.*}.root"

  # use runSim=true, runData=false
  root -b -l -q "macro/Fun4All_CaloTreeGen.C(0, \"$fileList\", \"$treeOutName\", true, false)"

else
    echo "[INFO] Detected data mode. runNumber = $runNumber"
    # Create a directory based on the run number
    outputDir="/sphenix/tg/tg01/bulk/jbennett/DirectPhotons/output/${runNumber}"
    mkdir -p ${outputDir}

    # Check if the file list is provided
    if [ ! -f "$fileList" ]; then
      echo "[ERROR] File list not provided or does not exist: $fileList"
      exit 1
    fi

    # Read the list of files from the text file
    mapfile -t rootFiles < "$fileList"

    # Construct the output filename using the first file in the list
    fileBaseName=$(basename "${rootFiles[0]}")
    treeOutName="${outputDir}/CaloTreeGen_${fileBaseName%.*}.root"

    # Run the ROOT macro with the list of input files
    root -b -l -q "macro/Fun4All_CaloTreeGen.C(0, \"$fileList\", \"$treeOutName\", false, true)"
fi
