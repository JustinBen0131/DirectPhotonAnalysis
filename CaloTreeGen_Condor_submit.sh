#!/usr/bin/env bash

# Ensure necessary directories exist
setup_directories() {
  mkdir -p /tmp/patsfan753_condor_logs
  mkdir -p /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/stdout
  mkdir -p /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/error
}

run_local_job() {
  local runNumber=$1

  if [ -z "$runNumber" ]; then
    runNumbers=($(cat RunCondorJobsOnTheseRuns.txt))
    runNumber="${runNumbers[0]}"
  fi

  echo "Running locally with run number: $runNumber"

  # Create a temporary paired list for the given run number
  local dst_list="dst_list/dst_calo_run2pp-000${runNumber}.list"
  local fileList="dst_list/dst_calo_run2pp-${runNumber}_input.list"

  mapfile -t files < "$dst_list"

  if [ ${#files[@]} -lt 2 ]; then
    echo "Insufficient files for pairing in $dst_list."
    exit 1
  fi

  # Create a text file containing two ROOT files
  echo -e "${files[0]}\n${files[1]}" > "$fileList"

  local outputDir="/sphenix/tg/tg01/bulk/jbennett/DirectPhotons/output/${runNumber}"
  mkdir -p ${outputDir}

  fileBaseName=$(basename "${files[0]}")
  treeOutName="${outputDir}/CaloTreeGen_${fileBaseName%.*}.root"

  root -b -l -q "macro/Fun4All_CaloTreeGen.C(0, \"$fileList\", \"$treeOutName\")"
}


submit_all_jobs() {
  echo "Submitting all jobs..."

  local runNumbers=($(cat RunCondorJobsOnTheseRuns.txt))
  
  # Ensure directories are set up
  setup_directories

  # Loop over each run number and create custom submission files
  for runNumber in "${runNumbers[@]}"; do
    # Read the list of filenames for this run number
    local dst_list="dst_list/dst_calo_run2pp-000${runNumber}.list"
    local paired_list="dst_list/dst_calo_run2pp-000${runNumber}_paired.list"
    
    if [ ! -f "$dst_list" ]; then
      echo "[ERROR] File not found: $dst_list"
      continue
    fi

    # Clear or create the paired list file
    > "$paired_list"

    # Read all filenames into an array
    mapfile -t filenames < "$dst_list"
    local total_files=${#filenames[@]}

    # Pair the files and write to the paired list
    local i=0
    while [ $i -lt $total_files ]; do
      local file1="${filenames[$i]}"
      i=$((i + 1))
      
      if [ $i -lt $total_files ]; then
        local file2="${filenames[$i]}"
        i=$((i + 1))
        echo -e "$file1\n$file2" > "input_files_${runNumber}_${i}.list"
        echo "input_files_${runNumber}_${i}.list" >> "$paired_list"
      else
        # Handle the odd file separately
        echo -e "$file1" > "input_files_${runNumber}_${i}.list"
        echo "input_files_${runNumber}_${i}.list" >> "$paired_list"
      fi
    done

    # Create the Condor submission file for the paired lists
    cat > isoCondor_${runNumber}.sub <<EOL
universe                = vanilla
executable              = CaloTreeGen_Condor.sh
arguments               = ${runNumber} \$(filename) \$(Cluster)
log                     = /tmp/patsfan753_condor_logs/job.\$(Cluster).\$(Process).log
output                  = /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/stdout/job.\$(Cluster).\$(Process).out
error                   = /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/error/job.\$(Cluster).\$(Process).err
request_memory          = 800MB
PeriodicHold            = (NumJobStarts>=1 && JobStatus == 1)
concurrency_limits      = CONCURRENCY_LIMIT_DEFAULT:100
queue filename from $paired_list
EOL

    # Submit the job
    condor_submit isoCondor_${runNumber}.sub
  done
}


# Main entry point for the script
main() {
  case "$1" in
    "local")
      runNumber=$2
      run_local_job "$runNumber"
      ;;
    *)
      submit_all_jobs
      ;;
  esac
}

# Call the main function with command-line arguments
main "$@"
