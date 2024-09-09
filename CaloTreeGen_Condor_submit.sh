#!/usr/bin/env bash

# Ensure necessary directories exist
setup_directories() {
  mkdir -p /tmp/patsfan753_condor_logs
  mkdir -p /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/stdout
  mkdir -p /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/error
}

# Submit a job to condor
submit_condor_job() {
  local runNumber=$1
  local filename=$2
  local fileBaseName=$(basename "$filename")

  # Create a submission file for each job
  cat > triggerCondor_single_${fileBaseName%.*}.sub <<EOL
universe                = vanilla
executable              = CaloTreeGen_Condor.sh
arguments               = ${runNumber} ${filename} \$(Cluster)
log                     = /tmp/patsfan753_condor_logs/job_single.\$(Cluster).\$(Process).log
output                  = /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/stdout/job_single.\$(Cluster).\$(Process).out
error                   = /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/error/job_single.\$(Cluster).\$(Process).err
request_memory          = 1.5GB
PeriodicHold            = (NumJobStarts>=1 && JobStatus == 1)
concurrency_limits      = CONCURRENCY_LIMIT_DEFAULT:100
queue
EOL

  # Submit the job
  condor_submit triggerCondor_single_${fileBaseName%.*}.sub
}

# Run the ROOT macro locally
run_local_job() {
  local runNumber=$1
  local filename=$(head -n 1 dst_list/dst_calo_run2pp-000${runNumber}.list)

  # Set the output directory
  local outputDir="/sphenix/tg/tg01/bulk/jbennett/DirectPhotons/output/${runNumber}"
  mkdir -p ${outputDir}

  # Construct the output filename with the base name
  local fileBaseName=$(basename "$filename")
  local treeOutName="${outputDir}/CaloTreeGen_${fileBaseName%.*}.root"

  # Run the ROOT macro locally
  root -b -l -q "macro/Fun4All_CaloTreeGen.C(0, \"$filename\", \"$treeOutName\")"
}

# Main function to submit the first 100 jobs for memory testing
test_condor_memory() {
  echo "Submitting the first 100 jobs..."

  local runNumbers=($(cat GoldenRunNumbers_afterRun46619_part4.txt))
  local runNumber=${runNumbers[0]}

  # Ensure directories are set up
  setup_directories

  # Extract first 100 lines or as many as exist if fewer
  local dst_list="dst_list/dst_calo_run2pp-000${runNumber}.list"
  
  if [ -f "$dst_list" ]; then
    total_lines=$(wc -l < "$dst_list")
    if [ "$total_lines" -lt 100 ]; then
      echo "Only $total_lines segments available, submitting all."
      job_limit=$total_lines
    else
      job_limit=100
    fi

    # Loop through the first `job_limit` lines of the dst list for the run number
    for filename in $(head -n $job_limit "$dst_list"); do
      if [ -n "$filename" ]; then
        submit_condor_job "$runNumber" "$filename"
      fi
    done

  else
    echo "File $dst_list does not exist. Exiting."
    exit 1
  fi
}

# Main function to submit all jobs for all run numbers
submit_all_jobs() {
  echo "Submitting all jobs..."

  local runNumbers=($(cat GoldenRunNumbers_afterRun46619_part2.txt))

  # Ensure directories are set up
  setup_directories

  # Loop over each run number and create custom submission files
  for runNumber in "${runNumbers[@]}"; do
    cat > triggerCondor_${runNumber}.sub <<EOL
universe                = vanilla
executable              = CaloTreeGen_Condor.sh
arguments               = ${runNumber} \$(filename) \$(Cluster)
log                     = /tmp/patsfan753_condor_logs/job.\$(Cluster).\$(Process).log
output                  = /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/stdout/job.\$(Cluster).\$(Process).out
error                   = /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/error/job.\$(Cluster).\$(Process).err
request_memory          = 1.5GB
PeriodicHold            = (NumJobStarts>=1 && JobStatus == 1)
concurrency_limits      = CONCURRENCY_LIMIT_DEFAULT:100
queue filename from dst_list/dst_calo_run2pp-000${runNumber}.list
EOL

    # Submit the job
    condor_submit triggerCondor_${runNumber}.sub
  done
}

# Main entry point for the script
main() {
  case "$1" in
    "local")
      echo "Running locally on the first file of the first run number..."
      runNumbers=($(cat GoldenRunNumbers_afterRun46619_part2.txt))
      run_local_job "${runNumbers[0]}"
      ;;
    "testCondorMemory")
      test_condor_memory
      ;;
    *)
      submit_all_jobs
      ;;
  esac
}

# Call the main function with command-line arguments
main "$@"
