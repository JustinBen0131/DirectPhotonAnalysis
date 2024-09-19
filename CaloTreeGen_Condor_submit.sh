#!/usr/bin/env bash

# Function to ensure necessary directories exist
setup_directories() {
  mkdir -p /tmp/patsfan753_condor_logs || { echo "Failed to create log directory"; exit 1; }
  mkdir -p /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/stdout || { echo "Failed to create stdout directory"; exit 1; }
  mkdir -p /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/error || { echo "Failed to create error directory"; exit 1; }
}

# Function to submit a job to Condor using local disk
submit_condor_job() {
  local runNumber=$1
  local filename=$2
  local fileBaseName=$(basename "$filename")

  # Create a Condor submission file for each job
  cat > isoCondor_single_${fileBaseName%.*}.sub <<EOL
universe                = vanilla
executable              = CaloTreeGen_Condor.sh
arguments               = ${runNumber} ${filename} \$(Cluster)
log                     = /tmp/patsfan753_condor_logs/job_single.\$(Cluster).\$(Process).log
output                  = /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/stdout/job_single.\$(Cluster).\$(Process).out
error                   = /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/error/job_single.\$(Cluster).\$(Process).err
request_memory          = 1.5GB
PeriodicHold            = (NumJobStarts>=1 && JobStatus == 1)
concurrency_limits      = CONCURRENCY_LIMIT_DEFAULT:100
should_transfer_files   = NO
queue
EOL

  # Submit the job to Condor
  condor_submit isoCondor_single_${fileBaseName%.*}.sub || { echo "Failed to submit job for ${filename}"; exit 1; }
}

# Function to run the ROOT macro in $_CONDOR_SCRATCH_DIR
run_condor_job() {
  local runNumber=$1
  local filename=$2

  # Set up environment
  source /opt/sphenix/core/bin/sphenix_setup.sh -n || { echo "Failed to source sphenix setup"; exit 1; }

  # Use local scratch directory for job execution
  if [[ -n "$_CONDOR_SCRATCH_DIR" && -d "$_CONDOR_SCRATCH_DIR" ]]; then
    cd "$_CONDOR_SCRATCH_DIR" || { echo "Failed to change directory to $_CONDOR_SCRATCH_DIR"; exit 1; }
    rsync -av "$(dirname "$0")/" . || { echo "Failed to rsync files to local scratch"; exit 1; }

    # Copy input file to local scratch
    cp "$filename" . || { echo "Failed to copy input file to local scratch"; exit 1; }

    # Define output file name
    local fileBaseName=$(basename "$filename")
    local treeOutName="CaloTreeGen_${fileBaseName%.*}.root"

    # Run the ROOT macro from local scratch
    root -b -l -q "macro/Fun4All_CaloTreeGen.C(0, \"$fileBaseName\", \"$treeOutName\")" || { echo "ROOT macro execution failed"; exit 1; }

    # Copy output back to central storage
    cp "$treeOutName" "/sphenix/tg/tg01/bulk/jbennett/DirectPhotons/output/${runNumber}/" || { echo "Failed to copy output file back to central storage"; exit 1; }

  else
    echo "condor scratch NOT set or not found"
    exit 1
  fi
}

# Function to submit the first 100 jobs for memory testing
test_condor_memory() {
  echo "Submitting the first 100 jobs..."

  # Read run numbers from file
  local runNumbers=($(cat GoldenRunNumbers_afterRun46619_part1.txt))
  local runNumber=${runNumbers[0]}

  # Ensure directories are set up
  setup_directories

  # Define the destination list file
  local dst_list="dst_list/dst_calo_run2pp-000${runNumber}.list"

  if [ -f "$dst_list" ]; then
    total_lines=$(wc -l < "$dst_list")
    job_limit=$(( total_lines < 100 ? total_lines : 100 ))

    echo "Submitting $job_limit jobs..."

    # Loop through the first `job_limit` lines of the destination list
    head -n $job_limit "$dst_list" | while read -r filename; do
      if [ -n "$filename" ]; then
        submit_condor_job "$runNumber" "$filename"
      fi
    done

  else
    echo "File $dst_list does not exist. Exiting."
    exit 1
  fi
}

# Function to submit all jobs for all run numbers
submit_all_jobs() {
  echo "Submitting all jobs..."

  # Read run numbers from file
  local runNumbers=($(cat GoldenRunNumbers_afterRun46619_part1.txt))

  # Ensure directories are set up
  setup_directories

  # Loop over each run number to create and submit Condor submission files
  for runNumber in "${runNumbers[@]}"; do
    cat > isoCondor_${runNumber}.sub <<EOL
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

    # Submit the job to Condor
    condor_submit isoCondor_${runNumber}.sub || { echo "Failed to submit job for run number ${runNumber}"; exit 1; }
  done
}

# Main entry point for the script
main() {
  case "$1" in
    "local")
      echo "Running locally on the first file of the first run number..."
      runNumbers=($(cat GoldenRunNumbers_afterRun46619_part2.txt))
      run_condor_job "${runNumbers[0]}" || { echo "Failed to run job locally"; exit 1; }
      ;;
    "test_condor_memory")
      test_condor_memory
      ;;
    *)
      submit_all_jobs
      ;;
  esac
}

# Call the main function with command-line arguments
main "$@"
