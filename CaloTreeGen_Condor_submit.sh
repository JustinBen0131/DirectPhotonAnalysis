#!/usr/bin/env bash

###########################################
# 1) Ensure necessary directories exist
###########################################
setup_directories() {
  mkdir -p /tmp/patsfan753_condor_logs
  mkdir -p /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/stdout
  mkdir -p /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/error
}

###########################################
# 2) Data: local
###########################################
run_local_job_data() {
  local runNumber="$1"

  source /opt/sphenix/core/bin/sphenix_setup.sh -n
  source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

  if [ -z "$runNumber" ]; then
    # If user didn't supply a run number, read from file
    runNumbers=($(cat RunCondorJobsOnTheseRuns.txt))
    runNumber="${runNumbers[0]}"
  fi

  echo "[INFO] Running locally with run number: $runNumber"

  local dst_list="dst_list/dst_calo_run2pp-000${runNumber}.list"
  local fileList="dst_list/dst_calo_run2pp-${runNumber}_input.list"

  mapfile -t files < "$dst_list"

  if [ ${#files[@]} -lt 1 ]; then
    echo "[ERROR] No files found in $dst_list."
    exit 1
  fi

  # We can pick the first file or first two files, etc.
  echo -e "${files[0]}" > "$fileList"

  # Prepare the output directory for this run
  local outputDir="/sphenix/tg/tg01/bulk/jbennett/DirectPhotons/output/${runNumber}"

  if [ -d "$outputDir" ]; then
    echo "[WARNING] Output directory already exists: $outputDir"
    echo "          Removing existing files to avoid mixing old/new data."
    rm -f "$outputDir"/*
  else
    mkdir -p "$outputDir"
  fi

  fileBaseName=$(basename "${files[0]}")
  treeOutName="${outputDir}/CaloTreeGen_${fileBaseName%.*}.root"

  echo "[INFO] Writing output to: $treeOutName"

  # Pass runData=true, runSim=false
  root -b -l -q "macro/Fun4All_CaloTreeGen.C(0, \"$fileList\", \"$treeOutName\", false, true)"
}


###########################################
# 3) Data: submit condor
###########################################
submit_all_jobs_data() {
  echo "[INFO] Submitting all data condor jobs..."

  setup_directories

  local runNumbers=($(cat RunCondorJobsOnTheseRuns.txt))
  for runNumber in "${runNumbers[@]}"; do

    local dst_list="dst_list/dst_calo_run2pp-000${runNumber}.list"
    local paired_list="dst_list/dst_calo_run2pp-000${runNumber}_paired.list"
    if [ ! -f "$dst_list" ]; then
      echo "[ERROR] File not found: $dst_list"
      continue
    fi

    # create or clear
    > "$paired_list"
    mapfile -t filenames < "$dst_list"
    local total_files=${#filenames[@]}

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
        # last odd file
        echo -e "$file1" > "input_files_${runNumber}_${i}.list"
        echo "input_files_${runNumber}_${i}.list" >> "$paired_list"
      fi
    done

    # Build a condor .sub
    cat > isoCondor_${runNumber}.sub <<EOL
universe                = vanilla
executable              = CaloTreeGen_Condor.sh
# pass runNumber and fileList, plus cluster
arguments               = ${runNumber} \$(filename) \$(Cluster)
log                     = /tmp/patsfan753_condor_logs/job.\$(Cluster).\$(Process).log
output                  = /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/stdout/job.\$(Cluster).\$(Process).out
error                   = /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/error/job.\$(Cluster).\$(Process).err
request_memory          = 800MB
PeriodicHold            = (NumJobStarts>=1 && JobStatus == 1)
concurrency_limits      = CONCURRENCY_LIMIT_DEFAULT:100
queue filename from $paired_list
EOL

    condor_submit isoCondor_${runNumber}.sub
  done
}


###########################################
# 4) SIM: local
###########################################
run_local_sim() {
  source /opt/sphenix/core/bin/sphenix_setup.sh -n
  source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

  local simList="../caloAnaSimListFile.list"
  if [ ! -f "$simList" ]; then
    echo "[ERROR] No sim list file found: $simList"
    exit 1
  fi

  mapfile -t simfiles < "$simList"
  if [ ${#simfiles[@]} -eq 0 ]; then
    echo "[ERROR] No simulation files in $simList!"
    exit 1
  fi

  local firstSimFile="${simfiles[0]}"
  echo "[INFO] Running local simulation for file: $firstSimFile"

  local outputDirSim="/sphenix/tg/tg01/bulk/jbennett/DirectPhotons/output/SimOut"

  if [ -d "$outputDirSim" ]; then
    echo "[WARNING] Simulation output directory already exists: $outputDirSim"
    echo "          Removing existing files to avoid mixing old/new data."
    rm -f "$outputDirSim"/*
  else
    mkdir -p "$outputDirSim"
  fi

  fileBaseName=$(basename "${firstSimFile}")
  local simOutName="${outputDirSim}/CaloTreeGenSim_${fileBaseName%.*}.root"

  echo "[INFO] Writing simulation output to: $simOutName"

  # runSim=true, runData=false
  root -b -l -q "macro/Fun4All_CaloTreeGen.C(0, \"$simList\", \"$simOutName\", true, false)"
}


###########################################
# 5) SIM: submit condor
###########################################
submit_all_jobs_sim() {
  echo "[INFO] Submitting all simulation condor jobs..."

  setup_directories

  local simList="dst_list/caloAnaSimListFile.list"
  if [ ! -f "$simList" ]; then
    echo "[ERROR] simList file not found: $simList"
    exit 1
  fi

  mapfile -t simfiles < "$simList"
  local total_sfiles=${#simfiles[@]}
  if [ $total_sfiles -eq 0 ]; then
    echo "[ERROR] No sim files in $simList."
    exit 1
  fi

  local simQueuedFileList="sim_segment_files.list"
  > "$simQueuedFileList"

  local i=0
  while [ $i -lt $total_sfiles ]; do
    local file="${simfiles[$i]}"
    i=$((i + 1))
    local simlistfile="sim_input_file_${i}.list"
    echo "$file" > "$simlistfile"
    echo "$simlistfile" >> "$simQueuedFileList"
  done

  cat > isoCondorSim.sub <<EOL
universe                = vanilla
executable              = CaloTreeGen_Condor.sh
arguments               = SIMMODE \$(filename) \$(Cluster)
log                     = /tmp/patsfan753_condor_logs/job.\$(Cluster).\$(Process).log
output                  = /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/stdout/job.\$(Cluster).\$(Process).out
error                   = /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/error/job.\$(Cluster).\$(Process).err
request_memory          = 800MB
PeriodicHold            = (NumJobStarts>=1 && JobStatus == 1)
concurrency_limits      = CONCURRENCY_LIMIT_DEFAULT:100
queue filename from $simQueuedFileList
EOL

  condor_submit isoCondorSim.sub
}


###########################################
# 6) MAIN: interpret arguments
###########################################
main() {
  # We want to let the user pass multiple arguments like 'localData localSim',
  # so let's parse them in a loop

  # Flags to track which steps to do
  local doLocalData="false"
  local doLocalSim="false"
  local doCondorData="false"
  local doCondorSim="false"

  if [ $# -eq 0 ]; then
    echo "[ERROR] No arguments specified!"
    echo " Usage: $0 [localData] [localSim] [submitData] [submitSim]"
    echo " E.g. $0 localData localSim   => run data & sim locally"
    echo " E.g. $0 submitData           => submit data to condor"
    echo " E.g. $0 submitData submitSim => submit both to condor"
    exit 1
  fi

  # parse each argument
  for arg in "$@"; do
    case "$arg" in
      localData)
        doLocalData="true"
        ;;
      localSim)
        doLocalSim="true"
        ;;
      submitData)
        doCondorData="true"
        ;;
      submitSim)
        doCondorSim="true"
        ;;
      *)
        echo "[WARNING] Unrecognized argument: $arg"
        ;;
    esac
  done

  # Now let's do them in an order, or in the same script you can do everything
  if [ "$doLocalData" = "true" ]; then
    run_local_job_data ""
  fi
  if [ "$doLocalSim" = "true" ]; then
    run_local_sim
  fi
  if [ "$doCondorData" = "true" ]; then
    submit_all_jobs_data
  fi
  if [ "$doCondorSim" = "true" ]; then
    submit_all_jobs_sim
  fi
}

main "$@"
