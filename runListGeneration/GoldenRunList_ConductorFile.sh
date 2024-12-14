#!/bin/bash

########################################
# Steps:
# 1. Extract initial run numbers from databases.
# 2. Apply Calorimeter QA filters.
# 3. Apply run duration and livetime cuts.
# 4. Identify runs with missing bad tower maps (without removing them).
# 5. Generate DST lists for the selected runs.
# 6. Summarize the statistics at each stage.
########################################

echo "========================================"
echo "Starting the Golden Run List Generation"
echo "========================================"

# Create necessary directories
echo "Creating necessary directories..."
mkdir -p FileLists/ dst_list/ list/
echo "Directories created."
echo "----------------------------------------"

# Set the working directory
workplace=$(pwd)
echo "Working directory is set to $workplace"
echo "----------------------------------------"

# Step 1: Run the Python script to get the initial run numbers
echo "Step 1: Running Python script to generate initial run lists..."

python_output=$(python3 << EOF
import pyodbc

def get_all_run_numbers(cursor):
    # Query to get run numbers with at least 1 million events and runnumber >= 47289
    query = """
    SELECT runnumber
    FROM datasets
    WHERE dsttype='DST_CALO_run2pp' and dataset = 'ana446_2024p007'
    GROUP BY runnumber
    HAVING SUM(events) >= 1000000 AND runnumber >= 47289;
    """
    cursor.execute(query)
    run_numbers = [row.runnumber for row in cursor.fetchall()]
    return run_numbers

def get_golden_run_numbers(detector, file_catalog_run_numbers, production_cursor):
    query = f"""
    SELECT runnumber
    FROM goodruns
    WHERE ({detector}_auto).runclass = 'GOLDEN'
    """
    production_cursor.execute(query)
    golden_run_numbers = {row.runnumber for row in production_cursor.fetchall()}
    return golden_run_numbers.intersection(set(file_catalog_run_numbers))

def main():
    # Connect to the FileCatalog database
    file_catalog_conn = pyodbc.connect("DSN=FileCatalog;UID=phnxrc;READONLY=True")
    file_catalog_cursor = file_catalog_conn.cursor()

    # Get unique run numbers with at least 1 million total events
    file_catalog_run_numbers = get_all_run_numbers(file_catalog_cursor)
    file_catalog_run_numbers.sort()
    with open('list/list_runnumber_all.txt', 'w') as f:
        for run_number in file_catalog_run_numbers:
            f.write(f"{run_number}\n")
    print(f"TOTAL_RUNS:{len(file_catalog_run_numbers)}")

    # Close the FileCatalog database connection
    file_catalog_conn.close()

    # Connect to the Production database
    production_conn = pyodbc.connect("DSN=Production_write")
    production_cursor = production_conn.cursor()

    # Filter 'GOLDEN' run numbers for each calorimeter
    detectors = ['emcal', 'ihcal', 'ohcal']
    golden_run_numbers = set(file_catalog_run_numbers)
    for detector in detectors:
        detector_golden_runs = get_golden_run_numbers(detector, file_catalog_run_numbers, production_cursor)
        golden_run_numbers = golden_run_numbers.intersection(detector_golden_runs)

    golden_run_numbers = sorted(golden_run_numbers)
    with open('list/Full_ppGoldenRunList.txt', 'w') as f:
        for run_number in golden_run_numbers:
            f.write(f"{run_number}\n")
    print(f"COMBINED_GOLDEN_RUNS:{len(golden_run_numbers)}")

    # Close the Production database connection
    production_conn.close()

if __name__ == "__main__":
    main()
EOF
)

echo "Python script execution completed."
echo "----------------------------------------"

# Parse the output to get the counts
total_runs=$(echo "$python_output" | grep 'TOTAL_RUNS' | cut -d':' -f2)
combined_golden_runs=$(echo "$python_output" | grep 'COMBINED_GOLDEN_RUNS' | cut -d':' -f2)

echo "Summary of initial counts:"
echo "Total runs after firmware fix and >1M events: $total_runs"
echo "Number of runs passing Calo QA (all three GOLDEN statuses): $combined_golden_runs"
echo "----------------------------------------"

# Step 2: Check that the initial golden run list exists
echo "Step 2: Checking for the initial golden run list..."

golden_run_list="list/Full_ppGoldenRunList.txt"
if [[ ! -f "$golden_run_list" ]]; then
    echo "[ERROR] The file $golden_run_list does not exist. Please check the path and try again."
    exit 1
fi

echo "Initial golden run list found: $golden_run_list"
echo "----------------------------------------"

# Old method function (incorrect) using SUM(raw)
# --------------------------------------------------
# This function aggregates event counts by summing the 'raw' column from the gl1_scalers table.
# We now know that this is not the correct measure of fully recorded events—it's just trigger-level counts.
# However, I keep this function to compare the old approach to the new correct one.
#
# Procedure:
# 1. Reads a list of runnumbers from a file.
# 2. Processes them in batches of 100 for efficiency.
# 3. For each batch, runs a query to sum up raw trigger counts (gl1_scalers).
# 4. Accumulates these into 'total_events'.
# 5. Outputs 'total_events' at the end (old, incorrect count).
get_total_events() {
    input_file=$1
    total_events=0
    batch_size=100
    run_numbers=()

    while IFS= read -r runnumber; do
        if [[ -n "$runnumber" ]]; then
            run_numbers+=("$runnumber")

            if [[ ${#run_numbers[@]} -ge $batch_size ]]; then
                run_list=$(printf ",%s" "${run_numbers[@]}")
                run_list=${run_list:1}

                query="SELECT SUM(raw) FROM gl1_scalers WHERE runnumber IN (${run_list});"
                result=$(psql -h sphnxdaqdbreplica -d daq -t -c "$query")
                events=$(echo "$result" | xargs)

                if [[ "$events" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
                    total_events=$(echo "$total_events + $events" | bc)
                fi

                run_numbers=()
            fi
        fi
    done < "$input_file"

    # Process leftover runs if less than batch_size
    if [[ ${#run_numbers[@]} -gt 0 ]]; then
        run_list=$(printf ",%s" "${run_numbers[@]}")
        run_list=${run_list:1}
        query="SELECT SUM(raw) FROM gl1_scalers WHERE runnumber IN (${run_list});"
        result=$(psql -h sphnxdaqdbreplica -d daq -t -c "$query")
        events=$(echo "$result" | xargs)
        if [[ "$events" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
            total_events=$(echo "$total_events + $events" | bc)
        fi
    fi

    echo "$total_events"
}

# New method function (correct) using filelist evt files: (lastevent - firstevent + 1)
# ------------------------------------------------------------------------------------
# This function calculates the actual number of fully recorded events by examining the raw GL1 .evt files.
# Instead of trigger-level counts, it uses 'filelist' (lastevent and firstevent) for each .evt segment.
#
# Summation of (lastevent - firstevent + 1) over all .evt segments for a run gives the true event count.
# It handles multiple segments automatically because we sum over all .evt files that match the pattern.
get_actual_events_from_evt() {
    input_file=$1
    total_events=0
    batch_size=100
    run_numbers=()

    while IFS= read -r runnumber; do
        if [[ -n "$runnumber" ]]; then
            run_numbers+=("$runnumber")
            if [[ ${#run_numbers[@]} -ge $batch_size ]]; then
                run_list=$(printf ",%s" "${run_numbers[@]}")
                run_list=${run_list:1}

                # Query to sum up actual recorded events from .evt files
                query="SELECT SUM(lastevent - firstevent + 1) FROM filelist WHERE runnumber IN (${run_list}) AND filename LIKE '%GL1_physics_gl1daq%.evt';"
                result=$(psql -h sphnxdaqdbreplica -d daq -t -c "$query")
                events=$(echo "$result" | xargs)

                if [[ "$events" =~ ^[0-9]+$ ]]; then
                    total_events=$(echo "$total_events + $events" | bc)
                fi
                run_numbers=()
            fi
        fi
    done < "$input_file"

    # Process leftover runs if less than batch_size
    if [[ ${#run_numbers[@]} -gt 0 ]]; then
        run_list=$(printf ",%s" "${run_numbers[@]}")
        run_list=${run_list:1}
        query="SELECT SUM(lastevent - firstevent + 1) FROM filelist WHERE runnumber IN (${run_list}) AND filename LIKE '%GL1_physics_gl1daq%.evt';"
        result=$(psql -h sphnxdaqdbreplica -d daq -t -c "$query")
        events=$(echo "$result" | xargs)
        if [[ "$events" =~ ^[0-9]+$ ]]; then
            total_events=$(echo "$total_events + $events" | bc)
        fi
    fi

    echo "$total_events"
}

#################################
# Version 1: Calo QA + run time > 5 mins + livetime > 80%
#################################

echo "----------------------------------------"
echo "Processing Version 1: Calo QA + run time > 5 mins + livetime > 80%"
echo "----------------------------------------"

# Step 3a: Apply run duration filter (>300 seconds)
input_file="list/Full_ppGoldenRunList.txt"  # Runs after Calo QA
output_file_duration_v1="list/list_runnumber_runtime_v1.txt"
> "$output_file_duration_v1"

total_runs_duration_v1=0
runs_dropped_runtime_v1=0

while IFS= read -r runnumber; do
    if [[ -z "$runnumber" ]]; then
        continue
    fi
    query="SELECT runnumber, EXTRACT(EPOCH FROM (ertimestamp - brtimestamp)) AS duration FROM run WHERE runnumber = ${runnumber};"
    result=$(psql -h sphnxdaqdbreplica -d daq -t -c "$query" | tr -d '[:space:]')
    duration=$(echo "$result" | awk -F '|' '{print $2}')

    if [[ $duration =~ ^[0-9]+([.][0-9]+)?$ ]]; then
        if (( $(echo "$duration > 300" | bc -l) )); then
            echo "$runnumber" >> "$output_file_duration_v1"
            total_runs_duration_v1=$((total_runs_duration_v1+1))
        else
            runs_dropped_runtime_v1=$((runs_dropped_runtime_v1+1))
        fi
    else
        runs_dropped_runtime_v1=$((runs_dropped_runtime_v1+1))
    fi
done < "$input_file"

echo "Total runs after run duration cut (>5 mins): $total_runs_duration_v1"
echo "Number of runs dropped due to run duration cut: $runs_dropped_runtime_v1"
echo "----------------------------------------"

# Step 4a: Apply live/raw ratio filter for trigger index 10 (>80%)
input_file="$output_file_duration_v1"
output_file_livetime_v1="list/list_runnumber_livetime_v1.txt"
bad_file_livetime_v1="list/list_runnumber_bad_livetime_v1.txt"
> "$output_file_livetime_v1"
> "$bad_file_livetime_v1"

total_runs_livetime_v1=0
runs_dropped_livetime_v1=0

while IFS= read -r runnumber; do
    if [[ -z "$runnumber" ]]; then
        continue
    fi
    index_to_check=10
    query="SELECT index, raw, live FROM gl1_scalers WHERE runnumber = ${runnumber} AND index = ${index_to_check};"
    result=$(psql -h sphnxdaqdbreplica -d daq -t -c "$query")

    index_pass=false
    while IFS='|' read -r index raw live; do
        index=$(echo "$index" | xargs)
        raw=$(echo "$raw" | xargs)
        live=$(echo "$live" | xargs)

        if [[ "$raw" =~ ^[0-9]+$ ]] && [[ "$live" =~ ^[0-9]+$ ]] && [ "$raw" -ne 0 ]; then
            ratio=$(echo "scale=2; $live / $raw * 100" | bc -l)
            if [ "$index" -eq "$index_to_check" ]; then
                if (( $(echo "$ratio >= 80" | bc -l) )); then
                    index_pass=true
                fi
            fi
        fi
    done <<< "$result"

    if [[ "$index_pass" == true ]]; then
        echo "$runnumber" >> "$output_file_livetime_v1"
        total_runs_livetime_v1=$((total_runs_livetime_v1+1))
    else
        echo "$runnumber" >> "$bad_file_livetime_v1"
        runs_dropped_livetime_v1=$((runs_dropped_livetime_v1+1))
    fi
done < "$input_file"

echo "Total runs after livetime cut (>80% for trigger index 10): $total_runs_livetime_v1"
echo "Number of runs dropped due to livetime cut: $runs_dropped_livetime_v1"
echo "----------------------------------------"

# Step 5a: Identify runs with missing bad tower maps (without removing them)
echo "Identifying runs with missing bad tower maps..."

input_file="$output_file_livetime_v1"
output_file_final_v1="FileLists/Full_ppGoldenRunList_Version1.txt"
bad_tower_runs_file="list/list_runs_missing_bad_tower_maps.txt"

# Copy the input file to the final output file
cp "$input_file" "$output_file_final_v1"

available_bad_tower_runs=$(find /cvmfs/sphenix.sdcc.bnl.gov/calibrations/sphnxpro/cdb/CEMC_BadTowerMap -name "*p0*" | cut -d '-' -f2 | cut -d c -f1 | sort | uniq)
echo "$available_bad_tower_runs" > list/available_bad_tower_runs.txt

grep -Fxvf list/available_bad_tower_runs.txt "$input_file" > "$bad_tower_runs_file"

total_runs_with_bad_tower=$(grep -Fxf list/available_bad_tower_runs.txt "$input_file" | wc -l)
total_runs_missing_bad_tower=$(wc -l < "$bad_tower_runs_file")

echo "Total runs with bad tower maps: $total_runs_with_bad_tower"
echo "Total runs missing bad tower maps: $total_runs_missing_bad_tower"
echo "List of runs missing bad tower maps saved to $bad_tower_runs_file"
echo "----------------------------------------"

rm list/available_bad_tower_runs.txt

cp "$output_file_final_v1" dst_list/Final_RunNumbers_After_All_Cuts.txt
echo "Final run numbers after all cuts have been copied to dst_list/Final_RunNumbers_After_All_Cuts.txt"
echo "----------------------------------------"

# Step 6: Create .list file for CreateDstList.pl command
echo "Creating .list file for CreateDstList.pl..."
cp "$output_file_final_v1" Full_ppGoldenRunList_Version1.list
echo ".list file created."
echo "----------------------------------------"

# Step 7: Remove existing list files in dst_list directory
echo "Removing existing list files in dst_list directory..."
rm -f dst_list/*.list
echo "Existing list files removed."
echo "----------------------------------------"

# Step 8: Create DST lists
echo "Changing to dst_list directory to generate DST lists..."
cd dst_list

echo "Running CreateDstList.pl to generate the DST list..."
CreateDstList.pl --build ana437 --cdb 2024p007 DST_CALO_run2pp --list ../Full_ppGoldenRunList_Version1.list
echo "DST list generated and saved to dst_list."
echo "----------------------------------------"

cd ..

# Calculate total events (old method)
total_events_initial=$(get_total_events 'list/list_runnumber_all.txt')
total_events_calo_qa=$(get_total_events 'list/Full_ppGoldenRunList.txt')
total_events_after_runtime=$(get_total_events 'list/list_runnumber_runtime_v1.txt')
total_events_after_livetime=$(get_total_events 'list/list_runnumber_livetime_v1.txt')
total_events_after_all_cuts=$(get_total_events "$output_file_final_v1")

# Calculate total events (new method)
actual_events_initial=$(get_actual_events_from_evt 'list/list_runnumber_all.txt')
actual_events_calo_qa=$(get_actual_events_from_evt 'list/Full_ppGoldenRunList.txt')
actual_events_after_runtime=$(get_actual_events_from_evt 'list/list_runnumber_runtime_v1.txt')
actual_events_after_livetime=$(get_actual_events_from_evt 'list/list_runnumber_livetime_v1.txt')
actual_events_after_all_cuts=$(get_actual_events_from_evt "$output_file_final_v1")

# Adjust the counts for runs
total_runs_after_all_cuts_v1=$(wc -l < "$output_file_final_v1")

# Compute percentages (old method)
percent_runs_calo_qa=$(echo "scale=2; $combined_golden_runs / $total_runs * 100" | bc)
percent_runs_after_runtime=$(echo "scale=2; $total_runs_duration_v1 / $total_runs * 100" | bc)
percent_runs_after_livetime=$(echo "scale=2; $total_runs_livetime_v1 / $total_runs * 100" | bc)
percent_runs_after_badtower=$(echo "scale=2; $total_runs_after_all_cuts_v1 / $total_runs * 100" | bc)

percent_events_calo_qa=$(echo "scale=2; $total_events_calo_qa / $total_events_initial * 100" | bc)
percent_events_after_runtime=$(echo "scale=2; $total_events_after_runtime / $total_events_initial * 100" | bc)
percent_events_after_livetime=$(echo "scale=2; $total_events_after_livetime / $total_events_initial * 100" | bc)
percent_events_after_all_cuts=$(echo "scale=2; $total_events_after_all_cuts / $total_events_initial * 100" | bc)

# Compute percentages (new method)
percent_actual_events_calo_qa=$(echo "scale=2; $actual_events_calo_qa / $actual_events_initial * 100" | bc)
percent_actual_events_after_runtime=$(echo "scale=2; $actual_events_after_runtime / $actual_events_initial * 100" | bc)
percent_actual_events_after_livetime=$(echo "scale=2; $actual_events_after_livetime / $actual_events_initial * 100" | bc)
percent_actual_events_after_all_cuts=$(echo "scale=2; $actual_events_after_all_cuts / $actual_events_initial * 100" | bc)

echo "========================================"
echo "Final Summary: Version 1 (Calo QA → Runtime → Livetime)"

# Print the table exactly as shown in the provided image:
# Columns:
# 1) Step
# 2) % of Initial Events (old method)
# 3) GL1 Raw Events (Previous approach)
# 4) Event Counts (Using raw GL1 '.evt' files)
# 5) Runs
# 6) % of Initial Runs

printf "%-50s | %-35s | %-35s | %-25s\n" \
"Stage" "GL1 Raw Events (Wrong Approach)" ".evt File Events (New Approach)" "Runs"
echo "--------------------------------------------------|-------------------------------------|-------------------------------------|-------------------------"

# 1) After firmware fix and >1M events
printf "%-50s | %-35s | %-35s | %-25s\n" \
"1) After firmware fix (47289) and >1M events" \
"${total_events_initial} (100%)" \
"${actual_events_initial} (100%)" \
"${total_runs} (100%)"

# 2) && Golden EMCal/HCal
printf "%-50s | %-35s | %-35s | %-25s\n" \
"2) && Golden EMCal/HCal" \
"${total_events_calo_qa} (${percent_events_calo_qa}%)" \
"${actual_events_calo_qa} (${percent_actual_events_calo_qa}%)" \
"${combined_golden_runs} (${percent_runs_calo_qa}%)"

# 3) && > 5 minutes
printf "%-50s | %-35s | %-35s | %-25s\n" \
"3) && > 5 minutes" \
"${total_events_after_runtime} (${percent_events_after_runtime}%)" \
"${actual_events_after_runtime} (${percent_actual_events_after_runtime}%)" \
"${total_runs_duration_v1} (${percent_runs_after_runtime}%)"

# 4) && MB livetime > 80%
printf "%-50s | %-35s | %-35s | %-25s\n" \
"4) && MB livetime > 80%" \
"${total_events_after_livetime} (${percent_events_after_livetime}%)" \
"${actual_events_after_livetime} (${percent_actual_events_after_livetime}%)" \
"${total_runs_livetime_v1} (${percent_runs_after_livetime}%)"

# 5) Final run list
printf "%-50s | %-35s | %-35s | %-25s\n" \
"5) Final run list" \
"${total_events_after_all_cuts} (${percent_events_after_all_cuts}%)" \
"${actual_events_after_all_cuts} (${percent_actual_events_after_all_cuts}%)" \
"${total_runs_after_all_cuts_v1} (${percent_runs_after_badtower}%)"

echo "================================================="
echo ""
echo "Previous Approach: Aggregates event counts by summing the 'raw' column from the gl1_scalers table"
echo "New Approach: Uses (lastevent - firstevent + 1) from the raw GL1 .evt files to get the actual recorded event count"
echo "========================================"

# Copy the missing bad tower maps file to the requested path
cp "$bad_tower_runs_file" /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/list_runs_missing_bad_tower_maps.txt

# Copy the runs that fail livetime cut to the requested path
cp "$bad_file_livetime_v1" /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/list_runnumber_bad_livetime_v1.txt

cp dst_list/Final_RunNumbers_After_All_Cuts.txt /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/FinalGoldenRunList_ana446_2024p007.txt

echo "Files for missing bad tower maps, livetime failures, and the final run list have been copied to:"
echo "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/"
echo "Done."
