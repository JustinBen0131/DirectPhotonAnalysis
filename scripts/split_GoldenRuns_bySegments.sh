#!/bin/bash

# Set paths
input_file="/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/GoldenRunNumbers_afterRun46619.txt"
dst_list_dir="/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/dst_list"
output_dir="/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp"
max_segments=10000  # The maximum number of segments per output file

# Initialize variables
current_segments=0
file_counter=1
current_output_file="${output_dir}/GoldenRunNumbers_afterRun46619_part${file_counter}.txt"

# Create the first output file
echo "Creating file: ${current_output_file}"
touch "${current_output_file}"

# Loop over each run number in the input file
while read -r runNumber; do
    # Define the path to the run's list file
    list_file="${dst_list_dir}/dst_calo_run2pp-000${runNumber}.list"

    # Count the number of segments (lines) in the list file
    if [ -f "${list_file}" ]; then
        segment_count=$(wc -l < "${list_file}")

        # Check if adding this run number will exceed the max segments
        if (( current_segments + segment_count > max_segments )); then
            # Start a new output file
            ((file_counter++))
            current_output_file="${output_dir}/GoldenRunNumbers_afterRun46619_part${file_counter}.txt"
            echo "Creating new file: ${current_output_file}"
            touch "${current_output_file}"

            # Reset segment count for the new file
            current_segments=0
        fi

        # Add the run number to the current output file
        echo "${runNumber}" >> "${current_output_file}"

        # Update the current segment count
        ((current_segments += segment_count))
    else
        echo "Warning: List file for run number ${runNumber} not found!"
    fi
done < "${input_file}"

echo "Processing complete."
