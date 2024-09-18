#!/bin/bash
echo "Creating dst_list directory..."
mkdir -p ../dst_list/
cd ../dst_list/
segment=0
echo "Processing GoldenRunNumbers_afterRun46619.txt..."

for i in $(cat ../GoldenRunNumbers_afterRun46619.txt); do
    echo "Processing run number: $i"
    > dst_calo_run2pp-000${i}.list
    low_range=$(( (i / 100)*100 ))
    high_range=$(( low_range + 100 ))

    echo "Looking for DST_CALO_run2pp files in range ${low_range}-${high_range}..."

    # Collecting DST_CALO_run2pp files only
    ls /sphenix/lustre01/sphnxpro/physics/slurp/caloy2fitting/ana430_2024p007/run_000${low_range}_000${high_range}/DST_CALO_run2pp_ana430_2024p007-000${i}* >> dst_calo_run2pp-000${i}.list 2>/dev/null
    ls /sphenix/lustre01/sphnxpro/physics/slurp/caloy2calib/ana430_2024p007/run_000${low_range}_000${high_range}/DST_CALO_run2pp_ana430_2024p007-000${i}* >> dst_calo_run2pp-000${i}.list 2>/dev/null
    ls /sphenix/lustre01/sphnxpro/physics/slurp/calophysics/ana430_2024p007/run_000${low_range}_000${high_range}/DST_CALO_run2pp_ana430_2024p007-000${i}* >> dst_calo_run2pp-000${i}.list 2>/dev/null

    line_count=$(wc -l < dst_calo_run2pp-000${i}.list)
    echo "Run $i has $line_count DST_CALO files."
    if [ $line_count -lt 100 ]; then
        echo "Warning: Run $i has less than 100 files"
    fi
    segment=$((segment + line_count))
done

echo "All run numbers processed. Results saved in dst_list directory."
echo "Total number of segments: $segment"
