#!/bin/bash
mkdir -p ../dst_list/
cd ../dst_list/
segment=0
for i in $(cat ../list/list_runnumber_golden.txt); do
    > dst_calo_run2pp-000${i}.list
    > dst_calofitting_run2pp-000${i}.list
    > dst_triggered_event_run2pp-000${i}.list
    low_range=$(( (i / 100)*100 ))
    high_range=$(( low_range + 100 ))
    ls /sphenix/lustre01/sphnxpro/physics/slurp/caloy2fitting/ana430_2024p007/run_000${low_range}_000${high_range}/DST_CALO_run2pp_ana430_2024p007-000${i}* >> dst_calo_run2pp-000${i}.list 2>/dev/null
    ls /sphenix/lustre01/sphnxpro/physics/slurp/caloy2calib/ana430_2024p007/run_000${low_range}_000${high_range}/DST_CALO_run2pp_ana430_2024p007-000${i}* >> dst_calo_run2pp-000${i}.list 2>/dev/null
    ls /sphenix/lustre01/sphnxpro/physics/slurp/calophysics/ana430_2024p007/run_000${low_range}_000${high_range}/DST_CALO_run2pp_ana430_2024p007-000${i}* >> dst_calo_run2pp-000${i}.list 2>/dev/null
    while IFS= read -r line; do
        end=$(echo "$line" | sed 's/.*ana430_2024p007//')
        fittingaddress1=$(ls /sphenix/lustre01/sphnxpro/physics/slurp/caloy2calib/ana430_2024p007/run_000${low_range}_000${high_range}/DST_CALOFITTING_run2pp_ana430_2024p007$end 2>/dev/null)
        fittingaddress2=$(ls /sphenix/lustre01/sphnxpro/physics/slurp/caloy2fitting/ana430_2024p007/run_000${low_range}_000${high_range}/DST_CALOFITTING_run2pp_ana430_2024p007$end 2>/dev/null)
        fittingaddress3=$(ls /sphenix/lustre01/sphnxpro/physics/slurp/calophysics/ana430_2024p007/run_000${low_range}_000${high_range}/DST_CALOFITTING_run2pp_ana430_2024p007$end 2>/dev/null)
        if [ -n "$fittingaddress1" ]; then
            echo "$fittingaddress1" >> dst_calofitting_run2pp-000${i}.list
        elif [ -n "$fittingaddress2" ]; then
            echo "$fittingaddress2" >> dst_calofitting_run2pp-000${i}.list
        elif [ -n "$fittingaddress3" ]; then
            echo "$fittingaddress3" >> dst_calofitting_run2pp-000${i}.list
        else
            echo "No file found for $end"
        fi
        triggeredaddress1=$(ls /sphenix/lustre01/sphnxpro/physics/slurp/caloy2calib/ana430_2024p007/run_000${low_range}_000${high_range}/DST_TRIGGERED_EVENT_run2pp_ana430_2024p007$end 2>/dev/null)
        triggeredaddress2=$(ls /sphenix/lustre01/sphnxpro/physics/slurp/caloy2fitting/ana430_2024p007/run_000${low_range}_000${high_range}/DST_TRIGGERED_EVENT_run2pp_ana430_2024p007$end 2>/dev/null)
        triggeredaddress3=$(ls /sphenix/lustre01/sphnxpro/physics/slurp/calophysics/ana430_2024p007/run_000${low_range}_000${high_range}/DST_TRIGGERED_EVENT_run2pp_ana430_2024p007$end 2>/dev/null)
        if [ -n "$triggeredaddress1" ]; then
            echo "$triggeredaddress1" >> dst_triggered_event_run2pp-000${i}.list
        elif [ -n "$triggeredaddress2" ]; then
            echo "$triggeredaddress2" >> dst_triggered_event_run2pp-000${i}.list
        elif [ -n "$triggeredaddress3" ]; then
            echo "$triggeredaddress3" >> dst_triggered_event_run2pp-000${i}.list
        else
            echo "No file found for $end"
        fi
    done < dst_calo_run2pp-000${i}.list
    line_count=$(wc -l < dst_calo_run2pp-000${i}.list)
    if [ $line_count -lt 100 ]; then
        echo "Run $i has less than 100 files"
    fi
    segment=$((segment + line_count))
done
echo "All runnumbers processed. Results saved in dst_list directory"
echo "Total number of segments: $segment
