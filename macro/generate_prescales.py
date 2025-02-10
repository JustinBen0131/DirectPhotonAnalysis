#!/usr/bin/env python3
import sys
import psycopg2
import os

# -----------------------------------------------------------------------------
# 1) Define the triggers you want to query
#    We'll map DB (long) name -> output filename.
# -----------------------------------------------------------------------------
OUTPUT_DIR = "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/preScaleFiles"

TRIGGERS = {
    "MBD N&S >= 1":               "MBD_NandS_geq_1.txt",
    "Photon 3 GeV + MBD NS >= 1": "Photon_3_GeV_plus_MBD_NS_geq_1.txt",
    "Photon 4 GeV + MBD NS >= 1": "Photon_4_GeV_plus_MBD_NS_geq_1.txt",
    "Photon 5 GeV + MBD NS >= 1": "Photon_5_GeV_plus_MBD_NS_geq_1.txt"
}


def get_prescale(runNumber: int, trigger_db_name: str, verbose: bool=False) -> float:
    """Query the DAQ DB for 'live' and 'scaled' for the given runNumber & trigger,
       return prescale = live/scaled if scaled>0, else -1."""
    if runNumber <= 0:
        if verbose:
            print(f"[get_prescale] Error: runNumber={runNumber} is invalid => returning -1.")
        return -1.0


    try:
        conn = psycopg2.connect(
            host="sphnxdaqdbreplica",
            port=5432,
            dbname="daq",
            user="phnxro",
            password=""
        )
    except Exception as e:
        print(f"[get_prescale] DB connection failed for run={runNumber}, error={e}")
        return -1.0

    # Build query
    query = f"""
    SELECT s.live, s.scaled
    FROM gl1_scalers s
    JOIN gl1_triggernames t
      ON (s.index = t.index
          AND s.runnumber BETWEEN t.runnumber AND t.runnumber_last)
    WHERE s.runnumber={runNumber}
      AND t.triggername='{trigger_db_name}'
    ORDER BY s.index;
    """

    scaleVal = -1.0
    try:
        with conn.cursor() as cur:
            cur.execute(query)
            row = cur.fetchone()
            if row is not None:
                liveVal, scaledVal = row[0], row[1]
                if scaledVal == 0.0:
                    # Trigger is OFF
                    if verbose:
                        print(f"\033[1m[WARNING]\033[0m Trigger OFF (scaled=0) for '{trigger_db_name}' run={runNumber}")
                    scaleVal = -1.0
                elif scaledVal is not None and scaledVal > 0.0:
                    scaleVal = float(liveVal / scaledVal)
                # else remain -1
    except Exception as e:
        print(f"[get_prescale] Query failed for run={runNumber}, trigger='{trigger_db_name}', error={e}")
    finally:
        conn.close()

    return scaleVal

# -----------------------------------------------------------------------------
# 3) Main script: usage => python generate_prescales.py <runlist.txt>
#    We'll read the runlist, then for each trigger, create an output file & write lines.
# -----------------------------------------------------------------------------
def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <runListFile>")
        sys.exit(1)

    runListFile = sys.argv[1]

    # Ensure OUTPUT_DIR exists
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR, exist_ok=True)

    # 3a) Read run numbers from file
    with open(runListFile, 'r') as f:
        runNumbers = [int(line.strip()) for line in f if line.strip()]

    # 3b) Open one output file per trigger for writing
    output_files = {}
    for dbName, outName in TRIGGERS.items():
        # Overwrite existing
        outPath = os.path.join(OUTPUT_DIR, outName)  # e.g. /path/preScaleFiles/MBD_NandS_geq_1.txt
        outfile = open(outPath, 'w')
        # Write a header
        outfile.write("# runNumber  preScale\n")
        output_files[dbName] = outfile

    # Track how many runs had "scaled=0" or no row => prescale=-1 for each trigger
    scaled_zero_counts = {dbName: 0 for dbName in TRIGGERS.keys()}

    # Also track sums of valid prescales & counts => so we can compute average
    prescale_sums   = {dbName: 0.0 for dbName in TRIGGERS.keys()}
    prescale_counts = {dbName: 0   for dbName in TRIGGERS.keys()}

    # 3c) For each run number in the input file, query each trigger => get prescale => write to appropriate file
    for runNo in runNumbers:
        for dbName, outName in TRIGGERS.items():
            prescale = get_prescale(runNo, dbName, verbose=False)
            # If prescale < 0 => interpret as "not found or scaled=0"
            if prescale < 0.0:
                scaled_zero_counts[dbName] += 1
                prescale_str = "-1"
            else:
                # A valid prescale
                prescale_str = f"{prescale:.3f}"
                prescale_sums[dbName]   += prescale
                prescale_counts[dbName] += 1

            # Write line: "runNumber  prescale"
            output_files[dbName].write(f"{runNo}  {prescale_str}\n")

    # 3d) close all output files
    for dbName, fobj in output_files.items():
        fobj.close()

    # 3e) Print summary
    total_runs = len(runNumbers)
    print("========================================")
    print("Prescale summary for each trigger:")
    for dbName in TRIGGERS.keys():
        outName = TRIGGERS[dbName]
        outPath = os.path.join(OUTPUT_DIR, outName)
        zeroCount = scaled_zero_counts[dbName]

        print(f"Trigger '{dbName}': output file => {outPath}")
        print(f"  => {zeroCount} runs had scaled=0 or not found (prescale=-1) out of {total_runs} total runs.")

        # Compute average among valid prescales
        valid_count = prescale_counts[dbName]
        if valid_count > 0:
            avg_prescale = prescale_sums[dbName] / valid_count
            print(f"  => Average valid prescale = {avg_prescale:.3f} over {valid_count} valid runs.")
        else:
            print("  => No valid prescales found for this trigger.")

        print("")

    print("Done.")


if __name__ == "__main__":
    main()
