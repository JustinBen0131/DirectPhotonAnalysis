#!/usr/bin/env python3

import pyodbc

def get_all_run_numbers(cursor):
    query = """
    SELECT runnumber
    FROM datasets
    WHERE filename LIKE 'DST_CALO_run2pp_ana437_2024p007-%' AND dataset = 'ana437_2024p007'
    GROUP BY runnumber
    HAVING SUM(events) >= 1000000;
    """
    cursor.execute(query)
    run_numbers = [row.runnumber for row in cursor.fetchall()]
    return run_numbers


def get_all_run_numbers_no_event_count(cursor):
    query = """
    SELECT DISTINCT runnumber
    FROM datasets
    WHERE runnumber > 46619;
    """
    cursor.execute(query
    all_runs_no_event_count = [row.runnumber for row in cursor.fetchall()]
    return all_runs_no_event_count



def get_good_run_numbers(production_cursor):
    query = """
    SELECT runnumber
    FROM goodruns
    """
    production_cursor.execute(query)
    good_run_numbers = {row.runnumber for row in production_cursor.fetchall()}
    return list(good_run_numbers)

def get_emcal_auto_golden_run_numbers(file_catalog_run_numbers, production_cursor):
    query = """
    SELECT runnumber
    FROM goodruns
    WHERE (emcal_auto).runclass = 'GOLDEN'
    """
    production_cursor.execute(query)
    emcal_auto_golden_run_numbers = {row.runnumber for row in production_cursor.fetchall()}
    return list(
        emcal_auto_golden_run_numbers.intersection(
            set(file_catalog_run_numbers)
        )
    )

def get_ihcal_auto_golden_run_numbers(file_catalog_run_numbers, production_cursor):
    query = """
    SELECT runnumber
    FROM goodruns
    WHERE (ihcal_auto).runclass = 'GOLDEN'
    """
    production_cursor.execute(query)
    ihcal_auto_golden_run_numbers = {row.runnumber for row in production_cursor.fetchall()}
    return list(
        ihcal_auto_golden_run_numbers.intersection(
            set(file_catalog_run_numbers)
        )
    )

def get_ohcal_auto_golden_run_numbers(file_catalog_run_numbers, production_cursor):
    query = """
    SELECT runnumber
    FROM goodruns
    WHERE (ohcal_auto).runclass = 'GOLDEN'
    """
    production_cursor.execute(query)
    ohcal_auto_golden_run_numbers = {row.runnumber for row in production_cursor.fetchall()}
    return list(
        ohcal_auto_golden_run_numbers.intersection(
            set(file_catalog_run_numbers)
        )
    )

def main():
    # Connect to the FileCatalog database
    file_catalog_conn = pyodbc.connect("DSN=FileCatalog;UID=phnxrc;READONLY=True")
    file_catalog_cursor = file_catalog_conn.cursor()

    all_runs_no_event_count = get_all_run_numbers_no_event_count(file_catalog_cursor)
    print(f"Total number of runs after 46619 (regardless of event count): {len(all_runs_no_event_count)}")


    # Get unique run numbers with at least 1 million total events
    file_catalog_run_numbers = get_all_run_numbers(file_catalog_cursor)
    file_catalog_run_numbers.sort()
    with open('FileLists/runNumbers_with_AtLeast_OneMilEventsAfterRun46619.txt', 'w') as f:
        for run_number in file_catalog_run_numbers:
            f.write(f"{run_number}\n")
    print(f"Number of all runs saved to list_runnumber_all.txt that have at least 1 million events: {len(file_catalog_run_numbers)}")

    # Close the FileCatalog database connection
    file_catalog_conn.close()

    # Connect to the Production database
    production_conn = pyodbc.connect("DSN=Production_write")
    production_cursor = production_conn.cursor()

    # Filter good run numbers
    good_run_numbers = get_good_run_numbers(production_cursor)
    runs_not_in_goodruns = set(file_catalog_run_numbers) - set(good_run_numbers)
    print(f"Number of runs not in the goodruns table: {len(runs_not_in_goodruns)}")

    # Filter 'GOLDEN' run numbers
    emcal_auto_golden_run_numbers = get_emcal_auto_golden_run_numbers(file_catalog_run_numbers, production_cursor)
    emcal_auto_golden_run_numbers.sort()
    with open('FileLists/GoldenEmcalRunList.txt', 'w') as f:
        for run_number in emcal_auto_golden_run_numbers:
            f.write(f"{run_number}\n")
    print(f"Number of EMCal GOLDEN runs saved to GoldenEmcalRunList.txt: {len(emcal_auto_golden_run_numbers)}")

    ihcal_auto_golden_run_numbers = get_ihcal_auto_golden_run_numbers(file_catalog_run_numbers, production_cursor)
    ihcal_auto_golden_run_numbers.sort()
    with open('FileLists/GoldenIHCalRunList.txt', 'w') as f:
        for run_number in ihcal_auto_golden_run_numbers:
            f.write(f"{run_number}\n")
    print(f"Number of IHCal GOLDEN runs saved to GoldenIHCalRunList.txt: {len(ihcal_auto_golden_run_numbers)}")

    ohcal_auto_golden_run_numbers = get_ohcal_auto_golden_run_numbers(file_catalog_run_numbers, production_cursor)
    ohcal_auto_golden_run_numbers.sort()
    with open('FileLists/GoldenOHCalRunList.txt', 'w') as f:
        for run_number in ohcal_auto_golden_run_numbers:
            f.write(f"{run_number}\n")
    print(f"Number of OHCal GOLDEN runs saved to GoldenOHCalRunList.txt: {len(ohcal_auto_golden_run_numbers)}")

    # Get the intersection of all sets
    golden_run_numbers = list(
        set(emcal_auto_golden_run_numbers).intersection(
            set(ihcal_auto_golden_run_numbers),
            set(ohcal_auto_golden_run_numbers),
            set(file_catalog_run_numbers)
        )
    )
    golden_run_numbers.sort()
    with open('FileLists/GoldenCalorimeterRunList.txt', 'w') as f:
        for run_number in golden_run_numbers:
            f.write(f"{run_number}\n")
    print(f"Number of Calo GOLDEN runs saved to list_runnumber_calo.txt: {len(golden_run_numbers)}")

    # Close the Production database connection
    production_conn.close()

if __name__ == "__main__":
    main()
