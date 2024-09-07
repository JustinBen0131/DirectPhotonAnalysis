#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <vector>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

// Define a map that stores the unique trigger mapping to run numbers
std::map<std::map<int, std::string>, std::vector<int>> uniqueTriggerMappings;
std::map<int, std::string> previousTriggerMap;


// Function to get trigger mapping for a specific run
std::map<int, std::string> getTriggerMappingForRun(int run) {
    std::map<int, std::string> triggerMap;

    // Connect to the database
    TSQLServer* db = TSQLServer::Connect("pgsql://sphnxdaqdbreplica:5432/daq", "phnxro", "");
    if (!db) {
        std::cerr << "Error: Failed to connect to the database." << std::endl;
        return triggerMap;
    }

    // Construct a query to get the trigger mapping for a specific run
    std::string query =
        // Generate a list of indices from 0 to 63
        "WITH all_indices AS ( "
        "SELECT generate_series(0, 63) AS index) "
    
        // Select the index and corresponding trigger name, using 'EMPTY' if no trigger name exists
        "SELECT ai.index, COALESCE(gt.triggername, 'EMPTY') AS triggername "
    
        // Join the generated indices with the gl1_triggernames table to find matching trigger names
        "FROM all_indices ai "
        "LEFT JOIN gl1_triggernames gt "
        "ON ai.index = gt.index "
    
        // Filter the results to get trigger names valid for the given run number
        "AND gt.runnumber <= " + std::to_string(run) + " "
        "AND gt.runnumber_last >= " + std::to_string(run) + " "
    
        // Sort the results by index
        "ORDER BY ai.index;";

    TSQLResult* res = db->Query(query.c_str());
    if (!res) {
        std::cerr << "Error: Failed to query trigger mapping for run " << run << std::endl;
        delete db;
        return triggerMap;
    }

    // Store the trigger mapping in the map
    TSQLRow* row;
    while ((row = res->Next())) {
        int index = std::stoi(row->GetField(0));
        std::string triggername = row->GetField(1);
        triggerMap[index] = triggername;  // Store in the run-specific trigger map
        delete row;
    }

    // Clean up
    delete res;
    db->Close();
    delete db;

    return triggerMap;
}

// Function to check for changes and update unique trigger mappings
void checkAndUpdateUniqueMappings(int run, const std::map<int, std::string>& triggerMap) {
    bool mappingExists = false;

    // Check if the current triggerMap matches any existing unique mappings
    for (auto& [uniqueMap, runList] : uniqueTriggerMappings) {
        if (uniqueMap == triggerMap) {
            // If the mapping already exists, add the run number to the corresponding run list
            runList.push_back(run);
            mappingExists = true;
            break;
        }
    }

    if (!mappingExists) {
        // If no matching mapping was found, create a new entry in the uniqueTriggerMappings
        uniqueTriggerMappings[triggerMap].push_back(run);
    }
}

// Function to print the saved mappings
void printUniqueMappings() {
    std::ofstream file1("runListTriggerLUTv1.txt");
    std::ofstream file2("runListTriggerLUTv2.txt");

    std::cout << "\nUnique Trigger Mappings:\n";
    int groupNumber = 1;
    for (const auto& [uniqueMap, runList] : uniqueTriggerMappings) {
        std::cout << "\nGroup " << groupNumber << ":\n";
        for (const auto& [index, triggername] : uniqueMap) {
            std::cout << "Index " << index << ": " << triggername << std::endl;
        }

        std::cout << "Runs: ";
        for (int run : runList) {
            std::cout << run << " ";
        }
        std::cout << "\n----------------------------------------\n";

        // Write the run list to the appropriate file
        std::ofstream* outFile = (groupNumber == 1) ? &file1 : &file2;
        for (int run : runList) {
            *outFile << run << "\n";
        }
        groupNumber++;
    }

    file1.close();
    file2.close();

    std::cout << "\nTrigger Map 1:\nIndex\tTrigger Name\n";
    for (const auto& [index, triggername] : uniqueTriggerMappings.begin()->first) {
        std::cout << index << "\t" << triggername << "\n";
    }
    std::cout << "\n" << uniqueTriggerMappings.begin()->second.size()
              << " run numbers have been added to 'runListTriggerLUTv1.txt' with a unique LUT.\n";

    std::cout << "\nTrigger Map 2:\nIndex\tTrigger Name\n";
    auto it = uniqueTriggerMappings.begin();
    std::advance(it, 1);
    for (const auto& [index, triggername] : it->first) {
        std::cout << index << "\t" << triggername << "\n";
    }
    std::cout << "\n" << it->second.size()
              << " run numbers have been added to 'runListTriggerLUTv2.txt' with a unique LUT.\n";
}

// Function to detect trigger index switches between runs
void detectTriggerSwitches(int run, const std::map<int, std::string>& currentTriggerMap) {
    if (!previousTriggerMap.empty()) {
        for (const auto& [index, triggername] : currentTriggerMap) {
            if (previousTriggerMap[index] != triggername) {
                std::cout << "TRIGGER INDICES SWITCHED AT RUN " << run << "\n";
                break;
            }
        }
    }
    previousTriggerMap = currentTriggerMap; // Update for next run
}


// Function to read run numbers from a file
std::vector<int> readRunNumbersFromFile(const std::string& filename) {
    std::vector<int> runNumbers;
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error: Failed to open file " << filename << std::endl;
        return runNumbers;
    }

    int runNumber;
    while (infile >> runNumber) {
        runNumbers.push_back(runNumber);
    }

    infile.close();
    return runNumbers;
}

// Main function for ROOT macro
void QueryGl1() {
    // Read run numbers from the file
    std::vector<int> runNumbers = readRunNumbersFromFile("../GoldenRunNumbers_afterRun44452.txt");

    // Process each run number from the file
    for (int run : runNumbers) {
        std::map<int, std::string> triggerMap = getTriggerMappingForRun(run);
        detectTriggerSwitches(run, triggerMap);  // Detect trigger index switches
        checkAndUpdateUniqueMappings(run, triggerMap);
    }

    // Print the saved mappings after processing all runs
    printUniqueMappings();
}
