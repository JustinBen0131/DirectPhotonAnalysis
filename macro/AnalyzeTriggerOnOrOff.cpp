#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <fstream>
#include <algorithm>

// ROOT headers
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

// Define the two mappings from trigger indices to names
std::unordered_map<int, std::string> triggerNameMap1 = {
    {0, "Clock"},
    {1, "ZDC_South"},
    {2, "ZDC_North"},
    {3, "ZDC_Coincidence"},
    {4, "HCAL_Singles"},
    {5, "HCAL_Coincidence"},
    {8, "MBD_S_geq_1"},
    {9, "MBD_N_geq_1"},
    {10, "MBD_NandS_geq_1"},
    {11, "MBD_NandS_geq_2"},
    {12, "MBD_NandS_geq_1_vtx_lessTehn_10_cm"},
    {13, "MBD_NandS_geq_1_vtx_lessThen_30_cm"},
    {14, "MBD_NandS_geq_1_vtx_lessThen_60_cm"},
    {15, "HCAL_Singles_plus_MBD_NS_geq_1"},
    {16, "Jet_4_GeV_plus_MBD_NS_geq_1"},
    {17, "Jet_6_GeV_plus_MBD_NS_geq_1"},
    {18, "Jet_8_GeV_plus_MBD_NS_geq_1"},
    {19, "Jet_10_GeV_plus_MBD_NS_geq_1"},
    {20, "Jet_4_GeV"},
    {21, "Jet_6_GeV"},
    {22, "Jet_8_GeV"},
    {23, "Jet_10_GeV"},
    {24, "Photon_1_GeV_plus_MBD_NS_geq_1"},
    {25, "Photon_2_GeV_plus_MBD_NS_geq_1"},
    {26, "Photon_3_GeV_plus_MBD_NS_geq_1"},
    {27, "Photon_4_GeV_plus_MBD_NS_geq_1"},
    {28, "Photon_1_GeV"},
    {29, "Photon_2_GeV"},
    {30, "Photon_3_GeV"},
    {31, "Photon_4_GeV"}
};

std::unordered_map<int, std::string> triggerNameMap2 = {
    {0, "Clock"},
    {1, "ZDC_South"},
    {2, "ZDC_North"},
    {3, "ZDC_Coincidence"},
    {4, "HCAL_Singles"},
    {5, "HCAL_Coincidence"},
    {8, "MBD_S_geq_1"},
    {9, "MBD_N_geq_1"},
    {10, "MBD_NandS_geq_1"},
    {11, "MBD_NandS_geq_2"},
    {12, "MBD_NandS_geq_1_vtx_lessThen_10_cm"},
    {13, "MBD_NandS_geq_1_vtx_lessThen_30_cm"},
    {14, "MBD_NandS_geq_1_vtx_lessThen_60_cm"},
    {15, "HCAL_Singles_plus_MBD_NS_geq_1"},
    {16, "Jet_6_GeV_plus_MBD_NS_geq_1"},
    {17, "Jet_8_GeV_plus_MBD_NS_geq_1"},
    {18, "Jet_10_GeV_plus_MBD_NS_geq_1"},
    {19, "Jet_12_GeV_plus_MBD_NS_geq_1"},
    {20, "Jet_6_GeV"},
    {21, "Jet_8_GeV"},
    {22, "Jet_10_GeV"},
    {23, "Jet_12_GeV"},
    {24, "Photon_2_GeV_plus_MBD_NS_geq_1"},
    {25, "Photon_3_GeV_plus_MBD_NS_geq_1"},
    {26, "Photon_4_GeV_plus_MBD_NS_geq_1"},
    {27, "Photon_5_GeV_plus_MBD_NS_geq_1"},
    {28, "Photon_2_GeV"},
    {29, "Photon_3_GeV"},
    {30, "Photon_4_GeV"},
    {31, "Photon_5_GeV"}
};


// Function to load run numbers from a text file into a vector
std::vector<int> loadRunNumbersFromFile(const std::string& filename) {
    std::vector<int> runNumbers;
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        try {
            int runNumber = std::stoi(line);
            runNumbers.push_back(runNumber);
        } catch (const std::invalid_argument&) {
            std::cerr << "Invalid run number in file: " << line << std::endl;
        }
    }
    file.close();
    return runNumbers;
}

// Function to get scale-down factors for a run
bool get_scaledowns(int runnumber, int scaledowns[]) {
    // Connect to the PostgreSQL server
    TSQLServer *db = TSQLServer::Connect("pgsql://sphnxdaqdbreplica:5432/daq", "phnxro", "");

    // Check if the connection was successful
    if (!db) {
        std::cerr << "Failed to connect to the database" << std::endl;
        return false; // Return false if the connection fails
    }

    // Variables to store query results
    TSQLRow *row;
    TSQLResult *res;
    char sql[1000]; // Buffer to hold SQL query strings

    // Loop over each of the 64 scaledown values
    for (int is = 0; is < 64; is++) {
        // Format the SQL query to retrieve the scaledown value for the given run number and bit position
        sprintf(sql, "SELECT scaledown%02d FROM gl1_scaledown WHERE runnumber = %d;", is, runnumber);

        // Execute the query
        res = db->Query(sql);

        // Check if the result is valid
        if (!res) {
            std::cerr << "Failed to execute query: " << sql << std::endl;
            delete db; // Close the database connection
            return false;
        }

        // Check the number of rows returned by the query
        int nrows = res->GetRowCount();
        if (nrows == 0) {
            std::cerr << "No scaledown data found for run number: " << runnumber << std::endl;
            delete res;
            delete db;
            return false;
        }

        // Get the first row (should only be one row)
        row = res->Next();
        if (!row) {
            std::cerr << "Failed to get row data for run number: " << runnumber << std::endl;
            delete res;
            delete db;
            return false;
        }

        // Get the scaledown value
        scaledowns[is] = std::stoi(row->GetField(0));

        // Clean up
        delete row;
        delete res;
    }

    delete db; // Close the database connection
    return true;
}

void analyzeTriggersAndWriteCSV(const std::vector<int>& runNumbers, const std::vector<int>& triggerIndices, const std::unordered_map<int, std::string>& triggerNameMap, const std::string& outputFilePath) {
    std::ofstream csvFile(outputFilePath);
    csvFile << "runNumber";

    for (int triggerIndex : triggerIndices) {
        csvFile << "," << triggerNameMap.at(triggerIndex);
    }
    csvFile << "\n";

    // Data structure to hold run numbers for each trigger bit
    struct TriggerRunData {
        std::vector<int> runsOn;
        std::vector<int> runsOff;
    };

    // Map from trigger index to TriggerRunData
    std::unordered_map<int, TriggerRunData> triggerDataMap;

    // Initialize the map with trigger indices
    for (int triggerIndex : triggerIndices) {
        triggerDataMap[triggerIndex] = TriggerRunData();
    }

    // Loop over each run number
    for (int runNumber : runNumbers) {
        int scaledowns[64] = {0};
        // Get scaledown factors for the run
        if (!get_scaledowns(runNumber, scaledowns)) {
            std::cerr << "Skipping run number: " << runNumber << " due to error in retrieving scaledowns." << std::endl;
            continue;
        }
        
        csvFile << runNumber;

        // Loop over trigger bits of interest
        for (int triggerIndex : triggerIndices) {
            int scaledown = scaledowns[triggerIndex];

            if (scaledown == -1) {
                // Trigger is OFF
                triggerDataMap[triggerIndex].runsOff.push_back(runNumber);
                csvFile << ",OFF";
            } else {
                // Trigger is ON
                triggerDataMap[triggerIndex].runsOn.push_back(runNumber);
                csvFile << ",ON";
            }
        }
        csvFile << "\n";
    }
    csvFile.close();
}

void summarizeTriggerCountsFromCSV(const std::string& csvFilePath, const std::vector<int>& triggerIndices, const std::unordered_map<int, std::string>& triggerNameMap) {
    std::ifstream file(csvFilePath);
    std::string line;
    
    // Skip the header line
    std::getline(file, line);
    
    // Map to store counts of 'ON' and 'OFF' for each trigger
    std::unordered_map<int, int> onCounts;
    std::unordered_map<int, int> offCounts;

    // Initialize counts for each trigger
    for (int triggerIndex : triggerIndices) {
        onCounts[triggerIndex] = 0;
        offCounts[triggerIndex] = 0;
    }

    // Read the CSV line by line
    while (std::getline(file, line)) {
        std::istringstream lineStream(line);
        std::string cell;
        int columnIndex = 0;
        int runNumber;
        
        // Parse the run number
        std::getline(lineStream, cell, ',');
        runNumber = std::stoi(cell);

        // Count 'ON'/'OFF' for each trigger
        for (int triggerIndex : triggerIndices) {
            std::getline(lineStream, cell, ',');
            if (cell == "ON") {
                onCounts[triggerIndex]++;
            } else if (cell == "OFF") {
                offCounts[triggerIndex]++;
            }
        }
    }
    file.close();

    // Output the summary counts to the terminal
    std::cout << "Summary of Trigger Counts:\n";
    for (int triggerIndex : triggerIndices) {
        std::cout << "Trigger " << triggerIndex << " (" << triggerNameMap.at(triggerIndex) << "):\n";
        std::cout << "  Runs with Trigger ON: " << onCounts[triggerIndex] << "\n";
        std::cout << "  Runs with Trigger OFF: " << offCounts[triggerIndex] << "\n";
        std::cout << "----------------------------------------\n";
    }
}

void analyzeTriggersAndWriteCombinedCSV(
    const std::vector<int>& runNumbers1,
    const std::vector<int>& runNumbers2,
    const std::unordered_map<int, std::string>& triggerNameMap1,
    const std::unordered_map<int, std::string>& triggerNameMap2,
    const std::string& outputFilePath) {

    // Combine the run numbers
    std::vector<int> allRunNumbers = runNumbers1;
    allRunNumbers.insert(allRunNumbers.end(), runNumbers2.begin(), runNumbers2.end());

    // Define the combined trigger list in the desired order
    std::vector<std::string> combinedTriggerNames = {
        "MBD_NandS_geq_1",
        "Photon_1_GeV_plus_MBD_NS_geq_1",
        "Photon_2_GeV_plus_MBD_NS_geq_1",
        "Photon_3_GeV_plus_MBD_NS_geq_1",
        "Photon_4_GeV_plus_MBD_NS_geq_1",
        "Photon_5_GeV_plus_MBD_NS_geq_1",
        "Photon_1_GeV",
        "Photon_2_GeV",
        "Photon_3_GeV",
        "Photon_4_GeV",
        "Photon_5_GeV"
    };

    // Map trigger names to indices for both trigger maps
    std::unordered_map<std::string, int> triggerNameToIndexMap1, triggerNameToIndexMap2;
    for (const auto& kv : triggerNameMap1) {
        triggerNameToIndexMap1[kv.second] = kv.first;
    }
    for (const auto& kv : triggerNameMap2) {
        triggerNameToIndexMap2[kv.second] = kv.first;
    }

    // Open the output CSV file
    std::ofstream csvFile(outputFilePath);

    // Write the header
    csvFile << "runNumber";
    for (const auto& triggerName : combinedTriggerNames) {
        csvFile << "," << triggerName;
    }
    csvFile << "\n";

    // Process each run
    for (int runNumber : allRunNumbers) {
        int scaledowns[64] = {0};

        // Get scaledown factors for the run
        if (!get_scaledowns(runNumber, scaledowns)) {
            std::cerr << "Skipping run number: " << runNumber << " due to error in retrieving scaledowns." << std::endl;
            continue;
        }

        csvFile << runNumber;

        // Determine which trigger map to use
        const std::unordered_map<std::string, int>* triggerNameToIndexMap;
        if (std::find(runNumbers1.begin(), runNumbers1.end(), runNumber) != runNumbers1.end()) {
            triggerNameToIndexMap = &triggerNameToIndexMap1;
        } else {
            triggerNameToIndexMap = &triggerNameToIndexMap2;
        }

        // For each trigger in the combined list
        for (const auto& triggerName : combinedTriggerNames) {
            auto it = triggerNameToIndexMap->find(triggerName);
            if (it != triggerNameToIndexMap->end()) {
                int triggerIndex = it->second;
                int scaledown = scaledowns[triggerIndex];
                if (scaledown == -1) {
                    csvFile << ",OFF";
                } else {
                    csvFile << ",ON";
                }
            } else {
                // Trigger not present in this trigger map, set to OFF
                csvFile << ",OFF";
            }
        }
        csvFile << "\n";
    }

    csvFile.close();
}

void analyzeCombinationsFromCSV(const std::string& csvFilePath) {
    // List of trigger names we're interested in
    std::vector<std::string> triggersOfInterest = {
        "MBD_NandS_geq_1",
        "Photon_1_GeV_plus_MBD_NS_geq_1",
        "Photon_2_GeV_plus_MBD_NS_geq_1",
        "Photon_3_GeV_plus_MBD_NS_geq_1",
        "Photon_4_GeV_plus_MBD_NS_geq_1",
        "Photon_5_GeV_plus_MBD_NS_geq_1"
    };

    // Map from trigger name to column index
    std::map<std::string, int> triggerToIndex;

    // List of runs
    std::vector<int> runNumbers;

    // Map from run number to trigger status
    std::map<int, std::map<std::string, std::string>> runTriggerStatus;

    // Open the CSV file
    std::ifstream file(csvFilePath);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << csvFilePath << std::endl;
        return;
    }

    // Read the header line
    std::string line;
    if (!std::getline(file, line)) {
        std::cerr << "Failed to read header from CSV file." << std::endl;
        return;
    }

    // Parse the header to get the indices of the triggers
    std::vector<std::string> headers;
    std::istringstream headerStream(line);
    std::string header;
    int colIndex = 0;
    while (std::getline(headerStream, header, ',')) {
        // Trim whitespace
        header.erase(0, header.find_first_not_of(" \t\r\n"));
        header.erase(header.find_last_not_of(" \t\r\n") + 1);

        headers.push_back(header);
        // If header is in triggersOfInterest, store its index
        if (std::find(triggersOfInterest.begin(), triggersOfInterest.end(), header) != triggersOfInterest.end()) {
            triggerToIndex[header] = colIndex;
        } else if (header == "runNumber") {
            triggerToIndex[header] = colIndex;
        }
        colIndex++;
    }

    // Read each line of the CSV
    while (std::getline(file, line)) {
        // Parse the line into cells
        std::vector<std::string> cells;
        std::istringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            // Trim whitespace and carriage returns
            cell.erase(0, cell.find_first_not_of(" \t\r\n"));
            cell.erase(cell.find_last_not_of(" \t\r\n") + 1);

            cells.push_back(cell);
        }

        // Ensure that the number of cells matches the number of headers
        if (cells.size() != headers.size()) {
            std::cerr << "Mismatch between number of cells and headers in line: " << line << std::endl;
            continue;
        }

        // Get run number
        int runNumber = std::stoi(cells[triggerToIndex["runNumber"]]);

        // Store the status of triggers of interest for this run
        for (const auto& trigger : triggersOfInterest) {
            if (triggerToIndex.find(trigger) != triggerToIndex.end()) {
                int idx = triggerToIndex[trigger];
                std::string status = cells[idx];
                // Trim whitespace
                status.erase(0, status.find_first_not_of(" \t\r\n"));
                status.erase(status.find_last_not_of(" \t\r\n") + 1);

                runTriggerStatus[runNumber][trigger] = status;
            } else {
                // If trigger not found in headers, assume 'OFF'
                runTriggerStatus[runNumber][trigger] = "OFF";
            }
        }
        // Keep track of run numbers
        runNumbers.push_back(runNumber);
    }

    file.close();

    // Generate all combinations of triggers of interest
    int n = triggersOfInterest.size();
    int totalCombinations = 1 << n; // 2^n combinations

    // Map from combination (set of triggers) to vector of run numbers
    std::map<std::set<std::string>, std::vector<int>> combinationToRuns;

    // For each combination
    for (int mask = 1; mask < totalCombinations; ++mask) {
        std::set<std::string> combination;
        // Build the combination based on the bits in mask
        for (int i = 0; i < n; ++i) {
            if (mask & (1 << i)) {
                combination.insert(triggersOfInterest[i]);
            }
        }

        // For each run, check if all triggers in the combination are 'ON'
        for (const auto& runNumber : runNumbers) {
            bool allOn = true;
            for (const auto& trigger : combination) {
                if (runTriggerStatus[runNumber][trigger] != "ON") {
                    allOn = false;
                    break;
                }
            }
            if (allOn) {
                combinationToRuns[combination].push_back(runNumber);
            }
        }
    }

    // Output the results
    std::cout << "\nSummary of Trigger Combinations:\n";
    // Sort combinations based on size and triggers
    std::vector<std::pair<std::set<std::string>, std::vector<int>>> sortedCombinations(
        combinationToRuns.begin(), combinationToRuns.end());

    std::sort(sortedCombinations.begin(), sortedCombinations.end(),
              [](const auto& a, const auto& b) {
                  // Compare based on the number of triggers active
                  int countA = a.first.size();
                  int countB = b.first.size();
                  if (countA != countB) {
                      return countA < countB;
                  } else {
                      // If same number of triggers, sort alphabetically
                      return a.first < b.first;
                  }
              });

    for (const auto& kv : sortedCombinations) {
        const std::set<std::string>& triggers = kv.first;
        const std::vector<int>& runs = kv.second;

        // Output the combination and counts
        std::cout << "Combination: ";
        for (const std::string& t : triggers) {
            std::cout << t << " ";
        }
        std::cout << "\n";
        std::cout << "Number of runs: " << runs.size() << "\n";
        std::cout << "Run numbers: ";
        for (size_t i = 0; i < runs.size(); ++i) {
            std::cout << runs[i];
            if ((i + 1) % 10 == 0) {
                std::cout << "\n             ";
            } else {
                std::cout << " ";
            }
        }
        std::cout << "\n-------------------------------------\n";
    }
}
void AnalyzeTriggerOnOrOff() {
    // First set: Using triggerNameMap1 and runListTriggerLUTv1.txt
    {
        // Path to the input text file containing run numbers
        std::string inputFile = "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/runListTriggerLUTv1.txt";

        // Load run numbers from file
        std::vector<int> runNumbers = loadRunNumbersFromFile(inputFile);

        if (runNumbers.empty()) {
            std::cerr << "No valid run numbers found in the file: " << inputFile << std::endl;
            return;
        }

        // Define trigger bits of interest
        std::vector<int> triggerBits = {10};
        for (int i = 24; i <= 31; ++i) {
            triggerBits.push_back(i);
        }

        std::string csvOutputFile = "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/triggerAnalysisLUTv1.csv";
        analyzeTriggersAndWriteCSV(runNumbers, triggerBits, triggerNameMap1, csvOutputFile);

        std::cout << "CSV file saved to: " << csvOutputFile << std::endl;
        summarizeTriggerCountsFromCSV(csvOutputFile, triggerBits, triggerNameMap1);
    }

    // Second set: Using triggerNameMap2 and runListTriggerLUTv2.txt
    {
        // Path to the input text file containing run numbers
        std::string inputFile = "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/runListTriggerLUTv2.txt";

        // Load run numbers from file
        std::vector<int> runNumbers = loadRunNumbersFromFile(inputFile);

        if (runNumbers.empty()) {
            std::cerr << "No valid run numbers found in the file: " << inputFile << std::endl;
            return;
        }

        // Define trigger bits of interest
        std::vector<int> triggerBits = {10};
        for (int i = 24; i <= 31; ++i) {
            triggerBits.push_back(i);
        }

        std::string csvOutputFile = "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/triggerAnalysisLUTv2.csv";
        analyzeTriggersAndWriteCSV(runNumbers, triggerBits, triggerNameMap2, csvOutputFile);

        std::cout << "CSV file saved to: " << csvOutputFile << std::endl;
        summarizeTriggerCountsFromCSV(csvOutputFile, triggerBits, triggerNameMap2);
    }
    // Combined processing
    {
         // Paths to the input text files containing run numbers
         std::string inputFile1 = "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/runListTriggerLUTv1.txt";
         std::string inputFile2 = "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/runListTriggerLUTv2.txt";

         // Load run numbers from files
         std::vector<int> runNumbers1 = loadRunNumbersFromFile(inputFile1);
         std::vector<int> runNumbers2 = loadRunNumbersFromFile(inputFile2);

         if (runNumbers1.empty() && runNumbers2.empty()) {
             std::cerr << "No valid run numbers found in the input files." << std::endl;
             return;
         }

         // Define the output file path for the combined CSV
         std::string combinedCsvOutputFile = "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/triggerAnalysisCombined.csv";

         // Generate the combined CSV
         analyzeTriggersAndWriteCombinedCSV(runNumbers1, runNumbers2, triggerNameMap1, triggerNameMap2, combinedCsvOutputFile);

         std::cout << "Combined CSV file saved to: " << combinedCsvOutputFile << std::endl;
        
         analyzeCombinationsFromCSV(combinedCsvOutputFile);
     }
}
