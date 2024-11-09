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

// Define the mapping from trigger indices to names
std::unordered_map<int, std::string> triggerNameMap = {
    {0, "Clock"},
    {1, "ZDC_South"},
    {2, "ZDC_North"},
    {3, "ZDC_Coincidence"},
    {4, "HCAL_Singles"},
    {5, "HCAL_Coincidence"},
    {8, "MBD_S_>=_1"},
    {9, "MBD_N_>=_1"},
    {10, "MBD_N&S_>=_1"},
    {11, "MBD_N&S_>=_2"},
    {12, "MBD_N&S_>=_1_vtx_<_10_cm"},
    {13, "MBD_N&S_>=_1_vtx_<_30_cm"},
    {14, "MBD_N&S_>=_1_vtx_<_60_cm"},
    {15, "HCAL_Singles_+_MBD_NS_>=_1"},
    {16, "Jet_4_GeV_+_MBD_NS_>=_1"},
    {17, "Jet_6_GeV_+_MBD_NS_>=_1"},
    {18, "Jet_8_GeV_+_MBD_NS_>=_1"},
    {19, "Jet_10_GeV_+_MBD_NS_>=_1"},
    {20, "Jet_4_GeV"},
    {21, "Jet_6_GeV"},
    {22, "Jet_8_GeV"},
    {23, "Jet_10_GeV"},
    {24, "Photon_1_GeV_+_MBD_NS_>=_1"},
    {25, "Photon_2_GeV_+_MBD_NS_>=_1"},
    {26, "Photon_3_GeV_+_MBD_NS_>=_1"},
    {27, "Photon_4_GeV_+_MBD_NS_>=_1"},
    {28, "Photon_1_GeV"},
    {29, "Photon_2_GeV"},
    {30, "Photon_3_GeV"},
    {31, "Photon_4_GeV"}
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

void analyzeTriggersAndWriteCSV(const std::vector<int>& runNumbers, const std::vector<int>& triggerIndices, const std::string& outputFilePath) {
    std::ofstream csvFile(outputFilePath);
    csvFile << "runNumber";

    for (int triggerIndex : triggerIndices) {
        csvFile << "," << triggerNameMap[triggerIndex];
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

    // Output the results
    for (int triggerIndex : triggerIndices) {
        std::string triggerName = triggerNameMap[triggerIndex];
        TriggerRunData& data = triggerDataMap[triggerIndex];

        std::cout << "Trigger Bit " << triggerIndex << " (" << triggerName << "):" << std::endl;

        // Output runs where trigger is ON
        std::cout << "  Runs with Trigger ON (" << data.runsOn.size() << " runs):" << std::endl;
        for (int runNumber : data.runsOn) {
            std::cout << "    " << runNumber << std::endl;
        }

        // Output runs where trigger is OFF
        std::cout << "  Runs with Trigger OFF (" << data.runsOff.size() << " runs):" << std::endl;
        for (int runNumber : data.runsOff) {
            std::cout << "    " << runNumber << std::endl;
        }

        std::cout << "----------------------------------------" << std::endl;
    }
}

void AnalyzeTriggerOnOrOff() {
    // Path to the input text file containing run numbers
    std::string inputFile = "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/runListTriggerLUTv1.txt";

    // Load run numbers from file
    std::vector<int> runNumbers = loadRunNumbersFromFile(inputFile);

    if (runNumbers.empty()) {
        std::cerr << "No valid run numbers found in the file: " << inputFile << std::endl;
        return 1;
    }

    // Define trigger bits of interest
    std::vector<int> triggerBits = {10};
    for (int i = 24; i <= 31; ++i) {
        triggerBits.push_back(i);
    }

    std::string csvOutputFile = "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/triggerAnalysisLUTv1.csv";
    analyzeTriggersAndWriteCSV(runNumbers, triggerBits, csvOutputFile);

    std::cout << "CSV file saved to: " << csvOutputFile << std::endl;
}

