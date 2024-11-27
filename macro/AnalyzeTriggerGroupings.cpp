#include "AnalyzeTriggerGroupings.h"
#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <memory>

// ANSI escape codes for colors
#define RESET   "\033[0m"
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define BOLD    "\033[1m"

bool enableFits = true; // Set to true if you want to enable the fits


std::map<std::tuple<
    std::string, // triggerGroupName (raw name)
    std::string, // triggerName (raw name)
    float,       // clusECore
    float,       // chi
    float,       // asymmetry
    float,       // pTMin
    float,       // pTMax
    float,       // isoEtMin
    float,       // isoEtMax
    std::string  // massWindowLabel
>, DataStructures::IsolatedPhotonLog> isolatedPhotonMap;

std::map<std::tuple<
    std::string, // triggerGroupName (raw name)
    std::string, // triggerName (raw name)
    float,       // clusECore
    float,       // chi
    float,       // asymmetry
    float,       // pTMin
    float,       // pTMax
    std::string  // massWindowLabel
>, DataStructures::TotalPhotonLog> totalPhotonMap;

std::map<std::tuple<
    std::string, // triggerGroupName (raw name)
    std::string, // triggerName (raw name)
    float,       // clusECore
    float,       // chi
    float,       // asymmetry
    float,       // pTMin
    float,       // pTMax
    std::string  // massWindowLabel
>, DataStructures::PtWeightingLog> pTweightingMap;

std::map<std::tuple<
    std::string, // TriggerGroupName
    std::string, // TriggerName
    float,       // ECore
    float,       // Chi
    float,       // Asymmetry
    float,       // pT Min
    float,       // pT Max
    float,       // isoMin
    float,       // isoMax
    std::string  // MassWindowLabel
>, DataStructures::IsolationData> dataMap_inMassWindow;

std::map<std::tuple<
    std::string, // TriggerGroupName
    std::string, // TriggerName
    float,       // ECore
    float,       // Chi
    float,       // Asymmetry
    float,       // pT Min
    float,       // pT Max
    float,       // isoMin
    float,       // isoMax
    std::string  // MassWindowLabel
>, DataStructures::IsolationData> dataMap_outsideMassWindow;

// Define GroupKey as a type alias for reuse across functions
using GroupKey = std::tuple<
    std::string, // TriggerGroupName
    std::string, // TriggerName
    float,       // ECore
    float,       // Chi
    float,       // Asymmetry
    std::string  // MassWindowLabel
>;

void estimateSigmoidParameters(TH1* ratioHist, double& amplitude, double& xOffset, double& slope) {
    // Estimate amplitude as the maximum value of the histogram
    amplitude = ratioHist->GetMaximum();

    // Find the x value where the histogram reaches 50% of the amplitude
    double halfMax = amplitude / 2.0;
    int nBins = ratioHist->GetNbinsX();
    double x50 = 0.0;
    bool found50 = false;
    
    for (int i = 1; i <= nBins; ++i) {
        double y = ratioHist->GetBinContent(i);
        if (y >= halfMax) {
            x50 = ratioHist->GetBinCenter(i);
            found50 = true;
            break;
        }
    }
    
    if (!found50) x50 = 1.0;  // Fallback if half-max point isn't found
    xOffset = x50;

    // Estimate slope based on the width over which the histogram rises from 20% to 80% of the amplitude
    double y20 = amplitude * 0.2;
    double y80 = amplitude * 0.8;
    double x20 = 0.0, x80 = 0.0;
    bool found20 = false, found80 = false;
    for (int i = 1; i <= nBins; ++i) {
        double y = ratioHist->GetBinContent(i);
        if (!found20 && y >= y20) {
            x20 = ratioHist->GetBinCenter(i);
            found20 = true;
        }
        if (!found80 && y >= y80) {
            x80 = ratioHist->GetBinCenter(i);
            found80 = true;
        }
        if (found20 && found80) {
            break;
        }
    }

    double deltaX = x80 - x20;
    slope = (deltaX > 0) ? 4.0 / deltaX : 0.1;  // Default to a gentle slope if deltaX is zero
}

std::map<std::set<std::string>, DataStructures::RunInfo> AnalyzeWhatTriggerGroupsAvailable(
    const std::string& csvFilePath,
    bool debugMode,
    const std::map<int, std::map<std::string, std::string>>& overrideTriggerStatus) {


    const std::vector<std::string>& allTriggers = TriggerConfig::allTriggers;
    const std::vector<std::string>& photonTriggers = TriggerConfig::photonTriggers;

    // Map from trigger name to column index in the CSV
    std::map<std::string, int> triggerToIndex;

    // Map from run number to set of triggers that are 'ON'
    std::map<int, std::set<std::string>> runToActiveTriggers;

    // Variables for debug output
    int totalRunsProcessed = 0;
    std::map<std::set<std::string>, std::vector<int>> initialCombinationToRuns;
    std::vector<std::pair<std::set<std::string>, DataStructures::RunInfo>> finalCombinations;

    // Open the CSV file
    std::ifstream file(csvFilePath);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << csvFilePath << std::endl;
        return {};
    }

    // Read the header line
    std::string line;
    if (!std::getline(file, line)) {
        std::cerr << "Failed to read header from CSV file." << std::endl;
        return {};
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
        // If header is in allTriggers or is 'runNumber', store its index
        if (std::find(allTriggers.begin(), allTriggers.end(), header) != allTriggers.end() || header == "runNumber") {
            triggerToIndex[header] = colIndex;
        }
        colIndex++;
    }

    if (triggerToIndex.find("runNumber") == triggerToIndex.end()) {
        std::cerr << "runNumber column not found in CSV header." << std::endl;
        return {};
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

        // Collect triggers that are 'ON' for this run
        std::set<std::string> activeTriggers;
        for (const auto& trigger : allTriggers) {
            if (triggerToIndex.find(trigger) != triggerToIndex.end()) {
                int idx = triggerToIndex[trigger];
                std::string status = cells[idx];
                // Trim whitespace
                status.erase(0, status.find_first_not_of(" \t\r\n"));
                status.erase(status.find_last_not_of(" \t\r\n") + 1);

                // Check if there is an override for this run and trigger
                auto runOverrideIt = overrideTriggerStatus.find(runNumber);
                if (runOverrideIt != overrideTriggerStatus.end()) {
                    const auto& triggerOverrides = runOverrideIt->second;
                    auto triggerOverrideIt = triggerOverrides.find(trigger);
                    if (triggerOverrideIt != triggerOverrides.end()) {
                        // Override exists, use it
                        status = triggerOverrideIt->second;
                    }
                }

                if (status == "ON") {
                    activeTriggers.insert(trigger);
                }
            }
        }
        // Store active triggers for this run
        runToActiveTriggers[runNumber] = activeTriggers;
    }

    file.close();

    // Total runs processed
    totalRunsProcessed = runToActiveTriggers.size();

    // Now generate all combinations of triggers we're interested in
    // Generate all subsets of photonTriggers
    std::vector<std::set<std::string>> triggerCombinations;
    int nPhotonTriggers = photonTriggers.size();
    int totalCombinations = 1 << nPhotonTriggers; // 2^n combinations

    for (int mask = 0; mask < totalCombinations; ++mask) {
        std::set<std::string> combination;
        combination.insert("MBD_NandS_geq_1"); // Always include MBD_NandS_geq_1
        for (int i = 0; i < nPhotonTriggers; ++i) {
            if (mask & (1 << i)) {
                combination.insert(photonTriggers[i]);
            }
        }
        triggerCombinations.push_back(combination);
    }

    // Map from combination to vector of run numbers
    std::map<std::set<std::string>, std::vector<int>> tempCombinationToRuns;

    // For each run, check which combinations it satisfies
    for (const auto& runEntry : runToActiveTriggers) {
        int runNumber = runEntry.first;
        const std::set<std::string>& activeTriggers = runEntry.second;

        // Proceed only if 'MBD_NandS_geq_1' is 'ON'
        if (activeTriggers.find("MBD_NandS_geq_1") != activeTriggers.end()) {
            // For each combination, check if it is satisfied
            for (const auto& combination : triggerCombinations) {
                // Check if all triggers in the combination are in activeTriggers
                bool satisfiesCombination = true;
                for (const auto& trigger : combination) {
                    if (activeTriggers.find(trigger) == activeTriggers.end()) {
                        satisfiesCombination = false;
                        break;
                    }
                }
                if (satisfiesCombination) {
                    // Add run number to this combination
                    tempCombinationToRuns[combination].push_back(runNumber);
                }
            }
        }
    }

    // Store initial combinations and their runs before splitting
    initialCombinationToRuns = tempCombinationToRuns;

    // Now, for each combination, split runs into runsBeforeFirmwareUpdate and runsAfterFirmwareUpdate
    struct TempRunInfo {
        std::vector<int> runsBeforeFirmwareUpdate;
        std::vector<int> runsAfterFirmwareUpdate;
    };
    std::map<std::set<std::string>, TempRunInfo> tempCombinationToRunInfo;

    // Collect combinations that include run 47289
    std::set<std::set<std::string>> combinationsWithRun47289;

    for (const auto& entry : tempCombinationToRuns) {
        const std::set<std::string>& combination = entry.first;
        const std::vector<int>& runs = entry.second;

        TempRunInfo runInfo;
        bool hasRun47289 = false;
        for (int runNumber : runs) {
            if (runNumber == 47289) {
                hasRun47289 = true;
            }
            if (runNumber < 47289) {
                runInfo.runsBeforeFirmwareUpdate.push_back(runNumber);
            } else {
                runInfo.runsAfterFirmwareUpdate.push_back(runNumber);
            }
        }

        if (hasRun47289) {
            combinationsWithRun47289.insert(combination);
        }

        tempCombinationToRunInfo[combination] = runInfo;
    }

    // Now, process runsBeforeFirmwareUpdate and runsAfterFirmwareUpdate separately

    // For runsBeforeFirmwareUpdate
    std::map<std::vector<int>, std::vector<std::set<std::string>>> runListToCombinationsBeforeFirmwareUpdate;

    for (const auto& entry : tempCombinationToRunInfo) {
        const std::set<std::string>& combination = entry.first;
        const std::vector<int>& runList = entry.second.runsBeforeFirmwareUpdate;
        if (runList.empty()) continue;
        std::vector<int> sortedRunList = runList;
        std::sort(sortedRunList.begin(), sortedRunList.end());

        runListToCombinationsBeforeFirmwareUpdate[sortedRunList].push_back(combination);
    }

    // For runsAfterFirmwareUpdate
    std::map<std::vector<int>, std::vector<std::set<std::string>>> runListToCombinationsAfterFirmwareUpdate;

    for (const auto& entry : tempCombinationToRunInfo) {
        const std::set<std::string>& combination = entry.first;
        const std::vector<int>& runList = entry.second.runsAfterFirmwareUpdate;
        if (runList.empty()) continue;
        std::vector<int> sortedRunList = runList;
        std::sort(sortedRunList.begin(), sortedRunList.end());

        runListToCombinationsAfterFirmwareUpdate[sortedRunList].push_back(combination);
    }

    // For each run list, find the largest combination(s) and remove subsets

    // Process before firmware update
    std::map<std::set<std::string>, std::vector<int>> filteredCombinationToRunsBeforeFirmwareUpdate;

    for (const auto& entry : runListToCombinationsBeforeFirmwareUpdate) {
        const std::vector<int>& runList = entry.first;
        const std::vector<std::set<std::string>>& combinations = entry.second;

        // Find the combination(s) with the maximum size (most triggers)
        size_t maxSize = 0;
        for (const auto& combo : combinations) {
            if (combo.size() > maxSize) {
                maxSize = combo.size();
            }
        }

        // Collect combinations with maximum size
        for (const auto& combo : combinations) {
            if (combo.size() == maxSize) {
                // Add to the filtered map
                filteredCombinationToRunsBeforeFirmwareUpdate[combo] = runList;
            }
        }
    }

    // Process after firmware update
    std::map<std::set<std::string>, std::vector<int>> filteredCombinationToRunsAfterFirmwareUpdate;

    for (const auto& entry : runListToCombinationsAfterFirmwareUpdate) {
        const std::vector<int>& runList = entry.first;
        const std::vector<std::set<std::string>>& combinations = entry.second;

        // Find the combination(s) with the maximum size (most triggers)
        size_t maxSize = 0;
        for (const auto& combo : combinations) {
            if (combo.size() > maxSize) {
                maxSize = combo.size();
            }
        }

        // Collect combinations with maximum size
        for (const auto& combo : combinations) {
            if (combo.size() == maxSize) {
                // Add to the filtered map
                filteredCombinationToRunsAfterFirmwareUpdate[combo] = runList;
            }
        }
    }

    // Now, assemble combinationToRuns
    std::map<std::set<std::string>, DataStructures::RunInfo> combinationToRuns;

    // Handle before firmware update
    for (const auto& entry : filteredCombinationToRunsBeforeFirmwareUpdate) {
        const std::set<std::string>& combination = entry.first;
        const std::vector<int>& runs = entry.second;

        DataStructures::RunInfo& runInfo = combinationToRuns[combination];
        runInfo.runsBeforeFirmwareUpdate = runs;
    }

    // Handle after firmware update
    for (const auto& entry : filteredCombinationToRunsAfterFirmwareUpdate) {
        const std::set<std::string>& combination = entry.first;
        const std::vector<int>& runs = entry.second;

        DataStructures::RunInfo& runInfo = combinationToRuns[combination];
        runInfo.runsAfterFirmwareUpdate = runs;
    }

    // Store final combinations for debug output
    finalCombinations.assign(combinationToRuns.begin(), combinationToRuns.end());

    // If debugMode is true, output the processing information and terminate
    if (debugMode) {
        // Output the detailed processing information
        std::cout << BOLD << BLUE << "\n===== Processing Summary =====\n" << RESET;

        // Output total number of runs processed
        std::cout << BOLD << "Total runs processed: " << totalRunsProcessed << RESET << "\n";

        // Output initial combinations before splitting
        std::cout << BOLD << "\nInitial Active Trigger Combinations (before splitting due to firmware update):\n" << RESET;
        // Prepare header
        std::cout << BOLD << std::left << std::setw(60) << "Combination" << std::right << std::setw(20) << "Number of Runs" << RESET << "\n";
        std::cout << std::string(80, '=') << "\n";

        for (const auto& entry : initialCombinationToRuns) {
            const std::set<std::string>& combination = entry.first;
            const std::vector<int>& runs = entry.second;

            // Build combination string
            std::string combinationStr;
            for (const auto& trigger : combination) {
                combinationStr += trigger + " ";
            }

            // Output combination and number of runs
            std::cout << std::left << std::setw(60) << combinationStr << std::right << std::setw(20) << runs.size() << "\n";
        }

        // Output combinations that had run 47289 and were split
        std::cout << BOLD << "\nCombinations that included run 47289 and were split due to firmware update:\n" << RESET;
        if (combinationsWithRun47289.empty()) {
            std::cout << "  None\n";
        } else {
            // Prepare header
            std::cout << BOLD << std::left << std::setw(60) << "Combination" << RESET << "\n";
            std::cout << std::string(60, '=') << "\n";

            for (const auto& combination : combinationsWithRun47289) {
                // Build combination string
                std::string combinationStr;
                for (const auto& trigger : combination) {
                    combinationStr += trigger + " ";
                }

                std::cout << std::left << std::setw(60) << combinationStr << "\n";
            }
        }

        // Output groups found with the same run numbers before firmware update
        std::cout << BOLD << "\nGroups with identical run numbers before firmware update:\n" << RESET;
        int groupIndex = 1;
        for (const auto& entry : runListToCombinationsBeforeFirmwareUpdate) {
            const std::vector<int>& runList = entry.first;
            const std::vector<std::set<std::string>>& combinations = entry.second;

            std::cout << YELLOW << BOLD << "\nGroup " << groupIndex++ << RESET << "\n";
            std::cout << "Run List (size " << runList.size() << "):\n";
            // Print run numbers in columns
            for (size_t i = 0; i < runList.size(); ++i) {
                std::cout << std::setw(8) << runList[i];
                if ((i + 1) % 10 == 0 || i == runList.size() - 1) {
                    std::cout << "\n";
                }
            }

            std::cout << "  Combinations:\n";
            for (const auto& combo : combinations) {
                // Build combination string
                std::string combinationStr;
                for (const auto& trigger : combo) {
                    combinationStr += trigger + " ";
                }
                std::cout << "    " << combinationStr << "\n";
            }
        }

        // Output groups found with the same run numbers after firmware update
        std::cout << BOLD << "\nGroups with identical run numbers after firmware update:\n" << RESET;
        groupIndex = 1;
        for (const auto& entry : runListToCombinationsAfterFirmwareUpdate) {
            const std::vector<int>& runList = entry.first;
            const std::vector<std::set<std::string>>& combinations = entry.second;

            std::cout << YELLOW << BOLD << "\nGroup " << groupIndex++ << RESET << "\n";
            std::cout << "Run List (size " << runList.size() << "):\n";
            // Print run numbers in columns
            for (size_t i = 0; i < runList.size(); ++i) {
                std::cout << std::setw(8) << runList[i];
                if ((i + 1) % 10 == 0 || i == runList.size() - 1) {
                    std::cout << "\n";
                }
            }

            std::cout << "  Combinations:\n";
            for (const auto& combo : combinations) {
                // Build combination string
                std::string combinationStr;
                for (const auto& trigger : combo) {
                    combinationStr += trigger + " ";
                }
                std::cout << "    " << combinationStr << "\n";
            }
        }

        // Output final combinations after splitting and filtering
        std::cout << BOLD << "\nFinal Active Trigger Combinations (after splitting and filtering):\n" << RESET;
        // Prepare header
        std::cout << BOLD << std::left << std::setw(60) << "Combination"
                  << std::right << std::setw(25) << "Runs Before Firmware Update"
                  << std::setw(25) << "Runs After Firmware Update" << RESET << "\n";
        std::cout << std::string(110, '=') << "\n";

        for (const auto& entry : finalCombinations) {
            const std::set<std::string>& combination = entry.first;
            const DataStructures::RunInfo& runInfo = entry.second;

            // Build combination string
            std::string combinationStr;
            for (const auto& trigger : combination) {
                combinationStr += trigger + " ";
            }

            std::cout << std::left << std::setw(60) << combinationStr;
            std::cout << std::right << std::setw(25) << runInfo.runsBeforeFirmwareUpdate.size();
            std::cout << std::setw(25) << runInfo.runsAfterFirmwareUpdate.size() << "\n";

            // Optionally, print run numbers
            if (!runInfo.runsBeforeFirmwareUpdate.empty()) {
                std::cout << "  Runs before firmware update (" << runInfo.runsBeforeFirmwareUpdate.size() << " runs):\n";
                for (size_t i = 0; i < runInfo.runsBeforeFirmwareUpdate.size(); ++i) {
                    std::cout << std::setw(8) << runInfo.runsBeforeFirmwareUpdate[i];
                    if ((i + 1) % 10 == 0 || i == runInfo.runsBeforeFirmwareUpdate.size() - 1) {
                        std::cout << "\n";
                    }
                }
            }
            if (!runInfo.runsAfterFirmwareUpdate.empty()) {
                std::cout << "  Runs after firmware update (" << runInfo.runsAfterFirmwareUpdate.size() << " runs):\n";
                for (size_t i = 0; i < runInfo.runsAfterFirmwareUpdate.size(); ++i) {
                    std::cout << std::setw(8) << runInfo.runsAfterFirmwareUpdate[i];
                    if ((i + 1) % 10 == 0 || i == runInfo.runsAfterFirmwareUpdate.size() - 1) {
                        std::cout << "\n";
                    }
                }
            }
        }

        std::cout << BOLD << BLUE << "\n===== End of Processing Summary =====\n" << RESET;

        // Terminate the program
        exit(0);
    }

    // Return the adjusted map
    return combinationToRuns;
}




void PrintSortedCombinations(const std::map<std::set<std::string>, DataStructures::RunInfo>& combinationToRuns) {
    std::vector<std::pair<std::set<std::string>, DataStructures::RunInfo>> sortedCombinations(
        combinationToRuns.begin(), combinationToRuns.end());
    
    // Sort combinations by the number of triggers in the set (descending)
    std::sort(sortedCombinations.begin(), sortedCombinations.end(),
              [](const auto& a, const auto& b) {
                  return a.first.size() > b.first.size();
              });

    for (const auto& entry : sortedCombinations) {
        const std::set<std::string>& combination = entry.first;
        const DataStructures::RunInfo& runInfo = entry.second;

        std::string combinationName;
        for (const auto& trigger : combination) {
            combinationName += trigger + " ";
        }

        std::cout << "Combination: " << combinationName << "\n";
        if (!runInfo.runsBeforeFirmwareUpdate.empty()) {
            std::cout << "  Runs before firmware update (" << runInfo.runsBeforeFirmwareUpdate.size() << " runs): ";
            for (int run : runInfo.runsBeforeFirmwareUpdate) {
                std::cout << run << " ";
            }
            std::cout << "\n";
        }
        if (!runInfo.runsAfterFirmwareUpdate.empty()) {
            std::cout << "  Runs after firmware update (" << runInfo.runsAfterFirmwareUpdate.size() << " runs): ";
            for (int run : runInfo.runsAfterFirmwareUpdate) {
                std::cout << run << " ";
            }
            std::cout << "\n";
        }
    }
}

void ProcessRunsForCombination(
    const std::string& combinationName,
    const std::vector<int>& runs,
    const std::set<std::string>& triggers,
    const std::string& outputDirectory,
    std::map<std::string, std::vector<int>>& combinationToValidRuns) {
    
    // Define the final output ROOT file path
    std::string finalRootFilePath = outputDirectory + "/" + combinationName + "_Combined.root";
    // Define the text file path to store valid runs
    std::string validRunsFilePath = outputDirectory + "/" + combinationName + "_ValidRuns.txt";
    
    // Check if both the final ROOT file and valid runs text file already exist to avoid overwriting
    bool rootFileExists = !gSystem->AccessPathName(finalRootFilePath.c_str());
    bool validRunsFileExists = !gSystem->AccessPathName(validRunsFilePath.c_str());
    
    if (rootFileExists && validRunsFileExists) {
        std::cout << "Final ROOT file and valid runs file already exist for combination: " << combinationName << ". Skipping merge." << std::endl;
        // Read valid runs from the text file
        std::vector<int> validRuns;
        std::ifstream validRunsFile(validRunsFilePath);
        if (validRunsFile.is_open()) {
            int runNumber;
            while (validRunsFile >> runNumber) {
                validRuns.push_back(runNumber);
            }
            validRunsFile.close();
            combinationToValidRuns[combinationName] = validRuns;
        } else {
            std::cerr << "Failed to open valid runs file: " << validRunsFilePath << std::endl;
        }
        return;
    }

    // Map to keep track of histograms (by trigger and histogram name)
    std::map<std::string, std::map<std::string, std::unique_ptr<TH1>>> mergedHistograms;
    
    // Collect histogram names
    std::set<std::string> histogramNames;
    
    // Vector to store valid run numbers for this combination
    std::vector<int> validRuns;
    
    // Iterate over each run number in the combination
    for (const auto& runNumber : runs) {
        // Define the run's ROOT file path
        std::stringstream ss;
        ss << outputDirectory << "/" << runNumber << "_HistOutput.root";
        std::string runRootFilePath = ss.str();
        
        std::cout << "Processing run: " << runNumber << ", file: " << runRootFilePath << std::endl;
        
        // Check if the run's ROOT file exists
        if (gSystem->AccessPathName(runRootFilePath.c_str())) {
            std::cerr << "Run ROOT file does not exist: " << runRootFilePath << ". Skipping this run." << std::endl;
            continue;
        }
        
        // Open the run's ROOT file
        std::cout << "Opening run file..." << std::endl;
        TFile runFile(runRootFilePath.c_str(), "READ");
        if (runFile.IsZombie() || !runFile.IsOpen()) {
            std::cerr << "Failed to open run ROOT file: " << runRootFilePath << ". Skipping this run." << std::endl;
            continue;
        }
        std::cout << "Run file opened successfully." << std::endl;
        
        bool runHasValidHistogram = false;
        
        
        // For each trigger in the combination
        for (const auto& trigger : triggers) {
            // Check if the trigger directory exists in the run file
            TDirectory* triggerDir = runFile.GetDirectory(trigger.c_str());
            if (!triggerDir) {
                std::cerr << "Trigger directory '" << trigger << "' not found in run " << runNumber << ". Skipping this trigger." << std::endl;
                continue;
            }
            
            // Get all histograms in the trigger directory
            TIter nextKey(triggerDir->GetListOfKeys());
            TKey* key;
            while ((key = (TKey*)nextKey())) {
                std::string className = key->GetClassName();
                if (className.find("TH1") != std::string::npos || className.find("TH2") != std::string::npos) {
                    TObject* obj = key->ReadObj();
                    TH1* hist = dynamic_cast<TH1*>(obj);
                    if (!hist) {
                        std::cerr << "Failed to read histogram '" << key->GetName() << "' in trigger '" << trigger << "' in run " << runNumber << std::endl;
                        delete obj;
                        continue;
                    }
                    std::string histName = hist->GetName();
                    
                    // Clone the histogram and set directory to nullptr
                    TH1* histClone = dynamic_cast<TH1*>(hist->Clone());
                    if (!histClone) {
                        std::cerr << "Failed to clone histogram: " << histName << " from run " << runNumber << std::endl;
                        delete obj;
                        continue;
                    }
                    histClone->SetDirectory(nullptr);
                    
                    // Check if we already have this histogram in mergedHistograms
                    auto& histMap = mergedHistograms[trigger];
                    auto it = histMap.find(histName);
                    if (it != histMap.end()) {
                        // Histogram already exists, add the new histogram to it
                        it->second->Add(histClone);
                        std::cout << "Added histogram '" << histName << "' from run " << runNumber << " to existing histogram in trigger '" << trigger << "'." << std::endl;
                        delete histClone;
                    } else {
                        // Add histClone to mergedHistograms
                        histMap[histName] = std::unique_ptr<TH1>(histClone);
                        std::cout << "Added histogram '" << histName << "' from run " << runNumber << " to merged histograms in trigger '" << trigger << "'." << std::endl;
                        // histClone is managed by unique_ptr in histMap
                    }
                    
                    // Clean up
                    delete obj;
                    
                    runHasValidHistogram = true;
                } else {
                    std::cout << "Skipping non-histogram object '" << key->GetName() << "' in trigger '" << trigger << "' in run " << runNumber << std::endl;
                }
            }
        }
        
        if (runHasValidHistogram) {
            // Add run number to validRuns
            validRuns.push_back(runNumber);
        }
        
        // Close the run file
        runFile.Close();
    }
    
    if (validRuns.empty()) {
        std::cout << "No valid runs found for combination: " << combinationName << ". Skipping this combination." << std::endl;
        return; // Replace continue with return
    }

    if (mergedHistograms.empty()) {
        std::cout << "No histograms were merged for combination: " << combinationName << ". Skipping writing output ROOT file." << std::endl;
        return; // Replace continue with return
    }

    // Write all merged histograms to the final ROOT file
    std::cout << "Writing merged histograms to final ROOT file: " << finalRootFilePath << std::endl;
    TFile finalFile(finalRootFilePath.c_str(), "RECREATE");
    if (finalFile.IsZombie() || !finalFile.IsOpen()) {
        std::cerr << "Failed to create final ROOT file: " << finalRootFilePath << std::endl;
        mergedHistograms.clear();
        return;
    }
    
    finalFile.cd();
    for (const auto& triggerHistPair : mergedHistograms) {
        const std::string& trigger = triggerHistPair.first;
        const auto& histMap = triggerHistPair.second;
        
        // Create the trigger directory
        TDirectory* triggerDir = finalFile.mkdir(trigger.c_str());
        if (!triggerDir) {
            std::cerr << "Failed to create directory for trigger: " << trigger << std::endl;
            continue;
        }
        
        triggerDir->cd();
        
        for (const auto& histPair : histMap) {
            const std::string& histName = histPair.first;
            TH1* hist = histPair.second.get();
            if (!hist) {
                std::cerr << "Null histogram encountered. Skipping." << std::endl;
                continue;
            }
            hist->Write();
            std::cout << "Histogram '" << histName << "' written to trigger directory '" << trigger << "' in final ROOT file." << std::endl;
        }
    }
    
    // Write and close the final ROOT file
    finalFile.Write();
    finalFile.Close();
    
    // Clean up merged histograms
    mergedHistograms.clear();
    
    std::cout << "Successfully created combined ROOT file: " << finalRootFilePath << std::endl;
    
    // Store the valid runs for this combination
    combinationToValidRuns[combinationName] = validRuns;
    
    // Write the valid runs to a text file
    std::ofstream validRunsFile(validRunsFilePath);
    if (validRunsFile.is_open()) {
        for (const auto& runNumber : validRuns) {
            validRunsFile << runNumber << "\n";
        }
        validRunsFile.close();
        std::cout << "Valid runs written to file: " << validRunsFilePath << std::endl;
    } else {
        std::cerr << "Failed to write valid runs to file: " << validRunsFilePath << std::endl;
    }
}


void ProcessAndMergeRootFiles(
    const std::map<std::set<std::string>, DataStructures::RunInfo>& combinationToRuns,
    const std::string& outputDirectory,
    std::map<std::string, std::vector<int>>& combinationToValidRuns) {
    
    std::cout << "Starting ProcessAndMergeRootFiles" << std::endl;
    // Disable automatic addition of histograms to directories
    TH1::AddDirectory(false);
    
    // Iterate over each trigger combination
    for (const auto& kv : combinationToRuns) {
        const std::set<std::string>& triggers = kv.first;
        const DataStructures::RunInfo& runInfo = kv.second;
        
        // Create a string to represent the combination for naming
        std::string baseCombinationName;
        for (const auto& trigger : triggers) {
            baseCombinationName += trigger + "_";
        }
        // Remove the trailing underscore
        if (!baseCombinationName.empty()) {
            baseCombinationName.pop_back();
        }

        // Check if this combination is exactly "MBD_NandS_geq_1"
        if (triggers.size() == 1 && triggers.count("MBD_NandS_geq_1") == 1) {
            // Combine runs before and after firmware update
            std::vector<int> allRuns = runInfo.runsBeforeFirmwareUpdate;
            allRuns.insert(allRuns.end(), runInfo.runsAfterFirmwareUpdate.begin(), runInfo.runsAfterFirmwareUpdate.end());

            std::string combinationName = baseCombinationName; // Do not append firmware update tags

            // Merge ROOT files for these runs
            ProcessRunsForCombination(combinationName, allRuns, triggers, outputDirectory, combinationToValidRuns);
        } else {
            // For other combinations, process runs before and after firmware update separately

            // Process runs before firmware update
            if (!runInfo.runsBeforeFirmwareUpdate.empty()) {
                std::string combinationName = baseCombinationName;
                if (!runInfo.runsAfterFirmwareUpdate.empty()) {
                    combinationName += "_beforeTriggerFirmwareUpdate";
                }
                const std::vector<int>& runs = runInfo.runsBeforeFirmwareUpdate;
                // Merge ROOT files for these runs
                ProcessRunsForCombination(combinationName, runs, triggers, outputDirectory, combinationToValidRuns);
            }
            
            // Process runs after firmware update
            if (!runInfo.runsAfterFirmwareUpdate.empty()) {
                std::string combinationName = baseCombinationName;
                if (!runInfo.runsBeforeFirmwareUpdate.empty()) {
                    combinationName += "_afterTriggerFirmwareUpdate";
                }
                const std::vector<int>& runs = runInfo.runsAfterFirmwareUpdate;
                // Merge ROOT files for these runs
                ProcessRunsForCombination(combinationName, runs, triggers, outputDirectory, combinationToValidRuns);
            }
        }
    }
    
    // Re-enable automatic addition of histograms to directories
    TH1::AddDirectory(true);
}




std::vector<std::string> ExtractTriggersFromFilename(const std::string& filename, const std::vector<std::string>& allTriggers) {
    std::string baseName = filename;
    std::string suffix = "_Combined.root";
    if (Utils::EndsWith(baseName, suffix)) {
        baseName = baseName.substr(0, baseName.length() - suffix.length());
    }

    // Strip firmware tag if present
    baseName = Utils::stripFirmwareTag(baseName);
    
    // Split baseName into tokens using underscores
    std::vector<std::string> filenameTokens;
    std::istringstream iss(baseName);
    std::string token;
    while (std::getline(iss, token, '_')) {
        filenameTokens.push_back(token);
    }

    // Tokenize each trigger
    std::vector<std::vector<std::string>> triggerTokensList;
    for (const auto& trigger : allTriggers) {
        std::vector<std::string> triggerTokens;
        std::istringstream triggerStream(trigger);
        std::string triggerToken;
        while (std::getline(triggerStream, triggerToken, '_')) {
            triggerTokens.push_back(triggerToken);
        }
        triggerTokensList.push_back(triggerTokens);
    }

    std::vector<std::string> triggersFound;

    size_t i = 0;
    while (i < filenameTokens.size()) {
        bool foundTrigger = false;
        // Try to match any trigger starting at position i
        for (size_t t = 0; t < allTriggers.size(); ++t) {
            const auto& triggerTokens = triggerTokensList[t];
            if (i + triggerTokens.size() <= filenameTokens.size()) {
                bool matches = true;
                for (size_t j = 0; j < triggerTokens.size(); ++j) {
                    if (filenameTokens[i + j] != triggerTokens[j]) {
                        matches = false;
                        break;
                    }
                }
                if (matches) {
                    // Found a trigger
                    triggersFound.push_back(allTriggers[t]);
                    i += triggerTokens.size();
                    foundTrigger = true;
                    break;
                }
            }
        }
        if (!foundTrigger) {
            // Move to next token
            ++i;
        }
    }

    return triggersFound;
}


TFitResultPtr PerformFitting(TH1* hPi0Mass, TF1*& totalFit, TF1*& gaussPi0Fit, TF1*& gaussEtaFit, TF1*& polyFit, double& fitStart, double& fitEnd) {
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);

    fitStart = 0.04;
    fitEnd = 0.9;

    // Define expected peak regions for pi0 and eta
    double pi0RangeMin = 0.1;
    double pi0RangeMax = 0.18;
    double etaRangeMin = 0.5;
    double etaRangeMax = 0.7;

    // Find maximum bin in pi0 peak region
    int binPi0Min = hPi0Mass->GetXaxis()->FindBin(pi0RangeMin);
    int binPi0Max = hPi0Mass->GetXaxis()->FindBin(pi0RangeMax);

    int maxBinPi0 = binPi0Min;
    double maxContentPi0 = 0.0;

    for (int bin = binPi0Min; bin <= binPi0Max; ++bin) {
        double content = hPi0Mass->GetBinContent(bin);
        if (content > maxContentPi0) {
            maxContentPi0 = content;
            maxBinPi0 = bin;
        }
    }

    double meanPi0Estimate = hPi0Mass->GetXaxis()->GetBinCenter(maxBinPi0);

    // Estimate background under pi0 peak
    double pi0BgLeftMin = 0.08;
    double pi0BgLeftMax = 0.1;
    double pi0BgRightMin = 0.18;
    double pi0BgRightMax = 0.2;

    int binPi0BgLeftMin = hPi0Mass->GetXaxis()->FindBin(pi0BgLeftMin);
    int binPi0BgLeftMax = hPi0Mass->GetXaxis()->FindBin(pi0BgLeftMax);
    int binPi0BgRightMin = hPi0Mass->GetXaxis()->FindBin(pi0BgRightMin);
    int binPi0BgRightMax = hPi0Mass->GetXaxis()->FindBin(pi0BgRightMax);

    double sumBgPi0 = 0.0;
    int nBinsBgPi0 = 0;

    for (int bin = binPi0BgLeftMin; bin <= binPi0BgLeftMax; ++bin) {
        sumBgPi0 += hPi0Mass->GetBinContent(bin);
        ++nBinsBgPi0;
    }
    for (int bin = binPi0BgRightMin; bin <= binPi0BgRightMax; ++bin) {
        sumBgPi0 += hPi0Mass->GetBinContent(bin);
        ++nBinsBgPi0;
    }

    double bgEstimatePi0 = (nBinsBgPi0 > 0) ? sumBgPi0 / nBinsBgPi0 : 0.0;
    double amplitudePi0Estimate = maxContentPi0 - bgEstimatePi0;
    amplitudePi0Estimate = std::max(amplitudePi0Estimate, 0.0); // Ensure non-negative amplitude

    // Sigma estimate for pi0
    double sigmaPi0Estimate = 0.012;  // Initial estimate
    double sigmaPi0Min = 0.008;
    double sigmaPi0Max = 0.02;

    // Find maximum bin in eta peak region
    int binEtaMin = hPi0Mass->GetXaxis()->FindBin(etaRangeMin);
    int binEtaMax = hPi0Mass->GetXaxis()->FindBin(etaRangeMax);

    int maxBinEta = binEtaMin;
    double maxContentEta = 0.0;

    for (int bin = binEtaMin; bin <= binEtaMax; ++bin) {
        double content = hPi0Mass->GetBinContent(bin);
        if (content > maxContentEta) {
            maxContentEta = content;
            maxBinEta = bin;
        }
    }

    double meanEtaEstimate = hPi0Mass->GetXaxis()->GetBinCenter(maxBinEta);

    // Estimate background under eta peak
    double etaBgLeftMin = 0.4;
    double etaBgLeftMax = 0.5;
    double etaBgRightMin = 0.7;
    double etaBgRightMax = 0.8;

    int binEtaBgLeftMin = hPi0Mass->GetXaxis()->FindBin(etaBgLeftMin);
    int binEtaBgLeftMax = hPi0Mass->GetXaxis()->FindBin(etaBgLeftMax);
    int binEtaBgRightMin = hPi0Mass->GetXaxis()->FindBin(etaBgRightMin);
    int binEtaBgRightMax = hPi0Mass->GetXaxis()->FindBin(etaBgRightMax);

    double sumBgEta = 0.0;
    int nBinsBgEta = 0;

    for (int bin = binEtaBgLeftMin; bin <= binEtaBgLeftMax; ++bin) {
        sumBgEta += hPi0Mass->GetBinContent(bin);
        ++nBinsBgEta;
    }
    for (int bin = binEtaBgRightMin; bin <= binEtaBgRightMax; ++bin) {
        sumBgEta += hPi0Mass->GetBinContent(bin);
        ++nBinsBgEta;
    }

    double bgEstimateEta = (nBinsBgEta > 0) ? sumBgEta / nBinsBgEta : 0.0;
    double amplitudeEtaEstimate = maxContentEta - bgEstimateEta;
    amplitudeEtaEstimate = std::max(amplitudeEtaEstimate, 0.0); // Ensure non-negative amplitude

    // Sigma estimate for eta
    double sigmaEtaEstimate = 0.03;  // Initial estimate
    double sigmaEtaMin = 0.02;
    double sigmaEtaMax = 0.05;

    // Define the totalFit function as two Gaussians (pi0 and eta) plus a fourth-order polynomial (pol4)
    totalFit = new TF1("totalFit", "gaus(0) + gaus(3) + pol4(6)", fitStart, fitEnd);
    totalFit->SetLineColor(kRed);

    // Set initial parameters for pi0 Gaussian (parameters 0-2)
    totalFit->SetParameter(0, amplitudePi0Estimate);
    totalFit->SetParameter(1, meanPi0Estimate);
    totalFit->SetParameter(2, sigmaPi0Estimate);

    // Set parameter limits for pi0 Gaussian
    totalFit->SetParLimits(0, 0, amplitudePi0Estimate * 2); // Amplitude positive and up to twice the estimate
    totalFit->SetParLimits(1, meanPi0Estimate - 0.01, meanPi0Estimate + 0.01); // Mean within ±10 MeV
    totalFit->SetParLimits(2, sigmaPi0Min, sigmaPi0Max); // Sigma within reasonable range

    // Set initial parameters for eta Gaussian (parameters 3-5)
    totalFit->SetParameter(3, amplitudeEtaEstimate);
    totalFit->SetParameter(4, meanEtaEstimate);
    totalFit->SetParameter(5, sigmaEtaEstimate);

    // Set parameter limits for eta Gaussian
    totalFit->SetParLimits(3, 0, amplitudeEtaEstimate * 2); // Amplitude positive and up to twice the estimate
    totalFit->SetParLimits(4, meanEtaEstimate - 0.02, meanEtaEstimate + 0.02); // Mean within ±20 MeV
    totalFit->SetParLimits(5, sigmaEtaMin, sigmaEtaMax); // Sigma within reasonable range

    // Initial parameters for polynomial background (parameters 6-10)
    for (int i = 6; i < 11; ++i) {
        totalFit->SetParameter(i, 0); // Start with zero coefficients
    }

    // Perform the fit
    TFitResultPtr fitResult = hPi0Mass->Fit("totalFit", "SR");

    // Check fit status and provide feedback
    int fitStatus = fitResult;
    if (fitStatus != 0) {
        std::cerr << "Fit did not converge properly. Fit status: " << fitStatus << std::endl;
    }

    // Separate Gaussian fits for pi0 and eta
    gaussPi0Fit = new TF1("gaussPi0Fit", "gaus", fitStart, fitEnd);
    gaussPi0Fit->SetParameters(totalFit->GetParameter(0), totalFit->GetParameter(1), totalFit->GetParameter(2));
    gaussPi0Fit->SetLineColor(kBlue);
    gaussPi0Fit->SetLineStyle(2);  // Dashed line for pi0

    gaussEtaFit = new TF1("gaussEtaFit", "gaus", fitStart, fitEnd);
    gaussEtaFit->SetParameters(totalFit->GetParameter(3), totalFit->GetParameter(4), totalFit->GetParameter(5));
    gaussEtaFit->SetLineColor(kGreen);
    gaussEtaFit->SetLineStyle(2);  // Dotted line for eta

    // Polynomial for the background (pol4 has 5 parameters)
    polyFit = new TF1("polyFit", "pol4", fitStart, fitEnd);
    for (int i = 6; i < 11; ++i) {  // Parameters 6 to 10 correspond to the pol4 background
        polyFit->SetParameter(i - 6, totalFit->GetParameter(i));
    }
    polyFit->SetLineColor(kOrange + 7);
    polyFit->SetLineStyle(2);

    return fitResult;
}


void printHistogramData(const std::vector<DataStructures::HistogramData>& histogramDataVector,
                        const std::string& outputDirPath,
                        const std::string& combinationName) {
    // Construct the CSV file path
    std::string csvFilePath = outputDirPath + "/InvariantMassInformation_" + combinationName + ".csv";

    // Open a CSV file for writing to the specified path
    std::ofstream csvFile(csvFilePath);
    
    // Check if the file opened successfully
    if (!csvFile.is_open()) {
        std::cerr << "Error: Could not open the file for writing at " << csvFilePath << std::endl;
        return;
    }

    // Write the CSV headers with separate columns for Ecore, Chi2, Asym, etc.
    csvFile << "Trigger Name,Ecore,Chi2,Asym,pTMin,pTMax,meanPi0,meanPi0Error,sigmaPi0,sigmaPi0Error,meanEta,meanEtaError,sigmaEta,sigmaEtaError\n";
    
    // Organize the data by cut combinations
    std::map<std::string, std::vector<DataStructures::HistogramData>> dataMap;

    // First, organize the data
    for (const auto& data : histogramDataVector) {
        std::ostringstream cutCombinationStream;
        cutCombinationStream << "E" << data.cuts.clusECore
                             << "_Chi" << data.cuts.chi
                             << "_Asym" << data.cuts.asymmetry;
        std::string cutCombination = cutCombinationStream.str();
        
        dataMap[cutCombination].push_back(data);
    }

    // Print the organized data to console and write it to CSV
    for (const auto& cutPair : dataMap) {
        // Decompose the cut combination
        float Ecore = cutPair.second[0].cuts.clusECore;
        float Chi2 = cutPair.second[0].cuts.chi;
        float Asym = cutPair.second[0].cuts.asymmetry;

        std::cout << "Cut Combination: E" << Ecore << "_Chi" << Chi2 << "_Asym" << Asym << "\n";
        
        // Print table header to console
        std::cout << "    pTMin - pTMax | meanPi0 ± error | sigmaPi0 ± error | meanEta ± error | sigmaEta ± error\n";
        std::cout << "    ----------------------------------------------------------------------------------------------------------\n";
        
        for (const auto& data : cutPair.second) {
            std::string pTBin = (data.cuts.pTMin == -1) ? "No pT bin" : Utils::formatToThreeSigFigs(data.cuts.pTMin) + " - " + Utils::formatToThreeSigFigs(data.cuts.pTMax);
            
            // Print to console
            std::cout << "    " << std::setw(14) << pTBin << " | "
                      << Utils::formatToThreeSigFigs(data.meanPi0) << " ± " << Utils::formatToThreeSigFigs(data.meanPi0Error) << " | "
                      << Utils::formatToThreeSigFigs(data.sigmaPi0) << " ± " << Utils::formatToThreeSigFigs(data.sigmaPi0Error) << " | "
                      << Utils::formatToThreeSigFigs(data.meanEta) << " ± " << Utils::formatToThreeSigFigs(data.meanEtaError) << " | "
                      << Utils::formatToThreeSigFigs(data.sigmaEta) << " ± " << Utils::formatToThreeSigFigs(data.sigmaEtaError) << "\n";

            // Write to CSV file
            csvFile << data.cuts.triggerName << ","
                    << Ecore << ","
                    << Chi2 << ","
                    << Asym << ","
                    << data.cuts.pTMin << ","
                    << data.cuts.pTMax << ","
                    << data.meanPi0 << ","
                    << data.meanPi0Error << ","
                    << data.sigmaPi0 << ","
                    << data.sigmaPi0Error << ","
                    << data.meanEta << ","
                    << data.meanEtaError << ","
                    << data.sigmaEta << ","
                    << data.sigmaEtaError << "\n";
        }
        std::cout << "\n"; // Empty line between cut combinations
    }
    std::cout << "\n"; // Empty line after all data

    // Close the CSV file
    csvFile.close();
    std::cout << "CSV data written to " << csvFilePath << std::endl;
}

void DrawInvMassCanvasText(const DataStructures::HistogramData& data, const std::string& triggerGroupName) {
    // Reconstruct pT range label
    std::ostringstream ptRangeLabel;
    if (data.cuts.pTMin != -1 && data.cuts.pTMax != -1) {
        ptRangeLabel << Utils::formatToThreeSigFigs(data.cuts.pTMin) << " - " << Utils::formatToThreeSigFigs(data.cuts.pTMax) << " GeV";
    }

    // Create two TLatex objects for the formatted output
    TLatex labelText, valueText;
    labelText.SetNDC();
    labelText.SetTextSize(0.03);
    labelText.SetTextColor(kRed);       // Set text color to red
    labelText.SetTextFont(62);          // Bold font for labels

    valueText.SetNDC();
    valueText.SetTextSize(0.03);
    valueText.SetTextColor(kBlack);     // Default color for values
    valueText.SetTextFont(42);          // Normal font for values

    
    // Add the 'Active Trigger Group' label above the current trigger label
    labelText.DrawLatex(0.2, 0.9, "Active Trigger Group:");
    valueText.DrawLatex(0.32, 0.9, triggerGroupName.c_str());

    std::string readableTriggerName = data.cuts.triggerName;
    auto triggerNameIt = TriggerConfig::triggerNameMap.find(data.cuts.triggerName);
    if (triggerNameIt != TriggerConfig::triggerNameMap.end()) {
        readableTriggerName = triggerNameIt->second;
    }

    // First column: Trigger information
    labelText.DrawLatex(0.5, 0.85, "Trigger:");
    valueText.DrawLatex(0.62, 0.85, readableTriggerName.c_str());

    labelText.DrawLatex(0.5, 0.8, "clusECore:");
    valueText.DrawLatex(0.62, 0.8, Utils::formatToThreeSigFigs(data.cuts.clusECore).c_str());

    labelText.DrawLatex(0.5, 0.75, "Chi2:");
    valueText.DrawLatex(0.62, 0.75, Utils::formatToThreeSigFigs(data.cuts.chi).c_str());

    labelText.DrawLatex(0.5, 0.7, "Asymmetry:");
    valueText.DrawLatex(0.62, 0.7, Utils::formatToThreeSigFigs(data.cuts.asymmetry).c_str());

    // If pT range is available, add it to the legend
    if (!ptRangeLabel.str().empty()) {
        labelText.DrawLatex(0.45, 0.65, "pT Range:");
        valueText.DrawLatex(0.62, 0.65, ptRangeLabel.str().c_str());
    }

    // Second column: Pi0 and Eta information
    labelText.DrawLatex(0.6, 0.55, "#pi^{0}:");
    labelText.DrawLatex(0.8, 0.55, "#eta:");

    labelText.DrawLatex(0.50, 0.5, "Mass:");
    valueText.DrawLatex(0.60, 0.5, Utils::formatToThreeSigFigs(data.meanPi0).c_str());
    valueText.DrawLatex(0.80, 0.5, Utils::formatToThreeSigFigs(data.meanEta).c_str());

    labelText.DrawLatex(0.50, 0.45, "Sigma:");
    valueText.DrawLatex(0.60, 0.45, Utils::formatToThreeSigFigs(data.sigmaPi0).c_str());
    valueText.DrawLatex(0.80, 0.45, Utils::formatToThreeSigFigs(data.sigmaEta).c_str());

    labelText.DrawLatex(0.50, 0.4, "S/B Ratio:");
    valueText.DrawLatex(0.60, 0.4, Utils::formatToThreeSigFigs(data.signalToBackgroundPi0Ratio).c_str());
    valueText.DrawLatex(0.80, 0.4, Utils::formatToThreeSigFigs(data.signalToBackgroundEtaRatio).c_str());

    // Third row: Mass ratio
    labelText.DrawLatex(0.5, 0.35, "Mass Ratio (#eta/#pi^{0}):");
    valueText.DrawLatex(0.68, 0.35, Utils::formatToThreeSigFigs(data.massRatio).c_str());
}

void ProcessInvariantMassHistograms(TFile* inputFile, const std::string& plotDirectory,
                                    const std::vector<std::string>& triggers,
                                    const std::map<std::string, int>& triggerColorMap,
                                    const std::string& combinationName) {
    // Function to parse histogram names and extract cut values, trigger name, and optional pT range
    auto parseHistName = [](const std::string& histName) -> DataStructures::CutValues {
        DataStructures::CutValues cuts;

        // Regex pattern to parse the histogram name
        std::regex re("invMass(?:_noPtBins)?_E([0-9]+(?:point[0-9]*)?)_Chi([0-9]+(?:point[0-9]*)?)_Asym([0-9]+(?:point[0-9]*)?)(?:_pT_([0-9]+(?:point[0-9]*)?)to([0-9]+(?:point[0-9]*)?))?_(.+)");
        std::smatch match;

        // Lambda function to convert strings with 'point' to float values
        auto convert = [](const std::string& input) -> float {
            std::string temp = input;
            size_t pointPos = temp.find("point");
            if (pointPos != std::string::npos) {
                temp.replace(pointPos, 5, ".");
            }
            try {
                return std::stof(temp);
            } catch (const std::exception&) {
                return 0.0f;
            }
        };

        // Check if the regex matches the histogram name
        if (std::regex_search(histName, match, re)) {
            if (match.size() >= 5) {
                cuts.clusECore = convert(match[1].str());
                cuts.chi = convert(match[2].str());
                cuts.asymmetry = convert(match[3].str());

                // Optional pT bin range (check if it was captured)
                if (match[4].matched && match[5].matched) {
                    cuts.pTMin = convert(match[4].str());
                    cuts.pTMax = convert(match[5].str());
                }

                cuts.triggerName = match[6].str();

                // Diagnostic prints
                std::cout << "Parsed histogram: " << histName << std::endl;
                std::cout << "  clusECore: " << cuts.clusECore << ", Chi: " << cuts.chi
                          << ", Asymmetry: " << cuts.asymmetry << ", pTMin: " << cuts.pTMin
                          << ", pTMax: " << cuts.pTMax << ", Trigger Name: " << cuts.triggerName << std::endl;
            }
        } else {
            std::cerr << "Error: Failed to parse histogram name: " << histName << std::endl;
        }

        return cuts;
    };

    // Vector to store histogram data for further analysis
    std::vector<DataStructures::HistogramData> histogramDataVector;
    
    std::string triggerGroupName = Utils::getTriggerCombinationName(combinationName, TriggerCombinationNames::triggerCombinationNameMap);
    
    for (const auto& trigger : triggers) {
        // Get the trigger directory
        TDirectory* triggerDir = inputFile->GetDirectory(trigger.c_str());
        if (!triggerDir) {
            std::cerr << "Trigger directory '" << trigger << "' not found in the input file. Skipping." << std::endl;
            continue;
        }

        // Iterate over histograms in the ROOT file
        TIter nextKey(triggerDir->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*)nextKey())) {
            TObject* obj = key->ReadObj();
            if (obj->InheritsFrom(TH1::Class())) {
                TH1* hist = (TH1*)obj;
                std::string histName = hist->GetName();

                // Check if the histogram is an invariant mass histogram
                if (histName.find("invMass_") == 0) {
                    // Parse the histogram name
                    DataStructures::CutValues cuts = parseHistName(histName);

                    // Check if the triggerName is in the triggers vector
                    if (std::find(triggers.begin(), triggers.end(), cuts.triggerName) == triggers.end()) {
                        // Not a trigger in this combination, skip
                        delete obj;
                        continue;
                    }

                    // Process the histogram
                    // Set axis labels
                    hist->GetXaxis()->SetTitle("M_{#gamma#gamma} [GeV]");
                    hist->GetYaxis()->SetTitle("Counts");

                    // Create the output directory structure
                    std::ostringstream cutDirStream;
                    cutDirStream << plotDirectory << "/E" << Utils::formatToThreeSigFigs(cuts.clusECore)
                                 << "_Chi" << Utils::formatToThreeSigFigs(cuts.chi)
                                 << "_Asym" << Utils::formatToThreeSigFigs(cuts.asymmetry);
                    std::string cutDirPath = cutDirStream.str();
                    gSystem->mkdir(cutDirPath.c_str(), true);  // Create cut directory

                    // If there's a pT range, create the pT folder inside the cut-specific directory
                    std::string outputDirPath = cutDirPath;
                    if (cuts.pTMin != -1 && cuts.pTMax != -1) {
                        std::ostringstream ptDirStream;
                        ptDirStream << cutDirPath << "/pT_" << Utils::formatToThreeSigFigs(cuts.pTMin)
                                    << "_to_" << Utils::formatToThreeSigFigs(cuts.pTMax);
                        outputDirPath = ptDirStream.str();
                        gSystem->mkdir(outputDirPath.c_str(), true);  // Create pT directory
                    }

                    // Construct the output file path for the histogram PNG
                    std::string outputFilePath = outputDirPath + "/" + histName + ".png";

                    // Create a canvas and draw the histogram
                    TCanvas canvas;
                    TF1* totalFit = nullptr;
                    TF1* gaussPi0Fit = nullptr;
                    TF1* gaussEtaFit = nullptr;
                    TF1* polyFit = nullptr;
                    double fitStart, fitEnd;

                    // Perform the fit
                    TFitResultPtr fitResult = PerformFitting(hist, totalFit, gaussPi0Fit, gaussEtaFit, polyFit, fitStart, fitEnd);

                    // Draw the histogram and the fits
                    hist->Draw();
                    gaussPi0Fit->Draw("SAME");
                    gaussEtaFit->Draw("SAME");
                    polyFit->Draw("SAME");
                    totalFit->Draw("SAME");

                    // Annotate the results
                    double meanPi0 = totalFit->GetParameter(1);
                    double meanPi0error = totalFit->GetParError(1);
                    double sigmaPi0 = totalFit->GetParameter(2);
                    double sigmaPi0error = totalFit->GetParError(2);
                    double meanEta = totalFit->GetParameter(4);
                    double meanEtaError = totalFit->GetParError(4);
                    double sigmaEta = totalFit->GetParameter(5);
                    double sigmaEtaError = totalFit->GetParError(5);
                    
                    // Calculate fit resolution with safeguards
                    double pi0FitResolution = std::numeric_limits<double>::quiet_NaN();
                    double pi0FitResolutionError = std::numeric_limits<double>::quiet_NaN();

                    // Debugging output for pi0 values
                    std::cout << "Calculating pi0FitResolution:" << std::endl;
                    std::cout << "  meanPi0: " << meanPi0 << ", meanPi0Error: " << meanPi0error << std::endl;
                    std::cout << "  sigmaPi0: " << sigmaPi0 << ", sigmaPi0Error: " << sigmaPi0error << std::endl;

                    if (meanPi0 > 1e-6 && sigmaPi0 > 1e-6) {  // Avoid near-zero values
                        pi0FitResolution = sigmaPi0 / meanPi0;
                        pi0FitResolutionError = pi0FitResolution * sqrt(pow(sigmaPi0error / sigmaPi0, 2) + pow(meanPi0error / meanPi0, 2));
                        std::cout << "  Calculated pi0FitResolution: " << pi0FitResolution << std::endl;
                        std::cout << "  Calculated pi0FitResolutionError: " << pi0FitResolutionError << std::endl;
                    } else {
                        std::cout << "  Skipping pi0FitResolution calculation due to small or zero values." << std::endl;
                    }

                    // Calculate eta fit resolution with safeguards
                    double etaFitResolution = std::numeric_limits<double>::quiet_NaN();
                    double etaFitResolutionError = std::numeric_limits<double>::quiet_NaN();

                    // Debugging output for eta values
                    std::cout << "Calculating etaFitResolution:" << std::endl;
                    std::cout << "  meanEta: " << meanEta << ", meanEtaError: " << meanEtaError << std::endl;
                    std::cout << "  sigmaEta: " << sigmaEta << ", sigmaEtaError: " << sigmaEtaError << std::endl;

                    if (meanEta > 1e-6 && sigmaEta > 1e-6) {  // Avoid near-zero values
                        etaFitResolution = sigmaEta / meanEta;
                        etaFitResolutionError = etaFitResolution * sqrt(pow(sigmaEtaError / sigmaEta, 2) + pow(meanEtaError / meanEta, 2));
                        std::cout << "  Calculated etaFitResolution: " << etaFitResolution << std::endl;
                        std::cout << "  Calculated etaFitResolutionError: " << etaFitResolutionError << std::endl;
                    } else {
                        std::cout << "  Skipping etaFitResolution calculation due to small or zero values." << std::endl;
                    }

                    double massRatio = meanEta / meanPi0;
                    
                    // Calculate the relative errors
                    double relativeErrorMeanEta = meanEtaError / meanEta;
                    double relativeErrorMeanPi0 = meanPi0error / meanPi0;

                    // Calculate the relative error of massRatio
                    double relativeErrorMassRatio = sqrt(
                        pow(relativeErrorMeanEta, 2) +
                        pow(relativeErrorMeanPi0, 2)
                    );

                    // Compute the absolute error of massRatio
                    double massRatioError = massRatio * relativeErrorMassRatio;
                    
                    // Clone histograms for signal and background extraction (for pi0)
                    TH1F* hSignalPi0 = (TH1F*)hist->Clone("hSignalPi0");
                    TH1F* hBackgroundPi0 = (TH1F*)hist->Clone("hBackgroundPi0");

                    // Clone histograms for signal and background extraction (for eta)
                    TH1F* hSignalEta = (TH1F*)hist->Clone("hSignalEta");
                    TH1F* hBackgroundEta = (TH1F*)hist->Clone("hBackgroundEta");

                    // Define signal and background bins for both pi0 and eta
                    int firstBinPi0 = hist->FindBin(std::max(meanPi0 - 2 * sigmaPi0, fitStart));
                    int lastBinPi0 = hist->FindBin(std::min(meanPi0 + 2 * sigmaPi0, fitEnd));

                    int firstBinEta = hist->FindBin(std::max(meanEta - 2 * sigmaEta, fitStart));
                    int lastBinEta = hist->FindBin(std::min(meanEta + 2 * sigmaEta, fitEnd));

                    // Extract signal and background for pi0
                    double binCenter, binContent, bgContent, binError;
                    for (int i = firstBinPi0; i <= lastBinPi0; ++i) {
                        binCenter = hist->GetBinCenter(i);
                        binContent = hist->GetBinContent(i);
                        binError = hist->GetBinError(i);
                        bgContent = totalFit->Eval(binCenter) - gaussPi0Fit->Eval(binCenter); // background model for pi0

                        bgContent = std::max(bgContent, 0.0);  // Ensure background is non-negative
                        hSignalPi0->SetBinContent(i, binContent - bgContent);
                        hBackgroundPi0->SetBinContent(i, bgContent);
                        hSignalPi0->SetBinError(i, binError);  // Error for signal
                        hBackgroundPi0->SetBinError(i, sqrt(bgContent));  // Error for background (Poisson statistics)
                    }

                    // Extract signal and background for eta
                    for (int i = firstBinEta; i <= lastBinEta; ++i) {
                        binCenter = hist->GetBinCenter(i);
                        binContent = hist->GetBinContent(i);
                        binError = hist->GetBinError(i);
                        bgContent = totalFit->Eval(binCenter) - gaussEtaFit->Eval(binCenter); // background model for eta

                        bgContent = std::max(bgContent, 0.0);  // Ensure background is non-negative
                        hSignalEta->SetBinContent(i, binContent - bgContent);
                        hBackgroundEta->SetBinContent(i, bgContent);
                        hSignalEta->SetBinError(i, binError);  // Error for signal
                        hBackgroundEta->SetBinError(i, sqrt(bgContent));  // Error for background (Poisson statistics)
                    }
                    std::ostringstream ptRangeLabel;
                    if (cuts.pTMin != -1 && cuts.pTMax != -1) {
                        ptRangeLabel << "pT: " << Utils::formatToThreeSigFigs(cuts.pTMin) << " - " << Utils::formatToThreeSigFigs(cuts.pTMax) << " GeV";
                    }
                    // Calculate signal and background yields and their errors (for pi0)
                    double signalPi0Yield, signalPi0Error, backgroundPi0Yield, backgroundPi0Error;
                    signalPi0Yield = hSignalPi0->IntegralAndError(firstBinPi0, lastBinPi0, signalPi0Error, "");
                    backgroundPi0Yield = hBackgroundPi0->IntegralAndError(firstBinPi0, lastBinPi0, backgroundPi0Error, "");

                    // Calculate signal and background yields and their errors (for eta)
                    double signalEtaYield, signalEtaError, backgroundEtaYield, backgroundEtaError;
                    signalEtaYield = hSignalEta->IntegralAndError(firstBinEta, lastBinEta, signalEtaError, "");
                    backgroundEtaYield = hBackgroundEta->IntegralAndError(firstBinEta, lastBinEta, backgroundEtaError, "");

                    // Calculate signal-to-background ratio for pi0
                    double signalToBackgroundPi0Ratio = backgroundPi0Yield > 0 ? signalPi0Yield / backgroundPi0Yield : 0;
                    double signalToBackgroundPi0Error = signalToBackgroundPi0Ratio > 0 ? signalToBackgroundPi0Ratio * sqrt(pow(signalPi0Error / signalPi0Yield, 2) + pow(backgroundPi0Error / backgroundPi0Yield, 2)) : 0;

                    // Calculate signal-to-background ratio for eta
                    double signalToBackgroundEtaRatio = backgroundEtaYield > 0 ? signalEtaYield / backgroundEtaYield : 0;
                    double signalToBackgroundEtaError = signalToBackgroundEtaRatio > 0 ? signalToBackgroundEtaRatio * sqrt(pow(signalEtaError / signalEtaYield, 2) + pow(backgroundEtaError / backgroundEtaYield, 2)) : 0;

                    // Get the corresponding trigger name using the helper function
                    std::string triggerNameLabel = cuts.triggerName;
                    
                    // Store the data in the HistogramData structure
                    DataStructures::HistogramData data;
                    data.cuts = cuts;
                    data.histName = histName;
                    data.meanPi0 = meanPi0;
                    data.meanPi0Error = meanPi0error;
                    data.sigmaPi0 = sigmaPi0;
                    data.sigmaPi0Error = sigmaPi0error;
                    data.meanEta = meanEta;
                    data.meanEtaError = meanEtaError;
                    data.sigmaEta = sigmaEta;
                    data.sigmaEtaError = sigmaEtaError;
                    data.massRatio = massRatio;
                    data.massRatioError = massRatioError;
                    data.signalPi0Yield = signalPi0Yield;
                    data.signalPi0Error = signalPi0Error;
                    data.backgroundPi0Yield = backgroundPi0Yield;
                    data.backgroundPi0Error = backgroundPi0Error;
                    data.signalToBackgroundPi0Ratio = signalToBackgroundPi0Ratio;
                    data.signalToBackgroundPi0Error = signalToBackgroundPi0Error;
                    data.signalEtaYield = signalEtaYield;
                    data.signalEtaError = signalEtaError;
                    data.backgroundEtaYield = backgroundEtaYield;
                    data.backgroundEtaError = backgroundEtaError;
                    data.signalToBackgroundEtaRatio = signalToBackgroundEtaRatio;
                    data.signalToBackgroundEtaError = signalToBackgroundEtaError;
                    data.pi0FitResolution = pi0FitResolution;
                    data.pi0FitResolutionError = pi0FitResolutionError;
                    data.etaFitResolution = etaFitResolution;
                    data.etaFitResolutionError = etaFitResolutionError;

                    // Add the data to the vector
                    histogramDataVector.push_back(data);

                    // Call the DrawInvMassCanvasText function
                    DrawInvMassCanvasText(data, triggerGroupName);
 

                    // Save the canvas as a PNG file
                    canvas.SaveAs(outputFilePath.c_str());
                    std::cout << "Saved: " << outputFilePath << std::endl;

                    // Clean up
                    delete totalFit;
                    delete hSignalPi0;
                    delete hBackgroundPi0;
                    delete hSignalEta;
                    delete hBackgroundEta;
                }
            }
            delete obj;
        }
    }
    printHistogramData(histogramDataVector, plotDirectory, combinationName);
}

void generateMesonPlotVsPt(
    const std::vector<double>& pTCenters,
    const std::vector<double>& meanValues,
    const std::vector<double>& meanErrors,
    const std::vector<std::string>& triggersUsed,
    const std::string& yAxisLabel,
    const std::string& outputFilePath,
    const std::string& triggerCombinationName,
    const std::string& cutCombination,
    double clusECore,
    double chi,
    double asymmetry,
    const std::map<std::string, int>& triggerColorMap,
    const std::map<std::string, std::string>& triggerNameMap,
    const std::vector<std::pair<double, double>>& pT_bins,
    double pTExclusionMax,
    double yMin = std::numeric_limits<double>::quiet_NaN(),
    double yMax = std::numeric_limits<double>::quiet_NaN()) {

    if (pTCenters.empty()) {
        return; // Nothing to plot if pTCenters is empty
    }

    // Create canvas
    TCanvas canvas("canvas", "", 800, 600);

    // Prepare bin edges for variable bin widths
    std::vector<double> binEdges;
    for (const auto& bin : pT_bins) {
        if (bin.first >= pTExclusionMax) {
            break;
        }
        binEdges.push_back(bin.first);
    }
    // Add the upper edge of the last included bin
    binEdges.push_back(pT_bins[binEdges.size() - 1].second);

    int nBins = binEdges.size() - 1;
    double* binEdgesArray = &binEdges[0];

    // Create a dummy histogram to set up the axes
    TH1F* hFrame = new TH1F("hFrame", "", nBins, binEdgesArray);
    hFrame->SetStats(0);
    hFrame->GetXaxis()->SetTitle("Leading Cluster p_{T} [GeV]");
    hFrame->GetYaxis()->SetTitle(yAxisLabel.c_str());

    // Set y-axis range
    if (std::isnan(yMin) || std::isnan(yMax)) {
        yMin = *std::min_element(meanValues.begin(), meanValues.end());
        yMax = *std::max_element(meanValues.begin(), meanValues.end());
        double yMargin = 0.05 * (yMax - yMin);
        yMin -= yMargin;
        yMax += yMargin;
    }
    hFrame->GetYaxis()->SetRangeUser(yMin, yMax);

    // Remove x-axis labels and ticks
    hFrame->GetXaxis()->SetLabelOffset(999);
    hFrame->GetXaxis()->SetTickLength(0);

    // Draw the frame
    hFrame->Draw();

    // Organize data by trigger
    std::map<std::string, std::vector<double>> triggerToPtCenters;
    std::map<std::string, std::vector<double>> triggerToMeanValues;
    std::map<std::string, std::vector<double>> triggerToMeanErrors;

    for (size_t i = 0; i < pTCenters.size(); ++i) {
        const std::string& trigger = triggersUsed[i];
        triggerToPtCenters[trigger].push_back(pTCenters[i]);
        triggerToMeanValues[trigger].push_back(meanValues[i]);
        triggerToMeanErrors[trigger].push_back(meanErrors[i]);
    }

    // Create legend
    TLegend* legend = new TLegend(0.55, 0.67, 0.88, 0.87);
    legend->SetTextSize(0.03);

    // Keep track of graphs to delete later
    std::vector<TGraphErrors*> graphs;

    // For each trigger, create a TGraphErrors and add to the canvas
    for (const auto& triggerDataPair : triggerToPtCenters) {
        const std::string& trigger = triggerDataPair.first;
        const std::vector<double>& pts = triggerDataPair.second;
        const std::vector<double>& means = triggerToMeanValues.at(trigger);
        const std::vector<double>& errors = triggerToMeanErrors.at(trigger);

        TGraphErrors* graph = new TGraphErrors();

        for (size_t i = 0; i < pts.size(); ++i) {
            double pT = pts[i];
            double mean = means[i];
            double meanError = errors[i];

            // Exclude points beyond pTExclusionMax
            if (pT >= pTExclusionMax) {
                continue;
            }

            // Find the bin corresponding to this pT
            int binIndex = -1;
            for (size_t j = 0; j < pT_bins.size(); ++j) {
                if (pT_bins[j].first >= pTExclusionMax) {
                    break;
                }
                if (pT >= pT_bins[j].first && pT < pT_bins[j].second) {
                    binIndex = j;
                    break;
                }
            }
            if (binIndex == -1) {
                continue;
            }

            // Calculate the bin center
            double binLowEdge = pT_bins[binIndex].first;
            double binUpEdge = pT_bins[binIndex].second;
            double binCenter = (binLowEdge + binUpEdge) / 2.0;

            int pointIndex = graph->GetN();
            graph->SetPoint(pointIndex, binCenter, mean);
            graph->SetPointError(pointIndex, 0, meanError); // No x errors
        }

        if (graph->GetN() == 0) {
            delete graph;
            continue;
        }

        // Set marker style and color
        int markerStyle = 20;
        int markerColor = kBlack;

        auto it_color = triggerColorMap.find(trigger);
        if (it_color != triggerColorMap.end()) {
            markerColor = it_color->second;
        }
        graph->SetMarkerStyle(markerStyle);

        graph->SetMarkerSize(1.0);

        graph->SetLineWidth(2);
        graph->SetMarkerColor(markerColor);
        graph->SetLineColor(markerColor);

        // Draw the graph
        graph->Draw("P SAME");

        // Add to legend
        std::string displayTriggerName = trigger;
        auto it_name = triggerNameMap.find(trigger);
        if (it_name != triggerNameMap.end()) {
            displayTriggerName = it_name->second;
        }
        legend->AddEntry(graph, displayTriggerName.c_str(), "p");

        // Store the graph for cleanup
        graphs.push_back(graph);
    }

    // Draw legend
    legend->Draw();

    // Draw custom x-axis ticks and labels
    double xMin = binEdges.front();
    double xMax = binEdges.back();
    double yAxisMin = hFrame->GetMinimum();
    double yAxisMax = hFrame->GetMaximum();

    double tickSize = (yAxisMax - yAxisMin) * 0.02;
    double labelOffset = (yAxisMax - yAxisMin) * 0.05;
    TLatex latex;
    latex.SetTextSize(0.035);
    latex.SetTextAlign(22); // Center alignment

    // Draw x-axis line
    TLine xAxisLine(xMin, yAxisMin, xMax, yAxisMin);
    xAxisLine.Draw();

    // Draw ticks and labels at bin edges
    for (size_t i = 0; i < binEdges.size(); ++i) {
        double xPos = binEdges[i];
        double yPos = yAxisMin;

        // Draw tick
        TLine* tick = new TLine(xPos, yPos, xPos, yPos - tickSize);
        tick->Draw();

        // Get pT value for label
        double pTValue = binEdges[i];

        // **Format label to show one decimal place**
        std::ostringstream labelStream;
        labelStream << std::fixed << std::setprecision(1) << pTValue;
        std::string label = labelStream.str();

        // Draw label
        latex.DrawLatex(xPos, yPos - labelOffset, label.c_str());
    }

    // Redraw the axes to ensure labels are on top
    canvas.RedrawAxis();

    // Add trigger and cut information on the plot
    TLatex labelText;
    labelText.SetNDC();
    labelText.SetTextSize(0.032);

    TLatex valueText;
    valueText.SetNDC();
    valueText.SetTextSize(0.028);

    // Use the utility function to get the display trigger group name
    std::string displayTriggerGroupName = Utils::getTriggerCombinationName(triggerCombinationName, TriggerCombinationNames::triggerCombinationNameMap);

    labelText.DrawLatex(0.2, 0.9, "#font[62]{Active Trigger Group:}");
    valueText.DrawLatex(0.42, 0.9, displayTriggerGroupName.c_str());

    labelText.DrawLatex(0.2, 0.83, "#font[62]{ECore #geq}");
    std::ostringstream eCoreWithUnit;
    eCoreWithUnit << clusECore << "   GeV";
    valueText.DrawLatex(0.42, 0.83, eCoreWithUnit.str().c_str());

    labelText.DrawLatex(0.2, 0.76, "#font[62]{#chi^{2} <}");
    std::ostringstream chiStr;
    chiStr << chi;
    valueText.DrawLatex(0.42, 0.76, chiStr.str().c_str());

    labelText.DrawLatex(0.2, 0.69, "#font[62]{Asymmetry <}");
    std::ostringstream asymmetryStr;
    asymmetryStr << asymmetry;
    valueText.DrawLatex(0.42, 0.69, asymmetryStr.str().c_str());

    // Ensure the directory exists
    std::string outputDirPath = outputFilePath.substr(0, outputFilePath.find_last_of("/"));
    gSystem->mkdir(outputDirPath.c_str(), true);

    // Save plot
    canvas.SaveAs(outputFilePath.c_str());
    std::cout << "Saved plot: " << outputFilePath << std::endl;

    // Clean up
    for (auto graph : graphs) {
        delete graph;
    }
    delete legend;
    delete hFrame;
}

void readHistogramDataForInvariantMasses(const std::string& csvFilePath,
                       std::vector<DataStructures::HistogramData>& dataVector,
                       std::set<std::string>& triggersInDataFile) {
    // Open the CSV file
    std::ifstream csvFile(csvFilePath);
    if (!csvFile.is_open()) {
        std::cerr << "Error: Could not open CSV file " << csvFilePath << std::endl;
        return;
    }

    // Read header line
    std::string line;
    std::getline(csvFile, line); // Skip header

    // Read data lines
    while (std::getline(csvFile, line)) {
        std::istringstream iss(line);
        std::string token;

        DataStructures::HistogramData data;

        // Read fields
        std::getline(iss, token, ',');
        data.cuts.triggerName = token;

        std::getline(iss, token, ',');
        data.cuts.clusECore = std::stof(token);

        std::getline(iss, token, ',');
        data.cuts.chi = std::stof(token);

        std::getline(iss, token, ',');
        data.cuts.asymmetry = std::stof(token);

        std::getline(iss, token, ',');
        data.cuts.pTMin = std::stof(token);

        std::getline(iss, token, ',');
        data.cuts.pTMax = std::stof(token);

        std::getline(iss, token, ',');
        data.meanPi0 = std::stod(token);

        std::getline(iss, token, ',');
        data.meanPi0Error = std::stod(token);

        std::getline(iss, token, ',');
        data.sigmaPi0 = std::stod(token);

        std::getline(iss, token, ',');
        data.sigmaPi0Error = std::stod(token);

        std::getline(iss, token, ',');
        data.meanEta = std::stod(token);

        std::getline(iss, token, ',');
        data.meanEtaError = std::stod(token);

        std::getline(iss, token, ',');
        data.sigmaEta = std::stod(token);

        std::getline(iss, token, ',');
        data.sigmaEtaError = std::stod(token);

        triggersInDataFile.insert(data.cuts.triggerName);  // Collect trigger names
        dataVector.push_back(data);
    }

    csvFile.close();
}

DataStructures::CutCombinationData processCutCombination(
    const std::vector<DataStructures::HistogramData>& dataList,
    const std::vector<std::string>& triggers,
    const std::map<std::string, double>& triggerEfficiencyPoints,
    const std::map<std::string, double>& triggerThresholds,
    double pTExclusionMax) {
    DataStructures::CutCombinationData result;

    // Function to extract photon threshold from trigger name if not in the map
    auto extractPhotonThreshold = [](const std::string& triggerName) -> double {
        std::regex re("Photon_(\\d+)_GeV");
        std::smatch match;
        if (std::regex_search(triggerName, match, re)) {
            if (match.size() >= 2) {
                return std::stod(match[1]);
            }
        }
        return 0.0; // Default to 0 if parsing fails
    };

    // Initialize the cuts
    result.clusECore = dataList[0].cuts.clusECore;
    result.chi = dataList[0].cuts.chi;
    result.asymmetry = dataList[0].cuts.asymmetry;

    // Organize data per pT bin and trigger
    std::map<std::pair<double, double>, std::map<std::string, DataStructures::HistogramData>> dataPerPtBin;

    for (const auto& data : dataList) {
        if (data.cuts.pTMin == -1 || data.cuts.pTMax == -1) {
            continue; // Skip data without pT bins
        }

        // Exclude pT bins beyond the exclusion maximum
        if (data.cuts.pTMin >= pTExclusionMax) {
            continue;
        }

        std::pair<double, double> pTBin = std::make_pair(data.cuts.pTMin, data.cuts.pTMax);
        dataPerPtBin[pTBin][data.cuts.triggerName] = data;

        // Debugging output
        std::cout << BLUE << "Added data for pT bin [" << data.cuts.pTMin << ", " << data.cuts.pTMax << "], "
                  << "trigger: " << data.cuts.triggerName << RESET << std::endl;
    }

    // Process each pT bin
    for (const auto& ptBinData : dataPerPtBin) {
        double pTMin = ptBinData.first.first;
        double pTMax = ptBinData.first.second;
        double pTCenter = (pTMin + pTMax) / 2.0;

        const auto& triggerDataMap = ptBinData.second;

        // Decide which trigger to use for the main plot
        std::string triggerToUse = "MBD_NandS_geq_1";
        double maxPhotonThreshold = 0.0;

        // List to keep track of efficient triggers at this pT bin
        std::vector<std::string> efficientTriggers;

        for (const auto& photonTrigger : triggers) {
            if (photonTrigger != "MBD_NandS_geq_1") {
                auto it = triggerEfficiencyPoints.find(photonTrigger);
                if (it != triggerEfficiencyPoints.end()) {
                    double x99 = it->second;
                    // Trigger is efficient if x99 <= pTMax
                    if (x99 <= pTMax) {
                        efficientTriggers.push_back(photonTrigger);
                    }
                }
            }
        }

        // Select the trigger with the highest photon threshold among efficient triggers
        for (const auto& efficientTrigger : efficientTriggers) {
            double photonThreshold = 0.0;
            auto thresholdIt = triggerThresholds.find(efficientTrigger);
            if (thresholdIt != triggerThresholds.end()) {
                photonThreshold = thresholdIt->second;
            } else {
                // Extract threshold from trigger name if not in the map
                photonThreshold = extractPhotonThreshold(efficientTrigger);
            }

            if (photonThreshold > maxPhotonThreshold) {
                maxPhotonThreshold = photonThreshold;
                triggerToUse = efficientTrigger;
            } else if (photonThreshold == maxPhotonThreshold) {
                // If thresholds are equal, prefer the one with lower x99 (more efficient)
                double x99_current = triggerEfficiencyPoints.at(triggerToUse);
                double x99_candidate = triggerEfficiencyPoints.at(efficientTrigger);
                if (x99_candidate < x99_current) {
                    triggerToUse = efficientTrigger;
                }
            }
        }
        // Debugging output
        std::cout << CYAN << "Processing pT bin [" << pTMin << ", " << pTMax << "], pTCenter: " << pTCenter << RESET << std::endl;
        std::cout << CYAN << "Efficient triggers at this pT:" << RESET << std::endl;
        for (const auto& etrig : efficientTriggers) {
            std::cout << CYAN << "  - " << etrig << RESET << std::endl;
        }
        std::cout << CYAN << "Selected trigger to use: " << triggerToUse << RESET << std::endl;

        // For the main plot, check if data from the selected trigger is available
        auto dataIt = triggerDataMap.find(triggerToUse);
        if (dataIt != triggerDataMap.end()) {
            const DataStructures::HistogramData& selectedData = dataIt->second;

            // Validate meanPi0 and meanEta before adding
            bool validPi0 = selectedData.meanPi0 > 0.0 && selectedData.meanPi0 < 0.4;
            bool validEta = selectedData.meanEta > 0.3 && selectedData.meanEta < 1.0;

            // Debugging output
            std::cout << GREEN << "Data found for trigger " << triggerToUse << " at pT bin [" << pTMin << ", " << pTMax << "]" << RESET << std::endl;
            if (validPi0) {
                std::cout << GREEN << "  Valid Pi0 mean: " << selectedData.meanPi0 << RESET << std::endl;
            } else {
                std::cout << YELLOW << "  Invalid Pi0 mean: " << selectedData.meanPi0 << RESET << std::endl;
            }

            if (validEta) {
                std::cout << GREEN << "  Valid Eta mean: " << selectedData.meanEta << RESET << std::endl;
            } else {
                std::cout << YELLOW << "  Invalid Eta mean: " << selectedData.meanEta << RESET << std::endl;
            }

            if (validPi0) {
                result.pTCentersPi0.push_back(pTCenter);
                result.meanPi0Values.push_back(selectedData.meanPi0);
                result.meanPi0Errors.push_back(selectedData.meanPi0Error);
                result.triggersUsedPi0.push_back(triggerToUse); // Track the trigger used
            }

            if (validEta) {
                result.pTCentersEta.push_back(pTCenter);
                result.meanEtaValues.push_back(selectedData.meanEta);
                result.meanEtaErrors.push_back(selectedData.meanEtaError);
                result.triggersUsedEta.push_back(triggerToUse); // Track the trigger used
            }
        } else {
            std::cout << RED << "Warning: Data not found for pT bin [" << pTMin << ", " << pTMax << "] and trigger " << triggerToUse << RESET << std::endl;
            // Debugging output: List available triggers for this pT bin
            std::cout << YELLOW << "Available triggers for this pT bin:" << RESET << std::endl;
            for (const auto& tDataPair : triggerDataMap) {
                std::cout << YELLOW << "  - " << tDataPair.first << RESET << std::endl;
            }
        }

        // Collect triggers for overlay plots
        for (const auto& tDataPair : triggerDataMap) {
            result.triggersInData.insert(tDataPair.first);
        }

        // Collect data for overlay plots
        for (const auto& tDataPair : triggerDataMap) {
            const std::string& triggerName = tDataPair.first;
            const DataStructures::HistogramData& data = tDataPair.second;

            // Validate meanPi0 and meanEta before adding
            bool validPi0 = data.meanPi0 > 0.0 && data.meanPi0 < 0.4;
            bool validEta = data.meanEta > 0.3 && data.meanEta < 1.0;

            if (validPi0) {
                result.triggerToPtCentersPi0[triggerName].push_back(pTCenter);
                result.triggerToMeanPi0Values[triggerName].push_back(data.meanPi0);
                result.triggerToMeanPi0Errors[triggerName].push_back(data.meanPi0Error);
            }

            if (validEta) {
                result.triggerToPtCentersEta[triggerName].push_back(pTCenter);
                result.triggerToMeanEtaValues[triggerName].push_back(data.meanEta);
                result.triggerToMeanEtaErrors[triggerName].push_back(data.meanEtaError);
            }
        }
    }

    return result;
}

void generateOverlayInvariantMassPlot(
    const std::map<std::string, std::vector<double>>& triggerToPtCenters,
    const std::map<std::string, std::vector<double>>& triggerToMeanValues,
    const std::map<std::string, std::vector<double>>& triggerToMeanErrors,
    const std::set<std::string>& triggersInData,
    const std::vector<std::pair<double, double>>& pT_bins,
    const std::string& mesonName, // e.g., "#pi^{0}" or "#eta"
    const std::string& yAxisTitle, // e.g., "Mean #pi^{0} Mass [GeV]"
    const double yBufferFraction, // e.g., 0.1 for 10% buffer
    const double pTExclusionMax,
    const std::string& outputFilePath,
    const std::string& plotDirectory,
    const std::string& cutCombination,
    const DataStructures::CutCombinationData& cutData,
    const std::map<std::string, int>& triggerColorMap,
    const std::map<std::string, std::string>& triggerNameMap) {

    if (triggerToPtCenters.empty()) {
        std::cout << "No data to plot." << std::endl;
        return;
    }

    // Create canvas
    TCanvas canvas("canvas", "", 800, 600);

    // Prepare bin edges for variable bin widths
    std::vector<double> binEdges;
    for (const auto& bin : pT_bins) {
        if (bin.first >= pTExclusionMax) {
            break;
        }
        binEdges.push_back(bin.first);
    }
    // Add the upper edge of the last included bin
    binEdges.push_back(pT_bins[binEdges.size() - 1].second);

    int nBins = binEdges.size() - 1;
    double* binEdgesArray = &binEdges[0];

    // Create a dummy histogram to set up the axes
    TH1F* hFrame = new TH1F("hFrame", "", nBins, binEdgesArray);
    hFrame->SetStats(0);
    hFrame->GetXaxis()->SetTitle("Leading Cluster p_{T} [GeV]");
    hFrame->GetYaxis()->SetTitle(yAxisTitle.c_str());

    // Remove x-axis labels and ticks
    hFrame->GetXaxis()->SetLabelOffset(999);
    hFrame->GetXaxis()->SetTickLength(0);

    // Draw the frame
    hFrame->Draw();

    // Declare legend pointer outside the conditional blocks
    TLegend* legend = nullptr;

    if (mesonName == "#pi^{0}") {
        // Create legend with specific coordinates for #pi^{0}
        legend = new TLegend(0.2, 0.2, 0.4, 0.4);
        legend->SetTextSize(0.03);
    } else {
        // Create legend with specific coordinates for other meson names
        legend = new TLegend(0.2, 0.2, 0.4, 0.4);
        legend->SetTextSize(0.03);
    }

    // Keep track of graphs to delete later
    std::vector<TGraphErrors*> graphs;

    // Variables to adjust y-axis range
    double yMinData = std::numeric_limits<double>::max();
    double yMaxData = std::numeric_limits<double>::lowest();

    // Convert triggersInData to a vector for indexing
    std::vector<std::string> triggerNames(triggersInData.begin(), triggersInData.end());
    int numTriggers = triggerNames.size();

    /*
     ADJUST OFFSET VALUE AS NECESSARY
     */
    double offsetValue;
    if (mesonName == "#eta") {
        offsetValue = 0.098;
    } else {
        offsetValue = 0.08;
    }

    // Calculate offsets based on the number of triggers
    std::vector<double> offsets;
    if (numTriggers == 1) {
        offsets.push_back(0.0);
    } else if (numTriggers == 2) {
        offsets.push_back(-0.5 * offsetValue);
        offsets.push_back(+0.5 * offsetValue);
    } else if (numTriggers == 3) {
        offsets.push_back(-0.5 * offsetValue);
        offsets.push_back(0.0);
        offsets.push_back(+0.5 * offsetValue);
    } else if (numTriggers == 4) {
        offsets.push_back(-1.5 * offsetValue);
        offsets.push_back(-0.5 * offsetValue);
        offsets.push_back(+0.5 * offsetValue);
        offsets.push_back(+1.5 * offsetValue);
    } else {
        // For N > 4, distribute offsets symmetrically around zero
        int midIndex = numTriggers / 2;
        for (int i = 0; i < numTriggers; ++i) {
            double offset = (i - midIndex) * offsetValue;
            // Adjust for even number of triggers
            if (numTriggers % 2 == 0) {
                offset += offsetValue / 2.0;
            }
            offsets.push_back(offset);
        }
    }

    // Loop over triggers with indexing
    for (size_t triggerIndex = 0; triggerIndex < triggerNames.size(); ++triggerIndex) {
        const std::string& triggerName = triggerNames[triggerIndex];
        double xOffset = offsets[triggerIndex];

        // Get the data
        auto it_pT = triggerToPtCenters.find(triggerName);
        auto it_meanValues = triggerToMeanValues.find(triggerName);
        auto it_meanErrors = triggerToMeanErrors.find(triggerName);

        if (it_pT == triggerToPtCenters.end() || it_meanValues == triggerToMeanValues.end() || it_meanErrors == triggerToMeanErrors.end()) {
            continue;
        }

        const std::vector<double>& pTCenters = it_pT->second;
        const std::vector<double>& meanValues = it_meanValues->second;
        const std::vector<double>& meanErrors = it_meanErrors->second;

        if (pTCenters.empty()) {
            continue; // Skip if no data for this trigger
        }

        TGraphErrors* graph = new TGraphErrors();

        for (size_t i = 0; i < pTCenters.size(); ++i) {
            double pT = pTCenters[i];
            double mean = meanValues[i];
            double meanError = meanErrors[i];

            // Exclude points beyond pTExclusionMax
            if (pT >= pTExclusionMax) {
                continue;
            }

            // Find the bin corresponding to this pT
            int binIndex = -1;
            for (size_t j = 0; j < pT_bins.size(); ++j) {
                if (pT_bins[j].first >= pTExclusionMax) {
                    break;
                }
                if (pT >= pT_bins[j].first && pT < pT_bins[j].second) {
                    binIndex = j;
                    break;
                }
            }
            if (binIndex == -1) {
                continue;
            }

            // Calculate the bin center
            double binLowEdge = pT_bins[binIndex].first;
            double binUpEdge = pT_bins[binIndex].second;
            double binCenter = (binLowEdge + binUpEdge) / 2.0;

            // Apply xOffset
            double xPos = binCenter + xOffset;

            int pointIndex = graph->GetN();
            graph->SetPoint(pointIndex, xPos, mean);
            graph->SetPointError(pointIndex, 0, meanError); // No x errors

            // Update y-axis range
            if (mean - meanError < yMinData) yMinData = mean - meanError;
            if (mean + meanError > yMaxData) yMaxData = mean + meanError;
        }

        if (graph->GetN() == 0) {
            delete graph;
            continue;
        }

        // Set marker style and color
        int markerStyle = 20;
        int markerColor = kBlack;
        auto it_color = triggerColorMap.find(triggerName);
        if (it_color != triggerColorMap.end()) {
            markerColor = it_color->second;
        }
        graph->SetMarkerStyle(markerStyle);
        double markerSize;
        if (mesonName == "#eta") {
            markerSize = 0.72;
        } else {
            markerSize = 0.85;
        }
        graph->SetMarkerSize(markerSize); // Adjust marker size if needed
        graph->SetLineWidth(2);
        graph->SetMarkerColor(markerColor);
        graph->SetLineColor(markerColor);

        // Draw the graph
        graph->Draw("P SAME");

        // Add to legend
        std::string displayTriggerName = triggerName;
        if (triggerNameMap.find(triggerName) != triggerNameMap.end()) {
            displayTriggerName = triggerNameMap.at(triggerName);
        }
        legend->AddEntry(graph, displayTriggerName.c_str(), "p");

        // Store the graph for cleanup
        graphs.push_back(graph);
    }

    // Adjust y-axis range with buffer
    if (yMinData >= yMaxData) {
        // If no valid data points were found, set default y-axis range
        yMinData = 0.0;
        yMaxData = 1.0;
    }
    double yRange = yMaxData - yMinData;
    double yBuffer = yRange * yBufferFraction;
    double yMin = yMinData - yBuffer;
    double yMax = yMaxData + yBuffer;

    if (mesonName == "#eta") {
        hFrame->GetYaxis()->SetRangeUser(0.4, 0.75);
    } else {
        hFrame->GetYaxis()->SetRangeUser(yMin, yMax);
    }

    // Draw legend
    legend->Draw();

    // Draw custom x-axis ticks and labels
    double xMin = binEdges.front();
    double xMax = binEdges.back();
    double yAxisMin = hFrame->GetMinimum();
    double yAxisMax = hFrame->GetMaximum();

    double tickSize = (yAxisMax - yAxisMin) * 0.02; // Adjust as needed
    double labelOffset = (yAxisMax - yAxisMin) * 0.05; // Adjust as needed
    TLatex latex;
    latex.SetTextSize(0.035);
    latex.SetTextAlign(22); // Center alignment

    // Draw x-axis line
    TLine xAxisLine(xMin, yAxisMin, xMax, yAxisMin);
    xAxisLine.Draw();

    // Draw ticks and labels at bin edges
    for (size_t i = 0; i < binEdges.size(); ++i) {
        double xPos = binEdges[i];
        double yPos = yAxisMin;

        // Draw tick
        TLine* tick = new TLine(xPos, yPos, xPos, yPos - tickSize);
        tick->Draw();

        // Get pT value for label
        double pTValue = binEdges[i];

        // Format label to show one decimal place
        std::ostringstream labelStream;
        labelStream << std::fixed << std::setprecision(1) << pTValue;
        std::string label = labelStream.str();

        // Draw label
        latex.DrawLatex(xPos, yPos - labelOffset, label.c_str());
    }

    // Redraw the axes to ensure labels are on top
    canvas.RedrawAxis();

    // Add cut combination information to the canvas
    TLatex labelText;
    labelText.SetNDC();
    labelText.SetTextSize(0.032);

    TLatex valueText;
    valueText.SetNDC();
    valueText.SetTextSize(0.028);

    // Use cutData to get the values
    labelText.DrawLatex(0.2, 0.87, "#font[62]{ECore #geq}");
    std::ostringstream eCoreWithUnit;
    eCoreWithUnit << cutData.clusECore << "   GeV";
    valueText.DrawLatex(0.32, 0.87, eCoreWithUnit.str().c_str());

    labelText.DrawLatex(0.2, 0.82, "#font[62]{#chi^{2} <}");
    std::ostringstream chiStr;
    chiStr << cutData.chi;
    valueText.DrawLatex(0.27, 0.82, chiStr.str().c_str());

    labelText.DrawLatex(0.2, 0.77, "#font[62]{Asymmetry <}");
    std::ostringstream asymmetryStr;
    asymmetryStr << cutData.asymmetry;
    valueText.DrawLatex(0.37, 0.77, asymmetryStr.str().c_str());

    // Ensure the directory exists
    gSystem->mkdir((plotDirectory + "/" + cutCombination).c_str(), true);

    // Save plot
    canvas.SaveAs(outputFilePath.c_str());
    std::cout << "Saved overlay plot: " << outputFilePath << std::endl;

    // Clean up
    for (auto graph : graphs) {
        delete graph;
    }
    delete legend;
    delete hFrame;
}

void ProcessMesonMassVsPt(
    const std::string& plotDirectory,
    const std::string& combinationName,
    const std::vector<std::string>& triggers,
    const std::map<std::string, double>& triggerEfficiencyPoints,
    const std::vector<std::pair<double, double>>& pT_bins,
    double pTExclusionMax = std::numeric_limits<double>::infinity(),
    double yBufferFractionPi0 = 0.1,
    double yBufferFractionEta = 0.1) {
    // Construct the CSV file path
    std::string csvFilePath = plotDirectory + "/InvariantMassInformation_" + combinationName + ".csv";

    std::vector<DataStructures::HistogramData> dataVector;
    std::set<std::string> triggersInDataFile;

    // Read the CSV data
    readHistogramDataForInvariantMasses(csvFilePath, dataVector, triggersInDataFile);

    // Check if dataVector is empty
    if (dataVector.empty()) {
        std::cerr << "Error: No data read from CSV file." << std::endl;
        return;
    }

    // Organize data by cut combinations
    std::map<std::string, std::vector<DataStructures::HistogramData>> dataByCuts;

    for (const auto& data : dataVector) {
        std::ostringstream cutCombinationStream;
        cutCombinationStream << "E" << Utils::formatToThreeSigFigs(data.cuts.clusECore)
                             << "_Chi" << Utils::formatToThreeSigFigs(data.cuts.chi)
                             << "_Asym" << Utils::formatToThreeSigFigs(data.cuts.asymmetry);
        std::string cutCombination = cutCombinationStream.str();

        dataByCuts[cutCombination].push_back(data);
    }

    for (const auto& cutPair : dataByCuts) {
        const std::string& cutCombination = cutPair.first;
        const std::vector<DataStructures::HistogramData>& dataList = cutPair.second;

        // Process data for this cut combination
        DataStructures::CutCombinationData cutData = processCutCombination(
            dataList,
            triggers,
            triggerEfficiencyPoints,
            TriggerConfig::triggerThresholds,
            pTExclusionMax
        );

        // Generate π⁰ mass vs pT plot if we have valid data
        if (!cutData.pTCentersPi0.empty()) {
            std::string outputFilePathPi0 = plotDirectory + "/" + cutCombination + "/meanPi0_vs_pT.png";
            generateMesonPlotVsPt(
                cutData.pTCentersPi0,
                cutData.meanPi0Values,
                cutData.meanPi0Errors,
                cutData.triggersUsedPi0,
                "Mean #pi^{0} Mass [GeV]",
                outputFilePathPi0,
                combinationName,
                cutCombination,
                cutData.clusECore,
                cutData.chi,
                cutData.asymmetry,
                TriggerConfig::triggerColorMap,
                TriggerConfig::triggerNameMap,
                pT_bins,
                8.0,  // pTExclusionMax for π⁰ plots
                0.13, 0.35 // yMin and yMax
            );
        } else {
            std::cout << "No valid π⁰ data to plot for cut combination " << cutCombination << std::endl;
        }

        // Generate η mass vs pT plot if we have valid data
        if (!cutData.pTCentersEta.empty()) {
            std::string outputFilePathEta = plotDirectory + "/" + cutCombination + "/meanEta_vs_pT.png";
            generateMesonPlotVsPt(
                cutData.pTCentersEta,
                cutData.meanEtaValues,
                cutData.meanEtaErrors,
                cutData.triggersUsedEta,
                "Mean #eta Mass [GeV]",
                outputFilePathEta,
                combinationName,
                cutCombination,
                cutData.clusECore,
                cutData.chi,
                cutData.asymmetry,
                TriggerConfig::triggerColorMap,
                TriggerConfig::triggerNameMap,
                pT_bins,
                20.0, // pTExclusionMax for η plots
                0.45, 0.75 // yMin and yMax
            );
        } else {
            std::cout << "No valid η data to plot for cut combination " << cutCombination << std::endl;
        }
        // Proceed with the overlay plots using cutData.triggerToPtCentersPi0, etc.

        if (!cutData.triggerToPtCentersPi0.empty()) {
            double pTExclusionMaxPi0 = 8.0; // For π⁰
            std::string outputFilePath = plotDirectory + "/" + cutCombination + "/meanPi0_vs_pT_Overlay.png";

            generateOverlayInvariantMassPlot(
                cutData.triggerToPtCentersPi0,
                cutData.triggerToMeanPi0Values,
                cutData.triggerToMeanPi0Errors,
                cutData.triggersInData,
                pT_bins,
                "#pi^{0}",
                "Mean #pi^{0} Mass [GeV]",
                yBufferFractionPi0,
                pTExclusionMaxPi0,
                outputFilePath,
                plotDirectory,
                cutCombination,
                cutData,
                TriggerConfig::triggerColorMap,
                TriggerConfig::triggerNameMap
            );
        } else {
            std::cout << "No valid π⁰ data to plot for cut combination " << cutCombination << std::endl;
        }

        // Similarly for η meson plots

        if (!cutData.triggerToPtCentersEta.empty()) {
            double pTExclusionMaxEta = 20.0; // For η
            std::string outputFilePath = plotDirectory + "/" + cutCombination + "/meanEta_vs_pT_Overlay.png";

            generateOverlayInvariantMassPlot(
                cutData.triggerToPtCentersEta,
                cutData.triggerToMeanEtaValues,
                cutData.triggerToMeanEtaErrors,
                cutData.triggersInData,
                pT_bins,
                "#eta",
                "Mean #eta Mass [GeV]",
                yBufferFractionEta,
                pTExclusionMaxEta,
                outputFilePath,
                plotDirectory,
                cutCombination,
                cutData,
                TriggerConfig::triggerColorMap,
                TriggerConfig::triggerNameMap
            );
        } else {
            std::cout << "No valid η data to plot for cut combination " << cutCombination << std::endl;
        }
    }
}


void AddLabelsToCanvas_isoHistsWithCuts(
    const DataStructures::CutValues& cuts,
    const std::string& massWindowLabel,
    const std::string& triggerName,
    const std::string& triggerGroupName,
    float pTMin,
    float pTMax,
    const std::string& histType,
    bool hasIsoEtRange,
    float isoEtMin,
    float isoEtMax) {

    // Map triggerGroupName and triggerName to human-readable names
    std::string readableTriggerGroupName = Utils::getTriggerCombinationName(
        triggerGroupName, TriggerCombinationNames::triggerCombinationNameMap);

    std::string readableTriggerName = triggerName;
    auto triggerNameIt = TriggerConfig::triggerNameMap.find(triggerName);
    if (triggerNameIt != TriggerConfig::triggerNameMap.end()) {
        readableTriggerName = triggerNameIt->second;
    }

    TLatex labelText, valueText;
    labelText.SetNDC();
    labelText.SetTextSize(0.027);
    labelText.SetTextColor(kRed);
    labelText.SetTextFont(62);

    valueText.SetNDC();
    valueText.SetTextSize(0.027);
    valueText.SetTextColor(kBlack);
    valueText.SetTextFont(42);

    labelText.DrawLatex(0.58, 0.9, "Active Trigger Group:");
    valueText.DrawLatex(0.77, 0.9, readableTriggerGroupName.c_str());

    labelText.DrawLatex(0.58, 0.85, "Trigger:");
    valueText.DrawLatex(0.7, 0.85, readableTriggerName.c_str());

    labelText.DrawLatex(0.58, 0.8, "p_{T} Range:");
    std::string pTRangeText = Utils::formatToThreeSigFigs(pTMin) + " to " + Utils::formatToThreeSigFigs(pTMax) + " GeV";
    valueText.DrawLatex(0.7, 0.8, pTRangeText.c_str());

    if (hasIsoEtRange) {
        labelText.DrawLatex(0.58, 0.75, "isoE_{T} Range:");
        std::string isoEtRangeText = Utils::formatToThreeSigFigs(isoEtMin) + " to " + Utils::formatToThreeSigFigs(isoEtMax) + " GeV";
        valueText.DrawLatex(0.7, 0.75, isoEtRangeText.c_str());
    }

    labelText.DrawLatex(0.58, hasIsoEtRange ? 0.7 : 0.75, "Cuts:");
    std::string cutsText = "E > " + Utils::formatToThreeSigFigs(cuts.clusECore) +
                           " GeV, Chi^{2} < " + Utils::formatToThreeSigFigs(cuts.chi) +
                           ", Asym < " + Utils::formatToThreeSigFigs(cuts.asymmetry);
    valueText.DrawLatex(0.7, hasIsoEtRange ? 0.7 : 0.75, cutsText.c_str());

    if (!massWindowLabel.empty()) {
        labelText.DrawLatex(0.58, hasIsoEtRange ? 0.65 : 0.7, "Mass Window:");
        valueText.DrawLatex(0.7, hasIsoEtRange ? 0.65 : 0.7, massWindowLabel.c_str());
    }

    labelText.DrawLatex(0.58, hasIsoEtRange ? 0.6 : 0.65, "Histogram Type:");
    valueText.DrawLatex(0.77, hasIsoEtRange ? 0.6 : 0.65, histType.c_str());
}

void ProcessIsolationEnergyHistogramsWithCuts(
    TFile* inputFile,
    const std::string& plotDirectory,
    const std::vector<std::string>& triggers,
    const std::string& combinationName) {
    
    // Use the raw combinationName as triggerGroupName in map keys
    std::string triggerGroupName = combinationName;

    auto parseIsolationHistNameWithCuts = [](const std::string& histName)
        -> std::tuple<DataStructures::CutValues, std::string, float, float, std::string, std::string, bool, float, float> {
        DataStructures::CutValues cuts;
        std::string massWindowLabel;
        float pTMin = -1;
        float pTMax = -1;
        std::string triggerName;
        std::string histType;
        bool hasIsoEtRange = false;
        float isoEtMin = 0.0f;
        float isoEtMax = 0.0f;

        std::smatch match;

        // Lambda to convert 'point' notation to float, handling negative numbers
        auto convert = [](const std::string& input) -> float {
            std::string temp = input;
            size_t pointPos = temp.find("point");
            if (pointPos != std::string::npos) {
                temp.replace(pointPos, 5, ".");
            }
            try {
                return std::stof(temp);
            } catch (const std::exception&) {
                return 0.0f;
            }
        };

        std::regex histPattern(
            "(h[12]_cluster_iso_Et|allPhotonCount|ptPhoton|isolatedPhotonCount)"
            "_E(-?[0-9]+(?:point[0-9]*)?)"
            "_Chi(-?[0-9]+(?:point[0-9]*)?)"
            "_Asym(-?[0-9]+(?:point[0-9]*)?)"
            "(?:(_inMassWindow|_outsideMassWindow))?"
            "(?:_isoEt_(-?[0-9]+(?:point[0-9]*)?)to(-?[0-9]+(?:point[0-9]*)?))?"
            "_pT_(-?[0-9]+(?:point[0-9]*)?)to(-?[0-9]+(?:point[0-9]*)?)"
            "_([^ ]+)"
        );

        if (std::regex_match(histName, match, histPattern)) {
            if (match.size() >= 10) {
                histType = match[1].str();
                cuts.clusECore = convert(match[2].str());
                cuts.chi = convert(match[3].str());
                cuts.asymmetry = convert(match[4].str());
                massWindowLabel = match[5].str(); // may be empty

                // Remove leading underscore from massWindowLabel if present
                if (!massWindowLabel.empty() && massWindowLabel[0] == '_') {
                    massWindowLabel = massWindowLabel.substr(1);
                }
                
                std::string isoEtMinStr = match[6].str();
                std::string isoEtMaxStr = match[7].str();
                if (!isoEtMinStr.empty() && !isoEtMaxStr.empty()) {
                    isoEtMin = convert(isoEtMinStr);
                    isoEtMax = convert(isoEtMaxStr);
                    hasIsoEtRange = true;
                }
                pTMin = convert(match[8].str());
                pTMax = convert(match[9].str());
                triggerName = match[10].str();
            }
        } else {
            std::cerr << "Error: Failed to parse histogram name: " << histName << std::endl;
        }

        return std::make_tuple(cuts, massWindowLabel, pTMin, pTMax, triggerName, histType, hasIsoEtRange, isoEtMin, isoEtMax);
    };

    for (const auto& trigger : triggers) {
        TDirectory* triggerDir = inputFile->GetDirectory(trigger.c_str());
        if (!triggerDir) {
            std::cerr << "Trigger directory '" << trigger << "' not found. Skipping." << std::endl;
            continue;
        }

        TIter nextKey(triggerDir->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*)nextKey())) {
            TObject* obj = key->ReadObj();
            if (!obj) continue;

            std::string histName = obj->GetName();

            if (obj->InheritsFrom(TH1::Class())) {
                TH1* hist = dynamic_cast<TH1*>(obj);

                auto [cuts, massWindowLabel, pTMin, pTMax, triggerName, histType, hasIsoEtRange, isoEtMin, isoEtMax] = parseIsolationHistNameWithCuts(histName);

                if (triggerName != trigger) {
                    delete obj;
                    continue;
                }

                // Create output directories
                std::ostringstream cutDirStream;
                cutDirStream << plotDirectory << "/E" << Utils::formatToThreeSigFigs(cuts.clusECore)
                             << "_Chi" << Utils::formatToThreeSigFigs(cuts.chi)
                             << "_Asym" << Utils::formatToThreeSigFigs(cuts.asymmetry);
                std::string cutDirPath = cutDirStream.str();
                gSystem->mkdir(cutDirPath.c_str(), true);

                std::string isolationDir = cutDirPath + "/isolationEnergies";
                gSystem->mkdir(isolationDir.c_str(), true);

                std::ostringstream ptDirStream;
                ptDirStream << isolationDir << "/pT_" << Utils::formatToThreeSigFigs(pTMin)
                            << "_to_" << Utils::formatToThreeSigFigs(pTMax);
                std::string ptDirPath = ptDirStream.str();
                gSystem->mkdir(ptDirPath.c_str(), true);

                std::string outputDirPath = ptDirPath;

                if (hasIsoEtRange) {
                    std::ostringstream isoEtDirStream;
                    isoEtDirStream << ptDirPath << "/isoEt_" << Utils::formatToThreeSigFigs(isoEtMin)
                                   << "_to_" << Utils::formatToThreeSigFigs(isoEtMax);
                    std::string isoEtDirPath = isoEtDirStream.str();
                    gSystem->mkdir(isoEtDirPath.c_str(), true);
                    outputDirPath = isoEtDirPath;
                }

                // Output file path
                std::string outputFilePath = outputDirPath + "/" + histName + ".png";

                // Draw histogram
                TCanvas canvas("canvas", "Histogram Canvas", 800, 600);
                if (hist->InheritsFrom(TH2::Class())) {
                    hist->Draw("COLZ");
                    canvas.SetLogz();
                } else {
                    hist->SetStats(true);
                    hist->Draw("HIST");
                    canvas.SetLogy();
                }

                // Add labels
                AddLabelsToCanvas_isoHistsWithCuts(cuts, massWindowLabel, triggerName, triggerGroupName, pTMin, pTMax, histType, hasIsoEtRange, isoEtMin, isoEtMax);

                // Save canvas
                canvas.SaveAs(outputFilePath.c_str());
                std::cout << "Saved: " << outputFilePath << std::endl;
                
                
                // Common keys for maps
                auto totalKey = std::make_tuple(
                    triggerGroupName,  // raw name
                    triggerName,       // raw name
                    cuts.clusECore,
                    cuts.chi,
                    cuts.asymmetry,
                    pTMin,
                    pTMax,
                    massWindowLabel
                );

                if (histType == "isolatedPhotonCount" && hasIsoEtRange) {
                    // Key includes isoEtMin and isoEtMax
                    auto isoKey = std::make_tuple(
                        triggerGroupName,  // raw name
                        triggerName,       // raw name
                        cuts.clusECore,
                        cuts.chi,
                        cuts.asymmetry,
                        pTMin,
                        pTMax,
                        isoEtMin,
                        isoEtMax,
                        massWindowLabel
                    );

                    // Fill IsolatedPhotonLog
                    DataStructures::IsolatedPhotonLog isoLog;
                    isoLog.triggerGroupName = triggerGroupName; // raw name
                    isoLog.triggerName = triggerName;           // raw name
                    isoLog.clusECore = cuts.clusECore;
                    isoLog.chi = cuts.chi;
                    isoLog.asymmetry = cuts.asymmetry;
                    isoLog.pTMin = pTMin;
                    isoLog.pTMax = pTMax;
                    isoLog.isoMin = isoEtMin;
                    isoLog.isoMax = isoEtMax;
                    isoLog.isolatedEntries = static_cast<int>(hist->GetEntries());
                    isoLog.massWindowLabel = massWindowLabel;

                    // Store in the map
                    isolatedPhotonMap[isoKey] = isoLog;
                } else if (histType == "allPhotonCount") {
                    // Fill TotalPhotonLog
                    DataStructures::TotalPhotonLog totalLog;
                    totalLog.triggerGroupName = triggerGroupName; // raw name
                    totalLog.triggerName = triggerName;           // raw name
                    totalLog.clusECore = cuts.clusECore;
                    totalLog.chi = cuts.chi;
                    totalLog.asymmetry = cuts.asymmetry;
                    totalLog.pTMin = pTMin;
                    totalLog.pTMax = pTMax;
                    totalLog.totalEntries = static_cast<int>(hist->GetEntries());
                    totalLog.massWindowLabel = massWindowLabel;

                    totalPhotonMap[totalKey] = totalLog;
                } else if (histType == "ptPhoton") {
                    // Calculate weighted average pT
                    double weightedSum = 0.0;
                    double totalCounts = 0.0;
                    int nBins = hist->GetNbinsX();
                    for (int i = 1; i <= nBins; ++i) {
                        double binContent = hist->GetBinContent(i);
                        double binCenter = hist->GetBinCenter(i);
                        weightedSum += binContent * binCenter;
                        totalCounts += binContent;
                    }
                    double weightedAveragePt = (totalCounts > 0) ? (weightedSum / totalCounts) : 0.0;

                    // Fill PtWeightingLog
                    DataStructures::PtWeightingLog ptLog;
                    ptLog.triggerGroupName = triggerGroupName; // raw name
                    ptLog.triggerName = triggerName;           // raw name
                    ptLog.clusECore = cuts.clusECore;
                    ptLog.chi = cuts.chi;
                    ptLog.asymmetry = cuts.asymmetry;
                    ptLog.pTMin = pTMin;
                    ptLog.pTMax = pTMax;
                    ptLog.weightedAveragePt = weightedAveragePt;
                    ptLog.massWindowLabel = massWindowLabel;

                    pTweightingMap[totalKey] = ptLog;
                }
            }

            delete obj;
        }
    }
}


void WriteIsolationDataToCSV(const std::string& outputFilePath) {
    std::ofstream outFile(outputFilePath);
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open CSV file for writing: " << outputFilePath << std::endl;
        return;
    }

    // Write CSV header
    outFile << "TriggerGroupName,TriggerName,ECore,Chi,Asymmetry,pT Min,pT Max,isoMin,isoMax,"
            << "Isolated Counts,Total Counts,Isolated/Total,Statistical Error,Weighted pT,"
            << "Bin Width,Bin Center,Isolated Yield,Isolated Yield Error,MassWindowLabel\n";


    // Iterate through isolatedPhotonMap and correlate with totalPhotonMap and pTweightingMap
    for (const auto& isoEntry : isolatedPhotonMap) {
        auto isoKey = isoEntry.first;
        const DataStructures::IsolatedPhotonLog& isoLog = isoEntry.second;

        // Extract the common key for totalPhotonMap and pTweightingMap (excluding isoEtMin and isoEtMax)
        auto totalKey = std::make_tuple(
            std::get<0>(isoKey),  // triggerGroupName
            std::get<1>(isoKey),  // triggerName
            std::get<2>(isoKey),  // clusECore
            std::get<3>(isoKey),  // chi
            std::get<4>(isoKey),  // asymmetry
            std::get<5>(isoKey),  // pTMin
            std::get<6>(isoKey),  // pTMax
            std::get<9>(isoKey)   // massWindowLabel
        );

        // Find corresponding entries in totalPhotonMap and pTweightingMap
        auto totalEntry = totalPhotonMap.find(totalKey);
        auto pTweightingEntry = pTweightingMap.find(totalKey);

        // Ensure corresponding entries exist in both maps
        if (totalEntry != totalPhotonMap.end() && pTweightingEntry != pTweightingMap.end()) {
            const DataStructures::TotalPhotonLog& totalLog = totalEntry->second;
            const DataStructures::PtWeightingLog& pTLog = pTweightingEntry->second;

            // Map triggerGroupName and triggerName to human-readable names
            std::string readableTriggerGroupName = Utils::getTriggerCombinationName(
                isoLog.triggerGroupName, TriggerCombinationNames::triggerCombinationNameMap);

            std::string readableTriggerName = isoLog.triggerName;
            auto triggerNameIt = TriggerConfig::triggerNameMap.find(isoLog.triggerName);
            if (triggerNameIt != TriggerConfig::triggerNameMap.end()) {
                readableTriggerName = triggerNameIt->second;
            }

            // Calculate the ratio (isolated / total)
            double ratio = (totalLog.totalEntries > 0) ? static_cast<double>(isoLog.isolatedEntries) / totalLog.totalEntries : 0.0;

            // Calculate the statistical error
            double error = 0.0;
            if (totalLog.totalEntries > 0 && isoLog.isolatedEntries > 0) {
                double isolatedError = std::sqrt(isoLog.isolatedEntries);
                double totalError = std::sqrt(totalLog.totalEntries);
                error = ratio * std::sqrt(
                    (isolatedError / isoLog.isolatedEntries) * (isolatedError / isoLog.isolatedEntries) +
                    (totalError / totalLog.totalEntries) * (totalError / totalLog.totalEntries)
                );
            }
            
            // Calculate bin width and bin center
            double binWidth = isoLog.pTMax - isoLog.pTMin;
            double binCenter = (isoLog.pTMin + isoLog.pTMax) / 2.0;

            // Calculate isolated yield and its error
            double isolatedYield = (binWidth > 0) ? (static_cast<double>(isoLog.isolatedEntries) / binWidth) : 0.0;
            double isolatedYieldError = (isoLog.isolatedEntries > 0) ? (std::sqrt(static_cast<double>(isoLog.isolatedEntries)) / binWidth) : 0.0;

            // Write CSV row with additional columns
            outFile << isoLog.triggerGroupName << ","
                    << isoLog.triggerName << ","
                    << isoLog.clusECore << ","
                    << isoLog.chi << ","
                    << isoLog.asymmetry << ","
                    << isoLog.pTMin << ","
                    << isoLog.pTMax << ","
                    << isoLog.isoMin << ","
                    << isoLog.isoMax << ","
                    << isoLog.isolatedEntries << ","
                    << totalLog.totalEntries << ","
                    << ratio << ","
                    << error << ","
                    << pTLog.weightedAveragePt << ","
                    << binWidth << ","
                    << binCenter << ","
                    << isolatedYield << ","
                    << isolatedYieldError << ","
                    << isoLog.massWindowLabel << "\n";
        } else {
            std::cerr << "Warning: Corresponding total or pT weighting data not found for key.\n";
        }
    }

    outFile.close();
    std::cout << "CSV file successfully written to " << outputFilePath << "\n";
}


void readDataFromCSV(
    const std::string& filename,
    std::map<std::tuple<
        std::string, // TriggerGroupName
        std::string, // TriggerName
        float,       // ECore
        float,       // Chi
        float,       // Asymmetry
        float,       // pT Min
        float,       // pT Max
        float,       // isoMin
        float,       // isoMax
        std::string  // MassWindowLabel
    >, DataStructures::IsolationData>& dataMap_inMassWindow,
    std::map<std::tuple<
        std::string, // TriggerGroupName
        std::string, // TriggerName
        float,       // ECore
        float,       // Chi
        float,       // Asymmetry
        float,       // pT Min
        float,       // pT Max
        float,       // isoMin
        float,       // isoMax
        std::string  // MassWindowLabel
    >, DataStructures::IsolationData>& dataMap_outsideMassWindow) {

    std::ifstream file(filename);
    std::string line;
    if (!file.is_open()) {
        std::cerr << "\033[31m[ERROR]\033[0m Error opening file: " << filename << std::endl;
        return;
    }

    // Skip the header line
    std::getline(file, line);
    std::cout << "\033[33m[INFO]\033[0m Skipping header: " << line << std::endl;

    int lineNumber = 1; // Start from 1 for the header

    // Read CSV data
    while (std::getline(file, line)) {
        lineNumber++;
        std::stringstream ss(line);
        std::string token, massWindowLabel;
        std::string triggerGroupName, triggerName;
        float eCore = 0.0f, chi = 0.0f, asym = 0.0f, ptMin = 0.0f, ptMax = 0.0f, isoMin = 0.0f, isoMax = 0.0f;
        int isolatedCounts = 0, totalCounts = 0;
        double ratio = 0.0, error = 0.0, weightedPt = 0.0;
        double binWidth = 0.0, binCenter = 0.0, isolatedYield = 0.0, isolatedYieldError = 0.0;

        // Parse values from CSV
        try {
            std::getline(ss, triggerGroupName, ',');
            std::getline(ss, triggerName, ',');

            std::getline(ss, token, ',');
            eCore = std::stof(token);

            std::getline(ss, token, ',');
            chi = std::stof(token);

            std::getline(ss, token, ',');
            asym = std::stof(token);

            std::getline(ss, token, ',');
            ptMin = std::stof(token);

            std::getline(ss, token, ',');
            ptMax = std::stof(token);

            std::getline(ss, token, ',');
            isoMin = std::stof(token);

            std::getline(ss, token, ',');
            isoMax = std::stof(token);

            // Parse the "Isolated Counts" and "Total Counts" columns
            std::getline(ss, token, ',');
            isolatedCounts = std::stoi(token);

            std::getline(ss, token, ',');
            totalCounts = std::stoi(token);

            // Parse the ratio (Isolated/Total)
            std::getline(ss, token, ',');
            ratio = std::stod(token);

            // Parse the statistical error
            std::getline(ss, token, ',');
            error = std::stod(token);

            // Parse the weighted pT value
            std::getline(ss, token, ',');
            weightedPt = std::stod(token);

            // Parse the Bin Width
            std::getline(ss, token, ',');
            binWidth = std::stod(token);

            // Parse the Bin Center
            std::getline(ss, token, ',');
            binCenter = std::stod(token);

            // Parse the Isolated Yield
            std::getline(ss, token, ',');
            isolatedYield = std::stod(token);

            // Parse the Isolated Yield Error
            std::getline(ss, token, ',');
            isolatedYieldError = std::stod(token);

            // Parse the MassWindowLabel
            std::getline(ss, massWindowLabel, ',');

            // Create the key tuple
            auto key = std::make_tuple(
                triggerGroupName,
                triggerName,
                eCore,
                chi,
                asym,
                ptMin,
                ptMax,
                isoMin,
                isoMax,
                massWindowLabel
            );

            // Create an IsolationData struct to hold the data
            DataStructures::IsolationData isoData;
            isoData.isolatedCounts = isolatedCounts;
            isoData.totalCounts = totalCounts;
            isoData.ratio = ratio;
            isoData.error = error;
            isoData.weightedPt = weightedPt;
            isoData.binWidth = binWidth;
            isoData.binCenter = binCenter;
            isoData.isolatedYield = isolatedYield;
            isoData.isolatedYieldError = isolatedYieldError;

            // Debugging output
            std::cout << "\033[32m[DEBUG]\033[0m Line " << lineNumber << ": Read data - "
                      << "TriggerGroupName: " << triggerGroupName << ", "
                      << "TriggerName: " << triggerName << ", "
                      << "ECore: " << eCore << ", "
                      << "Chi: " << chi << ", "
                      << "Asymmetry: " << asym << ", "
                      << "pT Min: " << ptMin << ", "
                      << "pT Max: " << ptMax << ", "
                      << "isoMin: " << isoMin << ", "
                      << "isoMax: " << isoMax << ", "
                      << "Isolated Counts: " << isolatedCounts << ", "
                      << "Total Counts: " << totalCounts << ", "
                      << "Ratio: " << ratio << ", "
                      << "Error: " << error << ", "
                      << "Weighted pT: " << weightedPt << ", "
                      << "Bin Width: " << binWidth << ", "
                      << "Bin Center: " << binCenter << ", "
                      << "Isolated Yield: " << isolatedYield << ", "
                      << "Isolated Yield Error: " << isolatedYieldError << ", "
                      << "MassWindowLabel: " << massWindowLabel << std::endl;

            // Add data to the appropriate map based on MassWindowLabel
            if (massWindowLabel == "inMassWindow") {
                dataMap_inMassWindow[key] = isoData;
            } else if (massWindowLabel == "outsideMassWindow") {
                dataMap_outsideMassWindow[key] = isoData;
            } else {
                std::cerr << "\033[31m[WARNING]\033[0m Line " << lineNumber << ": Unknown MassWindowLabel '"
                          << massWindowLabel << "'. Skipping this entry." << std::endl;
            }
        } catch (const std::exception& e) {
            std::cerr << "\033[31m[ERROR]\033[0m Line " << lineNumber << ": Exception occurred while parsing line. "
                      << "Error: " << e.what() << ". Line content: " << line << std::endl;
        }
    }
    file.close();

    // Summary of data read
    std::cout << "\033[33m[INFO]\033[0m Finished reading CSV file." << std::endl;
    std::cout << "\033[33m[INFO]\033[0m Total entries in dataMap_inMassWindow: " << dataMap_inMassWindow.size() << std::endl;
    std::cout << "\033[33m[INFO]\033[0m Total entries in dataMap_outsideMassWindow: " << dataMap_outsideMassWindow.size() << std::endl;
}


void GeneratePerTriggerSpectraPlots(
    const std::map<std::tuple<
        std::string, // TriggerGroupName
        std::string, // TriggerName
        float,       // ECore
        float,       // Chi
        float,       // Asymmetry,
        float,       // pT Min
        float,       // pT Max
        float,       // isoMin
        float,       // isoMax
        std::string  // MassWindowLabel
    >, DataStructures::IsolationData>& dataMap_inMassWindow,
    const std::map<std::tuple<
        std::string, // TriggerGroupName
        std::string, // TriggerName
        float,       // ECore
        float,       // Chi
        float,       // Asymmetry,
        float,       // pT Min
        float,       // pT Max
        float,       // isoMin
        float,       // isoMax
        std::string  // MassWindowLabel
    >, DataStructures::IsolationData>& dataMap_outsideMassWindow,
    const std::string& basePlotDirectory,
    const std::map<std::string, std::string>& triggerCombinationNameMap,
    const std::map<std::string, std::string>& triggerNameMap,
    const std::vector<std::pair<float, float>>& exclusionRanges) {
    
    // Define a list of marker styles to cycle through for different isoEt ranges
    std::vector<int> markerStyles = {20, 21, 22, 23, 29, 34, 35, 36, 38, 39};
    size_t markerStyleCount = markerStyles.size();
    
    // Structure to hold unique group keys
    struct GroupKey {
        std::string triggerGroupName;
        std::string triggerName;
        float eCore;
        float chi;
        float asymmetry;
        
        bool operator<(const GroupKey& other) const {
            return std::tie(triggerGroupName, triggerName, eCore, chi, asymmetry) <
            std::tie(other.triggerGroupName, other.triggerName, other.eCore, other.chi, other.asymmetry);
        }
    };
    
    // Structure to hold grouped data entries
    struct GroupedDataEntry {
        std::pair<float, float> isoEtRange; // {isoMin, isoMax}
        float ptMin;
        float ptMax;
        DataStructures::IsolationData isoData_in;
    };
    
    // Organize dataMap_inMassWindow into groups
    std::map<GroupKey, std::vector<GroupedDataEntry>> groupedData;
    
    for (const auto& [key, isoData_in] : dataMap_inMassWindow) {
        // Extract group key (excluding MassWindowLabel)
        GroupKey groupKey = {
            std::get<0>(key), // triggerGroupName
            std::get<1>(key), // triggerName
            std::get<2>(key), // ECore
            std::get<3>(key), // Chi
            std::get<4>(key)  // Asymmetry
        };
        
        // Extract isoEt range and pT bins
        float isoMin = std::get<7>(key);
        float isoMax = std::get<8>(key);
        float ptMin = std::get<5>(key);
        float ptMax = std::get<6>(key);
        
        // Check if the current isoEt range is in the exclusion list
        std::pair<float, float> currentIsoEtRange = {isoMin, isoMax};
        if (std::find(exclusionRanges.begin(), exclusionRanges.end(), currentIsoEtRange) != exclusionRanges.end()) {
            continue;  // Skip excluded isoEt ranges
        }
        
        
        // Create a GroupedDataEntry
        GroupedDataEntry entry = {
            {isoMin, isoMax},
            ptMin,
            ptMax,
            isoData_in
        };
        
        // Append to groupedData
        groupedData[groupKey].emplace_back(entry);
    }
    
    // Iterate over each group to create plots
    for (const auto& [groupKey, isoEtDataVec] : groupedData) {
        const std::string& triggerGroupName = groupKey.triggerGroupName;
        const std::string& triggerName = groupKey.triggerName;
        float eCore = groupKey.eCore;
        float chi = groupKey.chi;
        float asym = groupKey.asymmetry;
        
        // Map to human-readable names
        std::string readableTriggerGroupName = Utils::getTriggerCombinationName(
                                                                                triggerGroupName, triggerCombinationNameMap);
        
        std::string readableTriggerName = triggerName;
        auto triggerNameIt = triggerNameMap.find(triggerName);
        if (triggerNameIt != triggerNameMap.end()) {
            readableTriggerName = triggerNameIt->second;
        }
        
        // Define output directory
        std::ostringstream dirStream;
        dirStream << basePlotDirectory << "/" << triggerGroupName
        << "/E" << Utils::formatToThreeSigFigs(eCore)
        << "_Chi" << Utils::formatToThreeSigFigs(chi)
        << "_Asym" << Utils::formatToThreeSigFigs(asym)
        << "/Spectra/Overlay";
        std::string dirPath = dirStream.str();
        gSystem->mkdir(dirPath.c_str(), true);
        
        // Create a TCanvas
        TCanvas* canvas = new TCanvas("canvas", "Isolated Photon Spectra Overlay", 800, 600);
        canvas->SetLogy();

        // Create a TMultiGraph to hold all TGraphErrors
        TMultiGraph* multiGraph = new TMultiGraph();
        
        // Create a legend
        TLegend* legend = new TLegend(0.38, 0.78, 0.7, 0.9); // Adjust positions as needed
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.024);
        
        // Assign marker styles to isoEt ranges
        // Create a map from isoEt range to marker style
        std::map<std::pair<float, float>, int> isoEtRangeMarkerMap;
        size_t currentMarkerIndex = 0;
        for (const auto& entry : isoEtDataVec) {
            const auto& isoEtRange = entry.isoEtRange;
            if (isoEtRangeMarkerMap.find(isoEtRange) == isoEtRangeMarkerMap.end()) {
                isoEtRangeMarkerMap[isoEtRange] = markerStyles[currentMarkerIndex % markerStyleCount];
                currentMarkerIndex++;
            }
        }
        // Track the global maximum Y value for scaling
        double globalMaxY = 0.0;
        // Iterate over each isoEt range within the group
        for (const auto& entry : isoEtDataVec) {
            const auto& isoEtRange = entry.isoEtRange;
            float isoMin = isoEtRange.first;
            float isoMax = isoEtRange.second;
            float ptMin = entry.ptMin;
            float ptMax = entry.ptMax;
            const DataStructures::IsolationData& isoData_in = entry.isoData_in;
            
            // Retrieve the assigned marker style
            int markerStyle = isoEtRangeMarkerMap[isoEtRange];
            
            // Construct the key for outsideMassWindow using the full tuple
            std::tuple<std::string, std::string, float, float, float, float, float, float, float, std::string> outsideKey =
            std::make_tuple(
                            triggerGroupName,
                            triggerName,
                            eCore,
                            chi,
                            asym,
                            ptMin,
                            ptMax,
                            isoMin,
                            isoMax,
                            "outsideMassWindow"
                            );
            
            // Find the corresponding outsideMassWindow data
            auto it_outside = dataMap_outsideMassWindow.find(outsideKey);
            if (it_outside == dataMap_outsideMassWindow.end()) {
                std::cerr << "\033[31m[WARNING]\033[0m No 'outsideMassWindow' data found for isoEt range ["
                << isoMin << ", " << isoMax << "] in group '" << readableTriggerName << "'. Skipping.\n";
                continue;
            }
            const DataStructures::IsolationData& isoData_out = it_outside->second;
            globalMaxY = std::max(globalMaxY, std::max(isoData_in.isolatedYield, isoData_out.isolatedYield));
            
            // Create TGraphErrors for inMassWindow (blue)
            TGraphErrors* graphIn = new TGraphErrors();
            graphIn->SetPoint(0, isoData_in.binCenter, isoData_in.isolatedYield);
            graphIn->SetPointError(0, 0.0, isoData_in.isolatedYieldError);
            graphIn->SetMarkerStyle(markerStyle);
            graphIn->SetMarkerColor(kBlue);
            graphIn->SetLineColor(kBlue);
            graphIn->SetLineWidth(2);
            graphIn->SetTitle("Isolated Photon Spectra");
            
            // Create TGraphErrors for outsideMassWindow (red)
            TGraphErrors* graphOut = new TGraphErrors();
            graphOut->SetPoint(0, isoData_out.binCenter, isoData_out.isolatedYield);
            graphOut->SetPointError(0, 0.0, isoData_out.isolatedYieldError);
            graphOut->SetMarkerStyle(markerStyle);
            graphOut->SetMarkerColor(kRed);
            graphOut->SetLineColor(kRed);
            graphOut->SetLineWidth(2);
            graphOut->SetTitle("Isolated Photon Spectra");
            
            // Add graphs to the multigraph
            multiGraph->Add(graphIn, "P");
            multiGraph->Add(graphOut, "P");
            
            // Create descriptive legend entries
            std::ostringstream legendEntryIn, legendEntryOut;
            legendEntryIn << "Isolated Photons from Meson Decay Yield  (" << isoMin << " #leq E_{T,iso} < " << isoMax << " GeV)";
            legendEntryOut << "Prompt Photon Candidate Yield (" << isoMin << " #leq E_{T,iso} < " << isoMax << " GeV)";
            
            // To avoid duplicate legend entries, check if already added
            bool entryExistsIn = false;
            bool entryExistsOut = false;
            int entryCount = legend->GetListOfPrimitives()->GetEntries();
            for (int i = 0; i < entryCount; ++i) {
                TLegendEntry* entry = dynamic_cast<TLegendEntry*>(legend->GetListOfPrimitives()->At(i));
                if (entry) {
                    std::string label = entry->GetLabel();
                    if (label == legendEntryIn.str()) {
                        entryExistsIn = true;
                    }
                    if (label == legendEntryOut.str()) {
                        entryExistsOut = true;
                    }
                }
            }
            if (!entryExistsIn) {
                legend->AddEntry(graphIn, legendEntryIn.str().c_str(), "p");
            }
            if (!entryExistsOut) {
                legend->AddEntry(graphOut, legendEntryOut.str().c_str(), "p");
            }
        }
        
        // Determine Y-axis range based on globalMaxY
        // Ensure that globalMaxY is positive
        if (globalMaxY <= 0.0) {
            std::cerr << "\033[31m[ERROR]\033[0m Maximum Y value is non-positive for group '"
                      << readableTriggerName << "'. Skipping plot.\n";
            // Clean up and continue to next group
            delete multiGraph;
            delete legend;
            delete canvas;
            continue;
        }
        double yMax = globalMaxY * 1.5; // Scale maximum Y by 1.5 for better visualization
        
        // Set axis titles and dynamic Y-axis range
        std::string plotTitle = ("Isolated Photon Spectra for Trigger: " + readableTriggerName);
        multiGraph->SetTitle(plotTitle.c_str());
        multiGraph->GetXaxis()->SetTitle("Cluster p_{T} [GeV]");
        multiGraph->GetYaxis()->SetTitle("Yield");
        
        // Set the Y-axis range with a minimum to avoid log scale issues
        multiGraph->SetMinimum(1e-1);               // Set a minimum above zero for log scale
        multiGraph->SetMaximum(yMax);               // Set maximum based on data
        
        // Draw the multigraph
        multiGraph->Draw("A");
        
        // Draw the legend
        legend->Draw();
        
        // Add labels using TLatex in the top-left corner
        TLatex labelText;
        labelText.SetNDC();
        labelText.SetTextSize(0.025);       // Adjusted text size
        labelText.SetTextColor(kBlack);    // Ensured text color is black for readability
        
        double xStart = 0.55; // Starting x-coordinate (left side)
        double yStartLabel = 0.72; // Starting y-coordinate
        double yStepLabel = 0.05;  // Vertical spacing between lines
        
        // Prepare label strings
        std::string triggerGroupLabel = "Trigger Group: " + readableTriggerGroupName;
        std::string triggerNameLabel = "Trigger: " + readableTriggerName;
        std::string eCoreLabel = "ECore > " + Utils::formatToThreeSigFigs(eCore) + " GeV";
        std::string chiLabel = "Chi2 < " + Utils::formatToThreeSigFigs(chi);
        std::string asymLabel = "Asymmetry < " + Utils::formatToThreeSigFigs(asym);
        
        // Draw labels
        labelText.DrawLatex(xStart, yStartLabel, triggerGroupLabel.c_str());
        labelText.DrawLatex(xStart, yStartLabel - yStepLabel, triggerNameLabel.c_str());
        labelText.DrawLatex(xStart, yStartLabel - 2 * yStepLabel, eCoreLabel.c_str());
        labelText.DrawLatex(xStart, yStartLabel - 3 * yStepLabel, chiLabel.c_str());
        labelText.DrawLatex(xStart, yStartLabel - 4 * yStepLabel, asymLabel.c_str());
        
        // Force canvas update before saving
        canvas->Modified();
        canvas->Update();
        
        // Save the canvas
        std::ostringstream outputFilePathStream;
        outputFilePathStream << dirPath << "/OverlaySpectra_" << readableTriggerName << ".png";
        std::string outputFilePath = outputFilePathStream.str();
        canvas->SaveAs(outputFilePath.c_str());
        std::cout << "\033[33m[INFO]\033[0m Saved overlay spectra plot to " << outputFilePath << std::endl;
        
        // Clean up
        delete multiGraph;
        delete legend;
        delete canvas;
        // Note: TGraphErrors objects are managed by ROOT's ownership model once added to TMultiGraph
    }
}

// Define SpectraGroupKey structure
struct SpectraGroupKey {
    std::string triggerGroupName;
    float eCore;
    float chi;
    float asymmetry;

    bool operator<(const SpectraGroupKey& other) const {
        return std::tie(triggerGroupName, eCore, chi, asymmetry) <
               std::tie(other.triggerGroupName, other.eCore, other.chi, other.asymmetry);
    }
};

struct CombinedSpectraData {
    float pTCenter;
    std::string triggerUsed;
    double isolatedYield_in;         // Changed from float to double
    double isolatedYieldError_in;    // Changed from float to double
    double isolatedYield_out;        // Changed from float to double
    double isolatedYieldError_out;   // Changed from float to double
};

void SortAndCombineSpectraData(
    const std::map<std::tuple<
        std::string, // TriggerGroupName
        std::string, // TriggerName
        float,       // ECore
        float,       // Chi
        float,       // Asymmetry
        float,       // pT Min
        float,       // pT Max
        float,       // isoMin
        float,       // isoMax
        std::string  // MassWindowLabel
    >, DataStructures::IsolationData>& dataMap_inMassWindow,
    const std::map<std::tuple<
        std::string, // TriggerGroupName
        std::string, // TriggerName
        float,       // ECore
        float,       // Chi
        float,       // Asymmetry
        float,       // pT Min
        float,       // pT Max
        float,       // isoMin
        float,       // isoMax
        std::string  // MassWindowLabel
    >, DataStructures::IsolationData>& dataMap_outsideMassWindow,
    const std::map<std::string, std::map<std::string, double>>& combinationToTriggerEfficiencyPoints,
    const std::vector<std::pair<float, float>>& exclusionRanges,
    const std::vector<std::pair<double, double>>& pT_bins,
    double pTExclusionMax,
    std::map<SpectraGroupKey, std::map<float, CombinedSpectraData>>& combinedSpectraDataMap) {

    // Function to extract photon threshold from trigger name
    auto extractPhotonThreshold = [](const std::string& triggerName) -> double {
        std::regex re("Photon_(\\d+)_GeV");
        std::smatch match;
        if (std::regex_search(triggerName, match, re)) {
            if (match.size() >= 2) {
                return std::stod(match[1]);
            }
        }
        return 0.0; // Default to 0 for MinBias or parsing failure
    };

    // Function to get x99 (efficiency threshold) for a trigger
    auto getX99 = [&](const std::string& triggerGroupName, const std::string& triggerName) -> double {
        auto groupEffIt = combinationToTriggerEfficiencyPoints.find(triggerGroupName);
        if (groupEffIt != combinationToTriggerEfficiencyPoints.end()) {
            auto effIt = groupEffIt->second.find(triggerName);
            if (effIt != groupEffIt->second.end()) {
                return effIt->second;
            }
        }
        return std::numeric_limits<double>::max(); // Assign max value if not found
    };

    // Map to store data per group and pT bin
    std::map<SpectraGroupKey, std::map<std::pair<float, float>, std::map<std::string, std::map<std::string, DataStructures::IsolationData>>>>
        dataPerGroup;

    // Helper function to process data maps
    auto processDataMap = [&](const auto& dataMap) {
        for (const auto& entry : dataMap) {
            const auto& key = entry.first;
            const auto& isoData = entry.second;

            std::string triggerGroupName = std::get<0>(key);
            std::string triggerName = std::get<1>(key);
            float eCore = std::get<2>(key);
            float chi = std::get<3>(key);
            float asymmetry = std::get<4>(key);
            float pTMin = std::get<5>(key);
            float pTMax = std::get<6>(key);
            float isoMin = std::get<7>(key);
            float isoMax = std::get<8>(key);
            std::string massWindowLabel = std::get<9>(key); // "inMassWindow" or "outsideMassWindow"

            // Exclude isoEtRanges if necessary
            std::pair<float, float> isoEtRange = {isoMin, isoMax};
            if (std::find(exclusionRanges.begin(), exclusionRanges.end(), isoEtRange) != exclusionRanges.end()) {
                continue; // Exclude this isoEtRange
            }

            SpectraGroupKey groupKey{triggerGroupName, eCore, chi, asymmetry};
            std::pair<float, float> pTBin{pTMin, pTMax};

            dataPerGroup[groupKey][pTBin][triggerName][massWindowLabel] = isoData;
        }
    };

    // Process both inMassWindow and outsideMassWindow data
    processDataMap(dataMap_inMassWindow);
    processDataMap(dataMap_outsideMassWindow);

    // Now, for each group and pT bin, apply trigger selection logic
    for (const auto& groupEntry : dataPerGroup) {
        const SpectraGroupKey& groupKey = groupEntry.first;
        const auto& pTBinMap = groupEntry.second;

        std::string triggerGroupName = groupKey.triggerGroupName;

        // Collect triggers in this group
        std::set<std::string> triggersInGroup;
        for (const auto& pTBinEntry : pTBinMap) {
            const auto& triggerMap = pTBinEntry.second;
            for (const auto& triggerEntry : triggerMap) {
                triggersInGroup.insert(triggerEntry.first);
            }
        }

        // Convert to vector and sort triggers
        std::vector<std::string> sortedTriggerList(triggersInGroup.begin(), triggersInGroup.end());

        // Build triggerInfoMap: triggerName -> (photonThreshold, x99)
        std::map<std::string, std::pair<double, double>> triggerInfoMap;
        for (const auto& triggerName : sortedTriggerList) {
            double photonThreshold = extractPhotonThreshold(triggerName);
            double x99 = getX99(triggerGroupName, triggerName);
            triggerInfoMap[triggerName] = std::make_pair(photonThreshold, x99);
        }

        // Sort triggers according to the logic in processCutCombination
        std::sort(sortedTriggerList.begin(), sortedTriggerList.end(),
            [&](const std::string& a, const std::string& b) {
                double photonThresholdA = triggerInfoMap[a].first;
                double photonThresholdB = triggerInfoMap[b].first;

                if (photonThresholdA != photonThresholdB) {
                    return photonThresholdA > photonThresholdB; // Descending photon threshold
                } else {
                    double x99A = triggerInfoMap[a].second;
                    double x99B = triggerInfoMap[b].second;
                    return x99A < x99B; // Ascending x99
                }
            }
        );

        // Now, for each pT bin, select trigger and collect data
        for (const auto& pTBinEntry : pTBinMap) {
            const std::pair<float, float>& pTBin = pTBinEntry.first;
            float pTMin = pTBin.first;
            float pTMax = pTBin.second;
            float pTCenter = (pTMin + pTMax) / 2.0;

            // Apply pTExclusionMax
            if (pTCenter >= pTExclusionMax) {
                continue; // Exclude this pT bin
            }

            const auto& triggerMap = pTBinEntry.second;

            // Identify efficient triggers
            std::vector<std::string> efficientTriggers;
            for (const auto& triggerName : sortedTriggerList) {
                double x99 = triggerInfoMap[triggerName].second;
                if (x99 <= pTMax) {
                    efficientTriggers.push_back(triggerName);
                }
            }

            // Select trigger to use
            std::string triggerToUse = "MBD_NandS_geq_1";
            if (!efficientTriggers.empty()) {
                triggerToUse = efficientTriggers.front();
            }

            // Check if data is available for triggerToUse
            auto triggerDataIt = triggerMap.find(triggerToUse);
            if (triggerDataIt != triggerMap.end()) {
                const auto& massWindowDataMap = triggerDataIt->second;

                // Check that data is available for both inMassWindow and outsideMassWindow
                if (massWindowDataMap.find("inMassWindow") != massWindowDataMap.end() &&
                    massWindowDataMap.find("outsideMassWindow") != massWindowDataMap.end()) {
                    const DataStructures::IsolationData& isoData_in = massWindowDataMap.at("inMassWindow");
                    const DataStructures::IsolationData& isoData_out = massWindowDataMap.at("outsideMassWindow");

                    // Prepare CombinedSpectraData
                    CombinedSpectraData combinedData;
                    combinedData.pTCenter = pTCenter;
                    combinedData.triggerUsed = triggerToUse;
                    combinedData.isolatedYield_in = isoData_in.isolatedYield;
                    combinedData.isolatedYieldError_in = isoData_in.isolatedYieldError;
                    combinedData.isolatedYield_out = isoData_out.isolatedYield;
                    combinedData.isolatedYieldError_out = isoData_out.isolatedYieldError;

                    // Add to combinedSpectraDataMap
                    combinedSpectraDataMap[groupKey][pTCenter] = combinedData;
                } else {
                    // Missing data for inMassWindow or outsideMassWindow
                    std::cerr << "[WARNING] Missing data for inMassWindow or outsideMassWindow for trigger '" << triggerToUse
                              << "' at pT bin [" << pTMin << ", " << pTMax << "]. Skipping this pT bin.\n";
                }
            } else {
                // No data for selected trigger
                std::cerr << "[WARNING] No data for trigger '" << triggerToUse << "' at pT bin [" << pTMin << ", " << pTMax << "]. Skipping this pT bin.\n";
            }
        }
    }
}


void GenerateCombinedSpectraPlots(
    const std::map<SpectraGroupKey, std::map<float, CombinedSpectraData>>& combinedSpectraDataMap,
    const std::string& basePlotDirectory,
    const std::map<std::string, std::string>& triggerCombinationNameMap,
    const std::vector<std::pair<double, double>>& pT_bins,
    double pTExclusionMax) {
    std::cout << "[INFO] Starting GenerateCombinedSpectraPlots function.\n";

    for (const auto& [spectraGroupKey, ptDataMap] : combinedSpectraDataMap) {
        const std::string& triggerGroupName = spectraGroupKey.triggerGroupName;
        float eCore = spectraGroupKey.eCore;
        float chi = spectraGroupKey.chi;
        float asym = spectraGroupKey.asymmetry;

        // Map to human-readable names
        std::string readableTriggerGroupName = Utils::getTriggerCombinationName(
            triggerGroupName, triggerCombinationNameMap);

        std::cout << "[INFO] Processing plot for Trigger Group: " << readableTriggerGroupName
                  << ", ECore > " << eCore << " GeV, Chi2 < " << chi << ", Asymmetry < " << asym << ".\n";

        // Define output directory
        std::ostringstream dirStream;
        dirStream << basePlotDirectory << "/" << triggerGroupName
                  << "/E" << Utils::formatToThreeSigFigs(eCore)
                  << "_Chi" << Utils::formatToThreeSigFigs(chi)
                  << "_Asym" << Utils::formatToThreeSigFigs(asym)
                  << "/Spectra/Overlay";
        std::string dirPath = dirStream.str();
        gSystem->mkdir(dirPath.c_str(), true);

        // Create canvas
        TCanvas* canvas = new TCanvas("canvas", "Combined Isolated Photon Spectra", 800, 600);
        canvas->SetLogy();

        // Map to store graphs for cleanup
        std::vector<TGraphErrors*> graphs_in;
        std::vector<TGraphErrors*> graphs_out;

        // Legend
        TLegend* legend = new TLegend(0.55, 0.67, 0.88, 0.87);
        legend->SetBorderSize(0);
        legend->SetTextSize(0.03);

        // Map to organize data per trigger
        std::map<std::string, std::vector<CombinedSpectraData>> dataPerTrigger;

        // Collect data per trigger
        for (const auto& dataEntry : ptDataMap) {
            const CombinedSpectraData& data = dataEntry.second;
            dataPerTrigger[data.triggerUsed].push_back(data);
        }

        // Keep track of global Y-axis max and min
        double globalMaxY = 0.0;
        double globalMinY = std::numeric_limits<double>::max();

        // For each trigger
        for (const auto& triggerEntry : dataPerTrigger) {
            const std::string& triggerName = triggerEntry.first;
            const std::vector<CombinedSpectraData>& dataList = triggerEntry.second;

            std::vector<double> pTValues;
            std::vector<double> isolatedYields_in;
            std::vector<double> isolatedYieldsError_in;
            std::vector<double> isolatedYields_out;
            std::vector<double> isolatedYieldsError_out;

            for (const auto& data : dataList) {
                // Exclude data points with pT >= pTExclusionMax
                if (data.pTCenter >= pTExclusionMax) {
                    continue;
                }

                pTValues.push_back(data.pTCenter);
                isolatedYields_in.push_back(data.isolatedYield_in);
                isolatedYieldsError_in.push_back(data.isolatedYieldError_in);
                isolatedYields_out.push_back(data.isolatedYield_out);
                isolatedYieldsError_out.push_back(data.isolatedYieldError_out);
                globalMaxY = std::max({globalMaxY, data.isolatedYield_in, data.isolatedYield_out});
                globalMinY = std::min({globalMinY, data.isolatedYield_in, data.isolatedYield_out});
            }

            // Skip if no valid data points
            if (pTValues.empty()) {
                continue;
            }

            // Create TGraphErrors for inMassWindow data
            TGraphErrors* graphIn = new TGraphErrors(pTValues.size(),
                                                     pTValues.data(),
                                                     isolatedYields_in.data(),
                                                     nullptr,
                                                     isolatedYieldsError_in.data());
            // Set marker style and color
            int markerStyleIn = 20; // Closed circle
            int markerColor = kBlack;
            auto it_color = TriggerConfig::triggerColorMap.find(triggerName);
            if (it_color != TriggerConfig::triggerColorMap.end()) {
                markerColor = it_color->second;
            }
            graphIn->SetMarkerStyle(markerStyleIn);
            graphIn->SetMarkerColor(markerColor);
            graphIn->SetLineColor(markerColor);
            graphIn->SetLineWidth(2);

            // Create TGraphErrors for outsideMassWindow data
            TGraphErrors* graphOut = new TGraphErrors(pTValues.size(),
                                                      pTValues.data(),
                                                      isolatedYields_out.data(),
                                                      nullptr,
                                                      isolatedYieldsError_out.data());
            // Set marker style and color
            int markerStyleOut = 24; // Open circle
            graphOut->SetMarkerStyle(markerStyleOut);
            graphOut->SetMarkerColor(markerColor);
            graphOut->SetLineColor(markerColor);
            graphOut->SetLineWidth(2);

            // Draw graphs
            if (graphs_in.empty()) {
                graphIn->Draw("AP");
            } else {
                graphIn->Draw("P SAME");
            }
            graphOut->Draw("P SAME");

            // Add entries to legend
            std::string readableTriggerName = Utils::getTriggerCombinationName(triggerName, TriggerConfig::triggerNameMap);
            legend->AddEntry(graphIn, (readableTriggerName + " (In-Mass Window)").c_str(), "p");
            legend->AddEntry(graphOut, (readableTriggerName + " (Outside-Mass Window)").c_str(), "p");

            // Store graphs for cleanup
            graphs_in.push_back(graphIn);
            graphs_out.push_back(graphOut);
        }

        // Set Y-axis range
        if (globalMaxY > 0.0 && globalMinY > 0.0) {
            double yMin = globalMinY * 0.8;
            double yMax = globalMaxY * 1.2;
            gPad->Update();
            gPad->SetLogy();
            gPad->GetUymax();
            gPad->GetUymin();
            gPad->Update();
            graphs_in[0]->GetHistogram()->SetMinimum(yMin);
            graphs_in[0]->GetHistogram()->SetMaximum(yMax);
        }

        // Draw legend
        legend->Draw();

        // Add labels using TLatex
        TLatex labelText;
        labelText.SetNDC();
        labelText.SetTextSize(0.028);
        labelText.SetTextColor(kBlack);

        double xStart = 0.2;
        double yStartLabel = 0.85;
        double yStepLabel = 0.05;

        // Prepare label strings
        std::ostringstream oss;
        oss << "#font[62]{Active Trigger Group:} " << readableTriggerGroupName;
        labelText.DrawLatex(xStart, yStartLabel, oss.str().c_str());

        oss.str("");
        oss << "#font[62]{ECore #geq} " << eCore << " GeV";
        labelText.DrawLatex(xStart, yStartLabel - yStepLabel, oss.str().c_str());

        oss.str("");
        oss << "#font[62]{#chi^{2} <} " << chi;
        labelText.DrawLatex(xStart, yStartLabel - 2 * yStepLabel, oss.str().c_str());

        oss.str("");
        oss << "#font[62]{Asymmetry <} " << asym;
        labelText.DrawLatex(xStart, yStartLabel - 3 * yStepLabel, oss.str().c_str());

        // Update canvas and save
        canvas->Modified();
        canvas->Update();

        // Save the canvas
        std::string outputFilePath = dirPath + "/CombinedOverlaySpectra.png";
        canvas->SaveAs(outputFilePath.c_str());
        std::cout << "[INFO] Saved combined overlay spectra plot to " << outputFilePath << std::endl;

        // Clean up
        for (auto graph : graphs_in) {
            delete graph;
        }
        for (auto graph : graphs_out) {
            delete graph;
        }
        delete legend;
        delete canvas;
    }

    std::cout << "[INFO] Finished GenerateCombinedSpectraPlots function.\n";
}


// Helper Function to Generate Per-Trigger Plots
void GeneratePerTriggerIsoPlots(
    const std::map<std::tuple<
        std::string, // TriggerGroupName
        std::string, // TriggerName
        float,       // ECore
        float,       // Chi
        float,       // Asymmetry
        std::string  // MassWindowLabel
    >, std::map<std::pair<float, float>, std::vector<DataStructures::IsolationDataWithPt>>>& groupedData,
    const std::string& basePlotDirectory,
    const std::vector<std::pair<float, float>>& isoEtRanges,
    const std::vector<int>& isoEtColors,
    const std::vector<double>& referencePTGamma,
    const std::vector<double>& referenceRatio,
    const std::vector<double>& referenceStatError,
    const std::vector<double>& referenceTwoPTGamma,
    const std::vector<double>& referenceTwoRatio,
    const std::vector<double>& referenceTwoStatError,
    const std::map<std::string, std::string>& triggerCombinationNameMap,
    const std::map<std::string, std::string>& triggerNameMap,
    bool drawRefA,
    bool drawRefB,
    const std::vector<std::pair<float, float>>& exclusionRanges) {
    int groupCounter = 0;  // Initialize groupCounter

    // Now, generate standard plots for each trigger group
    for (const auto& groupEntry : groupedData) {
        groupCounter++;
        const auto& groupKey = groupEntry.first;
        const auto& isoEtDataMap = groupEntry.second;
        
        // Unpack the group key
        std::string triggerGroupName = std::get<0>(groupKey);
        std::string triggerName = std::get<1>(groupKey);
        float eCore = std::get<2>(groupKey);
        float chi = std::get<3>(groupKey);
        float asym = std::get<4>(groupKey);
        std::string massWindowLabel = std::get<5>(groupKey);
        
        // Map triggerGroupName and triggerName to human-readable names
        std::string readableTriggerGroupName = Utils::getTriggerCombinationName(
            triggerGroupName, TriggerCombinationNames::triggerCombinationNameMap);
        
        std::string readableTriggerName = triggerName;
        auto triggerNameIt = TriggerConfig::triggerNameMap.find(triggerName);
        if (triggerNameIt != TriggerConfig::triggerNameMap.end()) {
            readableTriggerName = triggerNameIt->second;
        }

        // Debugging output to verify mapping
        std::cout << "\033[34m[INFO]\033[0m Processing group " << groupCounter << ": "
                  << "TriggerGroupName: " << triggerGroupName << ", "
                  << "ReadableTriggerGroupName: " << readableTriggerGroupName << ", "
                  << "TriggerName: " << triggerName << ", "
                  << "ReadableTriggerName: " << readableTriggerName << ", "
                  << "ECore: " << eCore << ", "
                  << "Chi: " << chi << ", "
                  << "Asymmetry: " << asym << ", "
                  << "MassWindowLabel: " << massWindowLabel << std::endl;
    
        // Check if isoEtDataMap is empty
        if (isoEtDataMap.empty()) {
            std::cerr << "\033[31m[WARNING]\033[0m Group " << groupCounter << " has no data. Skipping." << std::endl;
            continue;
        }
        
        // Create output directories
        std::ostringstream dirStream;
        dirStream << basePlotDirectory << "/" << triggerGroupName
                  << "/E" << Utils::formatToThreeSigFigs(eCore)
                  << "_Chi" << Utils::formatToThreeSigFigs(chi)
                  << "_Asym" << Utils::formatToThreeSigFigs(asym);
        std::string dirPath = dirStream.str();
        gSystem->mkdir(dirPath.c_str(), true);
        
        std::string isolationDir = dirPath + "/isolationEnergies";
        gSystem->mkdir(isolationDir.c_str(), true);
        
        // Create a folder for mass window type
        std::string massWindowDir = isolationDir + "/" + massWindowLabel;
        gSystem->mkdir(massWindowDir.c_str(), true);
        
        // Create a TCanvas
        TCanvas canvas("canvas", "Isolation Data", 800, 600);

        // Create a TMultiGraph to overlay graphs
        TMultiGraph* multiGraph = new TMultiGraph();
        
        // Create a legend in the top-left corner
        TLegend legend(0.7, 0.65, 0.9, 0.8); // Adjust as needed
        legend.SetBorderSize(0);
        legend.SetTextSize(0.025);

        // Create the reference graphs and add to legend
        TGraphErrors* refGraphOne = nullptr;
        TGraphErrors* refGraphTwo = nullptr;

        if (drawRefA) {
            refGraphOne = new TGraphErrors(ReferenceData::referencePTGamma.size(),
                                           ReferenceData::referencePTGamma.data(),
                                           ReferenceData::referenceRatio.data(),
                                           nullptr,
                                           ReferenceData::referenceStatError.data());
            refGraphOne->SetMarkerStyle(20);
            refGraphOne->SetMarkerColor(kOrange);
            refGraphOne->SetLineColor(kOrange);
            refGraphOne->SetLineWidth(2);
            legend.AddEntry(refGraphOne, "#font[62]{PHENIX 2003 pp Run:} #frac{Isolated Direct}{All Direct}", "p");
        }

        if (drawRefB) {
            refGraphTwo = new TGraphErrors(ReferenceData::referenceTwoPTGamma.size(),
                                           ReferenceData::referenceTwoPTGamma.data(),
                                           ReferenceData::referenceTwoRatio.data(),
                                           nullptr,
                                           ReferenceData::referenceTwoStatError.data());
            refGraphTwo->SetMarkerStyle(20);
            refGraphTwo->SetMarkerColor(kOrange);
            refGraphTwo->SetLineColor(kOrange);
            refGraphTwo->SetLineWidth(2);
            legend.AddEntry(refGraphTwo, "#font[62]{PHENIX 2003 pp Run:} #frac{Isolated #pi^{0} Decay}{All #pi^{0} Decay}", "p");
        }

        // Loop over isoEt ranges and collect data for multigraph
        for (size_t i = 0; i < isoEtRanges.size(); ++i) {
            const auto& isoEtRange = isoEtRanges[i];

            // Skip excluded isoEt ranges
            if (std::find(exclusionRanges.begin(), exclusionRanges.end(), isoEtRange) != exclusionRanges.end()) {
                continue;
            }

            // Check if this isoEtRange is in the data
            auto it = isoEtDataMap.find(isoEtRange);
            if (it == isoEtDataMap.end()) {
                continue; // No data for this isoEtRange
            }

            const auto& isoDataList = it->second;
            int color = isoEtColors[i];

            // Prepare data vectors for plotting
            std::vector<double> weightedPts;
            std::vector<double> ratios;
            std::vector<double> errors;
            
            for (const auto& isoData : isoDataList) {
                if (isoData.weightedPt > 0) {
                    weightedPts.push_back(isoData.weightedPt);
                    ratios.push_back(isoData.ratio);
                    errors.push_back(isoData.error);
                } else {
                    std::cerr << "\033[31m[WARNING]\033[0m Invalid weightedPt (" << isoData.weightedPt
                              << ") in group " << groupCounter << ". Skipping this data point." << std::endl;
                }
            }
            
            if (weightedPts.empty()) {
                std::cerr << "\033[31m[WARNING]\033[0m No valid data points with positive weightedPt for isoEt range ["
                          << isoEtRange.first << ", " << isoEtRange.second << "] in group "
                          << groupCounter << ". Skipping this isoEt range." << std::endl;
                continue;
            }
            
            // Sort the data by weightedPt for better plotting
            std::vector<size_t> indices(weightedPts.size());
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(), [&](size_t i1, size_t i2) {
                return weightedPts[i1] < weightedPts[i2];
            });
            
            std::vector<double> sortedWeightedPts, sortedRatios, sortedErrors;
            for (size_t idx : indices) {
                sortedWeightedPts.push_back(weightedPts[idx]);
                sortedRatios.push_back(ratios[idx]);
                sortedErrors.push_back(errors[idx]);
            }
            
            // Debugging output
            std::cout << "\033[32m[DEBUG]\033[0m Data points for isoEt range [" << isoEtRange.first << ", " << isoEtRange.second << "] in group " << groupCounter << ":" << std::endl;
            for (size_t j = 0; j < sortedWeightedPts.size(); ++j) { // Changed loop variable to 'j' to avoid shadowing
                std::cout << "\033[32m[DEBUG]\033[0m   weightedPt: " << sortedWeightedPts[j]
                          << ", ratio: " << sortedRatios[j]
                          << ", error: " << sortedErrors[j] << std::endl;
            }

            // Create a TGraphErrors
            TGraphErrors* graph = new TGraphErrors(sortedWeightedPts.size(),
                                                   sortedWeightedPts.data(),
                                                   sortedRatios.data(),
                                                   nullptr,
                                                   sortedErrors.data());

            graph->SetMarkerStyle(20);
            graph->SetMarkerColor(color);
            graph->SetLineColor(color);
            graph->SetLineWidth(2);

            // Add to the multigraph
            multiGraph->Add(graph, "P");

            // Add entry to legend with the desired prefix
            std::ostringstream legendEntry;
            legendEntry << "Run 24 sPHENIX pp: " << isoEtRange.first << " #leq E_{T, iso} < " << isoEtRange.second << " GeV";
            legend.AddEntry(graph, legendEntry.str().c_str(), "p");
        }
        
        // Set titles
        std::string plotTitle = "Isolation Ratio vs pT (" + massWindowLabel + ")";
        multiGraph->SetTitle(plotTitle.c_str());
        multiGraph->GetXaxis()->SetTitle("Weighted p_{T} of Leading Cluster [GeV]");
        // Define y-axis title based on mass window type
        const std::string yAxisTitle = (massWindowLabel == "inMassWindow") ?
            "#frac{Isolated Photons from #pi^{0}/#eta Decays}{All Photons from #pi^{0}/#eta Decays}" :
            "#frac{Isolated Prompt Photons}{All Prompt Photons}";
        
        multiGraph->GetYaxis()->SetTitle(yAxisTitle.c_str());
        multiGraph->GetYaxis()->SetRangeUser(0, 2.0);
        multiGraph->GetXaxis()->SetLimits(2.0, 25.0);
        
        // Draw multigraph
        multiGraph->Draw("A");
        
        // Draw a dashed line at y = 1
        TLine* line = new TLine(2, 1, 25, 1);
        line->SetLineStyle(2); // Dashed line
        line->Draw();

        // Draw the reference graphs
        if (drawRefA && refGraphOne) {
            refGraphOne->Draw("P SAME");
        }

        if (drawRefB && refGraphTwo) {
            refGraphTwo->Draw("P SAME");
        }

        // Draw legend
        legend.Draw();
        
        // Add labels using TLatex in the top-left corner
        TLatex labelText;
        labelText.SetNDC();
        labelText.SetTextSize(0.0215);       // Adjusted text size
        labelText.SetTextColor(kBlack);    // Ensured text color is black
        double xStart = 0.18; // Starting x-coordinate (left side)
        double yStartLabel = 0.9; // Starting y-coordinate
        double yStepLabel = 0.035;  // Vertical spacing between lines

        // Prepare label strings
        std::string triggerGroupLabel = "#bf{Trigger Group:} " + readableTriggerGroupName;
        std::string triggerNameLabel = "bf{Trigger:} " + readableTriggerName;
        std::string eCoreLabel = "bf{ECore >} " + Utils::formatToThreeSigFigs(eCore) + " GeV";
        std::string chiLabel = "bf{Chi2} < " + Utils::formatToThreeSigFigs(chi);
        std::string asymLabel = "bf{Asymmetry} < " + Utils::formatToThreeSigFigs(asym);
        std::string massWindowLabelStr = "bf{Mass Window:} " + massWindowLabel;
        std::string coneRadiusLabel = "bf{#Delta R_{cone}} < 0.3";

        // Draw labels
        labelText.DrawLatex(xStart, yStartLabel, triggerGroupLabel.c_str());
        labelText.DrawLatex(xStart, yStartLabel - yStepLabel, triggerNameLabel.c_str());
        labelText.DrawLatex(xStart, yStartLabel - 2 * yStepLabel, eCoreLabel.c_str());
        labelText.DrawLatex(xStart, yStartLabel - 3 * yStepLabel, chiLabel.c_str());
        labelText.DrawLatex(xStart, yStartLabel - 4 * yStepLabel, asymLabel.c_str());
        labelText.DrawLatex(xStart, yStartLabel - 5 * yStepLabel, massWindowLabelStr.c_str());
        labelText.DrawLatex(xStart, yStartLabel - 6 * yStepLabel, coneRadiusLabel.c_str());

        // Force canvas update before saving
        canvas.Modified();
        canvas.Update();

        // Save the canvas
        std::ostringstream outputFilePathStream;
        outputFilePathStream << massWindowDir << "/IsolationRatio_vs_pT_" << triggerName << ".png";
        std::string outputFilePath = outputFilePathStream.str();
        canvas.SaveAs(outputFilePath.c_str());
        std::cout << "\033[33m[INFO]\033[0m Saved plot to " << outputFilePath << std::endl;

        // Clean up
        delete multiGraph;
        delete line;
        if (refGraphOne) delete refGraphOne;
        if (refGraphTwo) delete refGraphTwo;
        // Note: graphs added to multiGraph are owned by it and will be deleted automatically
    }
}

void SortAndCombineTriggers(
    const std::map<GroupKey, std::map<std::pair<float, float>, std::vector<DataStructures::IsolationDataWithPt>>>& groupedData,
    const std::map<std::string, std::map<std::string, double>>& combinationToTriggerEfficiencyPoints,
    std::map<std::string, std::vector<std::string>>& sortedTriggersByGroupName,
    std::map<std::string, std::map<std::pair<float, float>,
    std::vector<DataStructures::IsolationDataWithPt>>>& combinedTriggerDataMap) {
    // Function to extract photon threshold from trigger name
    auto extractPhotonThreshold = [](const std::string& triggerName) -> double {
        std::regex re("Photon_(\\d+)_GeV");
        std::smatch match;
        if (std::regex_search(triggerName, match, re)) {
            if (match.size() >= 2) {
                return std::stod(match[1]);
            }
        }
        return 0.0; // Default to 0 for MinBias or parsing failure
    };

    // Function to get x99 (efficiency threshold) for a trigger
    auto getX99 = [&](const std::string& triggerGroupName, const std::string& triggerName) -> double {
        auto groupEffIt = combinationToTriggerEfficiencyPoints.find(triggerGroupName);
        if (groupEffIt != combinationToTriggerEfficiencyPoints.end()) {
            auto effIt = groupEffIt->second.find(triggerName);
            if (effIt != groupEffIt->second.end()) {
                return effIt->second;
            }
        }
        return std::numeric_limits<double>::max(); // Assign max value if not found
    };

    // Step 1: Populate sortedTriggersByGroupName from groupedData
    for (const auto& groupEntry : groupedData) {
        const std::string& triggerGroupName = std::get<0>(groupEntry.first);
        const std::string& triggerName = std::get<1>(groupEntry.first);

        // Initialize the group if it doesn't exist
        if (sortedTriggersByGroupName.find(triggerGroupName) == sortedTriggersByGroupName.end()) {
            sortedTriggersByGroupName[triggerGroupName] = {};
        }

        // Add the trigger to the group if it's not already present
        if (std::find(sortedTriggersByGroupName[triggerGroupName].begin(),
                      sortedTriggersByGroupName[triggerGroupName].end(),
                      triggerName) == sortedTriggersByGroupName[triggerGroupName].end()) {
            sortedTriggersByGroupName[triggerGroupName].push_back(triggerName);
        }
    }

    // Step 2: Sort triggers within each group
    for (auto& groupEntry : sortedTriggersByGroupName) {
        const std::string& triggerGroupName = groupEntry.first;
        auto& triggerList = groupEntry.second;

        // Build a map of trigger to photon threshold and x99
        std::map<std::string, std::pair<double, double>> triggerInfoMap; // triggerName -> (photonThreshold, x99)
        for (const std::string& triggerName : triggerList) {
            double photonThreshold = extractPhotonThreshold(triggerName);
            double x99 = getX99(triggerGroupName, triggerName);
            triggerInfoMap[triggerName] = std::make_pair(photonThreshold, x99);
        }

        // Now sort the triggers
        std::sort(triggerList.begin(), triggerList.end(),
            [&](const std::string& a, const std::string& b) -> bool {
                double photonThresholdA = triggerInfoMap[a].first;
                double photonThresholdB = triggerInfoMap[b].first;

                if (photonThresholdA != photonThresholdB) {
                    return photonThresholdA > photonThresholdB; // Descending photon threshold
                } else {
                    double x99A = triggerInfoMap[a].second;
                    double x99B = triggerInfoMap[b].second;
                    return x99A < x99B; // Ascending x99
                }
            }
        );

        // Debugging output
        std::cout << "Trigger Group Name: " << triggerGroupName << "\n";
        std::cout << "Sorted Trigger List: ";
        for (const auto& trigger : triggerList) {
            double photonThreshold = triggerInfoMap[trigger].first;
            double x99 = triggerInfoMap[trigger].second;
            std::cout << trigger << " (Photon Threshold: " << photonThreshold << ", x99: " << x99 << "), ";
        }
        std::cout << "\n";
    }

    // Now process each group and isoEtRange
    for (const auto& [triggerGroupName, sortedTriggerList] : sortedTriggersByGroupName) {
        std::cout << "[PROCESSING] Combining triggers for group: " << triggerGroupName << "\n";

        // Collect all isoEtRanges for this group
        std::set<std::pair<float, float>> allIsoEtRanges;
        for (const auto& [groupKey, isoEtMap] : groupedData) {
            if (std::get<0>(groupKey) == triggerGroupName) {
                for (const auto& [isoEtRange, isoDataList] : isoEtMap) {
                    allIsoEtRanges.emplace(isoEtRange);
                }
            }
        }

        if (allIsoEtRanges.empty()) {
            std::cerr << "[ERROR] No isoEt ranges found for group '" << triggerGroupName << "'\n";
            continue;
        }

        // Collect all pT bins
        std::set<std::pair<float, float>> allPtBins;
        for (const auto& [groupKey, isoEtMap] : groupedData) {
            if (std::get<0>(groupKey) == triggerGroupName) {
                for (const auto& [isoEtRange, isoDataList] : isoEtMap) {
                    for (const auto& isoData : isoDataList) {
                        allPtBins.emplace(std::make_pair(isoData.ptMin, isoData.ptMax));
                    }
                }
            }
        }

        // Iterate over each isoEtRange
        for (const auto& isoEtRange : allIsoEtRanges) {
            std::cout << "[INFO] Processing isoEtRange: [" << isoEtRange.first << ", " << isoEtRange.second << "]\n";

            std::vector<DataStructures::IsolationDataWithPt> selectedDataPoints;

            // For each pT bin
            for (const auto& ptBin : allPtBins) {
                double pTMin = ptBin.first;
                double pTMax = ptBin.second;
                double pTCenter = (pTMin + pTMax) / 2.0;

                bool triggerAssigned = false;
                DataStructures::IsolationDataWithPt selectedIsoData;

                // Identify efficient triggers for this pT bin
                std::vector<std::string> efficientTriggers;
                for (const auto& triggerName : sortedTriggerList) {
                    double x99 = getX99(triggerGroupName, triggerName);
                    if (x99 <= pTMax) {
                        efficientTriggers.push_back(triggerName);
                    }
                }

                // Select the trigger to use
                std::string triggerToUse = "MBD_NandS_geq_1";
                if (!efficientTriggers.empty()) {
                    triggerToUse = efficientTriggers.front(); // First efficient trigger from sorted list
                }

                // Now get data for triggerToUse
                // Find the groupKey for triggerToUse
                bool foundGroupKey = false;
                GroupKey currentGroupKey;
                for (const auto& [gk, isoEtMap] : groupedData) {
                    if (std::get<0>(gk) == triggerGroupName && std::get<1>(gk) == triggerToUse) {
                        currentGroupKey = gk;
                        foundGroupKey = true;
                        break;
                    }
                }
                if (!foundGroupKey) {
                    std::cerr << "[ERROR] Group key not found for trigger '" << triggerToUse << "'\n";
                    continue;
                }

                // Get the isoDataList for this isoEtRange
                auto isoIt = groupedData.at(currentGroupKey).find(isoEtRange);
                if (isoIt != groupedData.at(currentGroupKey).end()) {
                    // Find the isoData that matches this pT bin
                    auto dataIt = std::find_if(isoIt->second.begin(), isoIt->second.end(),
                        [&](const DataStructures::IsolationDataWithPt& id) {
                            return std::abs(id.ptMin - pTMin) < 1e-6 && std::abs(id.ptMax - pTMax) < 1e-6;
                        });
                    if (dataIt != isoIt->second.end()) {
                        selectedIsoData = *dataIt;
                        triggerAssigned = true;
                        std::cout << "[DEBUG] Assigned to trigger '" << triggerToUse << "' for pT bin [" << pTMin << ", " << pTMax << "]\n";
                    } else {
                        std::cout << "[WARNING] Data not found for trigger '" << triggerToUse << "' with pT bin [" << pTMin << ", " << pTMax << "]\n";
                    }
                }

                // If trigger not assigned, skip this pT bin
                if (triggerAssigned) {
                    selectedDataPoints.push_back(selectedIsoData);
                } else {
                    std::cout << "[WARNING] No suitable trigger found for pT bin [" << pTMin << ", " << pTMax << "] in isoEtRange [" << isoEtRange.first << ", " << isoEtRange.second << "]\n";
                }
            }

            // Assign selectedDataPoints to combinedTriggerDataMap
            combinedTriggerDataMap[triggerGroupName][isoEtRange] = selectedDataPoints;

            // Debugging output for combined data
            std::cout << "[DEBUG] Combined data points for group '" << triggerGroupName << "', isoEtRange [" << isoEtRange.first << ", " << isoEtRange.second << "]: " << selectedDataPoints.size() << " points\n";
        }
    }

    std::cout << "[INFO] Trigger sorting and combining completed.\n";
}


void GenerateCombinedRatioPlot(
    const std::map<std::string, std::map<std::pair<float, float>, std::vector<DataStructures::IsolationDataWithPt>>>& combinedTriggerDataMap,
    const std::map<GroupKey, std::map<std::pair<float, float>, std::vector<DataStructures::IsolationDataWithPt>>>& groupedData,
    const std::string& basePlotDirectory,
    const std::vector<std::pair<double, double>>& pT_bins,
    double pTExclusionMax,
    const std::map<std::string, std::map<std::string, double>>& combinationToTriggerEfficiencyPoints) {
    // Iterate over combinedTriggerDataMap to generate plots
    for (const auto& [triggerGroupName, isoEtMap] : combinedTriggerDataMap) {
        std::cout << "\033[34m[INFO]\033[0m Processing trigger group: \033[1m" << triggerGroupName << "\033[0m\n";

        for (const auto& [isoEtRange, dataPoints] : isoEtMap) {
            const float isoMin = isoEtRange.first;
            const float isoMax = isoEtRange.second;

            // Check if this isoEtRange has data
            if (dataPoints.empty()) {
                std::cerr << "\033[31m[WARNING]\033[0m No combined data for isoEtRange ["
                          << isoMin << ", " << isoMax << "]. Skipping combined plot." << std::endl;
                continue;
            }
            std::cout << "\033[32m[DEBUG]\033[0m Found " << dataPoints.size()
                      << " data points for isoEtRange [" << isoMin << ", " << isoMax << "].\n";

            // Find associated groupKey to extract cut values and massWindowLabel
            GroupKey correspondingGroupKey;
            bool foundGroupKey = false;
            for (const auto& [gk, isoMap] : groupedData) {
                if (std::get<0>(gk) == triggerGroupName && isoMap.find(isoEtRange) != isoMap.end()) {
                    correspondingGroupKey = gk;
                    foundGroupKey = true;
                    break;
                }
            }
            if (!foundGroupKey) {
                std::cerr << "\033[31m[ERROR]\033[0m Could not find corresponding groupKey for TriggerGroupName: "
                          << triggerGroupName << " and isoEtRange: [" << isoMin << ", " << isoMax << "]. Skipping plot." << std::endl;
                continue;
            }

            float eCore = std::get<2>(correspondingGroupKey);
            float chi = std::get<3>(correspondingGroupKey);
            float asym = std::get<4>(correspondingGroupKey);
            std::string massWindowLabel = std::get<5>(correspondingGroupKey);

            // Map triggerGroupName to human-readable name
            std::string readableTriggerGroupName = Utils::getTriggerCombinationName(
                triggerGroupName, TriggerCombinationNames::triggerCombinationNameMap);
            
            // Create output directories
            std::ostringstream dirStream;
            dirStream << basePlotDirectory << "/" << triggerGroupName
                      << "/E" << Utils::formatToThreeSigFigs(eCore)
                      << "_Chi" << Utils::formatToThreeSigFigs(chi)
                      << "_Asym" << Utils::formatToThreeSigFigs(asym);
            std::string dirPath = dirStream.str();
            gSystem->mkdir(dirPath.c_str(), true);
            
            std::string isolationDir = dirPath + "/isolationEnergies";
            gSystem->mkdir(isolationDir.c_str(), true);
            
            // Create a folder for mass window type
            std::string massWindowDir = isolationDir + "/" + massWindowLabel;
            gSystem->mkdir(massWindowDir.c_str(), true);

            // Create a TCanvas for the combined plot
            std::ostringstream canvasNameStream;
            canvasNameStream << "CombinedIsolationRatio_vs_pT_" << isoMin << "_" << isoMax;
            TCanvas combinedCanvas(canvasNameStream.str().c_str(), "Combined Trigger Data", 800, 600);
            combinedCanvas.cd();

            // Prepare bin edges for variable bin widths
            std::vector<double> binEdges;
            for (const auto& bin : pT_bins) {
                if (bin.first >= pTExclusionMax) {
                    break;
                }
                binEdges.push_back(bin.first);
            }
            // Add the upper edge of the last included bin
            if (!binEdges.empty()) {
                if (pT_bins[binEdges.size() - 1].second < pTExclusionMax) {
                    binEdges.push_back(pT_bins[binEdges.size() - 1].second);
                } else {
                    binEdges.push_back(pTExclusionMax);
                }
            } else {
                // No bins to plot
                std::cerr << "\033[31m[WARNING]\033[0m No pT bins to plot. Skipping plot.\n";
                continue;
            }

            int nBins = binEdges.size() - 1;
            double* binEdgesArray = binEdges.data();

            // Create a dummy histogram to set up the axes
            TH1F* hFrame = new TH1F("hFrame", "", nBins, binEdgesArray);
            hFrame->SetStats(0);
            hFrame->GetXaxis()->SetTitle("Leading Cluster p_{T} [GeV]");
            // Set y-axis title
            const std::string yAxisTitle = (massWindowLabel == "inMassWindow") ?
                "#frac{Isolated Photons from #pi^{0}/#eta Decays}{All Photons from #pi^{0}/#eta Decays}" :
                "#frac{Isolated Prompt Photons}{All Prompt Photons}";
            hFrame->GetYaxis()->SetTitle(yAxisTitle.c_str());
            hFrame->GetYaxis()->SetRangeUser(0, 2.0);

            // Remove x-axis labels and ticks
            hFrame->GetXaxis()->SetLabelOffset(999);
            hFrame->GetXaxis()->SetTickLength(0);

            // Draw the frame
            hFrame->Draw("AXIS");

            // Create a legend
            TLegend combinedLegend(0.55, 0.67, 0.88, 0.87);
            combinedLegend.SetBorderSize(0);
            combinedLegend.SetTextSize(0.03);

            // Group data points by triggerName for coloring
            std::map<std::string, std::vector<DataStructures::IsolationDataWithPt>> dataByTrigger;
            for (const auto& isoData : dataPoints) {
                dataByTrigger[isoData.triggerName].push_back(isoData);
            }

            // Keep track of graphs to delete later
            std::vector<TGraphErrors*> graphs;

            // For each trigger, create a TGraphErrors and add to the canvas
            for (const auto& [triggerName, triggerDataPoints] : dataByTrigger) {
                std::vector<double> ptCenters;
                std::vector<double> ratios;
                std::vector<double> errors;

                for (const auto& isoData : triggerDataPoints) {
                    double ptMin = isoData.ptMin;
                    double ptMax = isoData.ptMax;

                    // Find the pT bin that matches ptMin and ptMax
                    bool foundBin = false;
                    double ptCenter = 0.0;
                    for (const auto& pT_bin : pT_bins) {
                        if (std::abs(pT_bin.first - ptMin) < 1e-6 && std::abs(pT_bin.second - ptMax) < 1e-6) {
                            // Found matching pT bin
                            ptCenter = (pT_bin.first + pT_bin.second) / 2.0;
                            foundBin = true;
                            break;
                        }
                    }
                    if (!foundBin) {
                        std::cerr << "\033[31m[WARNING]\033[0m Could not find matching pT bin for ptMin: " << ptMin << ", ptMax: " << ptMax << ". Skipping data point.\n";
                        continue;
                    }

                    // Exclude data points where ptCenter >= pTExclusionMax
                    if (ptCenter >= pTExclusionMax) {
                        continue;
                    }

                    ptCenters.push_back(ptCenter);
                    ratios.push_back(isoData.ratio);
                    errors.push_back(isoData.error);
                }

                if (ptCenters.empty()) {
                    std::cerr << "\033[31m[WARNING]\033[0m No valid data points for trigger: " << triggerName
                              << ". Skipping.\n";
                    continue;
                }

                // Create a TGraphErrors for this trigger
                TGraphErrors* graph = new TGraphErrors(ptCenters.size(),
                                                      ptCenters.data(),
                                                      ratios.data(),
                                                      nullptr,
                                                      errors.data());

                // Set marker style and color
                int markerStyle = 20;
                int markerColor = kBlack;

                auto it_color = TriggerConfig::triggerColorMap.find(triggerName);
                if (it_color != TriggerConfig::triggerColorMap.end()) {
                    markerColor = it_color->second;
                }
                graph->SetMarkerStyle(markerStyle);

                graph->SetMarkerSize(1.0);

                graph->SetLineWidth(2);
                graph->SetMarkerColor(markerColor);
                graph->SetLineColor(markerColor);

                // Draw the graph
                graph->Draw("P SAME");

                // Add entry to legend
                std::string readableTriggerName = triggerName;
                auto triggerNameIt = TriggerConfig::triggerNameMap.find(triggerName);
                if (triggerNameIt != TriggerConfig::triggerNameMap.end()) {
                    readableTriggerName = triggerNameIt->second;
                }
                combinedLegend.AddEntry(graph, readableTriggerName.c_str(), "p");
                
                std::cout << "\033[32m[DEBUG]\033[0m Added data for trigger: " << triggerName
                          << " (" << readableTriggerName << ") with " << ptCenters.size() << " points.\n";

                // Store the graph for cleanup
                graphs.push_back(graph);
            }

            // Draw a dashed line at y = 1
            TLine* combinedLine = new TLine(binEdges.front(), 1, binEdges.back(), 1);
            combinedLine->SetLineStyle(2); // Dashed line
            combinedLine->Draw("SAME");

            // Draw the legend
            combinedLegend.Draw("SAME");

            // Draw custom x-axis ticks and labels
            double xMin = binEdges.front();
            double xMax = binEdges.back();
            double yAxisMin = hFrame->GetMinimum();
            double yAxisMax = hFrame->GetMaximum();

            double tickSize = (yAxisMax - yAxisMin) * 0.02;
            double labelOffset = (yAxisMax - yAxisMin) * 0.05;
            TLatex latex;
            latex.SetTextSize(0.035);
            latex.SetTextAlign(22); // Center alignment

            // Draw x-axis line
            TLine xAxisLine(xMin, yAxisMin, xMax, yAxisMin);
            xAxisLine.Draw("SAME");

            // Draw ticks and labels at bin edges
            for (size_t i = 0; i < binEdges.size(); ++i) {
                double xPos = binEdges[i];
                double yPos = yAxisMin;

                // Draw tick
                TLine* tick = new TLine(xPos, yPos, xPos, yPos - tickSize);
                tick->Draw("SAME");

                // Get pT value for label
                double pTValue = binEdges[i];

                // Format label to show one decimal place
                std::ostringstream labelStream;
                labelStream << std::fixed << std::setprecision(1) << pTValue;
                std::string label = labelStream.str();

                // Draw label
                latex.DrawLatex(xPos, yPos - labelOffset, label.c_str());
            }

            // Redraw the axes to ensure labels are on top
            combinedCanvas.RedrawAxis();

            // Add labels using TLatex in the top-left corner
            TLatex labelText;
            labelText.SetNDC();
            labelText.SetTextSize(0.032);       // Adjust text size as needed

            TLatex valueText;
            valueText.SetNDC();
            valueText.SetTextSize(0.028);

            double xStart = 0.2; // Starting x-coordinate (left side)
            double yStartLabel = 0.9; // Starting y-coordinate
            double yStepLabel = 0.07;  // Vertical spacing between lines

            // Prepare label strings
            labelText.DrawLatex(xStart, yStartLabel, "#font[62]{Active Trigger Group:}");
            valueText.DrawLatex(xStart + 0.22, yStartLabel, readableTriggerGroupName.c_str());

            labelText.DrawLatex(xStart, yStartLabel - yStepLabel, "#font[62]{ECore #geq}");
            std::ostringstream eCoreWithUnit;
            eCoreWithUnit << eCore << "   GeV";
            valueText.DrawLatex(xStart + 0.22, yStartLabel - yStepLabel, eCoreWithUnit.str().c_str());

            labelText.DrawLatex(xStart, yStartLabel - 2 * yStepLabel, "#font[62]{#chi^{2} <}");
            std::ostringstream chiStr;
            chiStr << chi;
            valueText.DrawLatex(xStart + 0.22, yStartLabel - 2 * yStepLabel, chiStr.str().c_str());

            labelText.DrawLatex(xStart, yStartLabel - 3 * yStepLabel, "#font[62]{Asymmetry <}");
            std::ostringstream asymmetryStr;
            asymmetryStr << asym;
            valueText.DrawLatex(xStart + 0.22, yStartLabel - 3 * yStepLabel, asymmetryStr.str().c_str());

            labelText.DrawLatex(xStart, yStartLabel - 4 * yStepLabel, "#font[62]{Mass Window:}");
            valueText.DrawLatex(xStart + 0.22, yStartLabel - 4 * yStepLabel, massWindowLabel.c_str());

            labelText.DrawLatex(xStart, yStartLabel - 5 * yStepLabel, "#font[62]{#Delta R_{cone} <}");
            valueText.DrawLatex(xStart + 0.22, yStartLabel - 5 * yStepLabel, "0.3");

            // Force canvas update before saving
            combinedCanvas.Modified();
            combinedCanvas.Update();

            // Save the combined canvas
            std::ostringstream combinedOutputPathStream;
            combinedOutputPathStream << massWindowDir << "/CombinedIsolationRatio_vs_pT_" << isoMin << "_" << isoMax << ".png";
            std::string combinedOutputPath = combinedOutputPathStream.str();
            combinedCanvas.SaveAs(combinedOutputPath.c_str());
            std::cout << "\033[33m[INFO]\033[0m Saved combined plot to " << combinedOutputPath << std::endl;

            // Clean up
            delete hFrame;
            delete combinedLine;
            for (auto graph : graphs) {
                delete graph;
            }
            // The tick lines are managed by ROOT and don't need explicit deletion
        }
    }

    std::cout << "\033[34m[INFO]\033[0m Trigger sorting and combining completed.\n";
}


void ProcessIsolationData(
    const std::map<std::tuple<
        std::string, // TriggerGroupName
        std::string, // TriggerName
        float,       // ECore
        float,       // Chi
        float,       // Asymmetry
        float,       // pT Min
        float,       // pT Max
        float,       // isoMin
        float,       // isoMax
        std::string  // MassWindowLabel
    >, DataStructures::IsolationData>& dataMap,
      const std::string& basePlotDirectory,
      const std::vector<std::pair<float, float>>& exclusionRanges,
      const std::map<std::string, std::map<std::string, double>>& combinationToTriggerEfficiencyPoints,
      bool drawRefA = false,
      bool drawRefB = false) {

    // -----------------------------
    // ** isoEtRange Setup **
    // -----------------------------
    const auto& isoEtRangeColorMap = TriggerConfig::isoEtRangeColorMap;

    std::vector<std::pair<float, float>> isoEtRanges;
    std::vector<int> isoEtColors;
    for (const auto& entry : isoEtRangeColorMap) {
        isoEtRanges.push_back(entry.first);
        isoEtColors.push_back(entry.second);
    }

    // Map to hold grouped data: GroupKey -> (isoEtRange -> vector of IsolationDataWithPt)
    std::map<GroupKey, std::map<std::pair<float, float>, std::vector<DataStructures::IsolationDataWithPt>>> groupedData;

    std::cout << "\033[33m[INFO]\033[0m Starting to process isolation data..." << std::endl;

    // Populate the groupedData map
    for (const auto& entry : dataMap) {
        const auto& key = entry.first;
        const auto& isoData = entry.second;
        
        // Unpack the key
        std::string triggerGroupName = std::get<0>(key);
        std::string triggerName = std::get<1>(key);
        float eCore = std::get<2>(key);
        float chi = std::get<3>(key);
        float asym = std::get<4>(key);
        float ptMin = std::get<5>(key);
        float ptMax = std::get<6>(key);
        float isoMin = std::get<7>(key);
        float isoMax = std::get<8>(key);
        std::string massWindowLabel = std::get<9>(key);
        
        // Check if the current isoEt range is in the exclusion list
        std::pair<float, float> isoEtRange = std::make_pair(isoMin, isoMax);
        if (std::find(exclusionRanges.begin(), exclusionRanges.end(), isoEtRange) != exclusionRanges.end()) {
            continue;  // Skip excluded isoEt ranges
        }
        
        // Create the group key (without isoMin and isoMax)
        GroupKey groupKey = std::make_tuple(
            triggerGroupName,
            triggerName,
            eCore,
            chi,
            asym,
            massWindowLabel
        );
        
        // Create a data structure that includes pT info and isoEt range
        DataStructures::IsolationDataWithPt isoDataWithPt;
        isoDataWithPt.ptMin = ptMin;
        isoDataWithPt.ptMax = ptMax;
        isoDataWithPt.weightedPt = isoData.weightedPt;
        isoDataWithPt.ratio = isoData.ratio;
        isoDataWithPt.error = isoData.error;
        isoDataWithPt.isoMin = isoMin;
        isoDataWithPt.isoMax = isoMax;
        isoDataWithPt.triggerName = triggerName;
        
        // Add to the grouped data
        groupedData[groupKey][isoEtRange].push_back(isoDataWithPt);

        // Debugging output
        std::cout << "\033[32m[DEBUG]\033[0m Grouping data: TriggerGroupName: " << triggerGroupName
                  << ", TriggerName: " << triggerName
                  << ", ECore: " << eCore
                  << ", Chi: " << chi
                  << ", Asymmetry: " << asym
                  << ", MassWindowLabel: " << massWindowLabel
                  << ", isoEtRange: [" << isoMin << ", " << isoMax << "]"
                  << ". Added pT range [" << ptMin << ", " << ptMax << "] with weightedPt: " << isoDataWithPt.weightedPt
                  << ", ratio: " << isoDataWithPt.ratio << ", error: " << isoDataWithPt.error << std::endl;
    }

    std::cout << "\033[33m[INFO]\033[0m Total number of groups to process: " << groupedData.size() << std::endl;
    
    // -----------------------------
    // ** Trigger Sorting **
    // -----------------------------
    // Map to hold TriggerGroupName -> Sorted list of TriggerNames
    std::map<std::string, std::vector<std::string>> sortedTriggersByGroupName;
    // Map to hold combined data: TriggerGroupName -> isoEtRange -> vector of selected IsolationDataWithPt
    std::map<std::string, std::map<std::pair<float, float>, std::vector<DataStructures::IsolationDataWithPt>>> combinedTriggerDataMap;
    
    SortAndCombineTriggers(
        groupedData,
        combinationToTriggerEfficiencyPoints,
        sortedTriggersByGroupName,
        combinedTriggerDataMap
    );


    // Pass the map to GenerateCombinedRatioPlot if needed
    GenerateCombinedRatioPlot(
        combinedTriggerDataMap,
        groupedData,
        basePlotDirectory,
        DataStructures::pT_bins,
        20.0,
        combinationToTriggerEfficiencyPoints);

    
    // Prepare trigger combination and name maps
    const std::map<std::string, std::string>& triggerCombinationNameMap = TriggerCombinationNames::triggerCombinationNameMap;
    const std::map<std::string, std::string>& triggerNameMap = TriggerConfig::triggerNameMap;
    GeneratePerTriggerIsoPlots(
        groupedData,
        basePlotDirectory,
        isoEtRanges,
        isoEtColors,
        ReferenceData::referencePTGamma,
        ReferenceData::referenceRatio,
        ReferenceData::referenceStatError,
        ReferenceData::referenceTwoPTGamma,
        ReferenceData::referenceTwoRatio,
        ReferenceData::referenceTwoStatError,
        triggerCombinationNameMap,
        triggerNameMap,
        drawRefA,
        drawRefB,
        exclusionRanges
    );
    GeneratePerTriggerSpectraPlots(
        dataMap_inMassWindow,
        dataMap_outsideMassWindow,
        basePlotDirectory,
        triggerCombinationNameMap,
        triggerNameMap,
        exclusionRanges
    );

    std::map<SpectraGroupKey, std::map<float, CombinedSpectraData>> combinedSpectraDataMap;
    SortAndCombineSpectraData(
        dataMap_inMassWindow,
        dataMap_outsideMassWindow,
        combinationToTriggerEfficiencyPoints,
        exclusionRanges,
        DataStructures::pT_bins,    // Pass pT_bins
        20.0,                        // Pass pTExclusionMax
        combinedSpectraDataMap
    );
    GenerateCombinedSpectraPlots(
        combinedSpectraDataMap,
        basePlotDirectory,
        TriggerCombinationNames::triggerCombinationNameMap,
        DataStructures::pT_bins,    // Pass pT_bins
        20.0                        // Pass pTExclusionMax
    );

    std::cout << "\033[33m[INFO]\033[0m Finished processing isolation data." << std::endl;
}

void PlotRunByRunHistograms(
    const std::string& outputDirectory,
    const std::string& plotDirectory,
    const std::vector<std::string>& triggers,
    const std::vector<int>& runNumbers,
    const std::map<std::string, int>& triggerColorMap,
    const std::map<std::string, std::string>& triggerNameMap,
    const std::string& firmwareStatus) {

    // Create the run-by-run overlays directory
    std::string runByRunDir = plotDirectory + "/runByRun8by8overlays";
    gSystem->mkdir(runByRunDir.c_str(), true);

    // Create the directory for individual run plots
    std::string runByRunIndividualDir = plotDirectory + "/runByRunIndividual";
    gSystem->mkdir(runByRunIndividualDir.c_str(), true);

    // Determine the grid size
    const int nColumns = 9;
    const int nRows = 5;
    const int runsPerPage = nColumns * nRows;

    size_t totalRuns = runNumbers.size();
    size_t totalPages = (totalRuns + runsPerPage - 1) / runsPerPage;

    // Loop over pages
    for (size_t pageIndex = 0; pageIndex < totalPages; ++pageIndex) {
        // Create a canvas with multiple pads for overlay images
        std::ostringstream canvasName;
        canvasName << "canvas_page_" << pageIndex;
        TCanvas* canvas = new TCanvas(canvasName.str().c_str(), "Run-by-Run Overlay Plot", 2400, 1500);
        canvas->Divide(nColumns, nRows);

        // Vectors to store histograms and legends for this page
        std::vector<TH1*> pageClonedHists;
        std::vector<TLegend*> pageLegends;

        // Loop over runs in this page
        for (int padIndex = 1; padIndex <= runsPerPage; ++padIndex) {
            size_t runIndex = pageIndex * runsPerPage + padIndex - 1;
            if (runIndex >= totalRuns) {
                break; // No more runs
            }

            int runNumber = runNumbers[runIndex];
            std::string runFileName = std::to_string(runNumber) + "_HistOutput.root";
            std::string runFilePath = outputDirectory + "/" + runFileName;

            // Open the ROOT file for the run
            TFile* runFile = TFile::Open(runFilePath.c_str(), "READ");
            if (!runFile || runFile->IsZombie()) {
                std::cerr << "Error: Could not open run file " << runFilePath << std::endl;
                continue;
            }

            // Create vectors to store cloned histograms for individual plots
            std::vector<TH1*> clonedHistsIndividual;

            // Create a legend for the overlay canvas pad
            TLegend* legend = new TLegend(0.45, 0.6, 0.9, 0.9);
            legend->SetTextSize(0.04);
            legend->SetBorderSize(0);

            // Store the legend for later deletion
            pageLegends.push_back(legend);

            bool firstDraw = true;

            // Create an individual canvas for this run
            std::ostringstream individualCanvasName;
            individualCanvasName << "individualCanvas_run_" << runNumber;
            TCanvas* individualCanvas = new TCanvas(individualCanvasName.str().c_str(), "Individual Run Plot", 800, 600);
            individualCanvas->cd();

            // Set log scale if desired
            individualCanvas->SetLogy();

            // Create a legend for the individual canvas
            TLegend* individualLegend = new TLegend(0.45, 0.6, 0.9, 0.9);
            individualLegend->SetTextSize(0.04);
            individualLegend->SetBorderSize(0);

            bool firstDrawIndividual = true;

            // Loop over triggers and plot histograms
            for (const auto& trigger : triggers) {
                // Get the trigger directory
                TDirectory* triggerDir = runFile->GetDirectory(trigger.c_str());
                if (!triggerDir) {
                    std::cerr << "Trigger directory '" << trigger << "' not found in run file " << runFileName << std::endl;
                    continue;
                }

                // Construct histogram name
                std::string histName = "h8by8TowerEnergySum_" + trigger;

                // Get histogram from the trigger directory
                TH1* hist = (TH1*)triggerDir->Get(histName.c_str());
                if (!hist) {
                    std::cerr << "Warning: Histogram " << histName << " not found in run file " << runFileName << std::endl;
                    continue;
                }

                // Clone histogram to avoid modifying original
                TH1* histClone = (TH1*)hist->Clone();
                histClone->SetDirectory(0); // Detach from file

                // Store the cloned histogram for later deletion
                pageClonedHists.push_back(histClone);

                // Clone histogram for individual canvas
                TH1* histCloneIndividual = (TH1*)histClone->Clone();
                histCloneIndividual->SetDirectory(0); // Detach from file

                // Store the cloned histogram for later deletion
                clonedHistsIndividual.push_back(histCloneIndividual);

                int color = kBlack; // Default color
                auto it = triggerColorMap.find(trigger);
                if (it != triggerColorMap.end()) {
                    color = it->second;
                }

                histClone->SetLineColor(color);
                histClone->SetLineWidth(2);

                histCloneIndividual->SetLineColor(color);
                histCloneIndividual->SetLineWidth(2);

                // Set titles (optional)
                histClone->SetTitle("");
                histClone->GetXaxis()->SetTitle("");
                histClone->GetYaxis()->SetTitle("");

                histCloneIndividual->SetTitle("");
                histCloneIndividual->GetXaxis()->SetTitle("");
                histCloneIndividual->GetYaxis()->SetTitle("");

                // Adjust axis labels and titles
                histClone->GetXaxis()->SetLabelSize(0.07);
                histClone->GetYaxis()->SetLabelSize(0.07);

                histCloneIndividual->GetXaxis()->SetLabelSize(0.07);
                histCloneIndividual->GetYaxis()->SetLabelSize(0.07);

                // Plot into individual canvas
                individualCanvas->cd();
                if (firstDrawIndividual) {
                    histCloneIndividual->Draw("HIST");
                    firstDrawIndividual = false;
                } else {
                    histCloneIndividual->Draw("HIST SAME");
                }

                // Add to individual legend
                std::string displayTriggerName = trigger;
                if (triggerNameMap.find(trigger) != triggerNameMap.end()) {
                    displayTriggerName = triggerNameMap.at(trigger);
                }
                individualLegend->AddEntry(histCloneIndividual, displayTriggerName.c_str(), "l");

                // Plot into pad of overlay canvas
                canvas->cd(padIndex);
                gPad->SetLogy();

                if (firstDraw) {
                    histClone->Draw("HIST");
                    firstDraw = false;
                } else {
                    histClone->Draw("HIST SAME");
                }

                // Add to legend
                legend->AddEntry(histClone, displayTriggerName.c_str(), "l");
            }

            // Draw the legend on individual canvas
            individualCanvas->cd();
            individualLegend->Draw();

            // Draw the run number on the individual plot
            TLatex runNumberTextIndividual;
            runNumberTextIndividual.SetNDC();
            runNumberTextIndividual.SetTextAlign(13);
            runNumberTextIndividual.SetTextSize(0.08);
            runNumberTextIndividual.SetTextColor(kBlack);

            std::ostringstream runNumberStrIndividual;
            runNumberStrIndividual << "Run " << runNumber;

            runNumberTextIndividual.DrawLatex(0.1, 0.85, runNumberStrIndividual.str().c_str());

            // Update and save the individual canvas
            individualCanvas->Modified();
            individualCanvas->Update();

            std::ostringstream individualOutputFileName;
            individualOutputFileName << runByRunIndividualDir << "/Run" << runNumber << ".png";
            individualCanvas->SaveAs(individualOutputFileName.str().c_str());
            std::cout << "Saved individual run plot to " << individualOutputFileName.str() << std::endl;

            // Clean up individual canvas and legend
            delete individualCanvas;
            delete individualLegend;

            // Clean up cloned histograms for individual canvas
            for (auto hist : clonedHistsIndividual) {
                delete hist;
            }

            // Draw the legend on the pad
            canvas->cd(padIndex);
            legend->Draw();

            // Draw the run number on the pad
            TLatex runNumberText;
            runNumberText.SetNDC();
            runNumberText.SetTextAlign(13);
            runNumberText.SetTextSize(0.08);
            runNumberText.SetTextColor(kBlack);

            std::ostringstream runNumberStr;
            runNumberStr << "Run " << runNumber;

            runNumberText.DrawLatex(0.1, 0.85, runNumberStr.str().c_str());

            // Close the run file
            runFile->Close();
            delete runFile;
        }

        // Draw the firmware status on the overlay canvas
        if (!firmwareStatus.empty()) {
            canvas->cd();
            TLatex firmwareStatusText;
            firmwareStatusText.SetNDC();
            firmwareStatusText.SetTextAlign(22); // Centered
            firmwareStatusText.SetTextSize(0.02);
            firmwareStatusText.SetTextColor(kBlack);
            firmwareStatusText.DrawLatex(0.5, 0.95, firmwareStatus.c_str());
        }

        // Update the canvas
        canvas->Modified();
        canvas->Update();

        // Save the canvas
        std::ostringstream outputFileName;
        outputFileName << runByRunDir << "/RunOverlay_Page" << pageIndex + 1 << ".png";
        canvas->SaveAs(outputFileName.str().c_str());
        std::cout << "Saved run-by-run overlay plot to " << outputFileName.str() << std::endl;

        // Clean up histograms and legends for this page
        for (auto hist : pageClonedHists) {
            delete hist;
        }

        for (auto legend : pageLegends) {
            delete legend;
        }

        // Clean up
        delete canvas;
    }
}

void PlotCombinedHistograms(
    const std::string& outputDirectory,
    const std::vector<std::string>& combinedRootFiles,
    const std::map<std::string, std::vector<int>>& combinationToValidRuns,
    const std::string& fitFunctionType = "sigmoid", bool fitOnly = false) {
    
    // List of all triggers
    const std::vector<std::string>& allTriggers = TriggerConfig::allTriggers;
    const std::vector<std::string>& photonTriggers = TriggerConfig::photonTriggers;
    const auto& triggerColorMap = TriggerConfig::triggerColorMap;
    const auto& triggerNameMap = TriggerConfig::triggerNameMap;
    std::map<std::string, double> triggerEfficiencyPoints;
    std::map<std::string, std::map<std::string, double>> combinationToTriggerEfficiencyPoints;

    // Base directory for plots
    std::string basePlotDirectory = "/Users/patsfan753/Desktop/DirectPhotonAna/Plots";

    // Create the base plot directory if it doesn't exist
    gSystem->mkdir(basePlotDirectory.c_str(), true);
    

    // Function to draw run numbers on canvas
    auto drawRunNumbersOnCanvas = [](const std::vector<int>& runNumbers, const std::string& firmwareStatus) {
        std::cout << "drawRunNumbersOnCanvas: Number of runs = " << runNumbers.size() << std::endl;

        // Set up TLatex for run numbers
        TLatex runNumbersLatex;
        runNumbersLatex.SetNDC();
        runNumbersLatex.SetTextAlign(13); // Align at top left
        runNumbersLatex.SetTextColor(kBlack);

        // Define configuration variables
        double xStart, xEnd, yStart, yEnd, textSize;
        double xSpacingFactor, ySpacingFactor;
        int numColumns;

        // Customize settings based on the number of runs
        if (runNumbers.size() == 692) {
            numColumns = 27;
            textSize = 0.012;
            xStart = 0.54; xEnd = 0.93;
            yStart = 0.9; yEnd = 0.45;
            xSpacingFactor = 0.8;
            ySpacingFactor = 0.8;
        } else if (runNumbers.size() == 71) {
            numColumns = 5;
            textSize = 0.022;
            xStart = 0.54; xEnd = 0.93;
            yStart = 0.9; yEnd = 0.45;
            xSpacingFactor = 0.8;
            ySpacingFactor = 0.8;
        } else if (runNumbers.size() == 125) {
            numColumns = 8;
            textSize = 0.022;
            xStart = 0.52; xEnd = 0.93;
            yStart = 0.9; yEnd = 0.45;
            xSpacingFactor = 1.0;
            ySpacingFactor = 1.0;
        } else if (runNumbers.size() == 54) {
            numColumns = 3;
            textSize = 0.025;
            xStart = 0.57; xEnd = 0.94;
            yStart = 0.9; yEnd = 0.42;
            xSpacingFactor = 0.8;
            ySpacingFactor = 0.92;
        } else if (runNumbers.size() == 126) {
            numColumns = 8;
            textSize = 0.022;
            xStart = 0.52; xEnd = 0.93;
            yStart = 0.9; yEnd = 0.45;
            xSpacingFactor = 1.0;
            ySpacingFactor = 1.0;
        } else if (runNumbers.size() == 146) {
            numColumns = 8;
            textSize = 0.018;
            xStart = 0.54; xEnd = 0.93;
            yStart = 0.9; yEnd = 0.5;
            xSpacingFactor = 1.0;
            ySpacingFactor = 1.0;
        } else {
            numColumns = (runNumbers.size() <= 20) ? 2 : 5;
            textSize = 0.03;
            xStart = 0.7; xEnd = 0.9;
            yStart = 0.85; yEnd = 0.3;
            xSpacingFactor = 1.0;
            ySpacingFactor = 1.0;
        }
        double headerTextSize = 0.03;
        runNumbersLatex.SetTextSize(headerTextSize);
        std::ostringstream headerText;
        headerText << runNumbers.size() << " runs";
        if (!firmwareStatus.empty()) {
            headerText << " (" << firmwareStatus << ")";
        }
        runNumbersLatex.DrawLatex(0.42, 0.9, headerText.str().c_str());

        // Skip plotting run numbers for specific size
        if (runNumbers.size() == 257 || runNumbers.size() == 251 || runNumbers.size() == 102 || runNumbers.size() == 101 || runNumbers.size() == 240 || runNumbers.size() == 142 || runNumbers.size() == 440 || runNumbers.size() == 443 || runNumbers.size() == 725 || runNumbers.size() == 100 || runNumbers.size() == 250 || runNumbers.size() == 141 || runNumbers.size() == 254) {
            std::cout << "[INFO] Skipping run number plotting for runNumbers.size() = 692." << std::endl;
            return;
        }

        runNumbersLatex.SetTextSize(textSize);

        int numRows = (runNumbers.size() + numColumns - 1) / numColumns;

        // Create a grid to arrange run numbers
        std::vector<std::vector<std::string>> grid(numRows, std::vector<std::string>(numColumns, ""));

        // Fill the grid with run numbers
        int runIndex = 0;
        for (int col = 0; col < numColumns; ++col) {
            for (int row = 0; row < numRows; ++row) {
                if (runIndex < runNumbers.size()) {
                    grid[row][col] = std::to_string(runNumbers[runIndex]);
                    ++runIndex;
                }
            }
        }

        // Calculate spacing based on coordinates and number of rows/columns
        double xSpacing = xSpacingFactor * (xEnd - xStart) / numColumns;
        double ySpacing = ySpacingFactor * (yStart - yEnd) / (numRows + 1); // +1 for spacing

        // Draw run numbers in the grid
        for (int row = 0; row < numRows; ++row) {
            double yPos = yStart - (row + 1) * ySpacing;
            for (int col = 0; col < numColumns; ++col) {
                if (!grid[row][col].empty()) {
                    double xPos = xStart + col * xSpacing;
                    runNumbersLatex.DrawLatex(xPos, yPos, grid[row][col].c_str());
                }
            }
        }
    };
    // Loop over combined ROOT files
    for (const auto& rootFileName : combinedRootFiles) {
        std::string rootFilePath = outputDirectory + "/" + rootFileName;

        // Extract triggers from the filename
        std::vector<std::string> triggers = ExtractTriggersFromFilename(rootFileName, allTriggers);

        // Check for firmware tag in the filename
        std::string firmwareTag;
        if (Utils::EndsWith(rootFileName, "_beforeTriggerFirmwareUpdate_Combined.root")) {
            firmwareTag = "_beforeTriggerFirmwareUpdate";
        } else if (Utils::EndsWith(rootFileName, "_afterTriggerFirmwareUpdate_Combined.root")) {
            firmwareTag = "_afterTriggerFirmwareUpdate";
        }
        
        // Extract firmware status from firmwareTag
        std::string firmwareStatus;
        if (firmwareTag == "_beforeTriggerFirmwareUpdate") {
            firmwareStatus = "#bf{Before firmware update at run 47289}";
        } else if (firmwareTag == "_afterTriggerFirmwareUpdate") {
            firmwareStatus = "#bf{After firmware update at run 47289}";
        } else {
            firmwareStatus = "";
        }

        std::cout << "Processing file: " << rootFileName << std::endl;
        std::cout << "Triggers found: ";
        for (const auto& trigger : triggers) {
            std::cout << trigger << " ";
        }
        std::cout << std::endl;

        // Reconstruct combinationName with firmware tag if present
        std::string combinationName;
        for (const auto& trigger : triggers) {
            combinationName += trigger + "_";
        }
        // Remove the trailing underscore
        if (!combinationName.empty()) {
            combinationName.pop_back();
        }
        // Append firmware tag
        combinationName += firmwareTag;

        // Sanitize combinationName for use as a directory name
        std::string sanitizedCombinationName = combinationName;
        std::replace(sanitizedCombinationName.begin(), sanitizedCombinationName.end(), '/', '_');
        std::replace(sanitizedCombinationName.begin(), sanitizedCombinationName.end(), ' ', '_');

        // Create the subdirectory for this trigger combination
        std::string plotDirectory = basePlotDirectory + "/" + sanitizedCombinationName;
        gSystem->mkdir(plotDirectory.c_str(), true);

        // Open the ROOT file
        TFile* inputFile = TFile::Open(rootFilePath.c_str(), "READ");
        if (!inputFile || inputFile->IsZombie()) {
            std::cerr << "Error: Could not open file " << rootFilePath << std::endl;
            continue;
        }
        // -------------------- Overlay 8x8 Tower NRG --------------------
        // Create a canvas
        TCanvas* canvas = new TCanvas("canvas", "Overlay Plot", 800, 600);
        TLegend* legend = new TLegend(0.5, 0.58, 0.8, 0.88);
        legend->SetTextSize(0.028);
        canvas->SetLogy();
        bool firstDraw = true;
        // Loop over triggers and plot histograms
        for (const auto& trigger : triggers) {
            TDirectory* triggerDir = inputFile->GetDirectory(trigger.c_str());
            if (!triggerDir) {
                std::cerr << "Trigger directory '" << trigger << "' not found in file " << rootFileName << std::endl;
                continue;
            }
            std::string histName = "h8by8TowerEnergySum_" + trigger;

            TH1* hist = (TH1*)triggerDir->Get(histName.c_str());
            if (!hist) {
                std::cerr << "Warning: Histogram " << histName << " not found in file " << rootFileName << std::endl;
                continue;
            }
            // Clone histogram to avoid modifying original
            TH1* histClone = (TH1*)hist->Clone();
            histClone->SetDirectory(0); // Detach from file

            int color = kBlack; // Default color
            auto it = TriggerConfig::triggerColorMap.find(trigger);
            if (it != TriggerConfig::triggerColorMap.end()) {
                color = it->second;
            }
            histClone->SetLineColor(color);
            histClone->SetLineWidth(2);

            histClone->SetTitle(("Overlay for " + combinationName).c_str());
            histClone->GetXaxis()->SetTitle("Maximum 8x8 EMCal Tower Energy Sum [GeV]");
            histClone->GetYaxis()->SetTitle("Events");

            if (firstDraw) {
                histClone->Draw("HIST");
                firstDraw = false;
            } else {
                histClone->Draw("HIST SAME");
            }
            std::string displayTriggerName = trigger;
            auto it_name = TriggerConfig::triggerNameMap.find(trigger);
            if (it_name != TriggerConfig::triggerNameMap.end()) {
                displayTriggerName = it_name->second;
            }

            // Add histogram to legend with display name
            legend->AddEntry(histClone, displayTriggerName.c_str(), "l");
        }

        legend->Draw();

        auto it = combinationToValidRuns.find(combinationName);
        if (it != combinationToValidRuns.end()) {
            const std::vector<int>& validRuns = it->second;
            // Draw run numbers on canvas using the drawRunNumbersOnCanvas function
            drawRunNumbersOnCanvas(validRuns, firmwareStatus);
        } else {
            // If no valid runs info, indicate it
            TLatex text;
            text.SetNDC();
            text.SetTextAlign(13); // Align at top left
            text.SetTextSize(0.03);
            text.DrawLatex(0.15, 0.6, "Run numbers not available");
        }


        // Update the canvas
        canvas->Modified();
        canvas->Update();

        // Save the canvas
        std::string outputFileName = plotDirectory + "/h8by8TowerEnergySum_Overlay.png";
        canvas->SaveAs(outputFileName.c_str());
        std::cout << "Saved plot to " << outputFileName << std::endl;

        // Clean up
        delete canvas;

        // -------------------- Turn-On Plot --------------------
        // Identify photon triggers in the combination
        std::vector<std::string> photonTriggersInCombination;
        for (const auto& trigger : triggers) {
            if (std::find(photonTriggers.begin(), photonTriggers.end(), trigger) != photonTriggers.end()) {
                photonTriggersInCombination.push_back(trigger);
            }
        }
        // Proceed only if there are photon triggers
        if (!photonTriggersInCombination.empty()) {
            TDirectory* minbiasDir = inputFile->GetDirectory("MBD_NandS_geq_1");
            if (!minbiasDir) {
                std::cerr << "Warning: Minbias directory 'MBD_NandS_geq_1' not found in file " << rootFileName << std::endl;
            } else {
                // Get the minbias histogram
                std::string minbiasHistName = "h8by8TowerEnergySum_MBD_NandS_geq_1";
                TH1* minbiasHist = (TH1*)minbiasDir->Get(minbiasHistName.c_str());
                if (!minbiasHist) {
                    std::cerr << "Warning: Minbias histogram " << minbiasHistName << " not found in file " << rootFileName << std::endl;
                } else {
                    // Create a canvas for the turn-on plot
                    TCanvas* canvasTurnOn = new TCanvas("canvasTurnOn", "Turn-On Plot", 800, 600);
                    TLegend* legendTurnOn = new TLegend(0.18, 0.72, 0.45, 0.9);
                    legendTurnOn->SetTextSize(0.028);
                    bool firstDrawTurnOn = true;
                    
                    // Loop over photon triggers and plot ratios
                    for (const auto& photonTrigger : photonTriggersInCombination) {
                        // Get the photon trigger histogram from its directory
                        TDirectory* photonDir = inputFile->GetDirectory(photonTrigger.c_str());
                        if (!photonDir) {
                            std::cerr << "Warning: Photon trigger directory '" << photonTrigger << "' not found in file " << rootFileName << std::endl;
                            continue;
                        }
                        // Get the photon trigger histogram
                        std::string photonHistName = "h8by8TowerEnergySum_" + photonTrigger;
                        TH1* photonHist = (TH1*)photonDir->Get(photonHistName.c_str());
                        if (!photonHist) {
                            std::cerr << "Warning: Histogram " << photonHistName << " not found in file " << rootFileName << std::endl;
                            continue;
                        }
                        
                        // Compute the ratio with proper error propagation
                        std::string ratioHistName = "ratio_" + photonTrigger;
                        TH1* ratioHist = (TH1*)photonHist->Clone(ratioHistName.c_str());
                        ratioHist->SetDirectory(0); // Detach from file
                        
                        // Perform the division with proper error propagation
                        ratioHist->Divide(photonHist, minbiasHist, 1.0, 1.0, "B"); // "B" for binomial errors
                        
                        int color = kBlack; // Default color
                        auto it = TriggerConfig::triggerColorMap.find(photonTrigger);
                        if (it != TriggerConfig::triggerColorMap.end()) {
                            color = it->second;
                        }

                        ratioHist->SetMarkerStyle(20); // Filled circle
                        ratioHist->SetMarkerColor(color);
                        ratioHist->SetLineColor(color);
                        
                        ratioHist->SetTitle(("Turn-On Curve for " + combinationName).c_str());
                        ratioHist->GetXaxis()->SetTitle("Maximum 8x8 Energy Sum (EMCal) [GeV]");
                        ratioHist->GetYaxis()->SetTitle("Ratio to Minbias");
                        ratioHist->GetYaxis()->SetRangeUser(0, 2.0);
                        
                        if (firstDrawTurnOn) {
                            ratioHist->Draw("E1"); // Draw with error bars
                            firstDrawTurnOn = false;
                        } else {
                            ratioHist->Draw("E1 SAME");
                        }
                        
                        // Draw a dashed line at y = 1
                        TLine* line = new TLine(ratioHist->GetXaxis()->GetXmin(), 1, ratioHist->GetXaxis()->GetXmax(), 1);
                        line->SetLineStyle(1); // Dashed line
                        line->SetLineColor(kBlack); // Black color
                        line->Draw("SAME");
                        
                        std::string displayPhotonTriggerName = photonTrigger;
                        auto it_name = TriggerConfig::triggerNameMap.find(photonTrigger);
                        if (it_name != TriggerConfig::triggerNameMap.end()) {
                            displayPhotonTriggerName = it_name->second;
                        }

                        
                        std::ostringstream legendEntry;
                        legendEntry << displayPhotonTriggerName;
                        
                        // Prepare the key using combinationName with firmware tag
                        std::pair<std::string, std::string> key = std::make_pair(combinationName, photonTrigger);

                        // Try to find the fit parameters for the combination and trigger
                        auto it_fitParams = TriggerConfig::triggerFitParameters.find(key);

                        if (it_fitParams == TriggerConfig::triggerFitParameters.end()) {
                            // If not found, try with combinationName without firmware tag
                            std::string combinationNameWithoutFirmware = Utils::stripFirmwareTag(combinationName);
                            key = std::make_pair(combinationNameWithoutFirmware, photonTrigger);
                            it_fitParams = TriggerConfig::triggerFitParameters.find(key);
                        }

                        if (it_fitParams == TriggerConfig::triggerFitParameters.end()) {
                            // If not found, try with empty combinationName
                            key = std::make_pair("", photonTrigger);
                            it_fitParams = TriggerConfig::triggerFitParameters.find(key);
                        }


                        if (enableFits && it_fitParams != TriggerConfig::triggerFitParameters.end()) {
                            DataStructures::FitParameters params = it_fitParams->second;

                            TF1* fitFunc = nullptr;
                            if (fitFunctionType == "sigmoid") {
                                fitFunc = Utils::sigmoidFit(("fit_" + photonTrigger).c_str(), 0.0, 20.0,
                                                            params.amplitudeEstimate, params.slopeEstimate, params.xOffsetEstimate,
                                                            params.amplitudeMin, params.amplitudeMax,
                                                            params.slopeMin, params.slopeMax,
                                                            params.xOffsetMin, params.xOffsetMax);
                            } else if (fitFunctionType == "erf") {
                                fitFunc = Utils::erfFit(("fit_" + photonTrigger).c_str(), 0.0, 20.0,
                                                        params.amplitudeEstimate, params.xOffsetEstimate, params.sigmaEstimate,
                                                        params.amplitudeMin, params.amplitudeMax,
                                                        params.xOffsetMin, params.xOffsetMax,
                                                        params.sigmaMin, params.sigmaMax);
                            } else {
                                // Default to sigmoid if unrecognized fitFunctionType
                                fitFunc = Utils::sigmoidFit(("fit_" + photonTrigger).c_str(), 0.0, 20.0,
                                                            params.amplitudeEstimate, params.slopeEstimate, params.xOffsetEstimate,
                                                            params.amplitudeMin, params.amplitudeMax,
                                                            params.slopeMin, params.slopeMax,
                                                            params.xOffsetMin, params.xOffsetMax);
                            }

                            fitFunc->SetLineColor(color);
                            ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
                            ratioHist->Fit(fitFunc, "R");
                            fitFunc->Draw("SAME");

                            // Retrieve the fit parameters regardless of convergence
                            double A = fitFunc->GetParameter(0);
                            double A_error = fitFunc->GetParError(0);

                            double x99 = 0;
                            double x99_error = 0;

                            if (fitFunctionType == "sigmoid") {
                                double k = fitFunc->GetParameter(1);
                                double x0 = fitFunc->GetParameter(2);
                                double k_error = fitFunc->GetParError(1);
                                double x0_error = fitFunc->GetParError(2);

                                // Calculate the 99% efficiency point for sigmoid
                                x99 = x0 + (std::log(99) / k);

                                // Propagate error for x99 using partial derivatives
                                x99_error = std::sqrt(
                                    (x0_error * x0_error) +
                                    (std::pow((std::log(99) / (k * k)), 2) * k_error * k_error)
                                );

                                // Debugging output
                                std::cout << "Sigmoid Fit for " << photonTrigger
                                          << ": A = " << A << " ± " << A_error
                                          << ", k = " << k << " ± " << k_error
                                          << ", x0 = " << x0 << " ± " << x0_error
                                          << ", x99 = " << x99 << " ± " << x99_error << " GeV" << std::endl;

                            } else if (fitFunctionType == "erf") {
                                double x0 = fitFunc->GetParameter(1);
                                double sigma = fitFunc->GetParameter(2);
                                double x0_error = fitFunc->GetParError(1);
                                double sigma_error = fitFunc->GetParError(2);

                                // Calculate the 99% efficiency point for erf
                                double erfInvValue = TMath::ErfInverse(0.98); // erfInv(0.98)
                                x99 = x0 + sqrt(2) * sigma * erfInvValue;

                                // Propagate error for x99
                                x99_error = std::sqrt(
                                    x0_error * x0_error +
                                    (sqrt(2) * erfInvValue * sigma_error) * (sqrt(2) * erfInvValue * sigma_error)
                                );

                                // Debugging output
                                std::cout << "Erf Fit for " << photonTrigger
                                          << ": A = " << A << " ± " << A_error
                                          << ", x0 = " << x0 << " ± " << x0_error
                                          << ", sigma = " << sigma << " ± " << sigma_error
                                          << ", x99 = " << x99 << " ± " << x99_error << " GeV" << std::endl;
                            }

                            // Store the 99% efficiency point
                            triggerEfficiencyPoints[photonTrigger] = x99;
                            combinationToTriggerEfficiencyPoints[combinationName][photonTrigger] = x99;
                            // Append fit parameters to the legend entry with errors
                            legendEntry << ", 99% efficiency = " << std::fixed << std::setprecision(2) << x99 << " GeV";

                            // Draw a dashed vertical line at x99
                            if (x99 > ratioHist->GetXaxis()->GetXmin() && x99 < ratioHist->GetXaxis()->GetXmax()) {
                                TLine* verticalLine = new TLine(x99, 0, x99, 1); // Draw line up to y = 1
                                verticalLine->SetLineStyle(2); // Dashed line
                                verticalLine->SetLineColor(color); // Use the color of the current trigger
                                verticalLine->SetLineWidth(5); // Set the line width
                                verticalLine->Draw("SAME");
                            }
                        } else {
                            std::cerr << "No fit parameters found for trigger " << photonTrigger << " in combination " << combinationName << std::endl;
                        }

                        // Add the updated legend entry with the fit parameters
                        legendTurnOn->AddEntry(ratioHist, legendEntry.str().c_str(), "p");
                        std::cout << "Legend Entry Added: " << legendEntry.str() << std::endl;
                    }
                    
                    // Draw the legend
                    legendTurnOn->Draw();
                    
                    
                    // Add a separate legend for the 99% efficiency line
                    TLegend* legendEfficiencyLine = new TLegend(0.18, 0.62, 0.38, 0.72);
                    legendEfficiencyLine->SetTextSize(0.03);
                    legendEfficiencyLine->SetBorderSize(0);
                    legendEfficiencyLine->SetFillStyle(0);
                    
                    // Create a dummy line for the 99% efficiency point
                    TLine* dummyLine = new TLine(0, 0, 0, 0); // The coordinates do not matter since it is a dummy :(
                    dummyLine->SetLineStyle(2); // Dashed line
                    dummyLine->SetLineColor(kGray + 1); // Gray color
                    dummyLine->SetLineWidth(2);
                    
                    // Add the dummy line to the legend
                    legendEfficiencyLine->AddEntry(dummyLine, "99% Efficiency Point", "l"); // "l" for line
                    legendEfficiencyLine->SetLineWidth(2);
                    legendEfficiencyLine->Draw();
                    
                    // Draw the firmware status on the overlay canvas
                    if (!firmwareStatus.empty()) {
                        canvasTurnOn->cd();
                        TLatex firmwareStatusText;
                        firmwareStatusText.SetNDC();
                        firmwareStatusText.SetTextAlign(22); // Centered
                        firmwareStatusText.SetTextSize(0.03);
                        firmwareStatusText.SetTextColor(kBlack);
                        firmwareStatusText.DrawLatex(0.5, 0.96, firmwareStatus.c_str());
                    }

                    
                    canvasTurnOn->Modified();
                    canvasTurnOn->Update();
                    
                    // Save the canvas
                    std::string outputTurnOnFileName = plotDirectory + "/h8by8TowerEnergySum_TurnOn.png";
                    canvasTurnOn->SaveAs(outputTurnOnFileName.c_str());
                    std::cout << "Saved turn-on plot to " << outputTurnOnFileName << std::endl;
                    
                    // Clean up
                    delete canvasTurnOn;
                }
            }
        }
        if (fitOnly) {
            // Close the ROOT file and continue to the next iteration
            inputFile->Close();
            delete inputFile;
            
            // Now, generate run-by-run overlays
            auto it_run = combinationToValidRuns.find(combinationName);
            // Call PlotRunByRunHistograms with firmwareStatus
            if (it_run != combinationToValidRuns.end()) {
                const std::vector<int>& validRuns = it_run->second;
                PlotRunByRunHistograms(outputDirectory, plotDirectory, triggers, validRuns, triggerColorMap, triggerNameMap, firmwareStatus);
            } else {
                std::cout << "No valid runs found for combination: " << combinationName << std::endl;
            }
            
            continue; // Skip to the next ROOT file
        }

        // The following code will only execute if fitOnly is false
        // -------------------- Process Invariant Mass Histograms --------------------
        ProcessInvariantMassHistograms(inputFile, plotDirectory, triggers, triggerColorMap, combinationName);

        // -------------------- Process Meson Mass vs Pt --------------------
        ProcessMesonMassVsPt(plotDirectory, combinationName, triggers, triggerEfficiencyPoints, DataStructures::pT_bins);

        // -------------------- Process Isolation Energy Histograms --------------------
        ProcessIsolationEnergyHistogramsWithCuts(inputFile, plotDirectory, triggers, combinationName);

        inputFile->Close();
        delete inputFile;
        
        // Now, generate run-by-run overlays
        auto it_run = combinationToValidRuns.find(combinationName);
        // Call PlotRunByRunHistograms with firmwareStatus
        if (it_run != combinationToValidRuns.end()) {
            const std::vector<int>& validRuns = it_run->second;
            PlotRunByRunHistograms(outputDirectory, plotDirectory, triggers, validRuns, triggerColorMap, triggerNameMap, firmwareStatus);
        } else {
            std::cout << "No valid runs found for combination: " << combinationName << std::endl;
        }
    }

    std::string csvOutputPath = "/Users/patsfan753/Desktop/isolation_data.csv";
    WriteIsolationDataToCSV(csvOutputPath);
    readDataFromCSV(csvOutputPath, dataMap_inMassWindow, dataMap_outsideMassWindow);

    std::vector<std::pair<float, float>> exclusionRanges = {
        {-100, 10},
        {-10, 0},
        {0, 10}
    };

    ProcessIsolationData(
        dataMap_inMassWindow,
        basePlotDirectory,
        exclusionRanges,
        combinationToTriggerEfficiencyPoints,
        false,
        true);

    ProcessIsolationData(
        dataMap_outsideMassWindow,
        basePlotDirectory,
        exclusionRanges,
        combinationToTriggerEfficiencyPoints,
        true,
        false);
    
}

void AddLabelsToCanvas(
    const std::string& triggerGroupName,
    const std::string& displayTriggerName,
    const std::string& pTMin,
    const std::string& pTMax,
    const std::string& subtractionType) {

    TLatex labelText, valueText;
    labelText.SetNDC();
    labelText.SetTextSize(0.027);
    labelText.SetTextColor(kRed);
    labelText.SetTextFont(62);

    valueText.SetNDC();
    valueText.SetTextSize(0.027);
    valueText.SetTextColor(kBlack);
    valueText.SetTextFont(42);

    labelText.DrawLatex(0.18, 0.9, "Active Trigger Group:");

    // Map triggerGroupName to a readable form using triggerCombinationNameMap
    std::string readableTriggerGroupName = triggerGroupName;
    auto it = TriggerCombinationNames::triggerCombinationNameMap.find(triggerGroupName);
    if (it != TriggerCombinationNames::triggerCombinationNameMap.end()) {
        readableTriggerGroupName = it->second;
    }

    valueText.DrawLatex(0.37, 0.9, readableTriggerGroupName.c_str());

    labelText.DrawLatex(0.18, 0.85, "Trigger:");
    valueText.DrawLatex(0.3, 0.85, displayTriggerName.c_str());

    if (!pTMin.empty() && !pTMax.empty()) {
        labelText.DrawLatex(0.18, 0.8, "pT Range:");
        valueText.DrawLatex(0.3, 0.8, (pTMin + " to " + pTMax + " GeV").c_str());
    } else {
        labelText.DrawLatex(0.18, 0.8, "pT Range:");
        valueText.DrawLatex(0.3, 0.8, "All pT");
    }

    labelText.DrawLatex(0.18, 0.75, "UE Subtracted:");
    valueText.DrawLatex(0.34, 0.75, subtractionType.c_str());
}



void Process1DHistogram(
    TH1* hist,
    const std::smatch& match,
    const std::string& isolationDir,
    const std::string& triggerGroupName,
    bool hasPtBins = true) {

    std::string subtractionType = match[1].str();
    std::string pTMin, pTMax, triggerName;

    if (hasPtBins) {
        pTMin = match[2].str();
        pTMax = match[3].str();
        triggerName = match[4].str();
    } else {
        triggerName = match[2].str();
    }

    // Map the trigger name to a human-readable format
    std::string displayTriggerName = TriggerConfig::triggerNameMap.count(triggerName) > 0
                                     ? TriggerConfig::triggerNameMap.at(triggerName)
                                     : triggerName;

    // Construct the output directory path
    std::string subDir = isolationDir + "/" + subtractionType;
    gSystem->mkdir(subDir.c_str(), true);

    std::string outputFilePath;

    if (hasPtBins) {
        std::string pTDir = subDir + "/pT_" + pTMin + "_to_" + pTMax;
        gSystem->mkdir(pTDir.c_str(), true);
        outputFilePath = pTDir + "/" + hist->GetName() + ".png";
    } else {
        outputFilePath = subDir + "/" + hist->GetName() + ".png";
    }

    // Draw the histogram on a canvas with log scale
    TCanvas canvas("canvas", "Histogram Canvas", 800, 600);
    canvas.SetLogy();
    hist->Draw("HIST");

    // Add TLatex labels
    AddLabelsToCanvas(triggerGroupName, displayTriggerName, pTMin, pTMax, subtractionType);
    canvas.SaveAs(outputFilePath.c_str());

    std::cout << "Saved 1D histogram to: " << outputFilePath << std::endl;
}


void Process2DHistogram(
    TH2* hist,
    const std::smatch& match,
    const std::string& isolationDir,
    const std::string& triggerGroupName,
    bool hasPtBins = true) {

    std::string subtractionType = match[1].str();
    std::string pTMin, pTMax, triggerName;

    if (hasPtBins) {
        pTMin = match[2].str();
        pTMax = match[3].str();
        triggerName = match[4].str();
    } else {
        triggerName = match[2].str();
    }

    // Map the trigger name to a human-readable format
    std::string displayTriggerName = TriggerConfig::triggerNameMap.count(triggerName) > 0
                                     ? TriggerConfig::triggerNameMap.at(triggerName)
                                     : triggerName;

    // Construct the output directory path
    std::string subDir = isolationDir + "/" + subtractionType;
    gSystem->mkdir(subDir.c_str(), true);

    std::string outputFilePath;
    std::string saveDir; // Variable to store the directory where projections will be saved

    if (hasPtBins) {
        std::string pTDir = subDir + "/pT_" + pTMin + "_to_" + pTMax;
        gSystem->mkdir(pTDir.c_str(), true);
        outputFilePath = pTDir + "/" + hist->GetName() + ".png";
        saveDir = pTDir; // Assign pTDir to saveDir for use in projections
    } else {
        outputFilePath = subDir + "/" + hist->GetName() + ".png";
        saveDir = subDir; // Assign subDir to saveDir for use in projections
    }

    double clusterETMin, clusterETMax;
    double specificClusterET, specificIsoET;
    double windowSize = 0.5; // Define a window size around the specific values

    if (hasPtBins) {
        // Convert pTMin and pTMax to double
        try {
            double pTMinVal = std::stod(pTMin);
            double pTMaxVal = std::stod(pTMax);
            clusterETMin = pTMinVal;
            clusterETMax = pTMaxVal;
        }
        catch (const std::invalid_argument& e) {
            std::cerr << "\033[31m[ERROR]\033[0m Invalid pT bin values: pTMin = '" << pTMin
                      << "', pTMax = '" << pTMax << "'. Skipping histogram '" << hist->GetName() << "'.\n";
            return;
        }

        // Define specificClusterET as the midpoint of the pT bin
        specificClusterET = (clusterETMin + clusterETMax) / 2.0;
        specificIsoET = 0.0; // Focus on 0 isolation energy for projections
    }
    else {
        // For histograms without pT bins, use the full x-axis range
        clusterETMin = hist->GetXaxis()->GetXmin();
        clusterETMax = hist->GetXaxis()->GetXmax();

        // Define specificClusterET as the median or a meaningful value within the range
        specificClusterET = (clusterETMin + clusterETMax) / 2.0;

        // Define specificIsoET based on histogram's y-axis or set to 0
        specificIsoET = 0.0; // Modify as needed
    }

    // Adjust x-axis range for the pT bin if applicable
    if (hasPtBins) {
        hist->GetXaxis()->SetRangeUser(clusterETMin, clusterETMax);
    }

    // Draw the 2D histogram on a canvas
    TCanvas canvas("canvas2D", "2D Histogram Canvas", 800, 600);
    canvas.SetLogz();
    hist->Draw("COLZ");

    // Add TLatex labels
    AddLabelsToCanvas(triggerGroupName, displayTriggerName, hasPtBins ? pTMin : "", hasPtBins ? pTMax : "", subtractionType);

    // Draw axis titles
    hist->GetXaxis()->SetTitle("Cluster E_{T} (GeV)");
    hist->GetYaxis()->SetTitle("Cluster E_{T}^{iso} (GeV)");
    hist->GetZaxis()->SetTitle("Counts");

    // ------------------------------
    // Add Projections for Specific Events
    // ------------------------------

    // Create lines on the 2D histogram to indicate the specific values
    TLine* lineX_lower = new TLine(specificClusterET - windowSize, hist->GetYaxis()->GetXmin(),
                                   specificClusterET - windowSize, hist->GetYaxis()->GetXmax());
    lineX_lower->SetLineColor(kRed);
    lineX_lower->SetLineStyle(2); // Dashed line
    lineX_lower->Draw("Same");

    TLine* lineX_upper = new TLine(specificClusterET + windowSize, hist->GetYaxis()->GetXmin(),
                                   specificClusterET + windowSize, hist->GetYaxis()->GetXmax());
    lineX_upper->SetLineColor(kRed);
    lineX_upper->SetLineStyle(2); // Dashed line
    lineX_upper->Draw("Same");

    TLine* lineY_lower = new TLine(hist->GetXaxis()->GetXmin(), specificIsoET - windowSize,
                                   hist->GetXaxis()->GetXmax(), specificIsoET - windowSize);
    lineY_lower->SetLineColor(kBlue);
    lineY_lower->SetLineStyle(2); // Dashed line
    lineY_lower->Draw("Same");

    TLine* lineY_upper = new TLine(hist->GetXaxis()->GetXmin(), specificIsoET + windowSize,
                                   hist->GetXaxis()->GetXmax(), specificIsoET + windowSize);
    lineY_upper->SetLineColor(kBlue);
    lineY_upper->SetLineStyle(2); // Dashed line
    lineY_upper->Draw("Same");

    // Extract a subrange around the specific Cluster E_T and Isolation E_T_iso
    int binX_low = hist->GetXaxis()->FindBin(specificClusterET - windowSize);
    int binX_high = hist->GetXaxis()->FindBin(specificClusterET + windowSize);
    int binY_low = hist->GetYaxis()->FindBin(specificIsoET - windowSize);
    int binY_high = hist->GetYaxis()->FindBin(specificIsoET + windowSize);

    // Ensure bin ranges are within histogram limits
    binX_low = std::max(1, binX_low);
    binX_high = std::min(hist->GetXaxis()->GetNbins(), binX_high);
    binY_low = std::max(1, binY_low);
    binY_high = std::min(hist->GetYaxis()->GetNbins(), binY_high);

    // Create a projection for Cluster E_T
    TH1D* projX = hist->ProjectionX("_px", binY_low, binY_high);
    projX->SetTitle(Form("Projection on Cluster E_{T} around %.1f GeV", specificClusterET));
    projX->SetLineColor(kRed);
    projX->SetLineWidth(2);

    // Create a projection for Isolation E_T_iso
    TH1D* projY = hist->ProjectionY("_py", binX_low, binX_high);
    projY->SetTitle(Form("Projection on Isolation E_{T}^{iso} around %.1f GeV", specificIsoET));
    projY->SetLineColor(kBlue);
    projY->SetLineWidth(2);

    // Create a separate canvas for projections
    TCanvas* projCanvas = new TCanvas("projCanvas", "Projections for Specific Event", 1600, 600); // Increased width for better visibility
    projCanvas->Divide(2,1);

    // Draw Cluster E_T projection
    projCanvas->cd(1);
    projX->Draw();

    // Draw Isolation E_T_iso projection
    projCanvas->cd(2);
    projY->Draw();

    // Save the projections
    std::ostringstream projOutputPathStream;
    projOutputPathStream << saveDir << "/" << hist->GetName()
                         << "_Projections.png";
    std::string projOutputFilePath = projOutputPathStream.str();
    projCanvas->SaveAs(projOutputFilePath.c_str());
    std::cout << "\033[34m[INFO]\033[0m Saved projections to: " << projOutputFilePath << std::endl;

    // Save the original 2D histogram
    canvas.SaveAs(outputFilePath.c_str());
    std::cout << "Saved 2D histogram to: " << outputFilePath << std::endl;

    // Clean up dynamically allocated objects
    delete lineX_lower;
    delete lineX_upper;
    delete lineY_lower;
    delete lineY_upper;
    delete projX;
    delete projY;
    delete projCanvas;
    // Note: Do NOT delete 'canvas' as it is stack-allocated
}



void ProcessIsolationEnergyHistograms(
    const std::string& rootFilePath,
    const std::string& triggerCombinationDir,
    const std::map<std::string, std::string>& triggerCombinationNameMap) {

    std::cout << "Opening ROOT file: " << rootFilePath << std::endl;
    TFile* inputFile = TFile::Open(rootFilePath.c_str(), "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open file " << rootFilePath << std::endl;
        return;
    }

    std::string baseOutputDir = "/Users/patsfan753/Desktop/DirectPhotonAna/Plots";
    std::string triggerGroupName = Utils::getTriggerCombinationName(triggerCombinationDir, triggerCombinationNameMap);
    std::string isolationDir = baseOutputDir + "/" + triggerCombinationDir + "/isolationEnergies";

    if (gSystem->AccessPathName(isolationDir.c_str())) {
        if (gSystem->mkdir(isolationDir.c_str(), true) != 0) {
            std::cerr << "Error: Could not create directory " << isolationDir << std::endl;
            inputFile->Close();
            delete inputFile;
            return;
        }
    }

    std::regex hist1DPattern("h1_isoEt_(subtracted|unsubtracted)_pT_([0-9]+)to([0-9]+)_(.+)");
    std::regex hist1DPatternNoPtBins("h1_isoEt_(subtracted|unsubtracted)_(.+)");
    std::regex hist2DPattern("h2_cluster_iso_Et_(subtracted|unsubtracted)_pT_([0-9]+)to([0-9]+)_(.+)");
    std::regex hist2DPatternNoPtBins("h2_cluster_iso_Et_(subtracted|unsubtracted)_(.+)");
    std::smatch match;
    
    
    // Get list of triggers (directories) in the ROOT file
    std::vector<std::string> triggers;
    TList* dirList = inputFile->GetListOfKeys();
    TIter nextDir(dirList);
    TKey* dirKey;
    while ((dirKey = (TKey*)nextDir())) {
        TObject* dirObj = dirKey->ReadObj();
        if (dirObj->InheritsFrom(TDirectory::Class())) {
            std::string dirName = dirObj->GetName();
            triggers.push_back(dirName);
        }
        delete dirObj; // Don't forget to delete the object
    }

    for (const auto& trigger : triggers) {
        // Access the trigger directory within the ROOT file
        TDirectory* triggerDir = inputFile->GetDirectory(trigger.c_str());
        if (!triggerDir) {
            std::cerr << "Warning: Trigger directory '" << trigger << "' not found in the input file. Skipping." << std::endl;
            continue;
        }
        TIter nextKey(triggerDir->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*)nextKey())) {
            TObject* obj = key->ReadObj();
            if (!obj) continue;

            std::string histName = obj->GetName();

            // Check for 2D histogram first
            if (obj->InheritsFrom(TH2::Class())) {
                TH2* hist2D = dynamic_cast<TH2*>(obj);
                if (hist2D) {
                    if (std::regex_match(histName, match, hist2DPattern)) {
                        std::cout << "Processing 2D Histogram with pT bins: " << histName << std::endl;
                        Process2DHistogram(hist2D, match, isolationDir, triggerGroupName);
                    } else if (std::regex_match(histName, match, hist2DPatternNoPtBins)) {
                        std::cout << "Processing 2D Histogram without pT bins: " << histName << std::endl;
                        Process2DHistogram(hist2D, match, isolationDir, triggerGroupName, /*hasPtBins=*/false);
                    }
                }
            }
            // Check for 1D histogram
            else if (obj->InheritsFrom(TH1::Class())) {
                TH1* hist1D = dynamic_cast<TH1*>(obj);
                if (hist1D) {
                    if (std::regex_match(histName, match, hist1DPattern)) {
                        std::cout << "Processing 1D Histogram with pT bins: " << histName << std::endl;
                        Process1DHistogram(hist1D, match, isolationDir, triggerGroupName);
                    } else if (std::regex_match(histName, match, hist1DPatternNoPtBins)) {
                        std::cout << "Processing 1D Histogram without pT bins: " << histName << std::endl;
                        Process1DHistogram(hist1D, match, isolationDir, triggerGroupName, /*hasPtBins=*/false);
                    }
                }
            }

            delete obj; // Clean up
        }
    }
    inputFile->Close();
    delete inputFile;
}


void ProcessAllIsolationEnergies(
    const std::string& outputDirectory,
    const std::vector<std::string>& combinedRootFiles,
    const std::map<std::string, std::string>& triggerCombinationNameMap) {

    for (const auto& rootFileName : combinedRootFiles) {
        std::string rootFilePath = outputDirectory + "/" + rootFileName;

        // Extract the trigger combination name from the ROOT file name
        std::string combinationName = rootFileName.substr(0, rootFileName.find("_Combined.root"));

        std::cout << "Processing isolation energies for: " << combinationName << std::endl;

        // Call the function to process isolation energy histograms
        ProcessIsolationEnergyHistograms(rootFilePath, combinationName, triggerCombinationNameMap);
    }
}
/*
 To process only a specific trigger combination, you can call AnalyzeTriggerGroupings with the desired combination name:
 */
// AnalyzeTriggerGroupings("MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1_afterTriggerFirmwareUpdate");
void AnalyzeTriggerGroupings(std::string specificCombinationName = "") {
    gROOT->LoadMacro("sPhenixStyle.C");
    SetsPhenixStyle();

    std::string csvFilePath = "/Users/patsfan753/Desktop/DirectPhotonAna/triggerAnalysisCombined.csv";
    std::string outputDirectory = "/Users/patsfan753/Desktop/DirectPhotonAna/output";
    
    
    bool debugMode = false;
    std::map<int, std::map<std::string, std::string>> overrideTriggerStatus; // Empty map
    overrideTriggerStatus[46697]["Photon_3_GeV_plus_MBD_NS_geq_1"] = "OFF";

    // Get the map of trigger combinations to run numbers
    std::map<std::set<std::string>, DataStructures::RunInfo> combinationToRuns = AnalyzeWhatTriggerGroupsAvailable(csvFilePath, debugMode, overrideTriggerStatus);

    // Now, loop through the map and print out the groupings and structure
    std::cout << "\nSummary of Trigger Combinations:\n";

    // Call the new function to sort and print the combinations
    PrintSortedCombinations(combinationToRuns);

    // Map to store valid runs for each combination
    std::map<std::string, std::vector<int>> combinationToValidRuns;

    // Check if all combined ROOT files and valid runs files already exist
    bool allOutputFilesExist = true;
    for (const auto& kv : combinationToRuns) {
        const std::set<std::string>& triggers = kv.first;
        const DataStructures::RunInfo& runInfo = kv.second;

        // Create a string to represent the combination for naming
        std::string baseCombinationName;
        for (const auto& trigger : triggers) {
            baseCombinationName += trigger + "_";
        }
        // Remove the trailing underscore
        if (!baseCombinationName.empty()) {
            baseCombinationName.pop_back();
        }

        // Check runs before firmware update
        if (!runInfo.runsBeforeFirmwareUpdate.empty()) {
            std::string combinationName = baseCombinationName;
            if (!runInfo.runsAfterFirmwareUpdate.empty()) {
                combinationName += "_beforeTriggerFirmwareUpdate";
            }
            std::string finalRootFilePath = outputDirectory + "/" + combinationName + "_Combined.root";
            std::string validRunsFilePath = outputDirectory + "/" + combinationName + "_ValidRuns.txt";

            bool rootFileExists = !gSystem->AccessPathName(finalRootFilePath.c_str());
            bool validRunsFileExists = !gSystem->AccessPathName(validRunsFilePath.c_str());

            if (!rootFileExists || !validRunsFileExists) {
                allOutputFilesExist = false;
                break;
            } else {
                // Read valid runs from the text file
                std::vector<int> validRuns;
                std::ifstream validRunsFile(validRunsFilePath);
                if (validRunsFile.is_open()) {
                    int runNumber;
                    while (validRunsFile >> runNumber) {
                        validRuns.push_back(runNumber);
                    }
                    validRunsFile.close();
                    combinationToValidRuns[combinationName] = validRuns;
                    std::cout << "Combination: " << combinationName << ", Number of valid runs: " << validRuns.size() << std::endl;
                } else {
                    std::cerr << "Failed to open valid runs file: " << validRunsFilePath << std::endl;
                    allOutputFilesExist = false;
                    break;
                }
            }
        }

        // Check runs after firmware update
        if (!runInfo.runsAfterFirmwareUpdate.empty()) {
            std::string combinationName = baseCombinationName;
            if (!runInfo.runsBeforeFirmwareUpdate.empty()) {
                combinationName += "_afterTriggerFirmwareUpdate";
            }
            std::string finalRootFilePath = outputDirectory + "/" + combinationName + "_Combined.root";
            std::string validRunsFilePath = outputDirectory + "/" + combinationName + "_ValidRuns.txt";

            bool rootFileExists = !gSystem->AccessPathName(finalRootFilePath.c_str());
            bool validRunsFileExists = !gSystem->AccessPathName(validRunsFilePath.c_str());

            if (!rootFileExists || !validRunsFileExists) {
                allOutputFilesExist = false;
                break;
            } else {
                // Read valid runs from the text file
                std::vector<int> validRuns;
                std::ifstream validRunsFile(validRunsFilePath);
                if (validRunsFile.is_open()) {
                    int runNumber;
                    while (validRunsFile >> runNumber) {
                        validRuns.push_back(runNumber);
                    }
                    validRunsFile.close();
                    combinationToValidRuns[combinationName] = validRuns;
                    std::cout << "Combination: " << combinationName << ", Number of valid runs: " << validRuns.size() << std::endl;
                } else {
                    std::cerr << "Failed to open valid runs file: " << validRunsFilePath << std::endl;
                    allOutputFilesExist = false;
                    break;
                }
            }
        }
    }

    if (!allOutputFilesExist) {
        // Some output files are missing, proceed to process and merge root files
        ProcessAndMergeRootFiles(combinationToRuns, outputDirectory, combinationToValidRuns);
    } else {
        std::cout << "All output ROOT files and valid runs files already exist. Skipping ProcessAndMergeRootFiles." << std::endl;
    }

    // Get list of combined ROOT files
    std::vector<std::string> combinedRootFiles;
    void* dirp = gSystem->OpenDirectory(outputDirectory.c_str());
    const char* entry;
    while ((entry = gSystem->GetDirEntry(dirp))) {
        std::string fileName = entry;
        if (Utils::EndsWith(fileName, "_Combined.root")) {
            // If specificCombinationName is set, only add matching files
            if (specificCombinationName.empty() || fileName.find(specificCombinationName + "_Combined.root") != std::string::npos) {
                combinedRootFiles.push_back(fileName);
            }
        }
    }
    gSystem->FreeDirectory(dirp);

    // Now plot the combined histograms
    PlotCombinedHistograms(outputDirectory, combinedRootFiles, combinationToValidRuns, "sigmoid", false);

    ProcessAllIsolationEnergies(outputDirectory, combinedRootFiles, TriggerCombinationNames::triggerCombinationNameMap);
}
