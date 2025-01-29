#include "AnalyzeTriggerGroupings.h"
#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <utility>
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

bool enableFits = false; // Set to true if you want to enable the fits


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
>, DataStructures::IsolationData> isolatedPhotonMap;

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


std::map<std::set<std::string>, DataStructures::RunInfo>
AnalyzeWhatTriggerGroupsAvailable(
    const std::string& csvFilePath,
    bool debugMode,
    const std::map<int, std::map<std::string, std::string>>& overrideTriggerStatus)
{
    using namespace std;

    // 1) Pull from TriggerConfig
    const vector<string>& allTriggers    = TriggerConfig::allTriggers;
    const vector<string>& photonTriggers = TriggerConfig::photonTriggers;

    // Jet triggers (we skip Jet_6 if desired)
    vector<string> jetTriggers = {
        "Jet_8_GeV_plus_MBD_NS_geq_1",
        "Jet_10_GeV_plus_MBD_NS_geq_1",
        "Jet_12_GeV_plus_MBD_NS_geq_1"
    };

    // 2) Data structures
    map<string, int> triggerToIndex;
    map<int, set<string>> runToActiveTriggers;

    int totalRunsProcessed = 0;

    // For debug printing
    map<set<string>, vector<int>> initialCombinationToRuns;
    vector<pair<set<string>, DataStructures::RunInfo>> finalCombinations;

    // 3) Open CSV + parse header
    ifstream file(csvFilePath);
    if (!file.is_open()) {
        cerr << "Failed to open file: " << csvFilePath << endl;
        return {};
    }

    string line;
    if (!getline(file, line)) {
        cerr << "Failed to read header from CSV file." << endl;
        return {};
    }

    vector<string> headers;
    {
        istringstream headerStream(line);
        string hdr;
        int colIdx = 0;
        while (getline(headerStream, hdr, ',')) {
            // Trim
            hdr.erase(0, hdr.find_first_not_of(" \t\r\n"));
            hdr.erase(hdr.find_last_not_of(" \t\r\n") + 1);

            headers.push_back(hdr);

            // If "runNumber" or recognized trigger, store col index
            if (hdr == "runNumber" ||
                find(allTriggers.begin(), allTriggers.end(), hdr) != allTriggers.end())
            {
                triggerToIndex[hdr] = colIdx;
            }
            colIdx++;
        }
    }
    if (triggerToIndex.find("runNumber") == triggerToIndex.end()) {
        cerr << "runNumber column not found in CSV header." << endl;
        return {};
    }

    // 4) Read CSV lines => fill runToActiveTriggers
    while (getline(file, line)) {
        vector<string> cells;
        {
            istringstream rowStream(line);
            string c;
            while (getline(rowStream, c, ',')) {
                c.erase(0, c.find_first_not_of(" \t\r\n"));
                c.erase(c.find_last_not_of(" \t\r\n") + 1);
                cells.push_back(c);
            }
        }
        if (cells.size() != headers.size()) {
            cerr << "Mismatch (#cells vs. #headers) in line: " << line << endl;
            continue;
        }

        int runNumber = stoi(cells[triggerToIndex["runNumber"]]);

        set<string> activeTriggers;
        for (const auto& trig : allTriggers) {
            auto itIdx = triggerToIndex.find(trig);
            if (itIdx == triggerToIndex.end()) {
                // e.g. triggers not in CSV or ones you skip
                continue;
            }
            int idx = itIdx->second;
            string status = cells[idx];

            // Trim
            status.erase(0, status.find_first_not_of(" \t\r\n"));
            status.erase(status.find_last_not_of(" \t\r\n") + 1);

            // Overrides
            auto runOver = overrideTriggerStatus.find(runNumber);
            if (runOver != overrideTriggerStatus.end()) {
                auto trigOver = runOver->second.find(trig);
                if (trigOver != runOver->second.end()) {
                    status = trigOver->second;
                }
            }
            if (status == "ON") {
                activeTriggers.insert(trig);
            }
        }
        // Use std::move to quiet the unqualified call warnings
        runToActiveTriggers[runNumber] = std::move(activeTriggers);
    }
    file.close();
    totalRunsProcessed = static_cast<int>(runToActiveTriggers.size());

    // 5) Build photon combos & jet combos, then unify
    // 5A) photon combos
    vector<set<string>> photonCombinations;
    {
        int nPhot = static_cast<int>(photonTriggers.size());
        int totalPhotonSubsets = (1 << nPhot);
        for (int mask = 0; mask < totalPhotonSubsets; ++mask) {
            set<string> combo;
            combo.insert("MBD_NandS_geq_1");
            for (int i = 0; i < nPhot; i++) {
                if (mask & (1 << i)) {
                    combo.insert(photonTriggers[i]);
                }
            }
            // Again, use std::move to avoid the warning
            photonCombinations.push_back(std::move(combo));
        }
    }

    // 5B) jet combos
    vector<set<string>> jetCombinations;
    {
        int nJet = static_cast<int>(jetTriggers.size());
        int totalJetSubsets = (1 << nJet);
        for (int mask = 0; mask < totalJetSubsets; ++mask) {
            set<string> combo;
            combo.insert("MBD_NandS_geq_1");
            for (int i = 0; i < nJet; i++) {
                if (mask & (1 << i)) {
                    combo.insert(jetTriggers[i]);
                }
            }
            jetCombinations.push_back(std::move(combo));
        }
    }

    // unify them
    vector<set<string>> triggerCombinations;
    triggerCombinations.reserve(photonCombinations.size() + jetCombinations.size());

    // push photon combos
    for (auto& pc : photonCombinations) {
        triggerCombinations.push_back(std::move(pc));
    }
    // push jet combos
    for (auto& jc : jetCombinations) {
        triggerCombinations.push_back(std::move(jc));
    }

    // 6) For each run => which combos are satisfied
    map<set<string>, vector<int>> tempCombinationToRuns;
    for (auto& kv : runToActiveTriggers) {
        int runNumber = kv.first;
        const set<string>& active = kv.second;

        // must have MBD
        if (active.find("MBD_NandS_geq_1") == active.end()) {
            continue;
        }
        for (auto& combo : triggerCombinations) {
            bool satisfies = true;
            for (auto& trig : combo) {
                if (active.find(trig) == active.end()) {
                    satisfies = false;
                    break;
                }
            }
            if (satisfies) {
                tempCombinationToRuns[combo].push_back(runNumber);
            }
        }
    }
    initialCombinationToRuns = tempCombinationToRuns;

    // 7) Split runs by firmware
    struct TempRunInfo {
        vector<int> runsBeforeFirmwareUpdate;
        vector<int> runsAfterFirmwareUpdate;
    };
    map<set<string>, TempRunInfo> tempComboRunInfo;
    set<set<string>> combosWithRun47289;

    for (auto& kv2 : tempCombinationToRuns) {
        const auto& combo = kv2.first;
        const auto& runs  = kv2.second;

        TempRunInfo runInfo;
        bool has47289 = false;
        for (int rn : runs) {
            if (rn == 47289) {
                has47289 = true;
            }
            if (rn < 47289) {
                runInfo.runsBeforeFirmwareUpdate.push_back(rn);
            } else {
                runInfo.runsAfterFirmwareUpdate.push_back(rn);
            }
        }
        if (has47289) {
            combosWithRun47289.insert(combo);
        }
        tempComboRunInfo[combo] = std::move(runInfo);
    }

    // 8) Group combos by identical run-lists => pick largest combos
    map<vector<int>, vector<set<string>>> runListToCombinationsBeforeFW;
    map<vector<int>, vector<set<string>>> runListToCombinationsAfterFW;

    for (auto& kv3 : tempComboRunInfo) {
        const auto& combo   = kv3.first;
        const auto& runInfo = kv3.second;

        if (!runInfo.runsBeforeFirmwareUpdate.empty()) {
            vector<int> tmp = runInfo.runsBeforeFirmwareUpdate;
            sort(tmp.begin(), tmp.end());
            runListToCombinationsBeforeFW[tmp].push_back(combo);
        }
        if (!runInfo.runsAfterFirmwareUpdate.empty()) {
            vector<int> tmp = runInfo.runsAfterFirmwareUpdate;
            sort(tmp.begin(), tmp.end());
            runListToCombinationsAfterFW[tmp].push_back(combo);
        }
    }

    map<set<string>, vector<int>> filteredCombinationToRunsBeforeFirmwareUpdate;
    {
        for (auto& kvB : runListToCombinationsBeforeFW) {
            const auto& runList = kvB.first;
            const auto& combos  = kvB.second;

            size_t maxSize = 0;
            for (auto& c : combos) {
                if (c.size() > maxSize) {
                    maxSize = c.size();
                }
            }
            for (auto& c : combos) {
                if (c.size() == maxSize) {
                    filteredCombinationToRunsBeforeFirmwareUpdate[c] = runList;
                }
            }
        }
    }

    map<set<string>, vector<int>> filteredCombinationToRunsAfterFirmwareUpdate;
    {
        for (auto& kvA : runListToCombinationsAfterFW) {
            const auto& runList = kvA.first;
            const auto& combos  = kvA.second;

            size_t maxSize = 0;
            for (auto& c : combos) {
                if (c.size() > maxSize) {
                    maxSize = c.size();
                }
            }
            for (auto& c : combos) {
                if (c.size() == maxSize) {
                    filteredCombinationToRunsAfterFirmwareUpdate[c] = runList;
                }
            }
        }
    }

    // 9) Build final map => combinationToRuns
    map<set<string>, DataStructures::RunInfo> combinationToRuns;

    // fill "before FW"
    for (auto& kvBB : filteredCombinationToRunsBeforeFirmwareUpdate) {
        const auto& combo  = kvBB.first;
        const auto& runVec = kvBB.second;

        DataStructures::RunInfo& runInfo = combinationToRuns[combo];
        runInfo.runsBeforeFirmwareUpdate = runVec;
    }
    // fill "after FW"
    for (auto& kvAA : filteredCombinationToRunsAfterFirmwareUpdate) {
        const auto& combo  = kvAA.first;
        const auto& runVec = kvAA.second;

        DataStructures::RunInfo& runInfo = combinationToRuns[combo];
        runInfo.runsAfterFirmwareUpdate = runVec;
    }

    // for debug
    for (auto& kvC : combinationToRuns) {
        finalCombinations.emplace_back(kvC.first, kvC.second);
    }

    // 10) debugMode printing
    if (debugMode) {
        cout << BOLD << BLUE << "\n===== Processing Summary =====\n" << RESET;
        cout << BOLD << "Total runs processed: " << totalRunsProcessed << RESET << "\n";

        // 1) initial combos
        cout << BOLD << "\nInitial Active Trigger Combinations (before splitting due to firmware update):\n" << RESET;
        cout << BOLD << left << setw(60) << "Combination" << right << setw(20) << "Number of Runs" << RESET << "\n";
        cout << string(80, '=') << "\n";
        for (auto& kv0 : initialCombinationToRuns) {
            const auto& combo = kv0.first;
            const auto& runs  = kv0.second;

            string comboStr;
            for (auto& trig : combo) {
                comboStr += trig + " ";
            }
            cout << left << setw(60) << comboStr << right << setw(20) << runs.size() << "\n";
        }

        // 2) combos that had run 47289
        cout << BOLD << "\nCombinations that included run 47289 and were split due to firmware update:\n" << RESET;
        if (combosWithRun47289.empty()) {
            cout << "  None\n";
        } else {
            cout << BOLD << left << setw(60) << "Combination" << RESET << "\n";
            cout << string(60, '=') << "\n";
            for (auto& combination : combosWithRun47289) {
                string comboStr;
                for (auto& trig : combination) {
                    comboStr += trig + " ";
                }
                cout << left << setw(60) << comboStr << "\n";
            }
        }

        // 3) groups w/ same run numbers (before FW)
        cout << BOLD << "\nGroups with identical run numbers before firmware update:\n" << RESET;
        {
            int groupIndex = 1;
            for (auto& kvB : runListToCombinationsBeforeFW) {
                const auto& runList = kvB.first;
                const auto& combos  = kvB.second;

                cout << YELLOW << BOLD << "\nGroup " << groupIndex++ << RESET << "\n";
                cout << "Run List (size " << runList.size() << "):\n";
                for (size_t i=0; i<runList.size(); ++i) {
                    cout << setw(8) << runList[i];
                    if((i+1)%10==0 || i==runList.size()-1) {
                        cout << "\n";
                    }
                }
                cout << "  Combinations:\n";
                for (auto& c : combos) {
                    string comboStr;
                    for (auto& t : c) {
                        comboStr += t + " ";
                    }
                    cout << "    " << comboStr << "\n";
                }
            }
        }

        // 4) groups w/ same run numbers (after FW)
        cout << BOLD << "\nGroups with identical run numbers after firmware update:\n" << RESET;
        {
            int groupIndex = 1;
            for (auto& kvA : runListToCombinationsAfterFW) {
                const auto& runList = kvA.first;
                const auto& combos  = kvA.second;

                cout << YELLOW << BOLD << "\nGroup " << groupIndex++ << RESET << "\n";
                cout << "Run List (size " << runList.size() << "):\n";
                for (size_t i=0; i<runList.size(); ++i) {
                    cout << setw(8) << runList[i];
                    if((i+1)%10==0 || i==runList.size()-1) {
                        cout << "\n";
                    }
                }
                cout << "  Combinations:\n";
                for (auto& c : combos) {
                    string comboStr;
                    for (auto& t : c) {
                        comboStr += t + " ";
                    }
                    cout << "    " << comboStr << "\n";
                }
            }
        }

        // 5) final combos
        cout << BOLD << "\nFinal Active Trigger Combinations (after splitting and filtering):\n" << RESET;
        cout << BOLD << left << setw(60) << "Combination"
             << right << setw(25) << "Runs Before Firmware Update"
             << setw(25) << "Runs After Firmware Update" << RESET << "\n";
        cout << string(110, '=') << "\n";

        for (auto& fc : finalCombinations) {
            const auto& combo   = fc.first;
            const auto& runInfo = fc.second;

            string comboStr;
            for (auto& trig : combo) {
                comboStr += trig + " ";
            }
            cout << left << setw(60) << comboStr;
            cout << right << setw(25) << runInfo.runsBeforeFirmwareUpdate.size();
            cout << setw(25) << runInfo.runsAfterFirmwareUpdate.size() << "\n";

            if (!runInfo.runsBeforeFirmwareUpdate.empty()) {
                cout << "  Runs before firmware update ("
                     << runInfo.runsBeforeFirmwareUpdate.size() << " runs):\n";
                for (size_t i=0; i<runInfo.runsBeforeFirmwareUpdate.size(); i++) {
                    cout << setw(8) << runInfo.runsBeforeFirmwareUpdate[i];
                    if ((i+1)%10==0 || i==runInfo.runsBeforeFirmwareUpdate.size()-1) {
                        cout << "\n";
                    }
                }
            }
            if (!runInfo.runsAfterFirmwareUpdate.empty()) {
                cout << "  Runs after firmware update ("
                     << runInfo.runsAfterFirmwareUpdate.size() << " runs):\n";
                for (size_t i=0; i<runInfo.runsAfterFirmwareUpdate.size(); i++) {
                    cout << setw(8) << runInfo.runsAfterFirmwareUpdate[i];
                    if ((i+1)%10==0 || i==runInfo.runsAfterFirmwareUpdate.size()-1) {
                        cout << "\n";
                    }
                }
            }
        }

        cout << BOLD << BLUE << "\n===== End of Processing Summary =====\n" << RESET;
        // Possibly exit(0);
    }

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
    std::map<std::string, std::vector<int>>& combinationToValidRuns)
{
    // Define the final output ROOT file path
    std::string finalRootFilePath = outputDirectory + "/" + combinationName + "_Combined.root";
    // Define the text file path to store valid runs
    std::string validRunsFilePath = outputDirectory + "/" + combinationName + "_ValidRuns.txt";

    // New file to track runs that have zero data in all histograms
    std::string zeroDataFilePath = outputDirectory + "/RunsWithNoData.txt";
    std::ofstream zeroDataFile(zeroDataFilePath, std::ios::app); // append mode

    bool rootFileExists = !gSystem->AccessPathName(finalRootFilePath.c_str());
    bool validRunsFileExists = !gSystem->AccessPathName(validRunsFilePath.c_str());

    // If both files already exist, skip
    if (rootFileExists && validRunsFileExists) {
        std::cout << "Final ROOT file and valid runs file already exist for combination: "
                  << combinationName << ". Skipping merge.\n";
        // Load valid runs from text file
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

    // Map: triggerName -> (histogramName -> unique_ptr<TH1>)
    std::map<std::string, std::map<std::string, std::unique_ptr<TH1>>> mergedHistograms;

    // Vector of run numbers that have at least some valid histograms
    std::vector<int> validRuns;

    // Loop over each run in this combination
    for (int runNumber : runs) {
        std::stringstream ss;
        ss << outputDirectory << "/" << runNumber << "_HistOutput.root";
        std::string runRootFilePath = ss.str();

        std::cout << "\nProcessing run: " << runNumber
                  << ", file: " << runRootFilePath << std::endl;

        if (gSystem->AccessPathName(runRootFilePath.c_str())) {
            std::cerr << "Run ROOT file does not exist: " << runRootFilePath
                      << ". Skipping run.\n";
            continue;
        }
        std::cout << "Opening run file...\n";
        TFile runFile(runRootFilePath.c_str(), "READ");
        if (runFile.IsZombie() || !runFile.IsOpen()) {
            std::cerr << "Failed to open run ROOT file: " << runRootFilePath
                      << ". Skipping run.\n";
            continue;
        }
        std::cout << "Run file opened successfully.\n";

        bool runHasValidHistogram = false;   // indicates we found at least one histogram
        bool runHasNonZeroEntries = false;   // indicates at least one histogram has >0 entries

        // Loop over each trigger in the combination
        for (const auto& trigger : triggers) {
            // Grab the trigger directory
            TDirectory* triggerDir = runFile.GetDirectory(trigger.c_str());
            if (!triggerDir) {
                std::cerr << "Trigger directory '" << trigger
                          << "' not found for run " << runNumber << ".\n";
                continue;
            }

            // Read keys in that directory
            TIter nextKey(triggerDir->GetListOfKeys());
            TKey* key;
            while ((key = (TKey*)nextKey())) {
                std::string className = key->GetClassName();
                TObject* obj = key->ReadObj();
                if (!obj) {
                    continue;
                }

                // *** CRITICAL CHECKS FOR SEGV AVOIDANCE ***
                // 1) If it's not a TH1 or TH2, skip
                // 2) If you only want to handle 1D vs 2D separately, do it here:
                //    e.g. skip TProfile
                if (!obj->InheritsFrom(TH1::Class())) {
                    std::cout << "[DEBUG] Skipping non-TH1 object: "
                              << key->GetName() << " of class " << className << "\n";
                    delete obj;
                    continue;
                }
                // Optionally skip TProfile:
                if (obj->InheritsFrom(TProfile::Class())) {
                    std::cout << "[DEBUG] Skipping TProfile: "
                              << key->GetName() << "\n";
                    delete obj;
                    continue;
                }

                // Now we can treat it as a TH1
                TH1* hist = dynamic_cast<TH1*>(obj);
                if (!hist) {
                    // This should not happen given the check, but just in case:
                    std::cerr << "[ERROR] dynamic_cast<TH1*> failed for "
                              << key->GetName() << "\n";
                    delete obj;
                    continue;
                }

                std::string histName = hist->GetName();
                // *** If you want to skip merging 1D with 2D, check dimensions ***
                // e.g. if you want only 1D:
                // if (hist->GetDimension() != 1) { skip or handle 2D separately }

                // Clone the histogram
                TH1* histClone = dynamic_cast<TH1*>(hist->Clone());
                if (!histClone) {
                    std::cerr << "[ERROR] Clone failed for histogram: "
                              << histName << "\n";
                    delete obj;
                    continue;
                }
                histClone->SetDirectory(nullptr);

                // Merging logic
                auto& histMap = mergedHistograms[trigger];
                auto it = histMap.find(histName);
                if (it != histMap.end()) {
                    // *** Binning check to avoid Add(...) mismatch
                    TH1* existing = it->second.get();
                    bool canMerge = (existing->GetDimension() == histClone->GetDimension());

                    // (Optional) further check: same number of bins, same edges, etc.
                    if (canMerge) {
                        // e.g. check if # bins match
                        if (existing->GetNbinsX() != histClone->GetNbinsX()) {
                            std::cerr << "[WARNING] Bin mismatch merging '"
                                      << histName << "' from run " << runNumber
                                      << " => skipping.\n";
                            delete histClone;
                            delete obj;
                            continue;
                        }
                        // If 2D, also check Y dimension, etc.
                    }
                    else {
                        std::cerr << "[WARNING] Dimension mismatch merging '"
                                  << histName << "' => skipping.\n";
                        delete histClone;
                        delete obj;
                        continue;
                    }

                    // Safe to merge
                    existing->Add(histClone);
                    delete histClone;
                    std::cout << "Added histogram '" << histName
                              << "' from run " << runNumber
                              << " to existing histogram in trigger '"
                              << trigger << "'.\n";
                } else {
                    // First time seeing this histogram name for this trigger
                    histMap[histName] = std::unique_ptr<TH1>(histClone);
                    std::cout << "Added histogram '" << histName
                              << "' from run " << runNumber
                              << " to merged histograms in trigger '"
                              << trigger << "'.\n";
                }

                // Mark that we had at least one valid histogram
                runHasValidHistogram = true;
                // Check if this histogram has >0 entries
                if (hist->GetEntries() > 0) {
                    runHasNonZeroEntries = true;
                }

                delete obj; // Done with original
            }
        } // end trigger loop

        // If at least one valid histogram was found => record this run
        if (runHasValidHistogram) {
            validRuns.push_back(runNumber);
        }

        // If no histogram in this run had non-zero entries, log it
        if (!runHasNonZeroEntries) {
            zeroDataFile << runNumber << "\n";
            std::cout << "[INFO] Run " << runNumber
                      << " had NO non-zero hist entries => appended to: "
                      << zeroDataFilePath << "\n";
        }

        runFile.Close();
    } // end runs loop

    // If no runs found
    if (validRuns.empty()) {
        std::cout << "[INFO] No valid runs found for combination: "
                  << combinationName << ". Skipping.\n";
        return;
    }

    // If no histograms aggregated
    if (mergedHistograms.empty()) {
        std::cout << "[INFO] No histograms were merged for combination: "
                  << combinationName << ". Skipping output.\n";
        return;
    }

    // Write merged histograms to final ROOT file
    std::cout << "[INFO] Writing merged histograms => "
              << finalRootFilePath << "\n";
    TFile finalFile(finalRootFilePath.c_str(), "RECREATE");
    if (finalFile.IsZombie() || !finalFile.IsOpen()) {
        std::cerr << "[ERROR] Could not create " << finalRootFilePath << "\n";
        mergedHistograms.clear();
        return;
    }
    finalFile.cd();

    // For each trigger => create directory => write histograms
    for (const auto& [trigger, histMap] : mergedHistograms) {
        TDirectory* trigDir = finalFile.mkdir(trigger.c_str());
        if (!trigDir) {
            std::cerr << "[ERROR] Could not mkdir for trigger " << trigger << "\n";
            continue;
        }
        trigDir->cd();

        for (const auto& [histName, histPtr] : histMap) {
            if (!histPtr) {
                std::cerr << "[WARNING] Null histogram? Skipping " << histName << "\n";
                continue;
            }
            histPtr->Write();
            std::cout << "[INFO] Wrote '" << histName
                      << "' in directory '" << trigger << "'\n";
        }
    }

    finalFile.Write();
    finalFile.Close();
    mergedHistograms.clear();

    std::cout << "[INFO] Successfully created combined ROOT file: "
              << finalRootFilePath << "\n";

    // Update valid runs text file
    combinationToValidRuns[combinationName] = validRuns;
    std::ofstream validRunsFile(validRunsFilePath);
    if (validRunsFile.is_open()) {
        for (int rn : validRuns) {
            validRunsFile << rn << "\n";
        }
        validRunsFile.close();
        std::cout << "[INFO] Valid runs => " << validRunsFilePath << "\n";
    } else {
        std::cerr << "[ERROR] Failed to write " << validRunsFilePath << "\n";
    }
    // zeroDataFile closes automatically when function ends
}


/**
 * @brief Builds a sorted combination name starting with "MBD_NandS_geq_1"
 *        then listing other triggers in ascending threshold order.
 *
 *        Examples:
 *            Input set: {"MBD_NandS_geq_1", "Photon_3_GeV_plus_MBD_NS_geq_1", "Photon_5_GeV_plus_MBD_NS_geq_1"}
 *            Output:    "MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_5_GeV_plus_MBD_NS_geq_1"
 *
 *            Input set: {"MBD_NandS_geq_1", "Jet_8_GeV_plus_MBD_NS_geq_1", "Jet_10_GeV_plus_MBD_NS_geq_1"}
 *            Output:    "MBD_NandS_geq_1_Jet_8_GeV_plus_MBD_NS_geq_1_Jet_10_GeV_plus_MBD_NS_geq_1"
 *
 *        If the set is just {"MBD_NandS_geq_1"}, returns "MBD_NandS_geq_1".
 *
 */
static std::string buildSortedCombinationName(const std::set<std::string>& triggers)
{
    // If the set is just MBD alone => quick return
    if (triggers.size() == 1 && triggers.count("MBD_NandS_geq_1") == 1) {
        return "MBD_NandS_geq_1";
    }

    // We'll gather all triggers except "MBD_NandS_geq_1" in a vector
    // Then parse their threshold and sort them ascending.
    std::vector<std::string> otherTriggers;
    otherTriggers.reserve(triggers.size() - 1);

    for (const auto& trig : triggers) {
        if (trig == "MBD_NandS_geq_1") {
            continue; // skip
        }
        otherTriggers.push_back(trig);
    }

    // A small helper to parse the threshold integer from e.g. "Photon_4_GeV_plus_MBD_NS_geq_1"
    // or "Jet_8_GeV_plus_MBD_NS_geq_1".
    // You can adjust the regex as needed or do manual substring search.
    auto extractThreshold = [&](const std::string& triggerName) -> int {
        // Regex approach: match e.g. "Photon_(\d+)_GeV_plus" or "Jet_(\d+)_GeV_plus"
        static const std::regex re(R"((?:Photon|Jet)_(\d+)_GeV_plus)");
        std::smatch match;
        if (std::regex_search(triggerName, match, re)) {
            // capture group #1 = the digits
            return std::stoi(match[1].str());
        }
        return 0; // fallback if something goes wrong
    };

    // Sort otherTriggers by the integer threshold
    std::sort(otherTriggers.begin(), otherTriggers.end(),
              [&](const std::string& A, const std::string& B) {
                  int thrA = extractThreshold(A);
                  int thrB = extractThreshold(B);
                  return (thrA < thrB);
              }
    );

    // Now build final name, always start with "MBD_NandS_geq_1"
    // then underscore, then sorted triggers.
    std::ostringstream oss;
    oss << "MBD_NandS_geq_1";  // always first
    for (const auto& trig : otherTriggers) {
        oss << "_" << trig;
    }

    return oss.str();
}

/**
 * @brief Merges the histograms for each unique trigger combination into a
 *        single combined ROOT file, labeling them with a sorted combination name.
 *
 *        If a combination is exactly "MBD_NandS_geq_1", we combine all runs
 *        before/after firmware into one file. Otherwise, we split runs
 *        into "beforeTriggerFirmwareUpdate" and "afterTriggerFirmwareUpdate"
 *        as usual.
 *
 * @param combinationToRuns        The map of (trigger set) -> (RunInfo with runsBefore/After)
 * @param outputDirectory          Where to output the combined ROOT files
 * @param combinationToValidRuns   [Output] This map is filled with the final valid run list for each combination
 */
void ProcessAndMergeRootFiles(
    const std::map<std::set<std::string>, DataStructures::RunInfo>& combinationToRuns,
    const std::string& outputDirectory,
    std::map<std::string, std::vector<int>>& combinationToValidRuns)
{
    std::cout << "Starting ProcessAndMergeRootFiles" << std::endl;

    // Disable automatic addition of histograms to directories
    TH1::AddDirectory(false);

    // -------------------------------------------------------------------------
    //  Iterate over each trigger combination
    // -------------------------------------------------------------------------
    for (const auto& kv : combinationToRuns) {
        const std::set<std::string>& triggers = kv.first;
        const DataStructures::RunInfo& runInfo = kv.second;

        // ---------------------------------------------------------------------
        // Build the sorted combination name: always MBD first, then ascending thresholds
        // ---------------------------------------------------------------------
        std::string baseCombinationName = buildSortedCombinationName(triggers);

        // ---------------------------------------------------------------------
        // If combination is exactly {"MBD_NandS_geq_1"} => combine runs (before/after)
        // ---------------------------------------------------------------------
        if (triggers.size() == 1 && triggers.count("MBD_NandS_geq_1") == 1) {
            // Merge runs before & after FW into one
            std::vector<int> allRuns = runInfo.runsBeforeFirmwareUpdate;
            allRuns.insert(allRuns.end(),
                           runInfo.runsAfterFirmwareUpdate.begin(),
                           runInfo.runsAfterFirmwareUpdate.end());

            // This name is simply "MBD_NandS_geq_1" in that scenario
            std::string combinationName = baseCombinationName;
            
            // Merge ROOT files for these runs
            ProcessRunsForCombination(
                combinationName, allRuns, triggers,
                outputDirectory, combinationToValidRuns
            );
        }
        // ---------------------------------------------------------------------
        // Otherwise, handle runsBeforeFirmwareUpdate & runsAfterFirmwareUpdate separately
        // ---------------------------------------------------------------------
        else {
            // 1) Runs before firmware update
            if (!runInfo.runsBeforeFirmwareUpdate.empty()) {
                std::string combinationName = baseCombinationName;
                if (!runInfo.runsAfterFirmwareUpdate.empty()) {
                    // If we also have runsAfter, append suffix
                    combinationName += "_beforeTriggerFirmwareUpdate";
                }
                const std::vector<int>& runs = runInfo.runsBeforeFirmwareUpdate;

                ProcessRunsForCombination(
                    combinationName, runs, triggers,
                    outputDirectory, combinationToValidRuns
                );
            }

            // 2) Runs after firmware update
            if (!runInfo.runsAfterFirmwareUpdate.empty()) {
                std::string combinationName = baseCombinationName;
                if (!runInfo.runsBeforeFirmwareUpdate.empty()) {
                    combinationName += "_afterTriggerFirmwareUpdate";
                }
                const std::vector<int>& runs = runInfo.runsAfterFirmwareUpdate;

                ProcessRunsForCombination(
                    combinationName, runs, triggers,
                    outputDirectory, combinationToValidRuns
                );
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
    labelText.SetTextSize(0.025);
    labelText.SetTextColor(kRed);       // Set text color to red
    labelText.SetTextFont(62);          // Bold font for labels

    valueText.SetNDC();
    valueText.SetTextSize(0.024);
    valueText.SetTextColor(kBlack);
    valueText.SetTextFont(42);

    
    // Add the 'Active Trigger Group' label above the current trigger label
    labelText.DrawLatex(0.34, 0.9, "Active Trigger Group:");
    valueText.DrawLatex(0.51, 0.9, triggerGroupName.c_str());

    std::string readableTriggerName = data.cuts.triggerName;
    auto triggerNameIt = TriggerConfig::triggerNameMap.find(data.cuts.triggerName);
    if (triggerNameIt != TriggerConfig::triggerNameMap.end()) {
        readableTriggerName = triggerNameIt->second;
    }

    // First column: Trigger information
    labelText.DrawLatex(0.5, 0.85, "Trigger:");
    valueText.DrawLatex(0.62, 0.85, readableTriggerName.c_str());

    labelText.DrawLatex(0.5, 0.8, "Ecore:");
    valueText.DrawLatex(0.62, 0.8, Utils::formatToThreeSigFigs(data.cuts.clusECore).c_str());

    labelText.DrawLatex(0.5, 0.75, "#chi^{2}:");
    valueText.DrawLatex(0.62, 0.75, Utils::formatToThreeSigFigs(data.cuts.chi).c_str());

    labelText.DrawLatex(0.5, 0.7, "Asymmetry:");
    valueText.DrawLatex(0.62, 0.7, Utils::formatToThreeSigFigs(data.cuts.asymmetry).c_str());

    // If pT range is available, add it to the legend
    if (!ptRangeLabel.str().empty()) {
        labelText.DrawLatex(0.5, 0.65, "pT Range:");
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

void readHistogramDataForInvariantMasses(
    const std::string& csvFilePath,
    std::vector<DataStructures::HistogramData>& dataVector,
    std::set<std::string>& triggersInDataFile)
{
    // Helper: safely convert string to float
    auto safeToFloat = [&](const std::string& s, bool& ok) -> float {
        // Trim whitespace
        std::string tmp = s;
        tmp.erase(0, tmp.find_first_not_of(" \t\r\n"));
        tmp.erase(tmp.find_last_not_of(" \t\r\n") + 1);

        try {
            float val = std::stof(tmp);
            ok = true;
            return val;
        }
        catch (...) {
            ok = false;
            return 0.0f; // fallback
        }
    };

    // Helper: safely convert string to double
    auto safeToDouble = [&](const std::string& s, bool& ok) -> double {
        // Trim whitespace
        std::string tmp = s;
        tmp.erase(0, tmp.find_first_not_of(" \t\r\n"));
        tmp.erase(tmp.find_last_not_of(" \t\r\n") + 1);

        try {
            double val = std::stod(tmp);
            ok = true;
            return val;
        }
        catch (...) {
            ok = false;
            return 0.0; // fallback
        }
    };

    // Attempt to open file
    std::ifstream csvFile(csvFilePath);
    if (!csvFile.is_open()) {
        std::cerr << "[ERROR] Could not open CSV file for reading: " << csvFilePath << std::endl;
        return;
    }

    // First line is the CSV header, so read and discard (or check)
    std::string line;
    if (!std::getline(csvFile, line)) {
        std::cerr << "[ERROR] CSV file appears empty or missing header => " << csvFilePath << std::endl;
        return;
    }
    std::cout << "[INFO] Skipping header line: " << line << std::endl;

    int lineNum = 1; // we've read the header, so data lines start at lineNum=2

    while (std::getline(csvFile, line)) {
        ++lineNum;
        if (line.empty()) {
            // skip empty lines
            continue;
        }

        // Split by comma into tokens
        std::vector<std::string> tokens;
        {
            std::istringstream iss(line);
            std::string cell;
            while (std::getline(iss, cell, ',')) {
                tokens.push_back(cell);
            }
        }

        // We expect exactly 14 columns
        if (tokens.size() < 14) {
            std::cerr << "[WARNING] line " << lineNum
                      << ": Only " << tokens.size() << " columns found, need 14. Skipping row.\n";
            continue;
        }

        // Prepare a new structure
        DataStructures::HistogramData data;
        bool ok = false; // used to track numeric parse success

        // 1) Trigger Name
        data.cuts.triggerName = tokens[0];
        triggersInDataFile.insert(tokens[0]); // remember we saw this trigger

        // 2) Ecore
        float eCoreVal = safeToFloat(tokens[1], ok);
        if (!ok) {
            std::cerr << "[WARNING] line " << lineNum
                      << ": Ecore parse failure => '" << tokens[1] << "'\n";
            continue; // skip or set default
        }
        data.cuts.clusECore = eCoreVal;

        // 3) Chi2
        float chiVal = safeToFloat(tokens[2], ok);
        if (!ok) {
            std::cerr << "[WARNING] line " << lineNum
                      << ": Chi parse failure => '" << tokens[2] << "'\n";
            continue;
        }
        data.cuts.chi = chiVal;

        // 4) Asym
        float asymVal = safeToFloat(tokens[3], ok);
        if (!ok) {
            std::cerr << "[WARNING] line " << lineNum
                      << ": Asym parse failure => '" << tokens[3] << "'\n";
            continue;
        }
        data.cuts.asymmetry = asymVal;

        // 5) pTMin
        float pTMinVal = safeToFloat(tokens[4], ok);
        if (!ok) {
            std::cerr << "[WARNING] line " << lineNum
                      << ": pTMin parse failure => '" << tokens[4] << "'\n";
            continue;
        }
        data.cuts.pTMin = pTMinVal;

        // 6) pTMax
        float pTMaxVal = safeToFloat(tokens[5], ok);
        if (!ok) {
            std::cerr << "[WARNING] line " << lineNum
                      << ": pTMax parse failure => '" << tokens[5] << "'\n";
            continue;
        }
        data.cuts.pTMax = pTMaxVal;

        // 7) meanPi0
        double meanPi0Val = safeToDouble(tokens[6], ok);
        if (!ok) {
            std::cerr << "[WARNING] line " << lineNum
                      << ": meanPi0 parse failure => '" << tokens[6] << "'\n";
            continue;
        }
        data.meanPi0 = meanPi0Val;

        // 8) meanPi0Error
        double meanPi0ErrVal = safeToDouble(tokens[7], ok);
        if (!ok) {
            std::cerr << "[WARNING] line " << lineNum
                      << ": meanPi0Error parse failure => '" << tokens[7] << "'\n";
            continue;
        }
        data.meanPi0Error = meanPi0ErrVal;

        // 9) sigmaPi0
        double sigmaPi0Val = safeToDouble(tokens[8], ok);
        if (!ok) {
            std::cerr << "[WARNING] line " << lineNum
                      << ": sigmaPi0 parse failure => '" << tokens[8] << "'\n";
            continue;
        }
        data.sigmaPi0 = sigmaPi0Val;

        // 10) sigmaPi0Error
        double sigmaPi0ErrVal = safeToDouble(tokens[9], ok);
        if (!ok) {
            std::cerr << "[WARNING] line " << lineNum
                      << ": sigmaPi0Error parse failure => '" << tokens[9] << "'\n";
            continue;
        }
        data.sigmaPi0Error = sigmaPi0ErrVal;

        // 11) meanEta
        double meanEtaVal = safeToDouble(tokens[10], ok);
        if (!ok) {
            std::cerr << "[WARNING] line " << lineNum
                      << ": meanEta parse failure => '" << tokens[10] << "'\n";
            continue;
        }
        data.meanEta = meanEtaVal;

        // 12) meanEtaError
        double meanEtaErrVal = safeToDouble(tokens[11], ok);
        if (!ok) {
            std::cerr << "[WARNING] line " << lineNum
                      << ": meanEtaError parse failure => '" << tokens[11] << "'\n";
            continue;
        }
        data.meanEtaError = meanEtaErrVal;

        // 13) sigmaEta
        double sigmaEtaVal = safeToDouble(tokens[12], ok);
        if (!ok) {
            std::cerr << "[WARNING] line " << lineNum
                      << ": sigmaEta parse failure => '" << tokens[12] << "'\n";
            continue;
        }
        data.sigmaEta = sigmaEtaVal;

        // 14) sigmaEtaError
        double sigmaEtaErrVal = safeToDouble(tokens[13], ok);
        if (!ok) {
            std::cerr << "[WARNING] line " << lineNum
                      << ": sigmaEtaError parse failure => '" << tokens[13] << "'\n";
            continue;
        }
        data.sigmaEtaError = sigmaEtaErrVal;

        // If we made it here, the line had 14 valid columns
        dataVector.push_back(data);
    } // end while lines

    csvFile.close();

    std::cout << "[INFO] Finished reading CSV => " << csvFilePath << "\n"
              << "       Read " << dataVector.size() << " valid rows.\n";
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
                result.sigmaPi0Values.push_back(selectedData.sigmaPi0);
                result.sigmaPi0Errors.push_back(selectedData.sigmaPi0Error);
                result.resolutionPi0Values.push_back(selectedData.pi0FitResolution);
                result.resolutionPi0Errors.push_back(selectedData.pi0FitResolutionError);
                result.signalToBackgroundPi0Ratios.push_back(selectedData.signalToBackgroundPi0Ratio);
                result.signalToBackgroundPi0Errors.push_back(selectedData.signalToBackgroundPi0Error);
                result.pi0YieldValues.push_back(selectedData.signalPi0Yield);
                result.pi0YieldErrors.push_back(selectedData.signalPi0Error);

                result.triggersUsedPi0.push_back(triggerToUse);
            }

            if (validEta) {
                result.pTCentersEta.push_back(pTCenter);
                result.meanEtaValues.push_back(selectedData.meanEta);
                result.meanEtaErrors.push_back(selectedData.meanEtaError);
                result.sigmaEtaValues.push_back(selectedData.sigmaEta);
                result.sigmaEtaErrors.push_back(selectedData.sigmaEtaError);
                result.resolutionEtaValues.push_back(selectedData.etaFitResolution);
                result.resolutionEtaErrors.push_back(selectedData.etaFitResolutionError);
                result.signalToBackgroundEtaRatios.push_back(selectedData.signalToBackgroundEtaRatio);
                result.signalToBackgroundEtaErrors.push_back(selectedData.signalToBackgroundEtaError);
                result.etaYieldValues.push_back(selectedData.signalEtaYield);
                result.etaYieldErrors.push_back(selectedData.signalEtaError);

                result.triggersUsedEta.push_back(triggerToUse);
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

        for (const auto& tDataPair : triggerDataMap) {
            const std::string& triggerName = tDataPair.first;
            const DataStructures::HistogramData& data = tDataPair.second;

            bool validPi0 = data.meanPi0 > 0.0 && data.meanPi0 < 0.4;
            bool validEta = data.meanEta > 0.3 && data.meanEta < 1.0;

            if (validPi0) {
                result.triggerToPtCentersPi0[triggerName].push_back(pTCenter);
                result.triggerToMeanPi0Values[triggerName].push_back(data.meanPi0);
                result.triggerToMeanPi0Errors[triggerName].push_back(data.meanPi0Error);
                result.triggerToSigmaPi0Values[triggerName].push_back(data.sigmaPi0);
                result.triggerToSigmaPi0Errors[triggerName].push_back(data.sigmaPi0Error);
                result.triggerToResolutionPi0Values[triggerName].push_back(data.pi0FitResolution);
                result.triggerToResolutionPi0Errors[triggerName].push_back(data.pi0FitResolutionError);
                result.triggerToSignalToBackgroundPi0Ratios[triggerName].push_back(data.signalToBackgroundPi0Ratio);
                result.triggerToSignalToBackgroundPi0Errors[triggerName].push_back(data.signalToBackgroundPi0Error);
                result.triggerToPi0YieldValues[triggerName].push_back(data.signalPi0Yield);
                result.triggerToPi0YieldErrors[triggerName].push_back(data.signalPi0Error);
            }

            if (validEta) {
                result.triggerToPtCentersEta[triggerName].push_back(pTCenter);
                result.triggerToMeanEtaValues[triggerName].push_back(data.meanEta);
                result.triggerToMeanEtaErrors[triggerName].push_back(data.meanEtaError);
                result.triggerToSigmaEtaValues[triggerName].push_back(data.sigmaEta);
                result.triggerToSigmaEtaErrors[triggerName].push_back(data.sigmaEtaError);
                result.triggerToResolutionEtaValues[triggerName].push_back(data.etaFitResolution);
                result.triggerToResolutionEtaErrors[triggerName].push_back(data.etaFitResolutionError);
                result.triggerToSignalToBackgroundEtaRatios[triggerName].push_back(data.signalToBackgroundEtaRatio);
                result.triggerToSignalToBackgroundEtaErrors[triggerName].push_back(data.signalToBackgroundEtaError);
                result.triggerToEtaYieldValues[triggerName].push_back(data.signalEtaYield);
                result.triggerToEtaYieldErrors[triggerName].push_back(data.signalEtaError);
            }
        }
    }

    return result;
}

void generateOverlayInvariantMassPlot(
    const std::map<std::string, std::vector<double>>& triggerToPtCenters,
    const std::map<std::string, std::vector<double>>& triggerToYValues,
    const std::map<std::string, std::vector<double>>& triggerToYErrors,
    const std::set<std::string>& triggersInData,
    const std::vector<std::pair<double, double>>& pT_bins,
    const std::string& mesonName,  // e.g., "#pi^{0}" or "#eta"
    const std::string& yAxisTitle, // e.g., "Mean #pi^{0} Mass [GeV]" or "Yield"
    const double yBufferFraction,
    const double pTExclusionMax,
    const std::string& outputFilePath,
    const std::string& plotDirectory,
    const std::string& cutCombination,
    const DataStructures::CutCombinationData& cutData,
    const std::map<std::string, int>& triggerColorMap,
    const std::map<std::string, std::string>& triggerNameMap)
{
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
    if (!binEdges.empty()) {
        binEdges.push_back(pT_bins[binEdges.size()-1].second);
    } else {
        // If no bins are valid, we can't plot
        std::cerr << "No valid pT bins for plotting." << std::endl;
        return;
    }

    int nBins = (int)binEdges.size() - 1;
    double* binEdgesArray = &binEdges[0];

    // Create a dummy histogram to set up the axes
    TH1F* hFrame = new TH1F("hFrame", "", nBins, binEdgesArray);
    hFrame->SetStats(0);
    hFrame->GetXaxis()->SetTitle("Leading Cluster p_{T} [GeV]");
    hFrame->GetYaxis()->SetTitle(yAxisTitle.c_str());

    // Remove x-axis labels and ticks (we will draw custom labels)
    hFrame->GetXaxis()->SetLabelOffset(999);
    hFrame->GetXaxis()->SetTickLength(0);

    // Draw the frame
    hFrame->Draw();

    // Configure legend position and style
    TLegend* legend = nullptr;
    if (mesonName == "#pi^{0}") {
        // Legend for pi0
        legend = new TLegend(0.2, 0.2, 0.4, 0.4);
        legend->SetTextSize(0.03);
    } else {
        // Legend for eta or others
        legend = new TLegend(0.2, 0.2, 0.4, 0.4);
        legend->SetTextSize(0.03);
    }

    std::vector<TGraphErrors*> graphs;

    // Variables to adjust y-axis range
    double yMinData = std::numeric_limits<double>::max();
    double yMaxData = std::numeric_limits<double>::lowest();

    // Determine offsets for multiple triggers
    std::vector<std::string> triggerNames(triggersInData.begin(), triggersInData.end());
    int numTriggers = (int)triggerNames.size();

    double offsetValue;
    if (mesonName == "#eta") {
        offsetValue = 0.098;
    } else {
        offsetValue = 0.08;
    }

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

    // Loop over triggers and plot
    for (size_t triggerIndex = 0; triggerIndex < triggerNames.size(); ++triggerIndex) {
        const std::string& triggerName = triggerNames[triggerIndex];
        double xOffset = offsets[triggerIndex];

        // Retrieve data for this trigger
        auto it_pT = triggerToPtCenters.find(triggerName);
        auto it_yValues = triggerToYValues.find(triggerName);
        auto it_yErrors = triggerToYErrors.find(triggerName);

        if (it_pT == triggerToPtCenters.end() || it_yValues == triggerToYValues.end() || it_yErrors == triggerToYErrors.end()) {
            continue;
        }

        const std::vector<double>& pTCenters = it_pT->second;
        const std::vector<double>& values = it_yValues->second;
        const std::vector<double>& errors = it_yErrors->second;

        if (pTCenters.empty()) {
            continue; // No data for this trigger
        }

        TGraphErrors* graph = new TGraphErrors();

        for (size_t i = 0; i < pTCenters.size(); ++i) {
            double pT = pTCenters[i];
            double val = values[i];
            double err = errors[i];

            // Exclude points beyond pTExclusionMax
            if (pT >= pTExclusionMax) {
                continue;
            }

            // Find the bin
            int binIndex = -1;
            for (size_t j = 0; j < pT_bins.size(); ++j) {
                if (pT_bins[j].first >= pTExclusionMax) {
                    break;
                }
                if (pT >= pT_bins[j].first && pT < pT_bins[j].second) {
                    binIndex = (int)j;
                    break;
                }
            }
            if (binIndex == -1) {
                continue;
            }

            double binLowEdge = pT_bins[binIndex].first;
            double binUpEdge = pT_bins[binIndex].second;
            double binCenter = (binLowEdge + binUpEdge) / 2.0;

            double xPos = binCenter + xOffset;
            int pointIndex = graph->GetN();
            graph->SetPoint(pointIndex, xPos, val);
            graph->SetPointError(pointIndex, 0, err); // No x errors

            // Update y-axis range
            if (val - err < yMinData) yMinData = val - err;
            if (val + err > yMaxData) yMaxData = val + err;
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
        double markerSize = (mesonName == "#eta") ? 0.72 : 0.85;
        graph->SetMarkerSize(markerSize);
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

        graphs.push_back(graph);
    }

    // Adjust y-axis range
    if (yMinData >= yMaxData) {
        // If no valid data points were found
        yMinData = 0.0;
        yMaxData = 1.0;
    }
    double yRange = yMaxData - yMinData;
    double yBuffer = yRange * yBufferFraction;
    double yMin = yMinData - yBuffer;
    double yMax = yMaxData + yBuffer;

    if (mesonName == "#eta") {
        // If you want to force a specific range for eta, adjust here:
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

    double tickSize = (yAxisMax - yAxisMin) * 0.02;
    double labelOffset = (yAxisMax - yAxisMin) * 0.05;
    TLatex latex;
    latex.SetTextSize(0.035);
    latex.SetTextAlign(22);

    // Draw x-axis line
    TLine xAxisLine(xMin, yAxisMin, xMax, yAxisMin);
    xAxisLine.Draw();

    for (size_t i = 0; i < binEdges.size(); ++i) {
        double xPos = binEdges[i];
        double yPos = yAxisMin;

        // Draw tick
        TLine* tick = new TLine(xPos, yPos, xPos, yPos - tickSize);
        tick->Draw();

        // Label
        double pTValue = binEdges[i];
        std::ostringstream labelStream;
        labelStream << std::fixed << std::setprecision(1) << pTValue;
        std::string label = labelStream.str();

        latex.DrawLatex(xPos, yPos - labelOffset, label.c_str());
    }

    canvas.RedrawAxis();

    // Add cut combination information
    TLatex labelText;
    labelText.SetNDC();
    labelText.SetTextSize(0.032);

    TLatex valueText;
    valueText.SetNDC();
    valueText.SetTextSize(0.028);

    labelText.DrawLatex(0.2, 0.87, "#font[62]{ECore #geq}");
    {
        std::ostringstream eCoreWithUnit;
        eCoreWithUnit << cutData.clusECore << "   GeV";
        valueText.DrawLatex(0.32, 0.87, eCoreWithUnit.str().c_str());
    }

    labelText.DrawLatex(0.2, 0.82, "#font[62]{#chi^{2} <}");
    {
        std::ostringstream chiStr;
        chiStr << cutData.chi;
        valueText.DrawLatex(0.27, 0.82, chiStr.str().c_str());
    }

    labelText.DrawLatex(0.2, 0.77, "#font[62]{Asymmetry <}");
    {
        std::ostringstream asymmetryStr;
        asymmetryStr << cutData.asymmetry;
        valueText.DrawLatex(0.37, 0.77, asymmetryStr.str().c_str());
    }

    // Ensure directory
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
                0.13, 0.22 // yMin and yMax
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
            
            std::string outputFilePathSigma = plotDirectory + "/" + cutCombination + "/sigmaPi0_vs_pT_Overlay.png";
            generateOverlayInvariantMassPlot(
                cutData.triggerToPtCentersPi0,
                cutData.triggerToSigmaPi0Values,
                cutData.triggerToSigmaPi0Errors,
                cutData.triggersInData,
                pT_bins,
                "#pi^{0}",
                "#sigma_{#pi 0}",
                0.1,  // yBufferFraction
                8.0,  // pTExclusionMax
                                             outputFilePathSigma,
                plotDirectory,
                cutCombination,
                cutData,
                TriggerConfig::triggerColorMap,
                TriggerConfig::triggerNameMap
            );
            
            
            std::string outputFilePathSb = plotDirectory + "/" + cutCombination + "/sbRatioPi0_vs_pT_Overlay.png";
            generateOverlayInvariantMassPlot(
                cutData.triggerToPtCentersPi0,
                cutData.triggerToSignalToBackgroundPi0Ratios,
                cutData.triggerToSignalToBackgroundPi0Errors,
                cutData.triggersInData,
                pT_bins,
                "#pi^{0}",
                "Signal-to-Background Ratio",
                0.1,  // yBufferFraction
                8.0,  // pTExclusionMax
                outputFilePathSb,
                plotDirectory,
                cutCombination,
                cutData,
                TriggerConfig::triggerColorMap,
                TriggerConfig::triggerNameMap
            );
            
            std::string outputFilePathResolution = plotDirectory + "/" + cutCombination + "/resolutionPi0_vs_pT_Overlay.png";
            generateOverlayInvariantMassPlot(
                cutData.triggerToPtCentersPi0,
                cutData.triggerToResolutionPi0Values,
                cutData.triggerToResolutionPi0Errors,
                cutData.triggersInData,
                pT_bins,
                "#pi^{0}",
                "#pi^{0} Fit Resolution",
                0.1,  // yBufferFraction
                8.0,  // pTExclusionMax
                                             outputFilePathResolution,
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
            
            std::string outputFilePathSigma = plotDirectory + "/" + cutCombination + "/sigmaEta_vs_pT_Overlay.png";
            generateOverlayInvariantMassPlot(
                cutData.triggerToPtCentersEta,
                cutData.triggerToSigmaEtaValues,
                cutData.triggerToSigmaEtaErrors,
                cutData.triggersInData,
                pT_bins,
                "#eta",
                "#sigma_{#eta}",
                0.1,  // yBufferFraction
                8.0,  // pTExclusionMax
                                             outputFilePathSigma,
                plotDirectory,
                cutCombination,
                cutData,
                TriggerConfig::triggerColorMap,
                TriggerConfig::triggerNameMap
            );
            
            
            std::string outputFilePathSb = plotDirectory + "/" + cutCombination + "/sbRatioEta_vs_pT_Overlay.png";
            generateOverlayInvariantMassPlot(
                cutData.triggerToPtCentersEta,
                cutData.triggerToSignalToBackgroundEtaRatios,
                cutData.triggerToSignalToBackgroundEtaErrors,
                cutData.triggersInData,
                pT_bins,
                "#eta",
                "Signal-to-Background Ratio",
                0.1,  // yBufferFraction
                8.0,  // pTExclusionMax
                outputFilePathSb,
                plotDirectory,
                cutCombination,
                cutData,
                TriggerConfig::triggerColorMap,
                TriggerConfig::triggerNameMap
            );
            
            std::string outputFilePathResolution = plotDirectory + "/" + cutCombination + "/resolutionEta_vs_pT_Overlay.png";
            generateOverlayInvariantMassPlot(
                cutData.triggerToPtCentersEta,
                cutData.triggerToResolutionEtaValues,
                cutData.triggerToResolutionEtaErrors,
                cutData.triggersInData,
                pT_bins,
                "#eta",
                "#eta Fit Resolution",
                0.1,  // yBufferFraction
                8.0,  // pTExclusionMax
                                             outputFilePathResolution,
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





std::map<std::tuple<
    std::string, // TriggerGroupName
    std::string, // TriggerName
    float,       // ECore
    float,       // Chi
    float,       // Asymmetry
    float,       // pTMin
    float,       // pTMax
    float,       // isoEtMin
    float,       // isoEtMax
    std::string  // showerCutLabel
>, DataStructures::IsolationData> dataMap_inMassWindow;

std::map<std::tuple<
    std::string,
    std::string,
    float, float, float,
    float, float,
    float, float,
    std::string
>, DataStructures::IsolationData> dataMap_outsideMassWindow;

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
    float isoEtMax,
    const std::string& showerCutLabel) {

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

///
/// ProcessIsolationEnergyHistogramsWithCuts
///
/// This function:
///   1) Loops over a given set of trigger names (directories) in a ROOT file
///   2) Finds relevant histograms (like "isolatedPhotonCount_...", "allPhotonCount_...", etc.)
///   3) Parses each histogram's name to extract cut parameters (ECore, Chi, Asym, pT range, isoEt range, etc.)
///   4) Plots and saves a PNG of each histogram into an organized directory structure
///   5) Fills global data structures for later analysis:
///         - dataMap_inMassWindow[...]  : IsolatedPhotonLog for "inMassWindow"
///         - dataMap_outsideMassWindow[: IsolatedPhotonLog for "outsideMassWindow"
///         - totalPhotonMap[...]        : TotalPhotonLog for all photons
///         - pTweightingMap[...]        : PtWeightingLog for average pT
///
/// **Global Maps Filled**:
///   1) `dataMap_inMassWindow` and `dataMap_outsideMassWindow`
///      - Key: A 10‐tuple =
///         ( triggerGroupName,
///           triggerName,
///           eCore, chi, asym,
///           pTMin, pTMax,
///           isoEtMin, isoEtMax,
///           showerCutLabel )
///      - Value: `DataStructures::IsolationData`
///         (counts how many “isolated” photons passed the given cuts)
///      - The difference is: if the histogram name says "inMassWindow", we store in `dataMap_inMassWindow`.
///        If it says "outsideMassWindow", we store in `dataMap_outsideMassWindow`.
///
///   2) `totalPhotonMap`
///      - Key: An 8‐tuple =
///         ( triggerGroupName,
///           triggerName,
///           eCore, chi, asym,
///           pTMin, pTMax,
///           massWindowLabel )
///      - Value: `DataStructures::TotalPhotonLog`
///         (the total photon counts for that set of cuts)
///
///   3) `pTweightingMap`
///      - Key: Same 8‐tuple as `totalPhotonMap`.
///      - Value: `DataStructures::PtWeightingLog` (stores “weightedAveragePt” for that bin)
///
///
/// \param inputFile       - The ROOT file containing directories named by each trigger
/// \param plotDirectory   - Path to folder(s) where we save the .png images
/// \param triggers        - A list of triggers (directory names) we look for in the ROOT file
/// \param combinationName - The "triggerGroupName" stored in the final keys (like "Photon_HighEt")
///
void ProcessIsolationEnergyHistogramsWithCuts(
    TFile* inputFile,
    const std::string& plotDirectory,
    const std::vector<std::string>& triggers,
    const std::string& combinationName)
{
    std::cout << "[DEBUG] Entering ProcessIsolationEnergyHistogramsWithCuts() for '"
              << combinationName << "' with " << triggers.size() << " triggers.\n";

    // The raw combinationName is used as the "triggerGroupName" in the final map keys.
    // For example: "Photon_5GeV_plus_MBD" or "Photon_10_GeV" etc.
    std::string triggerGroupName = combinationName;

    //--------------------------------------------------------------------------------
    // A safer parse function that returns false if the crucial capture groups are empty
    //--------------------------------------------------------------------------------
    // This lambda tries to interpret the histogram name of the form:
    //   "isolatedPhotonCount_E1_Chi3_Asym0point7_inMassWindow_isoEt_2to4_pT_10to12_withShowerShapeCuts_TrigName"
    // extracting all the relevant pieces (ECore, Chi, Asym, pT min, pT max, isoEt range, etc.).
    auto parseIsolationHistNameWithCuts = [&](const std::string& histName)
      -> std::pair<bool, std::tuple<
            DataStructures::CutValues,  // cuts (ECore, Chi, Asym)
            std::string,                // massWindowLabel
            float, float,               // pTMin, pTMax
            std::string,                // actual triggerName
            std::string,                // histType
            bool, float, float,         // hasIsoEtRange, isoEtMin, isoEtMax
            std::string                 // showerCutLabel => "withShowerShapeCuts"/"withoutShowerShapeCuts"/""
         >>
    {
        using namespace std;
        using namespace DataStructures;

        bool success = false;
        CutValues cuts;             // Will hold eCore, chi, asym.
        std::string massWindowLabel;
        float pTMin       = -1.f;   // e.g. pT_10
        float pTMax       = -1.f;   // e.g.  to 12
        bool  hasIsoRange = false;  // e.g. isoEt_2to4 => isoRange
        float isoEtMin    = 0.f;
        float isoEtMax    = 0.f;
        std::string triggerName;
        std::string histType;
        std::string showerCutLabel;

        // Helper to replace "point" with '.' and parse float carefully
        auto convert = [&](const std::string& inputStr) -> float
        {
            if (inputStr.empty()) return 0.f; // early return
            std::string tmp = inputStr;
            size_t pos = tmp.find("point");
            if (pos != std::string::npos)
            {
                // Replace the literal substring "point" with "."
                tmp.replace(pos, 5, ".");
            }
            try {
                return std::stof(tmp);
            }
            catch (...)
            {
                return 0.f;
            }
        };

        // ----------------------------------------------------------------------------
        // This REGEX carefully parses a histogram name like:
        //  isolatedPhotonCount_E1_Chi3_Asym0point7_outsideMassWindow_isoEt_-100to6_pT_10to12_withShowerShapeCuts_Photon_5_GeV_plus_MBD_NS_geq_1
        //
        // Breakdown of the capturing groups:
        //  (1) => histType: "isolatedPhotonCount" | "allPhotonCount" | "ptPhoton" | "h1_cluster_iso_Et" | "h2_cluster_iso_Et"
        //  (2) => E (like "1" or "1point5")
        //  (3) => Chi (like "3" or "4point5")
        //  (4) => Asym (like "0point7")
        //  (5) => optional "inMassWindow" or "outsideMassWindow"
        //  (6) => isoEtMin (like "-100")
        //  (7) => isoEtMax (like "6")
        //  (8) => pTMin (like "10")
        //  (9) => pTMax (like "12")
        //  (10) => optional "withShowerShapeCuts" or "withoutShowerShapeCuts"
        //  (11) => everything else at the end => triggerName
        // ----------------------------------------------------------------------------
        static std::regex histPattern(
            R"REGEX(^(isolatedPhotonCount|allPhotonCount|ptPhoton|h[12]_cluster_iso_Et)_E([-+]?\d+(?:point\d+)?)_Chi([-+]?\d+(?:point\d+)?)_Asym([-+]?\d+(?:point\d+)?)(?:_(inMassWindow|outsideMassWindow))?(?:_isoEt_([-+]?\d+(?:point\d+)?)to([-+]?\d+(?:point\d+)?))?_pT_([-+]?\d+(?:point\d+)?)to([-+]?\d+(?:point\d+)?)(?:_(withShowerShapeCuts|withoutShowerShapeCuts)_)?([^ ]+)$)REGEX"
        );

        std::smatch match;
        if (!std::regex_match(histName, match, histPattern)) {
            return { false, {} };
        }

        // The capturing groups are, in order:
        //   match[1]  => histType
        //   match[2]  => ECore
        //   match[3]  => Chi
        //   match[4]  => Asym
        //   match[5]  => massWindowLabel   (optionally inMassWindow / outsideMassWindow)
        //   match[6]  => isoEtMin
        //   match[7]  => isoEtMax
        //   match[8]  => pTMin
        //   match[9]  => pTMax
        //   match[10] => showerCutLabel    (optionally with/withoutShowerShapeCuts)
        //   match[11] => triggerName       (the rest of the string)
        //
        // match.size() should be 12 total. We do a debug print below.
        std::cout << "[PARSE DEBUG] histName: " << histName << "\n";
        for (size_t i = 0; i < match.size(); ++i) {
            std::cout << "    match[" << i << "]: '" << match[i].str() << "'\n";
        }

        if (match.size() < 12)
        {
            // Not all groups captured
            std::cerr << "[PARSE DEBUG] Not enough captures => " << match.size() << "\n";
            return { false, {} };
        }

        // If we get here, it means the regex matched and we have the required #groups
        success = true;

        histType         = match[1].str();
        cuts.clusECore   = convert(match[2].str());
        cuts.chi         = convert(match[3].str());
        cuts.asymmetry   = convert(match[4].str());
        massWindowLabel  = match[5].str(); // might be empty

        std::string isoEtMinStr = match[6].str();
        std::string isoEtMaxStr = match[7].str();
        if (!isoEtMinStr.empty() && !isoEtMaxStr.empty())
        {
            isoEtMin    = convert(isoEtMinStr);
            isoEtMax    = convert(isoEtMaxStr);
            hasIsoRange = true;
        }

        std::string pTMinStr = match[8].str();
        std::string pTMaxStr = match[9].str();
        if (pTMinStr.empty() || pTMaxStr.empty()) {
            // crucial pT range was not found => fail parse
            std::cerr << "[PARSE DEBUG] Missing pT ranges => pTMinStr='"
                      << pTMinStr << "' pTMaxStr='" << pTMaxStr << "'\n";
            return { false, {} };
        }
        pTMin = convert(pTMinStr);
        pTMax = convert(pTMaxStr);

        showerCutLabel = match[10].str();      // might be empty
        triggerName    = match[11].str();      // the final piece

        // If the final triggerName has a leading underscore, strip it
        if (!triggerName.empty() && triggerName.front() == '_') {
            triggerName.erase(0, 1);
        }

        // Debug info
        std::cout << "[PARSE DEBUG] => histType='" << histType << "'\n"
                  << "                E=" << cuts.clusECore << "\n"
                  << "                Chi=" << cuts.chi << "\n"
                  << "                Asym=" << cuts.asymmetry << "\n"
                  << "                massWindowLabel='" << massWindowLabel << "'\n"
                  << "                isoEtMin=" << isoEtMin << ", isoEtMax=" << isoEtMax
                  << " (hasIsoRange=" << std::boolalpha << hasIsoRange << ")\n"
                  << "                pTMin=" << pTMin << ", pTMax=" << pTMax << "\n"
                  << "                showerCutLabel='" << showerCutLabel << "'\n"
                  << "                triggerName='" << triggerName << "'\n";

        return {
            success,
            std::make_tuple(
                cuts,
                massWindowLabel,
                pTMin, pTMax,
                triggerName,
                histType,
                hasIsoRange, isoEtMin, isoEtMax,
                showerCutLabel
            )
        };
    };

    // For each trigger name in the input vector
    for (const auto& trig : triggers)
    {
        std::cout << "[DEBUG] Checking directory for trigger='" << trig << "'...\n";
        // Retrieve the TDirectory matching that trigger name
        TDirectory* trigDir = inputFile->GetDirectory(trig.c_str());
        if (!trigDir)
        {
            std::cerr << "[WARNING] Trigger directory '" << trig
                      << "' not found in file. Skipping.\n";
            continue;
        }

        // Iterate over all objects (keys) in that directory
        TIter nextKey(trigDir->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*)nextKey()))
        {
            std::cout << "[DEBUG] Found object key->GetName()='" << key->GetName()
                      << "', class='" << key->GetClassName() << "'\n";

            // Attempt to read the TObject from that key
            TObject* obj = key->ReadObj();
            if (!obj) {
                std::cout << "[DEBUG] key->ReadObj() returned null.\n";
                continue;
            }

            // We only care about TH1-based histograms (TH1, TH2, etc.)
            if (!obj->InheritsFrom(TH1::Class())) {
                std::cout << "[DEBUG] Object '" << obj->GetName()
                          << "' is not TH1-based. Skipping.\n";
                delete obj;
                continue;
            }
            // Specifically skip TProfile objects
            if (obj->InheritsFrom(TProfile::Class())) {
                std::cout << "[DEBUG] Skipping TProfile: " << obj->GetName() << "\n";
                delete obj;
                continue;
            }

            // Now we have a TH1 or TH2
            TH1* hist = dynamic_cast<TH1*>(obj);
            if (!hist) {
                std::cerr << "[ERROR] dynamic_cast<TH1*> failed for '"
                          << obj->GetName() << "'\n";
                delete obj;
                continue;
            }

            std::string histName = hist->GetName();

            // Attempt to parse the histogram name
            auto [parsedOk, dataTuple] = parseIsolationHistNameWithCuts(histName);
            if (!parsedOk)
            {
                std::cout << "[DEBUG] parseIsolationHistNameWithCuts() => FAIL for '"
                          << histName << "'. Skipping.\n";
                delete obj;
                continue;
            }

            // If parsing succeeded, destructure the result
            auto [cuts,
                  massWindowLabel,
                  pTMin, pTMax,
                  actualTriggerName,
                  histType,
                  hasIsoEtRange, isoEtMin, isoEtMax,
                  showerCutLabel] = dataTuple;

            // Sometimes, the histogram’s embedded “triggerName” may differ from
            // the directory name `trig`. If they mismatch, skip it to avoid confusion.
            if (actualTriggerName != trig) {
                std::cout << "[DEBUG] 'actualTriggerName' mismatch => skip.\n";
                delete obj;
                continue;
            }

            // ------------------------------------------------------------------
            // 4) Build the output directories for saving the .png
            // ------------------------------------------------------------------
            // We create subdirectories that reflect ECore, Chi, Asym, pT range, isoEt range, etc.
            std::ostringstream cutDirStr;
            cutDirStr << plotDirectory
                      << "/E"   << Utils::formatToThreeSigFigs(cuts.clusECore)
                      << "_Chi" << Utils::formatToThreeSigFigs(cuts.chi)
                      << "_Asym"<< Utils::formatToThreeSigFigs(cuts.asymmetry);
            std::string cutDirPath = cutDirStr.str();
            gSystem->mkdir(cutDirPath.c_str(), true);

            std::string isoDirPath = cutDirPath + "/isolationEnergies";
            gSystem->mkdir(isoDirPath.c_str(), true);

            // pT subdirectory
            std::ostringstream ptDirStr;
            ptDirStr << isoDirPath << "/pT_"
                     << Utils::formatToThreeSigFigs(pTMin)
                     << "_to_"
                     << Utils::formatToThreeSigFigs(pTMax);
            std::string ptDirPath = ptDirStr.str();
            gSystem->mkdir(ptDirPath.c_str(), true);

            // If there's an isoEt range, go one level deeper
            std::string finalDirPath = ptDirPath;
            if (hasIsoEtRange)
            {
                std::ostringstream isoEtSubStr;
                isoEtSubStr << ptDirPath
                            << "/isoEt_"
                            << Utils::formatToThreeSigFigs(isoEtMin)
                            << "_to_"
                            << Utils::formatToThreeSigFigs(isoEtMax);
                finalDirPath = isoEtSubStr.str();
                gSystem->mkdir(finalDirPath.c_str(), true);
            }

            // 5) The PNG file path we’ll save
            std::string outPngPath = finalDirPath + "/" + histName + ".png";

            // 6) Draw + Save the histogram
            //    Make a canvas, draw the hist, possibly set log scale, and write to disk
            TCanvas canvas("canvas","Histogram Canvas",800,600);
            if (hist->InheritsFrom(TH2::Class()))
            {
                hist->Draw("COLZ");
                canvas.SetLogz();
            }
            else
            {
                hist->SetStats(true);
                hist->Draw("HIST");
                canvas.SetLogy();
            }

            // 7) Add custom labels (if you have a helper function for that).
            //    This is presumably to add text about eCore, chi, asym, etc. onto the canvas.
            AddLabelsToCanvas_isoHistsWithCuts(
                cuts,
                massWindowLabel,
                trig,                 // the actual trigger
                triggerGroupName,     // combinationName
                pTMin, pTMax,
                histType,
                hasIsoEtRange,
                isoEtMin, isoEtMax,
                showerCutLabel
            );

            // Save the canvas to outPngPath
            canvas.SaveAs(outPngPath.c_str());
            std::cout << "[INFO] Saved => " << outPngPath << "\n";

            // -------------------------------------------------------------
            //  Fill the data structures if relevant
            // -------------------------------------------------------------
            // (A) Build "totalKey" ignoring isoEt range
            //     This is the 8-tuple key for totalPhotonMap and pTweightingMap:
            //       (groupName, trig, eCore, chi, asym, pTMin, pTMax, massWindowLabel)
            auto totalKey = std::make_tuple(
                triggerGroupName,
                trig, // raw name
                cuts.clusECore,
                cuts.chi,
                cuts.asymmetry,
                pTMin,
                pTMax,
                massWindowLabel
            );

            // (B) Based on histType => fill your global maps
            //     1) "isolatedPhotonCount" => dataMap_inMassWindow / dataMap_outsideMassWindow
            //     2) "allPhotonCount" => totalPhotonMap
            //     3) "ptPhoton"       => pTweightingMap
            if ((histType == "isolatedPhotonCount") && hasIsoEtRange)
            {
                // Build an isoKey that includes showerCutLabel, isoEt, etc.
                // Key is a 10-tuple:
                //   ( trigGroupName, trig, eCore, chi, asym, pTMin, pTMax, isoEtMin, isoEtMax, showerCutLabel )
                auto isoKey = std::make_tuple(
                    triggerGroupName,
                    trig,
                    cuts.clusECore,
                    cuts.chi,
                    cuts.asymmetry,
                    pTMin,
                    pTMax,
                    isoEtMin,
                    isoEtMax,
                    showerCutLabel  // <--- crucial
                );

                // Build the IsolatedPhotonLog to store # of isolated photons
                DataStructures::IsolationData isoData;
                isoData.triggerGroupName  = triggerGroupName;
                isoData.triggerName       = trig;
                isoData.clusECore         = cuts.clusECore;
                isoData.chi               = cuts.chi;
                isoData.asymmetry         = cuts.asymmetry;
                isoData.pTMin             = pTMin;
                isoData.pTMax             = pTMax;
                isoData.isoMin            = isoEtMin;
                isoData.isoMax            = isoEtMax;
                isoData.isolatedCounts   = static_cast<int>(hist->GetEntries());
                isoData.massWindowLabel   = massWindowLabel;
                isoData.showerCutLabel = showerCutLabel;

                // Now store it in the appropriate map based on massWindowLabel
                if (massWindowLabel == "inMassWindow") {
                    // dataMap_inMassWindow => key => isoKey
                    dataMap_inMassWindow[isoKey] = isoData;
                } else if (massWindowLabel == "outsideMassWindow") {
                    // dataMap_outsideMassWindow => key => isoKey
                    dataMap_outsideMassWindow[isoKey] = isoData;
                } else {
                    std::cerr << "[WARNING] Unknown massWindowLabel='"
                              << massWindowLabel << "' => skip.\n";
                }
            }
            else if (histType == "allPhotonCount")
            {
                // This histogram indicates total # of photons (not necessarily isolated).
                // We store that in totalPhotonMap with an 8-tuple key.
                auto totalKey = std::make_tuple(
                    triggerGroupName,
                    trig,
                    cuts.clusECore,
                    cuts.chi,
                    cuts.asymmetry,
                    pTMin,
                    pTMax,
                    massWindowLabel
                );
                DataStructures::TotalPhotonLog totalLog;
                totalLog.triggerGroupName = triggerGroupName;
                totalLog.triggerName      = trig;
                totalLog.clusECore        = cuts.clusECore;
                totalLog.chi              = cuts.chi;
                totalLog.asymmetry        = cuts.asymmetry;
                totalLog.pTMin            = pTMin;
                totalLog.pTMax            = pTMax;
                totalLog.totalEntries     = static_cast<int>(hist->GetEntries());
                totalLog.massWindowLabel  = massWindowLabel;

                totalPhotonMap[totalKey] = totalLog;
            }
            else if (histType == "ptPhoton")
            {
                // This histogram typically accumulates the pT distribution of photons
                // We compute a weighted average pT and store it in pTweightingMap.
                double wSum=0, countVal=0;
                for(int i=1;i<=hist->GetNbinsX();i++){
                    double bc=hist->GetBinContent(i);
                    double xCtr=hist->GetBinCenter(i);
                    wSum += bc*xCtr;
                    countVal+=bc;
                }
                double wAvgPt= (countVal>0? wSum/countVal : 0);

                auto totalKey = std::make_tuple(
                    triggerGroupName,
                    trig,
                    cuts.clusECore,
                    cuts.chi,
                    cuts.asymmetry,
                    pTMin,
                    pTMax,
                    massWindowLabel
                );
                DataStructures::PtWeightingLog ptLog;
                ptLog.triggerGroupName   = triggerGroupName;
                ptLog.triggerName        = trig;
                ptLog.clusECore          = cuts.clusECore;
                ptLog.chi                = cuts.chi;
                ptLog.asymmetry          = cuts.asymmetry;
                ptLog.pTMin              = pTMin;
                ptLog.pTMax              = pTMax;
                ptLog.weightedAveragePt  = wAvgPt;
                ptLog.massWindowLabel    = massWindowLabel;

                pTweightingMap[totalKey] = ptLog;
            }

            // Clean up the histogram object after we’ve processed it
            delete obj;
        } // end while
    } // end for triggers
}

void WriteIsolationDataToCSV(
    const std::string& outputFilePath,
    const std::map<
        std::tuple<
            std::string, // TriggerGroupName
            std::string, // TriggerName
            float,       // ECore
            float,       // Chi
            float,       // Asymmetry
            float,       // pT Min
            float,       // pT Max
            float,       // isoMin
            float,       // isoMax
            std::string  // showerCutLabel
        >,
        DataStructures::IsolationData
    >& dataMap_inMassWindow,
    const std::map<
        std::tuple<
            std::string,
            std::string,
            float, float, float,
            float, float,
            float, float,
            std::string
        >,
        DataStructures::IsolationData
    >& dataMap_outsideMassWindow
)
{
    // Try opening the CSV file
    std::ofstream outFile(outputFilePath);
    if (!outFile.is_open()) {
        std::ostringstream err;
        err << "[ERROR] Could not open CSV for writing: " << outputFilePath;
        std::cerr << err.str() << "\n";
        throw std::runtime_error(err.str());
    }

    // Write the CSV header
    outFile
      << "TriggerGroupName,TriggerName,ECore,Chi,Asymmetry,pT Min,pT Max,"
      << "isoMin,isoMax,ShowerShapeStatus,"
      << "Isolated Counts,Total Counts,Isolated/Total,Statistical Error,"
      << "Weighted pT,Bin Width,Bin Center,"
      << "Isolated Yield,Isolated Yield Error,MassWindowLabel\n";

    // ------------------------------------------------------------
    // 1) Write rows from dataMap_inMassWindow
    // ------------------------------------------------------------
    for (const auto& kv : dataMap_inMassWindow) {
        // Key is the 10‐tuple
        const auto& isoKey  = kv.first;
        const auto& isoData = kv.second; // final "IsolationData"

        // Unpack the key
        std::string trigGroup  = std::get<0>(isoKey);
        std::string trigName   = std::get<1>(isoKey);
        float eCore            = std::get<2>(isoKey);
        float chi              = std::get<3>(isoKey);
        float asym             = std::get<4>(isoKey);
        float pTMin            = std::get<5>(isoKey);
        float pTMax            = std::get<6>(isoKey);
        float isoMin           = std::get<7>(isoKey);
        float isoMax           = std::get<8>(isoKey);
        std::string showerCut  = std::get<9>(isoKey);

        // This entire map => "inMassWindow"
        std::string massWindowLabel = "inMassWindow";

        // Grab data from isoData
        int    isoCounts    = isoData.isolatedCounts;
        int    totalCounts  = isoData.totalCounts;
        double ratioVal     = isoData.ratio;
        double errVal       = isoData.error;
        double weightedPt   = isoData.weightedPt;
        double binWidth     = isoData.binWidth;
        double binCenter    = isoData.binCenter;
        double isoYield     = isoData.isolatedYield;
        double isoYieldErr  = isoData.isolatedYieldError;

        // Optional: do a quick sanity check
        if (totalCounts < isoCounts) {
            std::cerr << "[WARNING] inMassWindow: totalCounts < isoCounts for ("
                      << trigGroup << "," << trigName << ") => logic error?\n";
        }
        if (ratioVal < 0.0) {
            std::cerr << "[WARNING] inMassWindow: ratio < 0.0 => " << ratioVal
                      << " for (" << trigGroup << "," << trigName << ").\n";
        }

        // Write the row to CSV
        outFile
          << trigGroup << ","
          << trigName  << ","
          << eCore     << ","
          << chi       << ","
          << asym      << ","
          << pTMin     << ","
          << pTMax     << ","
          << isoMin    << ","
          << isoMax    << ","
          << showerCut << ","   // showerCut
          << isoCounts << ","
          << totalCounts << ","
          << ratioVal << ","
          << errVal << ","
          << weightedPt << ","
          << binWidth << ","
          << binCenter << ","
          << isoYield << ","
          << isoYieldErr << ","
          << massWindowLabel
          << "\n";
    }

    // ------------------------------------------------------------
    // 2) Write rows from dataMap_outsideMassWindow
    // ------------------------------------------------------------
    for (const auto& kv : dataMap_outsideMassWindow) {
        // Key is the same 10‐tuple
        const auto& isoKey  = kv.first;
        const auto& isoData = kv.second;

        // Unpack key
        std::string trigGroup = std::get<0>(isoKey);
        std::string trigName  = std::get<1>(isoKey);
        float eCore           = std::get<2>(isoKey);
        float chi             = std::get<3>(isoKey);
        float asym            = std::get<4>(isoKey);
        float pTMin           = std::get<5>(isoKey);
        float pTMax           = std::get<6>(isoKey);
        float isoMin          = std::get<7>(isoKey);
        float isoMax          = std::get<8>(isoKey);
        std::string showerCut = std::get<9>(isoKey);

        // This entire map => "outsideMassWindow"
        std::string massWindowLabel = "outsideMassWindow";

        // Grab data
        int isoCounts     = isoData.isolatedCounts;
        int totalCounts   = isoData.totalCounts;
        double ratioVal   = isoData.ratio;
        double errVal     = isoData.error;
        double weightedPt = isoData.weightedPt;
        double binWidth   = isoData.binWidth;
        double binCenter  = isoData.binCenter;
        double isoYield   = isoData.isolatedYield;
        double isoYieldErr= isoData.isolatedYieldError;

        // Quick sanity checks
        if (totalCounts < isoCounts) {
            std::cerr << "[WARNING] outsideMassWindow: totalCounts < isoCounts for ("
                      << trigGroup << "," << trigName << ").\n";
        }
        if (ratioVal < 0.0) {
            std::cerr << "[WARNING] outsideMassWindow: ratio < 0.0 => " << ratioVal
                      << " for (" << trigGroup << "," << trigName << ").\n";
        }

        // Write row
        outFile
          << trigGroup << ","
          << trigName  << ","
          << eCore     << ","
          << chi       << ","
          << asym      << ","
          << pTMin     << ","
          << pTMax     << ","
          << isoMin    << ","
          << isoMax    << ","
          << showerCut << ","
          << isoCounts << ","
          << totalCounts << ","
          << ratioVal << ","
          << errVal << ","
          << weightedPt << ","
          << binWidth << ","
          << binCenter << ","
          << isoYield << ","
          << isoYieldErr << ","
          << massWindowLabel
          << "\n";
    }

    // Done writing
    outFile.close();
    std::cout << "[INFO] CSV file written => " << outputFilePath << "\n";
}


void readDataFromCSV(
    const std::string& filename,
    // The 10-element key for inMassMap
    std::map<
        std::tuple<
            std::string, // TriggerGroupName
            std::string, // TriggerName
            float,       // ECore
            float,       // Chi
            float,       // Asymmetry
            float,       // pT Min
            float,       // pT Max
            float,       // isoMin
            float,       // isoMax
            std::string  // showerCutLabel
        >,
        DataStructures::IsolationData
    >& dataMap_inMassWindow,
    // The 10-element key for outMassMap
    std::map<
        std::tuple<
            std::string,
            std::string,
            float, float, float,
            float, float,
            float, float,
            std::string
        >,
        DataStructures::IsolationData
    >& dataMap_outsideMassWindow
)
{
    // Attempt to open
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::ostringstream err;
        err << "[ERROR] Cannot open CSV => " << filename;
        std::cerr << err.str() << "\n";
        throw std::runtime_error(err.str());
    }

    // Skip header
    std::string header;
    if(!std::getline(fin, header)) {
        std::ostringstream err;
        err << "[ERROR] File is empty => " << filename;
        std::cerr << err.str() << "\n";
        throw std::runtime_error(err.str());
    }
    std::cout << "[INFO] Reading CSV from " << filename << ", skipping header:\n   "
              << header << "\n";

    int lineNum = 1;
    std::string line;
    while (std::getline(fin, line)) {
        lineNum++;
        if(line.empty()) continue;

        // Split by comma
        std::vector<std::string> tokens;
        {
            std::istringstream iss(line);
            std::string cell;
            while(std::getline(iss, cell, ',')) {
                tokens.push_back(cell);
            }
        }

        // We expect exactly 20 columns
        if (tokens.size() < 20) {
            std::ostringstream warn;
            warn << "[WARNING] line " << lineNum << ": only "
                 << tokens.size() << " columns => need 20. Skipping.";
            std::cerr << warn.str() << "\n";
            continue;
        }

        // Helper lambdas
        auto trim = [&](std::string & s){
            while(!s.empty() && (s.front()==' ' || s.front()=='\t')) s.erase(0,1);
            while(!s.empty() && (s.back()==' '  || s.back()=='\t'))  s.pop_back();
        };
        auto safeFloat=[&](const std::string&s){
            std::string t = s; trim(t);
            try { return std::stof(t); }
            catch(...){
                std::ostringstream err;
                err << "[CRITICAL] line " << lineNum << " => cannot parse float: '"
                    << t << "'. Abort.";
                std::cerr << err.str() << "\n";
                throw std::runtime_error(err.str());
            }
        };
        auto safeDouble=[&](const std::string&s){
            std::string t = s; trim(t);
            try { return std::stod(t); }
            catch(...){
                std::ostringstream err;
                err << "[CRITICAL] line " << lineNum << " => cannot parse double: '"
                    << t << "'. Abort.";
                std::cerr << err.str() << "\n";
                throw std::runtime_error(err.str());
            }
        };
        auto safeInt=[&](const std::string&s){
            std::string t = s; trim(t);
            try { return std::stoi(t); }
            catch(...){
                std::ostringstream err;
                err << "[CRITICAL] line " << lineNum << " => cannot parse int: '"
                    << t << "'. Abort.";
                std::cerr << err.str() << "\n";
                throw std::runtime_error(err.str());
            }
        };

        // Now parse the columns in order:
        // (1) TriggerGroupName
        std::string trigGroupName = tokens[0]; trim(trigGroupName);
        // (2) TriggerName
        std::string trigName = tokens[1]; trim(trigName);

        // (3) ECore
        float eCore = safeFloat(tokens[2]);
        // (4) Chi
        float chi   = safeFloat(tokens[3]);
        // (5) Asym
        float asym  = safeFloat(tokens[4]);
        // (6) pTMin
        float pTMin = safeFloat(tokens[5]);
        // (7) pTMax
        float pTMax = safeFloat(tokens[6]);
        // (8) isoMin
        float isoMin = safeFloat(tokens[7]);
        // (9) isoMax
        float isoMax = safeFloat(tokens[8]);

        // (10) showerCutLabel
        std::string showerCutLabel = tokens[9];
        trim(showerCutLabel);

        // (11) IsolatedCounts
        int isoCounts = safeInt(tokens[10]);
        // (12) TotalCounts
        int totCounts = safeInt(tokens[11]);

        // (13) ratio
        double ratioVal = safeDouble(tokens[12]);
        // (14) StatError
        double statErr  = safeDouble(tokens[13]);
        // (15) WeightedPt
        double wPt      = safeDouble(tokens[14]);
        // (16) BinWidth
        double binWidth = safeDouble(tokens[15]);
        // (17) BinCenter
        double binCenter= safeDouble(tokens[16]);
        // (18) IsolatedYield
        double isoYield = safeDouble(tokens[17]);
        // (19) IsolatedYieldError
        double isoYieldErr = safeDouble(tokens[18]);
        // (20) massWindowLabel
        std::string massWindowLabel = tokens[19];
        trim(massWindowLabel);

        // Construct an IsolationData struct
        DataStructures::IsolationData isoData;
        isoData.isolatedCounts      = isoCounts;
        isoData.totalCounts         = totCounts;
        isoData.ratio               = ratioVal;
        isoData.error               = statErr;
        isoData.weightedPt          = wPt;
        isoData.binWidth            = binWidth;
        isoData.binCenter           = binCenter;
        isoData.isolatedYield       = isoYield;
        isoData.isolatedYieldError  = isoYieldErr;
        isoData.showerCutLabel      = showerCutLabel;
        isoData.massWindowLabel     = massWindowLabel;

        // Build the 10-tuple key (excluding massWindowLabel)
        auto key = std::make_tuple(
            trigGroupName,
            trigName,
            eCore,
            chi,
            asym,
            pTMin,
            pTMax,
            isoMin,
            isoMax,
            showerCutLabel
        );

        // Decide which map to store in (inMass vs. outMass)
        if (massWindowLabel == "inMassWindow") {
            dataMap_inMassWindow[key] = isoData;
        }
        else if (massWindowLabel == "outsideMassWindow") {
            dataMap_outsideMassWindow[key] = isoData;
        }
        else {
            std::ostringstream warn;
            warn << "[WARNING] line " << lineNum << ": Unknown MassWindowLabel='"
                 << massWindowLabel << "'. Skipping line.\n";
            std::cerr << warn.str();
            continue;
        }

        // Debug print
        std::cout << "[DEBUG] line " << lineNum << ": trigGroup=" << trigGroupName
                  << ", trig=" << trigName
                  << ", E=" << eCore << ", Chi=" << chi << ", Asym=" << asym
                  << ", pT=[" << pTMin << "," << pTMax << "], iso=[" << isoMin << "," << isoMax << "]"
                  << ", showerCut='" << showerCutLabel
                  << "', massWin='" << massWindowLabel << "'\n";
    }

    fin.close();
    std::cout << "[INFO] Finished reading CSV => " << filename << "\n"
              << "       dataMap_inMassWindow size=" << dataMap_inMassWindow.size() << "\n"
              << "       dataMap_outsideMassWindow size=" << dataMap_outsideMassWindow.size() << "\n";
}

/**
 * \struct SpectraGroupKey
 * A key for grouping triggers by (triggerGroupName, eCore, chi, asymmetry).
 *
 * Used to cluster data that share these 4 parameters, ignoring massWindowLabel
 * and ignoring showerCutLabel for the top-level grouping.
 */
struct SpectraGroupKey {
    std::string triggerGroupName;
    float eCore;
    float chi;
    float asymmetry;

    // Ensure sorting so we can use it in std::map
    bool operator<(const SpectraGroupKey& other) const {
        return std::tie(triggerGroupName, eCore, chi, asymmetry) <
               std::tie(other.triggerGroupName, other.eCore, other.chi, other.asymmetry);
    }
};

/**
 * \function GeneratePerTriggerSpectraPlots
 *
 * This function creates yield‐vs‐pT plots for isolated photons,
 * overlaying in‐mass‐window vs. outside‐mass‐window. For each top‐level group
 * (defined by [triggerGroupName, eCore, chi, asymmetry]), we further split
 * by (triggerName) and showerCutLabel (with/withoutShowerShapeCuts).
 *
 * We read from two global-like maps:
 *   1) dataMap_inMassWindow
 *   2) dataMap_outsideMassWindow
 * each keyed by a 10‐tuple
 *   ( trigGroupName, trigName, eCore, chi, asym, pTMin, pTMax, isoMin, isoMax, massWindowLabel ).
 *
 * Inside each 10‐tuple's IsolatedPhotonLog, we have the important fields:
 *   - pTMin, pTMax => used to compute bin center if needed
 *   - showerCutLabel => "withShowerShapeCuts" or "withoutShowerShapeCuts"
 *   - isolatedYield, isolatedYieldError => we treat them as y‐values
 *   - binCenter => can be read from isoData if you prefer
 *
 * We produce TGraphs for inMassWindow vs outsideMassWindow yields,
 * color-coded, on a log‐y canvas, with pT bins from your pT_bins vector.
 *
 * \param dataMap_inMassWindow   A map of IsolatedPhotonLog for "inMassWindow"
 * \param dataMap_outsideMassWindow  A map of IsolatedPhotonLog for "outsideMassWindow"
 * \param basePlotDirectory  Where to save the .png plots
 * \param triggerCombinationNameMap   Convert raw triggerGroupName to human-readable
 * \param triggerNameMap     Convert raw triggerName to human-readable
 * \param exclusionRanges    For ignoring certain isoEt ranges
 * \param pT_bins           pT bin edges (vector of <low, high>)
 * \param pTExclusionMax    Exclude any data where bin center >= this
 */
void GeneratePerTriggerSpectraPlots(
    const std::map<
        std::tuple<
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
        >,
        DataStructures::IsolationData
    >& dataMap_inMassWindow,

    const std::map<
        std::tuple<
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
        >,
        DataStructures::IsolationData
    >& dataMap_outsideMassWindow,

    const std::string& basePlotDirectory,
    const std::map<std::string, std::string>& triggerCombinationNameMap,
    const std::map<std::string, std::string>& triggerNameMap,
    const std::vector<std::pair<float, float>>& exclusionRanges,
    const std::vector<std::pair<double, double>>& pT_bins,
    double pTExclusionMax
)
{
    std::cout << "\n[INFO] Entering GeneratePerTriggerSpectraPlots()...\n"
              << "       dataMap_inMassWindow size=" << dataMap_inMassWindow.size() << "\n"
              << "       dataMap_outsideMassWindow size=" << dataMap_outsideMassWindow.size()
              << "\n\n";

    //
    // STEP 1) We define a simpler grouping key ignoring massWindowLabel
    //         and ignoring showerCutLabel. We'll store them in a bigMap
    //         structure for plotting.
    //
    struct SimpleGroupKey {
        std::string trigGroup;
        std::string trigName;
        float eCore;
        float chi;
        float asym;

        bool operator<(const SimpleGroupKey& o) const {
            return std::tie(trigGroup, trigName, eCore, chi, asym)
                 < std::tie(o.trigGroup, o.trigName, o.eCore, o.chi, o.asym);
        }
    };

    // bigMap[ SimpleGroupKey ][ showerCutLabel ][ massWindowLabel ][ isoEtRange => IsolatedPhotonLog ]
    // isoEtRange is a float-pair {isoMin, isoMax}
    // massWindowLabel is "inMassWindow" or "outsideMassWindow"
    using IsoKey   = std::pair<float,float>;
    using IsoMap   = std::map<IsoKey, DataStructures::IsolationData>;
    using MWmap    = std::map<std::string, IsoMap>;        // massWindow => IsoMap
    using ShowerMap= std::map<std::string, MWmap>;         // showerCut => MWmap
    std::map<SimpleGroupKey, ShowerMap> bigMap;

    // A lambda that “consumes” data from either dataMap_inMassWindow or dataMap_outsideMassWindow
    auto fillInMassData = [&](const auto& dataMap) {
        for (const auto& kv : dataMap) {
            const auto& key   = kv.first; // 10‐tuple (trigGroupName, trigName, eCore, etc., massWindowLabel)
            const auto& isoLog= kv.second; // IsolatedPhotonLog

            // Unpack the key
            std::string trigGroup = std::get<0>(key);
            std::string trigName  = std::get<1>(key);
            float eCore           = std::get<2>(key);
            float chi             = std::get<3>(key);
            float asym            = std::get<4>(key);
            float pTmin           = std::get<5>(key);
            float pTmax           = std::get<6>(key);
            float isoMin          = std::get<7>(key);
            float isoMax          = std::get<8>(key);
            std::string massWin   = std::get<9>(key); // "inMassWindow" or "outsideMassWindow"

            // Possibly skip isoEtRange if it’s in exclusionRanges
            std::pair<float,float> isoEtRange(isoMin, isoMax);
            if (std::find(exclusionRanges.begin(), exclusionRanges.end(), isoEtRange) != exclusionRanges.end()) {
                std::cout << "[DEBUG] Skipping isoEtRange={" << isoMin << "," << isoMax
                          << "} due to exclusionRanges.\n";
                continue;
            }

            // Build the “simple” group key ignoring massWindow and ignoring showerCut
            SimpleGroupKey groupK{
                trigGroup,
                trigName,
                eCore,
                chi,
                asym
            };

            // We do need showerCut from isoLog:
            std::string showerCutLabel = isoLog.showerCutLabel;
            if (showerCutLabel.empty()) {
                // If empty => we can treat as "withoutShowerShapeCuts" or
                // throw an error
                std::cerr << "[WARNING] Missing showerCutLabel => defaulting to 'unknownShowerCut'.\n";
                showerCutLabel = "unknownShowerCut";
            }

            // Insert into bigMap
            bigMap[groupK][showerCutLabel][massWin][ isoEtRange ] = isoLog;
        }
    };

    // Fill bigMap from both inMassWindow & outsideMassWindow data
    fillInMassData(dataMap_inMassWindow);
    fillInMassData(dataMap_outsideMassWindow);

    //
    // STEP 2) For each SimpleGroupKey, we have showerCut => "withShowerShapeCuts" or "withoutShowerShapeCuts".
    //         We'll produce a TCanvas overlaying “inMassWindow” vs “outsideMassWindow” yields.
    //
    int groupCount = 0;
    for (const auto& groupEntry : bigMap) {
        groupCount++;
        const auto& sKey      = groupEntry.first;   // (trigGroup, trigName, eCore, chi, asym)
        const auto& showerMap = groupEntry.second;  // showerCut => { massWindow => isoEtRange => isoLog }

        // Unpack
        std::string trigGroupName = sKey.trigGroup;
        std::string trigName      = sKey.trigName;
        float eCore               = sKey.eCore;
        float chi                 = sKey.chi;
        float asym                = sKey.asym;

        // Convert raw group/trigger names => human-readable
        std::string readableGroup = trigGroupName;
        if (auto it = triggerCombinationNameMap.find(trigGroupName); it != triggerCombinationNameMap.end()) {
            readableGroup = it->second;
        }
        std::string readableTrig = trigName;
        if (auto it2 = triggerNameMap.find(trigName); it2 != triggerNameMap.end()) {
            readableTrig = it2->second;
        }

        // Print debug
        std::cout << "[INFO] Processing group #" << groupCount
                  << ": trigGroup='" << trigGroupName
                  << "', trigName='" << trigName
                  << "', eCore=" << eCore
                  << ", chi=" << chi
                  << ", asym=" << asym
                  << "\n";

        // Prepare an output directory
        std::ostringstream dirStream;
        dirStream << basePlotDirectory << "/" << trigGroupName
                  << "/E"   << Utils::formatToThreeSigFigs(eCore)
                  << "_Chi" << Utils::formatToThreeSigFigs(chi)
                  << "_Asym"<< Utils::formatToThreeSigFigs(asym)
                  << "/Spectra_byShowerCut";
        std::string baseDir = dirStream.str();
        if (gSystem->mkdir(baseDir.c_str(), true) != 0) {
            // If mkdir fails in an unexpected way (except EEXIST) => throw
            // We'll ignore if it already exists
        }

        // We'll specifically handle each showerCut: withShowerShapeCuts / withoutShowerShapeCuts
        for (auto desiredCutLabel : {std::string("withShowerShapeCuts"), std::string("withoutShowerShapeCuts")}) {

            // see if we have data
            auto cutIter = showerMap.find(desiredCutLabel);
            if (cutIter == showerMap.end()) {
                std::cout << "[WARNING] No data found for showerCutLabel='"
                          << desiredCutLabel << "' => skipping.\n";
                continue;
            }

            // massWinMap => massWindowLabel => isoEtRange => isoLog
            const auto& massWinMap = cutIter->second;

            // We'll gather data for inMassWindow => vectors (ptCenter_in, yield_in, yerr_in)
            // and outsideMassWindow => vectors (ptCenter_out, yield_out, yerr_out)
            std::vector<double> ptCenter_in,  yield_in,  yerr_in;
            std::vector<double> ptCenter_out, yield_out, yerr_out;

            // Helper lambda to fill the above vectors by scanning isoEtRanges
            // We rely on isoLog.binCenter, isoLog.isolatedYield, isoLog.isolatedYieldError
            auto fillVectors = [&](const std::string& mwLabel,
                                   std::vector<double>& xVec,
                                   std::vector<double>& yVec,
                                   std::vector<double>& errVec)
            {
                auto itMassWin = massWinMap.find(mwLabel);
                if (itMassWin == massWinMap.end()) {
                    // no data for that massWindow => do nothing
                    return;
                }
                const auto& isoEtMap = itMassWin->second; // isoEtRange => isoLog

                for (const auto& isoRangeEntry : isoEtMap) {
                    const auto& isoRange = isoRangeEntry.first; // (isoMin, isoMax)
                    const auto& isoLog   = isoRangeEntry.second;
                    
                    // If isoLog doesn't have a binCenter or yield, that's a problem
                    double ptC  = isoLog.binCenter; // we assume you've stored the "binCenter" in isoLog
                    if (ptC >= pTExclusionMax) {
                        std::cout << "[DEBUG] ptC=" << ptC << " >= pTExclusionMax="
                                  << pTExclusionMax << ", skipping.\n";
                        continue;
                    }

                    double isoY = isoLog.isolatedYield;
                    double isoErr= isoLog.isolatedYieldError;

                    // If negative => problem => throw
                    if (isoY < 0.0) {
                        std::ostringstream err;
                        err << "[CRITICAL] Negative yield => " << isoY
                            << " for massWin=" << mwLabel
                            << ", isoRange=[" << isoRange.first << "," << isoRange.second << "]";
                        std::cerr << err.str() << "\n";
                        throw std::runtime_error(err.str());
                    }

                    xVec.push_back(ptC);
                    yVec.push_back(isoY);
                    errVec.push_back(isoErr);
                }
            };

            // Fill for inMassWindow, outsideMassWindow
            fillVectors("inMassWindow",      ptCenter_in,  yield_in,  yerr_in);
            fillVectors("outsideMassWindow", ptCenter_out, yield_out, yerr_out);

            // If both sets are empty => skip
            if (ptCenter_in.empty() && ptCenter_out.empty()) {
                std::cout << "[WARNING] No points found for showerCut='"
                          << desiredCutLabel << "' => skip.\n";
                continue;
            }

            // Make a TCanvas for the yield overlay
            std::ostringstream cname;
            cname << "canvas_" << desiredCutLabel;
            TCanvas* c1 = new TCanvas(cname.str().c_str(), "Isolated Photon Spectra", 800, 600);
            c1->SetLogy();

            // Build variable bin edges from pT_bins up to pTExclusionMax
            std::vector<double> binEdges;
            binEdges.reserve(pT_bins.size());
            for (auto& bin : pT_bins) {
                if (bin.first >= pTExclusionMax) break;
                binEdges.push_back(bin.first);
            }
            if (!binEdges.empty()) {
                double lastHigh = pT_bins[binEdges.size()-1].second;
                if (lastHigh < pTExclusionMax) {
                    binEdges.push_back(lastHigh);
                } else {
                    binEdges.push_back(pTExclusionMax);
                }
            }

            if (binEdges.size() < 2) {
                std::ostringstream err;
                err << "[CRITICAL] No valid pT bin edges. binEdges.size()="
                    << binEdges.size()
                    << " => cannot plot => abort.\n";
                std::cerr << err.str();
                delete c1;
                throw std::runtime_error(err.str());
            }

            // Make a dummy histogram to define the axis
            int nBins = binEdges.size() - 1;
            TH1F* hFrame = new TH1F("hFrame", "", nBins, &binEdges[0]);
            hFrame->SetStats(false);
            hFrame->GetXaxis()->SetTitle("Leading Cluster p_{T} [GeV]");
            hFrame->GetYaxis()->SetTitle("Isolated Photon Yield");
            hFrame->GetXaxis()->SetLabelOffset(999); // hide numeric labels
            hFrame->GetXaxis()->SetTickLength(0);
            hFrame->Draw("AXIS");

            // Build a TLegend
            TLegend* legend = new TLegend(0.55, 0.7, 0.88, 0.85);
            legend->SetBorderSize(0);
            legend->SetTextSize(0.025);

            // Build TGraphErrors for inMass
            TGraphErrors* grIn = nullptr;
            if (!ptCenter_in.empty()) {
                grIn = new TGraphErrors((int)ptCenter_in.size(),
                                        ptCenter_in.data(),
                                        yield_in.data(),
                                        nullptr,
                                        yerr_in.data());
                grIn->SetMarkerStyle(20);
                grIn->SetMarkerColor(kRed);
                grIn->SetLineColor(kRed);
                grIn->Draw("P SAME");
                legend->AddEntry(grIn, "In-Mass Window", "p");
            }

            // Build TGraphErrors for outMass
            TGraphErrors* grOut = nullptr;
            if (!ptCenter_out.empty()) {
                grOut = new TGraphErrors((int)ptCenter_out.size(),
                                         ptCenter_out.data(),
                                         yield_out.data(),
                                         nullptr,
                                         yerr_out.data());
                grOut->SetMarkerStyle(24);
                grOut->SetMarkerColor(kBlue);
                grOut->SetLineColor(kBlue);
                grOut->Draw("P SAME");
                legend->AddEntry(grOut, "Outside-Mass Window", "p");
            }

            // compute global min/max among the yields
            double gMin = 1e30, gMax=1e-30;
            auto scanMinMax = [&](const std::vector<double>& vals){
                for (auto v : vals) {
                    if (v > 0.0) {
                        gMax = std::max(gMax, v);
                        gMin = std::min(gMin, v);
                    }
                }
            };
            scanMinMax(yield_in);
            scanMinMax(yield_out);

            if (gMax < 1e-25) {
                // no positive yields => throw
                std::ostringstream err;
                err << "[CRITICAL] No positive yields for '"
                    << desiredCutLabel << "' => cannot set log scale.\n";
                std::cerr << err.str();
                delete c1;
                throw std::runtime_error(err.str());
            }
            if (gMin < 1e-25) gMin = 1e-25; // avoid zero in log scale
            hFrame->GetYaxis()->SetRangeUser(gMin * 0.5, gMax * 2.0);

            // draw legend
            legend->Draw();

            // Add textual labels
            TLatex label;
            label.SetNDC();
            label.SetTextSize(0.025);
            double xL = 0.15, yL = 0.88;

            {
                std::ostringstream oss;
                oss << "#bf{Trig Group:} " << readableGroup;
                label.DrawLatex(xL, yL, oss.str().c_str());
            }
            {
                std::ostringstream oss;
                oss << "#bf{Trigger:} " << readableTrig;
                label.DrawLatex(xL, yL - 0.04, oss.str().c_str());
            }
            {
                std::ostringstream oss;
                oss << "#bf{ECore >=} " << eCore << " GeV, #chi^{2}<" << chi << ", Asym<" << asym;
                label.DrawLatex(xL, yL - 0.08, oss.str().c_str());
            }
            {
                std::ostringstream oss;
                oss << "#bf{ShowerCut:} " << desiredCutLabel;
                label.DrawLatex(xL, yL - 0.12, oss.str().c_str());
            }

            // Build the custom x-axis
            double yLow = hFrame->GetMinimum();
            double yHigh= hFrame->GetMaximum();
            double tickSize = (std::log10(yHigh) - std::log10(yLow)) * 0.02;
            double labelOffset = (std::log10(yHigh) - std::log10(yLow)) * 0.05;

            TLine xAxisLine(binEdges.front(), yLow, binEdges.back(), yLow);
            xAxisLine.Draw("SAME");

            TLatex xLabel;
            xLabel.SetTextSize(0.035);
            xLabel.SetTextAlign(22);

            for (size_t i=0; i<binEdges.size(); i++) {
                double xPos = binEdges[i];
                TLine* tick = new TLine(xPos, yLow, xPos, yLow / std::pow(10.0, tickSize));
                tick->Draw("SAME");

                std::ostringstream lb;
                lb << std::fixed << std::setprecision(1) << binEdges[i];
                xLabel.DrawLatex(xPos, yLow / std::pow(10.0, labelOffset), lb.str().c_str());
            }

            c1->RedrawAxis();
            c1->Modified();
            c1->Update();

            // Save
            std::ostringstream outName;
            outName << baseDir << "/Spectra_"
                    << readableTrig << "_"
                    << desiredCutLabel << ".png";
            std::string outPath = outName.str();
            c1->SaveAs(outPath.c_str());
            std::cout << "[INFO] Saved => " << outPath << "\n";

            // cleanup
            delete grIn;
            delete grOut;
            delete legend;
            delete hFrame;
            delete c1;
        } // end for each showerLabel
    } // end for each SimpleGroupKey

    std::cout << "[INFO] Finished GeneratePerTriggerSpectraPlots() successfully.\n";
}



/**
 * \brief producePlotForSubset
 *
 * This function creates a ratio plot (Isolated / All photons) vs. pT
 * for a specific subset of data. That subset is identified by:
 *   - A range of isoEt (isoEtMin, isoEtMax)
 *   - A massWindowLabel ("inMassWindow" or "outsideMassWindow")
 *   - A showerCutLabel ("withShowerShapeCuts" or "withoutShowerShapeCuts")
 *
 * The actual data is passed in via a map:
 *       dataSubset[ isoEtRange ] => vector<IsolationDataWithPt>
 * Each `IsolationDataWithPt` should already have `ratio` and `error` members populated,
 * so we can directly plot them. We also show optional PHENIX reference data overlays
 * if the user sets drawRefA or drawRefB to true.
 *
 * \param dataSubset         A map of (isoEtMin, isoEtMax) => a vector of IsolationDataWithPt
 * \param subfolder          The directory path where the output PNG is saved
 * \param showerCutLabel     "withShowerShapeCuts" or "withoutShowerShapeCuts"
 * \param triggerName        The raw trigger name (e.g. "Photon_5GeV_plus_MBD")
 * \param readableTriggerName A more descriptive name (e.g. "Photon5 + MBD")
 * \param readableTriggerGroupName For labeling which group these triggers belong to
 * \param eCore              The minimum cluster ECore cut
 * \param chi                The chi^2 cut
 * \param asym               The asymmetry cut
 * \param massWindowLabel    "inMassWindow" or "outsideMassWindow"
 * \param exclusionRanges    If a certain isoEtRange is in exclusionRanges, we skip it
 * \param pT_bins           The pT bins used to derive a bin center
 * \param pTExclusionMax     We ignore data where the pT center is >= this threshold
 * \param drawRefA           Whether to overlay a reference data set (A)
 * \param drawRefB           Whether to overlay a second reference data set (B)
 * \param referencePTGamma   pT points for the reference dataset A
 * \param referenceRatio     ratio values for reference dataset A
 * \param referenceStatError ratio errors for reference dataset A
 * \param referenceTwoPTGamma pT points for the reference dataset B
 * \param referenceTwoRatio   ratio values for reference dataset B
 * \param referenceTwoStatError ratio errors for reference dataset B
 */
void producePlotForSubset(
    const std::map<std::pair<float, float>,
                   std::vector<DataStructures::IsolationDataWithPt>>& dataSubset,
    const std::string& subfolder,
    const std::string& showerCutLabel,
    const std::string& triggerName,
    const std::string& readableTriggerName,
    const std::string& readableTriggerGroupName,
    float eCore,
    float chi,
    float asym,
    const std::string& massWindowLabel,
    const std::vector<std::pair<float, float>>& exclusionRanges,
    const std::vector<std::pair<double, double>>& pT_bins,
    double pTExclusionMax,
    bool drawRefA,
    bool drawRefB,
    const std::vector<double>& referencePTGamma,
    const std::vector<double>& referenceRatio,
    const std::vector<double>& referenceStatError,
    const std::vector<double>& referenceTwoPTGamma,
    const std::vector<double>& referenceTwoRatio,
    const std::vector<double>& referenceTwoStatError
)
{
    // Basic debug print
    std::cout << "[DEBUG] producePlotForSubset(): subfolder='" << subfolder
              << "', showerCutLabel='" << showerCutLabel
              << "', dataSubset.size()=" << dataSubset.size() << "\n";

    //
    // 1) Early exit if subset is empty
    //
    if (dataSubset.empty()) {
        std::cout << "   [INFO] dataSubset is empty => skipping plot for '"
                  << showerCutLabel << "'.\n\n";
        return;
    }

    //
    // 2) Create a TCanvas for drawing the ratio plot
    //
    TCanvas canvas("canvas", "Isolated/All Photon Ratio", 800, 600);

    //
    // 3) Build pT bin edges from your pT_bins vector
    //
    std::vector<double> binEdges;
    binEdges.reserve(pT_bins.size());

    // Fill binEdges until we exceed pTExclusionMax
    for (const auto& bin : pT_bins) {
        if (bin.first >= pTExclusionMax) {
            std::cout << "[DEBUG] bin.first=" << bin.first
                      << " >= pTExclusionMax=" << pTExclusionMax
                      << " => break from bin loop.\n";
            break;
        }
        binEdges.push_back(bin.first);
    }

    // Add one last edge if possible
    if (!binEdges.empty()) {
        double lastHigh = pT_bins[ binEdges.size() - 1 ].second;
        if (lastHigh < pTExclusionMax) {
            binEdges.push_back(lastHigh);
        } else {
            binEdges.push_back(pTExclusionMax);
        }
    }

    if (binEdges.size() < 2) {
        std::cerr << "[WARNING] No valid pT bins => cannot plot => skipping '"
                  << showerCutLabel << "'.\n";
        return;
    }

    int nBins = static_cast<int>(binEdges.size()) - 1;
    if (nBins <= 0) {
        std::cerr << "[ERROR] nBins <= 0 => no valid histogram bins => skipping '"
                  << showerCutLabel << "'.\n";
        return;
    }

    //
    // 4) Build a dummy histogram to set up the axis
    //
    TH1F* hFrame = new TH1F("hFrame", "", nBins, binEdges.data());
    hFrame->SetStats(false);

    // Determine Y-axis label from massWindow
    std::string yTitle;
    if (massWindowLabel == "inMassWindow") {
        yTitle = "#frac{Isolated Photons from #pi^{0}/#eta}{All Photons from #pi^{0}/#eta}";
    } else {
        yTitle = "#frac{Isolated Prompt Photons}{All Prompt Photons}";
    }
    hFrame->GetYaxis()->SetTitle(yTitle.c_str());
    hFrame->GetYaxis()->SetRangeUser(0, 2.0);

    // X-axis labeling
    hFrame->GetXaxis()->SetTitle("Leading Cluster p_{T} [GeV]");
    hFrame->GetXaxis()->SetLabelOffset(999); // hide numeric labels
    hFrame->GetXaxis()->SetTickLength(0);

    hFrame->Draw("AXIS");

    //
    // 5) Prepare a TLegend
    //
    TLegend legend(0.45, 0.65, 0.88, 0.85);
    legend.SetBorderSize(0);
    legend.SetTextSize(0.025);

    //
    // 6) We'll collect (pTcenter, ratio, error) in arrays for TGraph
    //
    std::vector<double> ptCenters;
    std::vector<double> ratios;
    std::vector<double> errors;
    ptCenters.reserve(128);
    ratios.reserve(128);
    errors.reserve(128);

    //
    // 7) Loop over each isoEtRange => gather ratio points if they're not excluded
    //
    for (const auto& kv : dataSubset) {
        const auto& isoEtRange = kv.first;   // (isoMin, isoMax)
        const auto& isoDataVec = kv.second;  // vector<IsolationDataWithPt>

        // If isoEtRange is in exclusionRanges => skip
        if (std::find(exclusionRanges.begin(), exclusionRanges.end(), isoEtRange)
            != exclusionRanges.end())
        {
            std::cout << "[DEBUG] isoEtRange={" << isoEtRange.first << ","
                      << isoEtRange.second << "} is in exclusionRanges => skip.\n";
            continue;
        }

        if (isoDataVec.empty()) {
            std::cout << "[DEBUG] isoDataVec empty for isoEtRange => skip.\n";
            continue;
        }

        // For each IsolationDataWithPt => find a pT bin center
        for (const auto& isoDataWPt : isoDataVec) {
            double ptMin = isoDataWPt.ptMin;
            double ptMax = isoDataWPt.ptMax;

            // Find which pT bin (ptMin, ptMax) belongs to in pT_bins
            bool foundBin = false;
            double ptCenter = 0.0;
            for (const auto& bin : pT_bins) {
                if (std::fabs(bin.first - ptMin) < 1e-6 &&
                    std::fabs(bin.second - ptMax) < 1e-6)
                {
                    ptCenter = 0.5 * (bin.first + bin.second);
                    foundBin = true;
                    break;
                }
            }
            if (!foundBin) {
                std::cout << "[DEBUG] pT bin not found => ("
                          << ptMin << "," << ptMax << ") => skip.\n";
                continue;
            }

            if (ptCenter >= pTExclusionMax) {
                std::cout << "[DEBUG] pTcenter=" << ptCenter
                          << " >= pTExclusionMax=" << pTExclusionMax
                          << " => skip.\n";
                continue;
            }

            // Push ratio & error from isoDataWPt
            ptCenters.push_back(ptCenter);
            ratios.push_back(isoDataWPt.ratio);
            errors.push_back(isoDataWPt.error);
        }
    }

    // If no points => throw an error
    if (ptCenters.empty()) {
        std::ostringstream oss;
        oss << "[CRITICAL ERROR] No valid ratio points found for subset='"
            << showerCutLabel << "'. Check logic.\n";
        std::cerr << oss.str();
        throw std::runtime_error(oss.str());
    }

    //
    // 8) Build a TGraphErrors for ratio vs. pT
    //
    TGraphErrors* mainGraph = new TGraphErrors(
        static_cast<int>(ptCenters.size()),
        ptCenters.data(),
        ratios.data(),
        nullptr,
        errors.data()
    );
    // Stylize
    mainGraph->SetMarkerStyle(20);
    mainGraph->SetMarkerColor(kBlack);
    mainGraph->SetLineColor(kBlack);
    mainGraph->SetLineWidth(2);
    mainGraph->Draw("P SAME");

    // Add a legend entry for the subset
    {
        std::ostringstream label;
        label << "#bf{Run24:} " << readableTriggerName << " [" << showerCutLabel << "]";
        legend.AddEntry(mainGraph, label.str().c_str(), "p");
    }

    //
    // 9) Draw a dashed horizontal line at y=1
    //
    TLine* unityLine = new TLine(binEdges.front(), 1.0, binEdges.back(), 1.0);
    unityLine->SetLineStyle(2);
    unityLine->Draw("SAME");

    //
    // 10) Optionally overlay reference data sets
    //
    TGraphErrors* refGraphA = nullptr;
    TGraphErrors* refGraphB = nullptr;

    // (A) Reference A (PHENIX direct photon ratio, for example)
    if (drawRefA) {
        std::vector<double> refX, refY, refErr;
        refX.reserve(referencePTGamma.size());
        refY.reserve(referencePTGamma.size());
        refErr.reserve(referencePTGamma.size());

        for (size_t i=0; i<referencePTGamma.size(); i++) {
            double ptVal = referencePTGamma[i];
            if (ptVal >= pTExclusionMax) continue;

            // find which bin this ptVal belongs to
            bool foundRefBin = false;
            double binCenter = 0.0;
            for (auto& b : pT_bins) {
                if (ptVal >= b.first && ptVal < b.second) {
                    binCenter = 0.5*(b.first + b.second);
                    foundRefBin = true;
                    break;
                }
            }
            if (!foundRefBin) continue;

            // shift horizontally a bit for clarity
            binCenter += 0.25;

            refX.push_back(binCenter);
            refY.push_back(referenceRatio[i]);
            refErr.push_back(referenceStatError[i]);
        }

        if (!refX.empty()) {
            refGraphA = new TGraphErrors(
                static_cast<int>(refX.size()),
                refX.data(),
                refY.data(),
                nullptr,
                refErr.data()
            );
            refGraphA->SetMarkerStyle(24);
            refGraphA->SetMarkerColor(kRed);
            refGraphA->SetLineColor(kRed);
            refGraphA->SetLineWidth(2);
            refGraphA->Draw("P SAME");
            legend.AddEntry(refGraphA, "#font[62]{PHENIX 2003 pp:} Isolated Direct / All Direct", "p");
        }
    }

    // (B) Reference B (PHENIX pi0 ratio, for example)
    if (drawRefB) {
        std::vector<double> ref2X, ref2Y, ref2Err;
        ref2X.reserve(referenceTwoPTGamma.size());
        ref2Y.reserve(referenceTwoPTGamma.size());
        ref2Err.reserve(referenceTwoPTGamma.size());

        for (size_t i=0; i<referenceTwoPTGamma.size(); i++) {
            double ptVal = referenceTwoPTGamma[i];
            if (ptVal >= pTExclusionMax) continue;

            bool foundRefBin = false;
            double binCenter = 0.0;
            for (auto& b : pT_bins) {
                if (ptVal >= b.first && ptVal < b.second) {
                    binCenter = 0.5*(b.first + b.second);
                    foundRefBin = true;
                    break;
                }
            }
            if (!foundRefBin) continue;

            // shift left for clarity
            binCenter -= 0.1;

            ref2X.push_back(binCenter);
            ref2Y.push_back(referenceTwoRatio[i]);
            ref2Err.push_back(referenceTwoStatError[i]);
        }

        if (!ref2X.empty()) {
            refGraphB = new TGraphErrors(
                static_cast<int>(ref2X.size()),
                ref2X.data(),
                ref2Y.data(),
                nullptr,
                ref2Err.data()
            );
            refGraphB->SetMarkerStyle(25);
            refGraphB->SetMarkerColor(kBlue);
            refGraphB->SetLineColor(kBlue);
            refGraphB->SetLineWidth(2);
            refGraphB->Draw("P SAME");
            legend.AddEntry(refGraphB, "#font[62]{PHENIX 2003 pp:} Isolated #pi^{0} / All #pi^{0}", "p");
        }
    }

    //
    // 11) Draw the legend
    //
    legend.Draw();

    //
    // 12) Manually add X-axis ticks & labels
    //
    double xMin = binEdges.front();
    double xMax = binEdges.back();
    double yMin = hFrame->GetMinimum();
    double yMax = hFrame->GetMaximum();

    double tickSize   = 0.02*(yMax - yMin);
    double labelOffset= 0.05*(yMax - yMin);

    // Draw x-axis line
    TLine xAxisLine(xMin, yMin, xMax, yMin);
    xAxisLine.Draw("SAME");

    // Label each bin edge
    TLatex latex;
    latex.SetTextSize(0.035);
    latex.SetTextAlign(22);

    for (size_t i=0; i<binEdges.size(); i++) {
        double xPos = binEdges[i];
        double yPos = yMin;
        TLine* tLine = new TLine(xPos, yPos, xPos, yPos - tickSize);
        tLine->Draw("SAME");

        std::ostringstream lb;
        lb << std::fixed << std::setprecision(1) << binEdges[i];
        latex.DrawLatex(xPos, yPos - labelOffset, lb.str().c_str());
    }

    canvas.RedrawAxis();

    //
    // 13) Add textual labels in the top-left corner
    //
    TLatex labelText;
    labelText.SetNDC();
    labelText.SetTextSize(0.025);
    labelText.SetTextColor(kBlack);

    double xL = 0.20;
    double yL = 0.90;
    double yStep = 0.045;

    {
        std::ostringstream oss;
        oss << "#font[62]{Trigger Group:} " << readableTriggerGroupName;
        labelText.DrawLatex(xL, yL, oss.str().c_str());
    }
    {
        std::ostringstream oss;
        oss << "#font[62]{Trigger:} " << readableTriggerName
            << " [" << showerCutLabel << "]";
        labelText.DrawLatex(xL, yL - yStep, oss.str().c_str());
    }
    {
        std::ostringstream oss;
        oss << "#font[62]{ECore #geq} " << eCore << " GeV";
        labelText.DrawLatex(xL, yL - 2*yStep, oss.str().c_str());
    }
    {
        std::ostringstream oss;
        oss << "#font[62]{#chi^{2} <} " << chi;
        labelText.DrawLatex(xL, yL - 3*yStep, oss.str().c_str());
    }
    {
        std::ostringstream oss;
        oss << "#font[62]{Asymmetry <} " << asym;
        labelText.DrawLatex(xL, yL - 4*yStep, oss.str().c_str());
    }
    {
        std::ostringstream oss;
        oss << "#font[62]{Mass Window:} " << massWindowLabel;
        labelText.DrawLatex(xL, yL - 5*yStep, oss.str().c_str());
    }
    {
        std::ostringstream oss;
        oss << "#font[62]{#Delta R_{cone} <} 0.3";
        labelText.DrawLatex(xL, yL - 6*yStep, oss.str().c_str());
    }
    {
        std::ostringstream oss;
        oss << "#font[62]{E_{T, iso} <} 6 GeV";
        labelText.DrawLatex(xL, yL - 7*yStep, oss.str().c_str());
    }

    canvas.Modified();
    canvas.Update();

    //
    // 14) Save the resulting plot
    //
    std::ostringstream outname;
    outname << subfolder
            << "/IsolationRatio_vs_pT_"
            << triggerName
            << "_" << showerCutLabel << ".png";
    std::string outPath = outname.str();
    canvas.SaveAs(outPath.c_str());
    std::cout << "[INFO] Saved ratio plot => " << outPath << "\n";

    //
    // 15) Cleanup
    //
    delete unityLine;
    delete hFrame;
    if (refGraphA) delete refGraphA;
    if (refGraphB) delete refGraphB;
    // mainGraph is owned by ROOT, so no manual delete
}


/**
 * \brief GeneratePerTriggerIsoPlots
 *
 * This function iterates a “grouped” data structure:
 *   groupedData[ GroupKey ] => ( isoEtRange => vector<IsolationDataWithPt> )
 * where `GroupKey` is typically:
 *   ( TriggerGroupName, TriggerName, Ecore, chi, Asym, massWindowLabel ).
 *
 * We then partition each group’s isoEt data by "withShowerShapeCuts" vs. "withoutShowerShapeCuts",
 * storing them into two separate subset maps. Finally, we call `producePlotForSubset()`
 * on each subset to produce the ratio (Isolated/All) vs. pT plots.
 *
 * The final directory structure is something like:
 *    <basePlotDir>/<TriggerGroupName>/E[..]_Chi[..]_Asym[..]/isolationEnergies/<massWindowLabel>
 *      +-> withShowerShapeCuts/IsolationRatio_vs_pT_...png
 *      +-> withoutShowerShapeCuts/IsolationRatio_vs_pT_...png
 *
 * \param groupedData   The aggregated isoEt data by group
 * \param basePlotDirectory  Where to store the outputs
 * \param isoEtRanges   (Optional) known isoEtRanges (unused here, except for color references)
 * \param isoEtColors   (Optional) a matching color scheme for isoEtRanges
 * \param referencePTGamma, referenceRatio, referenceStatError   => For optional overlay (A)
 * \param referenceTwoPTGamma, referenceTwoRatio, referenceTwoStatError => For optional overlay (B)
 * \param triggerCombinationNameMap For converting raw triggerGroupName => human-readable
 * \param triggerNameMap            For converting raw triggerName => human-readable
 * \param drawRefA, drawRefB   Whether to overlay the reference data sets
 * \param exclusionRanges      If an isoEt range is in this set, skip it
 * \param pT_bins             pT bin edges
 * \param pTExclusionMax      skip points >= this pT
 */
void GeneratePerTriggerIsoPlots(
    const std::map<GroupKey,
                   std::map<std::pair<float, float>,
                            std::vector<DataStructures::IsolationDataWithPt>>>& groupedData,
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
    const std::vector<std::pair<float, float>>& exclusionRanges,
    const std::vector<std::pair<double, double>>& pT_bins,
    double pTExclusionMax
)
{
    std::cout << "[DEBUG] Entering GeneratePerTriggerIsoPlots()...\n"
              << "        groupedData.size() = " << groupedData.size()
              << "\n\n";

    int groupCounter = 0;

    // For each (GroupKey => isoEtData)
    for (const auto& groupEntry : groupedData) {
        groupCounter++;

        // The groupKey typically is (triggerGroupName, triggerName, eCore, chi, asym, massWindow)
        const auto& groupKey  = groupEntry.first;
        const auto& isoEtData = groupEntry.second;  // isoEtRange => vector<IsolationDataWithPt>

        // Unpack groupKey
        std::string triggerGroupName = std::get<0>(groupKey);
        std::string triggerName      = std::get<1>(groupKey);
        float eCore                  = std::get<2>(groupKey);
        float chi                    = std::get<3>(groupKey);
        float asym                   = std::get<4>(groupKey);
        std::string massWindowLabel  = std::get<5>(groupKey);

        // Convert to readable strings
        std::string readableGroup = triggerGroupName;
        if (auto it = triggerCombinationNameMap.find(triggerGroupName); it != triggerCombinationNameMap.end()) {
            readableGroup = it->second;
        }
        std::string readableTrig = triggerName;
        if (auto it2 = triggerNameMap.find(triggerName); it2 != triggerNameMap.end()) {
            readableTrig = it2->second;
        }

        // Debug info
        std::cout << "[DEBUG] ---------------------------------------------------------\n";
        std::cout << "[DEBUG] Processing group # " << groupCounter << " / " << groupedData.size()
                  << ":\n"
                  << "        TriggerGroupName = '" << triggerGroupName << "'\n"
                  << "        ReadableGroupName= '" << readableGroup << "'\n"
                  << "        TriggerName      = '" << triggerName << "'\n"
                  << "        ReadableTrigName = '" << readableTrig << "'\n"
                  << "        ECore=" << eCore << ", Chi=" << chi << ", Asym=" << asym << "\n"
                  << "        MassWindowLabel='" << massWindowLabel << "'\n"
                  << "        isoEtData.size()= " << isoEtData.size() << "\n\n";

        if (isoEtData.empty()) {
            std::cerr << "[WARNING] Group # " << groupCounter
                      << " => isoEtData is empty => skipping.\n\n";
            continue;
        }

        // Build directories
        std::ostringstream dirPathStream;
        dirPathStream << basePlotDirectory << "/"
                      << triggerGroupName
                      << "/E"   << Utils::formatToThreeSigFigs(eCore)
                      << "_Chi" << Utils::formatToThreeSigFigs(chi)
                      << "_Asym"<< Utils::formatToThreeSigFigs(asym);

        std::string dirPath = dirPathStream.str();
        gSystem->mkdir(dirPath.c_str(), true);

        std::string isoDir = dirPath + "/isolationEnergies";
        gSystem->mkdir(isoDir.c_str(), true);

        std::string massWinDir = isoDir + "/" + massWindowLabel;
        gSystem->mkdir(massWinDir.c_str(), true);

        // We'll produce subfolders for "withShowerShapeCuts" vs. "withoutShowerShapeCuts"
        std::string withShowerFolder    = massWinDir + "/withShowerShapeCuts";
        std::string withoutShowerFolder = massWinDir + "/withoutShowerShapeCuts";
        gSystem->mkdir(withShowerFolder.c_str(), true);
        gSystem->mkdir(withoutShowerFolder.c_str(), true);

        // Next, partition isoEtData by showerCutLabel
        // We'll define:
        //   withShowerMap[ isoEtRange ] => vector<IsolationDataWithPt>
        //   withoutShowerMap[ isoEtRange ] => vector<IsolationDataWithPt>
        std::map<std::pair<float, float>, std::vector<DataStructures::IsolationDataWithPt>> withShowerMap;
        std::map<std::pair<float, float>, std::vector<DataStructures::IsolationDataWithPt>> withoutShowerMap;

        // Loop over each isoEtRange
        for (const auto& isoRangeEntry : isoEtData) {
            const auto& isoRange   = isoRangeEntry.first;   // (isoMin, isoMax)
            const auto& isoDataVec = isoRangeEntry.second;  // vector<IsolationDataWithPt>

            // Split them
            std::vector<DataStructures::IsolationDataWithPt> tmpWith, tmpWithout;
            tmpWith.reserve(isoDataVec.size());
            tmpWithout.reserve(isoDataVec.size());

            for (const auto& isoDataWPt : isoDataVec) {
                const std::string& cutLabel = isoDataWPt.isoData.showerCutLabel;
                if (cutLabel == "withShowerShapeCuts") {
                    tmpWith.push_back(isoDataWPt);
                } else if (cutLabel == "withoutShowerShapeCuts") {
                    tmpWithout.push_back(isoDataWPt);
                } else {
                    // Unknown showerCut => skip or log
                    std::cerr << "[WARNING] Unknown showerCutLabel='"
                              << cutLabel << "' => skipping.\n";
                }
            }

            // Insert into final maps if not empty
            if (!tmpWith.empty()) {
                withShowerMap[ isoRange ] = tmpWith;
            }
            if (!tmpWithout.empty()) {
                withoutShowerMap[ isoRange ] = tmpWithout;
            }
        }

        // Now produce ratio plots for each subset
        //  (1) withShowerShapeCuts
        if (!withShowerMap.empty()) {
            producePlotForSubset(
                withShowerMap,
                withShowerFolder,
                "withShowerShapeCuts",
                triggerName,
                readableTrig,
                readableGroup,
                eCore,
                chi,
                asym,
                massWindowLabel,
                exclusionRanges,
                pT_bins,
                pTExclusionMax,
                drawRefA,
                drawRefB,
                referencePTGamma,
                referenceRatio,
                referenceStatError,
                referenceTwoPTGamma,
                referenceTwoRatio,
                referenceTwoStatError
            );
        } else {
            std::cerr << "[WARNING] 'withShowerShapeCuts' map is empty => skipping.\n";
        }

        // (2) withoutShowerShapeCuts
        if (!withoutShowerMap.empty()) {
            producePlotForSubset(
                withoutShowerMap,
                withoutShowerFolder,
                "withoutShowerShapeCuts",
                triggerName,
                readableTrig,
                readableGroup,
                eCore,
                chi,
                asym,
                massWindowLabel,
                exclusionRanges,
                pT_bins,
                pTExclusionMax,
                drawRefA,
                drawRefB,
                referencePTGamma,
                referenceRatio,
                referenceStatError,
                referenceTwoPTGamma,
                referenceTwoRatio,
                referenceTwoStatError
            );
        } else {
            std::cerr << "[WARNING] 'withoutShowerShapeCuts' map is empty => skipping.\n";
        }

        std::cout << "[DEBUG] Finished group # " << groupCounter << ".\n\n";
    }

    std::cout << "[DEBUG] Finished GeneratePerTriggerIsoPlots() for all groups.\n\n";
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
            labelText.SetTextSize(0.024);       // Adjust text size as needed

            TLatex valueText;
            valueText.SetNDC();
            valueText.SetTextSize(0.024);

            double xStart = 0.2; // Starting x-coordinate (left side)
            double yStartLabel = 0.9; // Starting y-coordinate
            double yStepLabel = 0.04;  // Vertical spacing between lines

            // Prepare label strings
            labelText.DrawLatex(xStart, yStartLabel, "#font[62]{Active Trigger Group:}");
            valueText.DrawLatex(xStart + 0.2, yStartLabel, readableTriggerGroupName.c_str());

            labelText.DrawLatex(xStart, yStartLabel - yStepLabel, "#font[62]{ECore #geq}");
            std::ostringstream eCoreWithUnit;
            eCoreWithUnit << eCore << "   GeV";
            valueText.DrawLatex(xStart + 0.15, yStartLabel - yStepLabel, eCoreWithUnit.str().c_str());

            labelText.DrawLatex(xStart, yStartLabel - 2 * yStepLabel, "#font[62]{#chi^{2} <}");
            std::ostringstream chiStr;
            chiStr << chi;
            valueText.DrawLatex(xStart + 0.15, yStartLabel - 2 * yStepLabel, chiStr.str().c_str());

            labelText.DrawLatex(xStart, yStartLabel - 3 * yStepLabel, "#font[62]{Asymmetry <}");
            std::ostringstream asymmetryStr;
            asymmetryStr << asym;
            valueText.DrawLatex(xStart + 0.15, yStartLabel - 3 * yStepLabel, asymmetryStr.str().c_str());

            labelText.DrawLatex(xStart, yStartLabel - 4 * yStepLabel, "#font[62]{Mass Window:}");
            valueText.DrawLatex(xStart + 0.15, yStartLabel - 4 * yStepLabel, massWindowLabel.c_str());

            labelText.DrawLatex(xStart, yStartLabel - 5 * yStepLabel, "#font[62]{#Delta R_{cone} <}");
            valueText.DrawLatex(xStart + 0.15, yStartLabel - 5 * yStepLabel, "0.3");
            
            labelText.DrawLatex(xStart, yStartLabel - 6 * yStepLabel, "#font[62]{E_{T, iso} <}");
            valueText.DrawLatex(xStart + 0.15, yStartLabel - 6 * yStepLabel, "6 GeV");

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


// -----------------------------------------------------------------------------
//   1) PrepareDataForIsolationPurity
// -----------------------------------------------------------------------------
void PrepareDataForIsolationPurity(
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

    // We now group by an 8-element key that ends with showerCutLabel:
    std::map<std::tuple<
        std::string, // TriggerGroupName
        std::string, // TriggerName
        float,       // ECore
        float,       // Chi
        float,       // Asym
        float,       // isoMin
        float,       // isoMax
        std::string  // showerCutLabel
    >,
    std::vector<DataStructures::IsolationDataWithPt>>& groupedData
)
{
    for (const auto& [key, isoData] : dataMap)
    {
        // Unpack the old key
        const std::string& triggerGroupName = std::get<0>(key);
        const std::string& triggerName      = std::get<1>(key);
        float eCore                         = std::get<2>(key);
        float chi                           = std::get<3>(key);
        float asym                          = std::get<4>(key);
        float ptMin                         = std::get<5>(key);
        float ptMax                         = std::get<6>(key);
        float isoMin                        = std::get<7>(key);
        float isoMax                        = std::get<8>(key);
        // We ignore massWindowLabel here, because we want to combine them.

        // The actual showerCutLabel comes from the IsolationData itself:
        const std::string& showerCutLabel = isoData.showerCutLabel;

        // Build the new grouping key with 8 fields:
        auto newGroupKey = std::make_tuple(
            triggerGroupName,
            triggerName,
            eCore,
            chi,
            asym,
            isoMin,
            isoMax,
            showerCutLabel // new
        );

        // Prepare the data-with-pt object
        DataStructures::IsolationDataWithPt isoDataWithPt;
        isoDataWithPt.ptMin   = ptMin;
        isoDataWithPt.ptMax   = ptMax;
        isoDataWithPt.isoData = isoData;

        groupedData[newGroupKey].push_back(isoDataWithPt);
    }
}


void GenerateIsolationPurityPlots(
    // Notice the group key now has *eight* elements, with showerCutLabel at index 7:
    const std::map<std::tuple<
        std::string, // TriggerGroupName
        std::string, // TriggerName
        float,       // ECore
        float,       // Chi
        float,       // Asymmetry
        float,       // isoMin
        float,       // isoMax
        std::string  // showerCutLabel => "withShowerShapeCuts" / "withoutShowerShapeCuts"
    >,
    std::vector<DataStructures::IsolationDataWithPt>>& groupedData,
    const std::string& basePlotDirectory,
    const std::map<std::string, std::string>& triggerCombinationNameMap,
    const std::map<std::string, std::string>& triggerNameMap,
    const std::vector<std::pair<double, double>>& pT_bins,
    double pTExclusionMax
) {
    std::cout << "[INFO] Generating Isolation Purity Plots (split by showerCutLabel)...\n";

    int plotCounter = 0;

    // -------------------------------------------------------------------------
    // Iterate over each group
    // -------------------------------------------------------------------------
    for (const auto& [groupKey, isoDataList] : groupedData)
    {
        plotCounter++;

        // Unpack the 8-element key
        std::string triggerGroupName  = std::get<0>(groupKey);
        std::string triggerName       = std::get<1>(groupKey);
        float eCore                   = std::get<2>(groupKey);
        float chi                     = std::get<3>(groupKey);
        float asym                    = std::get<4>(groupKey);
        float isoMin                  = std::get<5>(groupKey);
        float isoMax                  = std::get<6>(groupKey);
        std::string showerCutLabel    = std::get<7>(groupKey);

        // Convert to more readable strings (if available)
        std::string readableTriggerGroupName =
            Utils::getTriggerCombinationName(triggerGroupName, triggerCombinationNameMap);

        std::string readableTriggerName = triggerName;
        auto it_nameMap = triggerNameMap.find(triggerName);
        if (it_nameMap != triggerNameMap.end()) {
            readableTriggerName = it_nameMap->second;
        }

        std::cout << "\n[DEBUG] Processing Purity Plot #" << plotCounter << ":\n"
                  << "         GroupName='"       << readableTriggerGroupName << "'\n"
                  << "         Trigger='"          << readableTriggerName      << "'\n"
                  << "         ECore="             << eCore
                  << ", Chi="                    << chi
                  << ", Asymmetry="             << asym << "\n"
                  << "         isoMin="            << isoMin
                  << ", isoMax="                << isoMax << "\n"
                  << "         showerCutLabel='"   << showerCutLabel << "'\n"
                  << "         #IsoDataPoints="    << isoDataList.size() << "\n";

        // ---------------------------------------------------------------------
        // Create the base directory:
        //   <basePlotDirectory>/<triggerGroupName>/E..._Chi..._Asym.../isolationEnergies
        // Then add a subdirectory for the showerCutLabel
        // ---------------------------------------------------------------------
        std::ostringstream dirStream;
        dirStream << basePlotDirectory << "/" << triggerGroupName
                  << "/E"   << Utils::formatToThreeSigFigs(eCore)
                  << "_Chi" << Utils::formatToThreeSigFigs(chi)
                  << "_Asym"<< Utils::formatToThreeSigFigs(asym)
                  << "/isolationEnergies";
        std::string mainDirPath = dirStream.str();
        gSystem->mkdir(mainDirPath.c_str(), true);

        // Next level: withShowerShapeCuts / withoutShowerShapeCuts folder
        std::string showerCutDir = mainDirPath + "/" + showerCutLabel;
        gSystem->mkdir(showerCutDir.c_str(), true);

        // ---------------------------------------------------------------------
        // Prepare a TCanvas
        // ---------------------------------------------------------------------
        TCanvas canvas("canvas","Isolation Purity",800,600);

        // ---------------------------------------------------------------------
        // Prepare pT bin edges
        // ---------------------------------------------------------------------
        std::vector<double> binEdges;
        for (const auto& bin : pT_bins) {
            if (bin.first >= pTExclusionMax) {
                break;
            }
            binEdges.push_back(bin.first);
        }
        // Add the upper edge if possible
        if (!binEdges.empty()) {
            double lastSecond = pT_bins[binEdges.size()-1].second;
            if (lastSecond < pTExclusionMax) {
                binEdges.push_back(lastSecond);
            } else {
                binEdges.push_back(pTExclusionMax);
            }
        } else {
            std::cerr << "[WARNING] No valid pT bins => skipping.\n";
            continue;
        }

        int nBins = static_cast<int>(binEdges.size()) - 1;
        if (nBins <= 0) {
            std::cerr << "[ERROR] Not enough bins => skipping.\n";
            continue;
        }
        double* binArray = &binEdges[0];

        // ---------------------------------------------------------------------
        // Create dummy hist for axis
        // ---------------------------------------------------------------------
        TH1F* hFrame = new TH1F("hFrame","",nBins, binArray);
        hFrame->SetStats(0);
        hFrame->GetXaxis()->SetTitle("Leading Cluster p_{T} [GeV]");
        hFrame->GetYaxis()->SetTitle("#frac{N_{isolated prompt photons}}{N_{all isolated photons}}");
        hFrame->GetYaxis()->SetRangeUser(0,1.2);

        // Hide default numeric x-axis
        hFrame->GetXaxis()->SetLabelOffset(999);
        hFrame->GetXaxis()->SetTickLength(0);
        hFrame->Draw("AXIS");

        // ---------------------------------------------------------------------
        // A TLegend
        // ---------------------------------------------------------------------
        TLegend legend(0.2,0.75,0.4,0.85);
        legend.SetBorderSize(0);
        legend.SetTextSize(0.026);

        // ---------------------------------------------------------------------
        // We'll collect purity for each pT bin by counting:
        //   - outsideMassWindow => numerator
        //   - total => denominator
        // and skip zero denominators
        // ---------------------------------------------------------------------
        std::map<double, std::pair<int,int>> ptCounts;
        // Key= pT center, Value= <outsideMassCounts, totalIsolatedCounts>

        // Loop over isoDataList
        for (auto& isoDataWithPt : isoDataList) {
            double pTmin = isoDataWithPt.ptMin;
            double pTmax = isoDataWithPt.ptMax;
            double pTCenter = 0.5*(pTmin+pTmax);

            // skip if > pTExclusionMax
            if (pTCenter >= pTExclusionMax) {
                continue;
            }

            // retrieve massWindow label from isoData:
            const auto& isoData = isoDataWithPt.isoData;
            // increment total
            ptCounts[pTCenter].second += isoData.isolatedCounts;

            // if outsideMassWindow => numerator
            if (isoData.massWindowLabel == "outsideMassWindow") {
                ptCounts[pTCenter].first += isoData.isolatedCounts;
            }
        }

        // Now compute purity & fill vectors
        std::vector<double> vCenters;
        std::vector<double> vPurities;
        std::vector<double> vErrors;

        for (auto& kv : ptCounts) {
            double pTCenter = kv.first;
            int outsideCt   = kv.second.first;
            int totalCt     = kv.second.second;

            if (totalCt == 0) {
                continue;
            }

            double purity = static_cast<double>(outsideCt) / totalCt;
            double err    = std::sqrt( purity*(1.0 - purity) / static_cast<double>(totalCt) );

            vCenters.push_back(pTCenter);
            vPurities.push_back(purity);
            vErrors.push_back(err);

            std::cout << "[DEBUG] pTCenter=" << pTCenter
                      << ", outside=" << outsideCt
                      << ", total=" << totalCt
                      << " => purity=" << purity
                      << " +/- " << err << "\n";
        }

        if (vCenters.empty()) {
            std::cerr << "[WARNING] No valid data => skipping.\n";
            delete hFrame;
            continue;
        }

        // ---------------------------------------------------------------------
        // Build TGraphErrors
        // ---------------------------------------------------------------------
        TGraphErrors* gPurity = new TGraphErrors(
            static_cast<int>(vCenters.size()),
            &vCenters[0],
            &vPurities[0],
            nullptr,
            &vErrors[0]
        );

        // figure out color from the triggerName (optional)
        int markerColor = kBlack;
        auto itColor = TriggerConfig::triggerColorMap.find(triggerName);
        if (itColor != TriggerConfig::triggerColorMap.end()) {
            markerColor = itColor->second;
        }
        gPurity->SetMarkerStyle(20);
        gPurity->SetMarkerColor(markerColor);
        gPurity->SetLineColor(markerColor);
        gPurity->SetLineWidth(2);
        gPurity->Draw("P SAME");

        legend.AddEntry(gPurity, readableTriggerName.c_str(), "p");
        legend.Draw();

        // ---------------------------------------------------------------------
        // Add a dashed line at y=1
        // ---------------------------------------------------------------------
        double xMin = binEdges.front();
        double xMax = binEdges.back();
        TLine lineY1(xMin, 1.0, xMax, 1.0);
        lineY1.SetLineStyle(2);
        lineY1.Draw("SAME");

        // ---------------------------------------------------------------------
        // Manual x-axis ticks
        // ---------------------------------------------------------------------
        double yMin = hFrame->GetMinimum();
        double yMax = hFrame->GetMaximum();
        double tickSize    = (yMax - yMin)*0.02;
        double labelOffset = (yMax - yMin)*0.05;

        TLatex latex;
        latex.SetTextSize(0.035);
        latex.SetTextAlign(22);

        TLine xAxisLine(xMin, yMin, xMax, yMin);
        xAxisLine.Draw("SAME");

        for (size_t i = 0; i < binEdges.size(); i++) {
            double bx = binEdges[i];
            TLine* tick = new TLine(bx, yMin, bx, yMin - tickSize);
            tick->Draw("SAME");

            std::ostringstream lb;
            lb << std::fixed << std::setprecision(1) << bx;
            latex.DrawLatex(bx, yMin - labelOffset, lb.str().c_str());
        }

        // ---------------------------------------------------------------------
        // Add textual labels (left side)
        // ---------------------------------------------------------------------
        canvas.RedrawAxis();

        TLatex lbl;
        lbl.SetNDC();
        lbl.SetTextSize(0.024);
        lbl.SetTextColor(kBlack);

        double xLbl = 0.18;
        double yLbl = 0.45;
        double dy   = 0.045;

        {
            std::ostringstream oss;
            oss << "#font[62]{TriggerGroup:} " << readableTriggerGroupName;
            lbl.DrawLatex(xLbl, yLbl, oss.str().c_str());
        }
        {
            std::ostringstream oss;
            oss << "#font[62]{Trigger:} " << readableTriggerName;
            lbl.DrawLatex(xLbl, yLbl - dy, oss.str().c_str());
        }
        {
            std::ostringstream oss;
            oss << "#font[62]{ECore #geq} " << eCore << " GeV";
            lbl.DrawLatex(xLbl, yLbl - 2*dy, oss.str().c_str());
        }
        {
            std::ostringstream oss;
            oss << "#font[62]{#chi^{2} <} " << chi;
            lbl.DrawLatex(xLbl, yLbl - 3*dy, oss.str().c_str());
        }
        {
            std::ostringstream oss;
            oss << "#font[62]{Asymmetry <} " << asym;
            lbl.DrawLatex(xLbl, yLbl - 4*dy, oss.str().c_str());
        }
        {
            std::ostringstream oss;
            oss << "#font[62]{#Delta R_{cone} <} 0.3";
            lbl.DrawLatex(xLbl, yLbl - 5*dy, oss.str().c_str());
        }
        {
            std::ostringstream oss;
            oss << "#font[62]{E_{T, iso} <} 6 GeV";
            lbl.DrawLatex(xLbl, yLbl - 6*dy, oss.str().c_str());
        }

        // ---------------------------------------------------------------------
        // Show the showerCutLabel
        // ---------------------------------------------------------------------
        {
            std::ostringstream oss;
            oss << "#font[62]{ShowerCut:} " << showerCutLabel;
            lbl.DrawLatex(xLbl, yLbl - 7*dy, oss.str().c_str());
        }

        // ---------------------------------------------------------------------
        // ** ADD EXTRA LINES if we have "withShowerShapeCuts" **
        // ---------------------------------------------------------------------
        if (showerCutLabel == "withShowerShapeCuts")
        {
            // Example of how to put them on 2 lines for neatness
            double yCuts = yLbl - 8*dy;  // Start a bit lower
            {
                std::ostringstream oss;
                // E_{3x7}/E_{7x7} < 0.9, E_{T}^{HCal, 3x3}/ E_{T, cluster} < 0.1
                oss << "#font[62]{E_{3x7}/E_{7x7} < 0.9,  E_{T}^{HCal,3x3}/E_{T, cluster} < 0.1}";
                lbl.DrawLatex(xLbl, yCuts, oss.str().c_str());
            }
            {
                std::ostringstream oss;
                // E_{1x1}/E_{cluster} < 0.98, <w_{72}> > 0.75
                oss << "#font[62]{E_{1x1}/E_{cluster} < 0.98,   <w_{72}> > 0.75}";
                lbl.DrawLatex(xLbl, yCuts - dy, oss.str().c_str());
            }
        }

        canvas.Modified();
        canvas.Update();

        // ---------------------------------------------------------------------
        // Save the plot
        // ---------------------------------------------------------------------
        std::ostringstream outName;
        outName << showerCutDir
                << "/IsolationPurity_vs_pT_"
                << triggerName << ".png";
        std::string outPath = outName.str();
        canvas.SaveAs(outPath.c_str());
        std::cout << "[INFO] Saved => " << outPath << "\n";

        // Cleanup
        delete hFrame;
        // The TCanvas destructor will clean up TLines, TGraph, etc.
    }

    std::cout << "[INFO] Done generating isolation purity plots, with separate subfolders.\n";
}

// Define CombinedPurityData structure
struct CombinedPurityData {
    float pTCenter;
    double purity;
    double purityError;
};

// Function to sort and combine purity data
void SortAndCombinePurityData(
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
    std::map<SpectraGroupKey, std::map<float, CombinedPurityData>>& combinedPurityDataMap
) {
    // Map to store counts per groupKey and pT bin
    std::map<SpectraGroupKey, std::map<std::pair<float, float>, std::pair<int, int>>> dataPerGroup;

    // Helper function to process data maps
    auto processDataMap = [&](const auto& dataMap, bool isOutsideMassWindow) {
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
            // MassWindowLabel is not used here

            // Exclude isoEtRanges if necessary
            std::pair<float, float> isoEtRange = {isoMin, isoMax};
            if (std::find(exclusionRanges.begin(), exclusionRanges.end(), isoEtRange) != exclusionRanges.end()) {
                continue; // Exclude this isoEtRange
            }

            SpectraGroupKey groupKey{triggerGroupName, eCore, chi, asymmetry};
            std::pair<float, float> pTBin{pTMin, pTMax};

            // Initialize counts if not already present
            auto& counts = dataPerGroup[groupKey][pTBin];
            if (isOutsideMassWindow) {
                counts.first += isoData.isolatedCounts; // outsideCounts
            }
            counts.second += isoData.isolatedCounts;    // totalIsolatedCounts
        }
    };

    // Process both inMassWindow and outsideMassWindow data
    processDataMap(dataMap_inMassWindow, false);
    processDataMap(dataMap_outsideMassWindow, true);

    // Now, for each group and pT bin, calculate purity
    for (const auto& groupEntry : dataPerGroup) {
        const SpectraGroupKey& groupKey = groupEntry.first;
        const auto& pTBinMap = groupEntry.second;

        for (const auto& pTBinEntry : pTBinMap) {
            const std::pair<float, float>& pTBin = pTBinEntry.first;
            float ptMin = pTBin.first;
            float ptMax = pTBin.second;
            float pTCenter = (ptMin + ptMax) / 2.0;

            // Apply pTExclusionMax
            if (pTCenter >= pTExclusionMax) {
                continue; // Exclude this pT bin
            }

            const auto& counts = pTBinEntry.second;

            int outsideCounts = counts.first;
            int totalIsolatedCounts = counts.second;

            if (totalIsolatedCounts == 0) {
                std::cerr << "[WARNING] Total isolated counts is zero for pT center: " << pTCenter << " GeV. Skipping purity calculation.\n";
                continue;
            }

            double purity = static_cast<double>(outsideCounts) / totalIsolatedCounts;
            double error = std::sqrt(purity * (1 - purity) / totalIsolatedCounts); // Binomial error

            // Prepare CombinedPurityData
            CombinedPurityData combinedData;
            combinedData.pTCenter = pTCenter;
            combinedData.purity = purity;
            combinedData.purityError = error;

            // Add to combinedPurityDataMap
            combinedPurityDataMap[groupKey][pTCenter] = combinedData;
        }
    }
}

void GenerateCombinedPurityPlots(
    const std::map<SpectraGroupKey, std::map<float, CombinedPurityData>>& combinedPurityDataMap,
    const std::string& basePlotDirectory,
    const std::map<std::string, std::string>& triggerCombinationNameMap,
    const std::vector<std::pair<double, double>>& pT_bins,
    double pTExclusionMax
) {
    std::cout << "[INFO] Starting GenerateCombinedPurityPlots function.\n";

    for (const auto& [spectraGroupKey, ptDataMap] : combinedPurityDataMap) {
        const std::string& triggerGroupName = spectraGroupKey.triggerGroupName;
        float eCore = spectraGroupKey.eCore;
        float chi = spectraGroupKey.chi;
        float asym = spectraGroupKey.asymmetry;

        // Map to human-readable names
        std::string readableTriggerGroupName = Utils::getTriggerCombinationName(
            triggerGroupName, triggerCombinationNameMap);

        std::cout << "[INFO] Processing plot for Trigger Group: " << readableTriggerGroupName
                  << ", ECore ≥ " << eCore << " GeV, Chi² < " << chi << ", Asymmetry < " << asym << ".\n";

        // Define output directory
        std::ostringstream dirStream;
        dirStream << basePlotDirectory << "/" << triggerGroupName
                  << "/E" << Utils::formatToThreeSigFigs(eCore)
                  << "_Chi" << Utils::formatToThreeSigFigs(chi)
                  << "_Asym" << Utils::formatToThreeSigFigs(asym)
                  << "/isolationEnergies";
        std::string dirPath = dirStream.str();
        gSystem->mkdir(dirPath.c_str(), true);

        // Create canvas
        TCanvas* canvas = new TCanvas("canvas", "Combined Isolation Purity", 800, 600);

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
            binEdges.push_back(std::min(pT_bins[binEdges.size() - 1].second, pTExclusionMax));
        } else {
            // No bins to plot
            std::cerr << "[WARNING] No pT bins to plot. Skipping plot.\n";
            continue;
        }

        int nBins = binEdges.size() - 1;
        double* binEdgesArray = binEdges.data();

        // Create a dummy histogram to set up the axes
        TH1F* hFrame = new TH1F("hFrame", "", nBins, binEdgesArray);
        hFrame->SetStats(0);
        hFrame->GetXaxis()->SetTitle("Leading Cluster p_{T} [GeV]");
        hFrame->GetYaxis()->SetTitle("#frac{N_{isolated prompt photons}}{N_{all isolated photons}}");
        hFrame->GetYaxis()->SetRangeUser(0, 1.2);

        // Remove x-axis labels and ticks
        hFrame->GetXaxis()->SetLabelOffset(999);
        hFrame->GetXaxis()->SetTickLength(0);

        // Draw the frame
        hFrame->Draw("AXIS");

        // Prepare vectors to collect data
        std::vector<double> ptCenters;
        std::vector<double> purities;
        std::vector<double> errors;

        for (const auto& [pTCenter, data] : ptDataMap) {
            // Exclude data points with pT ≥ pTExclusionMax
            if (pTCenter >= pTExclusionMax) {
                continue;
            }

            ptCenters.push_back(data.pTCenter);
            purities.push_back(data.purity);
            errors.push_back(data.purityError);
        }

        if (ptCenters.empty()) {
            std::cerr << "[WARNING] No valid data points for Trigger Group: " << readableTriggerGroupName << ". Skipping.\n";
            delete hFrame;
            delete canvas;
            continue;
        }

        // Create a TGraphErrors
        TGraphErrors* graph = new TGraphErrors(ptCenters.size(),
                                               ptCenters.data(),
                                               purities.data(),
                                               nullptr,
                                               errors.data());

        graph->SetMarkerStyle(20);
        graph->SetMarkerColor(kBlack);
        graph->SetLineColor(kBlack);
        graph->SetLineWidth(2);

        // Draw the graph
        graph->Draw("P SAME");

        // Draw custom x-axis ticks and labels
        double xMin = binEdges.front();
        double xMax = binEdges.back();
        double yAxisMin = hFrame->GetMinimum();
        double yAxisMax = hFrame->GetMaximum();

        TLine* yLine = new TLine(xMin, 1.0, xMax, 1.0);
        yLine->SetLineColor(kBlack);
        yLine->SetLineStyle(2); // Dashed
        yLine->SetLineWidth(2); // Optional: set line width for better visibility
        yLine->Draw("SAME");

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
        canvas->RedrawAxis();

        // Add labels using TLatex
        TLatex labelText;
        labelText.SetNDC();
        labelText.SetTextSize(0.023);
        labelText.SetTextColor(kBlack);

        double xStart = 0.2;
        double yStartLabel = 0.5;
        double yStepLabel = 0.025;

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
        std::string outputFilePath = dirPath + "/CombinedIsolationPurity.png";
        canvas->SaveAs(outputFilePath.c_str());
        std::cout << "[INFO] Saved combined isolation purity plot to " << outputFilePath << std::endl;

        // Clean up
        delete hFrame;
        delete canvas;
        delete graph; // Don't forget to delete the graph
        delete yLine;
    }

    std::cout << "[INFO] Finished GenerateCombinedPurityPlots function.\n";
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
    
    using ExtendedGroupKey = std::tuple<
        std::string, // TriggerGroupName
        std::string, // TriggerName
        float,       // ECore
        float,       // Chi
        float,       // Asymmetry
        std::string, // MassWindowLabel
        std::string  // ShowerCutLabel
    >;

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

    std::cout << "[INFO] Starting to process isolation data...\n";

    // -------------------------------------------------
    // 1) Populate the aggregator from the input dataMap
    // -------------------------------------------------
    for (const auto& entry : dataMap)
    {
        // key is the 10-tuple with ECore, chi, asym, pT, iso, massWindow
        const auto& key     = entry.first;
        const auto& isoData = entry.second; // includes showerCutLabel

        // Unpack
        std::string triggerGroupName = std::get<0>(key);
        std::string triggerName      = std::get<1>(key);
        float eCore                  = std::get<2>(key);
        float chi                    = std::get<3>(key);
        float asym                  = std::get<4>(key);
        float ptMin                 = std::get<5>(key);
        float ptMax                 = std::get<6>(key);
        float isoMin                = std::get<7>(key);
        float isoMax                = std::get<8>(key);
        std::string massWindowLabel = std::get<9>(key);

        // Possibly skip if isoEtRange is in 'exclusionRanges'
        std::pair<float,float> isoEtRange = {isoMin, isoMax};
        if (std::find(exclusionRanges.begin(), exclusionRanges.end(), isoEtRange)
            != exclusionRanges.end())
        {
            continue; // skip
        }

        // Retrieve the showerCutLabel from the isoData
        const std::string& showerCutLabel = isoData.showerCutLabel;

        GroupKey groupKey = std::make_tuple(
            triggerGroupName,
            triggerName,
            eCore,
            chi,
            asym,
            massWindowLabel
        );

        // Build your IsolationDataWithPt struct
        DataStructures::IsolationDataWithPt isoDataWithPt;
        isoDataWithPt.ptMin      = ptMin;
        isoDataWithPt.ptMax      = ptMax;
        isoDataWithPt.isoMin     = isoMin;
        isoDataWithPt.isoMax     = isoMax;
        isoDataWithPt.triggerName= triggerName;
        isoDataWithPt.weightedPt = isoData.weightedPt;
        isoDataWithPt.ratio      = isoData.ratio;
        isoDataWithPt.error      = isoData.error;
        isoDataWithPt.isoData    = isoData;   // holds isolatedCounts, totalCounts, etc.

        // Insert into aggregator
        groupedData[groupKey][isoEtRange].push_back(isoDataWithPt);
    }

    
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
         exclusionRanges,
         DataStructures::pT_bins,    // Added pT_bins
         20.0                        // Added pTExclusionMax
     );
    GeneratePerTriggerSpectraPlots(
        dataMap_inMassWindow,
        dataMap_outsideMassWindow,
        basePlotDirectory,
        triggerCombinationNameMap,
        triggerNameMap,
        exclusionRanges,
        DataStructures::pT_bins, // Ensure this is defined and accessible
        20.0                     // Adjust pTExclusionMax as needed
    );
    
    // Prepare data for isolation purity plots
    std::map<std::tuple<
        std::string, // TriggerGroupName
        std::string, // TriggerName
        float,       // ECore
        float,       // Chi
        float,       // Asymmetry
        float,       // isoMin
        float,       // isoMax
        std::string  // MassWindowLabel (not used)
    >, std::vector<DataStructures::IsolationDataWithPt>> groupedDataForPurity;

    // Combine data from both mass windows
    PrepareDataForIsolationPurity(dataMap_inMassWindow, groupedDataForPurity);
    PrepareDataForIsolationPurity(dataMap_outsideMassWindow, groupedDataForPurity);

    // Generate isolation purity plots
    GenerateIsolationPurityPlots(
        groupedDataForPurity,
        basePlotDirectory,
        triggerCombinationNameMap,
        triggerNameMap,
        DataStructures::pT_bins,
        20.0
    );

    // Sort and combine purity data
    std::map<SpectraGroupKey, std::map<float, CombinedPurityData>> combinedPurityDataMap;
    SortAndCombinePurityData(
        dataMap_inMassWindow,
        dataMap_outsideMassWindow,
        combinationToTriggerEfficiencyPoints,
        exclusionRanges,
        DataStructures::pT_bins,    // Pass pT_bins
        20.0,                        // Pass pTExclusionMax
        combinedPurityDataMap
    );

    // Generate combined purity plots
    GenerateCombinedPurityPlots(
        combinedPurityDataMap,
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
    const std::string& fitFunctionType = "sigmoid",
    bool fitOnly = true,
    const std::string& histogramType = "maxEcore")
{
    // --------------------------------------------------------------------------
    // 1) Determine histogram prefix based on histogramType (for photons/minbias)
    // --------------------------------------------------------------------------
    std::cout << "[DEBUG] PlotCombinedHistograms() starting...\n";
    std::cout << "       # of combinedRootFiles = " << combinedRootFiles.size() << "\n";
    std::cout << "       histogramType = " << histogramType << "\n";

    std::string histPrefix;
    std::string xAxisTitle;
    std::string xAxisTitleTurnOn;

    if (histogramType == "maxEcore")
    {
        histPrefix        = "h_maxEnergyClus_";
        xAxisTitle        = "Maximum Cluster Energy [GeV]";
        xAxisTitleTurnOn  = "Maximum Cluster Energy [GeV]";
    }
    else
    {
        // Default to 8x8 Tower Energy hist
        histPrefix        = "h8by8TowerEnergySum_";
        xAxisTitle        = "Maximum 8x8 EMCal Tower Energy Sum [GeV]";
        xAxisTitleTurnOn  = "Maximum 8x8 Energy Sum [GeV]";
    }

    // Jet triggers always use this, ignoring histogramType:
    const std::string jetHistPrefix = "h_jet_energy_";

    // --------------------------------------------------------------------------
    // 2) Remove trailing underscore from histPrefix for file naming
    // --------------------------------------------------------------------------
    std::string histPrefixForFile = histPrefix;
    if (!histPrefixForFile.empty() && histPrefixForFile.back() == '_') {
        histPrefixForFile.pop_back();
    }

    // --------------------------------------------------------------------------
    // 3) Pull in triggers, color maps, etc. from TriggerConfig
    // --------------------------------------------------------------------------
    const std::vector<std::string>& allTriggers    = TriggerConfig::allTriggers;
    const std::vector<std::string>& photonTriggers = TriggerConfig::photonTriggers;
    const auto& triggerColorMap   = TriggerConfig::triggerColorMap;
    const auto& triggerNameMap    = TriggerConfig::triggerNameMap;

    // For storing 99%-efficiency points from the fits
    std::map<std::string, double> triggerEfficiencyPoints;
    std::map<std::string, std::map<std::string, double>> combinationToTriggerEfficiencyPoints;

    // --------------------------------------------------------------------------
    // 4) Base directory for plots
    // --------------------------------------------------------------------------
    std::string basePlotDirectory = "/Users/patsfan753/Desktop/DirectPhotonAna/Plots";
    gSystem->mkdir(basePlotDirectory.c_str(), true);

    // --------------------------------------------------------------------------
    // 5) Function to draw run numbers on canvas
    // --------------------------------------------------------------------------
    auto drawRunNumbersOnCanvas = [](const std::vector<int>& runNumbers, const std::string& firmwareStatus)
    {
        std::cout << "[DEBUG] drawRunNumbersOnCanvas: # of runs = " << runNumbers.size() << "\n";

        TLatex runNumbersLatex;
        runNumbersLatex.SetNDC();
        runNumbersLatex.SetTextAlign(13); // top-left
        runNumbersLatex.SetTextColor(kBlack);

        double xStart, xEnd, yStart, yEnd, textSize;
        double xSpacingFactor, ySpacingFactor;
        int numColumns;

        // Example-based customizing
        if (runNumbers.size() == 692) {
            numColumns = 27;
            textSize   = 0.012;
            xStart     = 0.54; xEnd = 0.93;
            yStart     = 0.9;  yEnd = 0.45;
            xSpacingFactor = 0.8;
            ySpacingFactor = 0.8;
        }
        else if (runNumbers.size() == 71) {
            numColumns = 5;
            textSize   = 0.022;
            xStart     = 0.54; xEnd = 0.93;
            yStart     = 0.9;  yEnd = 0.45;
            xSpacingFactor = 0.8;
            ySpacingFactor = 0.8;
        }
        else if (runNumbers.size() == 125) {
            numColumns = 8;
            textSize   = 0.022;
            xStart     = 0.52; xEnd = 0.93;
            yStart     = 0.9;  yEnd = 0.45;
            xSpacingFactor = 1.0;
            ySpacingFactor = 1.0;
        }
        else if (runNumbers.size() == 54) {
            numColumns = 3;
            textSize   = 0.025;
            xStart     = 0.57; xEnd = 0.94;
            yStart     = 0.9;  yEnd = 0.42;
            xSpacingFactor = 0.8;
            ySpacingFactor = 0.92;
        }
        else if (runNumbers.size() == 126) {
            numColumns = 8;
            textSize   = 0.022;
            xStart     = 0.52; xEnd = 0.93;
            yStart     = 0.9;  yEnd = 0.45;
            xSpacingFactor = 1.0;
            ySpacingFactor = 1.0;
        }
        else if (runNumbers.size() == 146) {
            numColumns = 8;
            textSize   = 0.018;
            xStart     = 0.54; xEnd = 0.93;
            yStart     = 0.9;  yEnd = 0.5;
            xSpacingFactor = 1.0;
            ySpacingFactor = 1.0;
        }
        else {
            numColumns = (runNumbers.size() <= 20) ? 2 : 5;
            textSize   = 0.03;
            xStart     = 0.7;  xEnd = 0.9;
            yStart     = 0.85; yEnd = 0.3;
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

        // Some sets we skip
        if (runNumbers.size() == 646 || runNumbers.size() == 251 || runNumbers.size() == 102 ||
            runNumbers.size() == 101 || runNumbers.size() == 240 || runNumbers.size() == 142 ||
            runNumbers.size() == 440 || runNumbers.size() == 443 || runNumbers.size() == 725 ||
            runNumbers.size() == 100 || runNumbers.size() == 250 || runNumbers.size() == 141 ||
            runNumbers.size() == 254 || runNumbers.size() == 641 || runNumbers.size() == 291 ||
            runNumbers.size() == 326 || runNumbers.size() == 723 || runNumbers.size() == 382 ||
            runNumbers.size() == 347 || runNumbers.size() == 88  || runNumbers.size() == 644) {
            std::cout << "[INFO] Skipping run number plotting for runNumbers.size() = 692.\n";
            return;
        }

        runNumbersLatex.SetTextSize(textSize);

        int numRows = (runNumbers.size() + numColumns - 1) / numColumns;
        std::vector<std::vector<std::string>> grid(numRows, std::vector<std::string>(numColumns, ""));
        int runIndex = 0;
        for (int col = 0; col < numColumns; ++col) {
            for (int row = 0; row < numRows; ++row) {
                if (runIndex < (int)runNumbers.size()) {
                    grid[row][col] = std::to_string(runNumbers[runIndex]);
                    ++runIndex;
                }
            }
        }

        double xSpacing = xSpacingFactor * (xEnd - xStart) / numColumns;
        double ySpacing = ySpacingFactor * (yStart - yEnd) / (numRows + 1);

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

    // --------------------------------------------------------------------------
    // 6) Loop over each combined ROOT file
    // --------------------------------------------------------------------------
    for (const auto& rootFileName : combinedRootFiles)
    {
        std::string rootFilePath = outputDirectory + "/" + rootFileName;
        std::cout << "\n[DEBUG] ------------------------------------------------------\n";
        std::cout << "[DEBUG] Processing combined ROOT file: " << rootFilePath << "\n";

        // Extract triggers
        std::vector<std::string> triggers = ExtractTriggersFromFilename(rootFileName, allTriggers);
        std::cout << "[DEBUG] Extracted triggers => [ ";
        for (auto& t : triggers) std::cout << t << " ";
        std::cout << "]\n";

        // Firmware tag
        std::string firmwareTag;
        if (Utils::EndsWith(rootFileName, "_beforeTriggerFirmwareUpdate_Combined.root")) {
            firmwareTag = "_beforeTriggerFirmwareUpdate";
        }
        else if (Utils::EndsWith(rootFileName, "_afterTriggerFirmwareUpdate_Combined.root")) {
            firmwareTag = "_afterTriggerFirmwareUpdate";
        }

        std::string firmwareStatus;
        if (firmwareTag == "_beforeTriggerFirmwareUpdate") {
            firmwareStatus = "#bf{Before firmware update at run 47289}";
        }
        else if (firmwareTag == "_afterTriggerFirmwareUpdate") {
            firmwareStatus = "#bf{After firmware update at run 47289}";
        }

        // Build combinationName
        std::string combinationName;
        for (const auto& trig : triggers) {
            combinationName += trig + "_";
        }
        if (!combinationName.empty()) combinationName.pop_back();
        combinationName += firmwareTag;

        // Sanitize
        std::string sanitizedCombinationName = combinationName;
        std::replace(sanitizedCombinationName.begin(), sanitizedCombinationName.end(), '/', '_');
        std::replace(sanitizedCombinationName.begin(), sanitizedCombinationName.end(), ' ', '_');

        // Make subdir for plots
        std::string plotDirectory = basePlotDirectory + "/" + sanitizedCombinationName;
        gSystem->mkdir(plotDirectory.c_str(), true);

        // ----------------------------------------------------------------------
        // Open combined root file
        // ----------------------------------------------------------------------
        TFile* inputFile = TFile::Open(rootFilePath.c_str(), "READ");
        if (!inputFile || inputFile->IsZombie()) {
            std::cerr << "[ERROR] Could not open file => " << rootFilePath << "\n";
            if (inputFile) { delete inputFile; }
            continue;
        }
        std::cout << "[DEBUG] Successfully opened => " << rootFilePath << "\n";

        // ----------------------------------------------------------------------
        // 7) PHOTON overlay with histPrefix
        // ----------------------------------------------------------------------
        TCanvas* canvas = new TCanvas("canvas", "Overlay Plot", 800, 600);
        TLegend* legend = new TLegend(0.5, 0.58, 0.8, 0.88);
        legend->SetTextSize(0.028);
        canvas->SetLogy();
        bool firstDraw = true;

        for (const auto& trigger : triggers)
        {
            // Exclude Jet triggers from this overlay
            if (trigger == "Jet_8_GeV_plus_MBD_NS_geq_1" ||
                trigger == "Jet_10_GeV_plus_MBD_NS_geq_1" ||
                trigger == "Jet_12_GeV_plus_MBD_NS_geq_1")
            {
                continue;
            }

            TDirectory* triggerDir = inputFile->GetDirectory(trigger.c_str());
            if (!triggerDir) {
                std::cerr << "[WARN] Directory '" << trigger
                          << "' not found => " << rootFileName << "\n";
                continue;
            }

            std::string histName = histPrefix + trigger;
            TH1* hist = (TH1*)triggerDir->Get(histName.c_str());
            if (!hist) {
                std::cerr << "[WARN] Hist => " << histName
                          << " not found => " << rootFileName << "\n";
                continue;
            }

            TH1* histClone = (TH1*)hist->Clone();
            histClone->SetDirectory(0);

            int color = kBlack;
            auto itC = triggerColorMap.find(trigger);
            if (itC != triggerColorMap.end()) color = itC->second;

            histClone->SetLineColor(color);
            histClone->SetLineWidth(2);
            histClone->SetTitle(("Overlay for " + combinationName).c_str());
            histClone->GetXaxis()->SetTitle(xAxisTitle.c_str());
            histClone->GetYaxis()->SetTitle("Prescaled Counts");

            if (firstDraw) {
                histClone->Draw("HIST");
                firstDraw = false;
            }
            else {
                histClone->Draw("HIST SAME");
            }

            // Legend
            std::string displayTriggerName = trigger;
            auto it_name = triggerNameMap.find(trigger);
            if (it_name != triggerNameMap.end()) {
                displayTriggerName = it_name->second;
            }
            legend->AddEntry(histClone, displayTriggerName.c_str(), "l");
        }

        legend->Draw();

        auto itRunNumbers = combinationToValidRuns.find(combinationName);
        if (itRunNumbers != combinationToValidRuns.end()) {
            drawRunNumbersOnCanvas(itRunNumbers->second, firmwareStatus);
        }
        else {
            TLatex text;
            text.SetNDC();
            text.SetTextAlign(13);
            text.SetTextSize(0.03);
            text.DrawLatex(0.15, 0.6, "Run numbers not available");
        }

        canvas->Modified();
        canvas->Update();
        std::string outFileName = plotDirectory + "/" + histPrefixForFile + "_Overlay.png";
        canvas->SaveAs(outFileName.c_str());
        std::cout << "[INFO] Saved overlay => " << outFileName << "\n";
        delete canvas;

        // ----------------------------------------------------------------------
        // 8) Photon turn-on
        // ----------------------------------------------------------------------
        std::vector<std::string> photonTriggersInCombination;
        for (const auto& trig : triggers) {
            if (std::find(photonTriggers.begin(), photonTriggers.end(), trig) != photonTriggers.end()) {
                photonTriggersInCombination.push_back(trig);
            }
        }
        if (!photonTriggersInCombination.empty()) {
            TDirectory* minbiasDir = inputFile->GetDirectory("MBD_NandS_geq_1");
            if (!minbiasDir) {
                std::cerr << "[WARN] MBD_NandS_geq_1 dir missing => " << rootFileName << "\n";
            }
            else {
                std::string minbiasHistName = histPrefix + "MBD_NandS_geq_1";
                TH1* minbiasHist = (TH1*)minbiasDir->Get(minbiasHistName.c_str());
                if (!minbiasHist) {
                    std::cerr << "[WARN] minbias hist => " << minbiasHistName
                              << " not found => " << rootFileName << "\n";
                }
                else {
                    TCanvas* cPhotonTurnOn = new TCanvas("cPhotonTurnOn", "Photon Turn-On", 800,600);
                    TLegend* legPhotonTurnOn = new TLegend(0.18, 0.72, 0.45, 0.9);
                    legPhotonTurnOn->SetTextSize(0.028);
                    bool firstDrawTurnOn = true;

                    for (const auto& photonTrigger : photonTriggersInCombination) {
                        TDirectory* pDir = inputFile->GetDirectory(photonTrigger.c_str());
                        if (!pDir) {
                            std::cerr << "[WARN] Photon trig dir => " << photonTrigger
                                      << " not found => " << rootFileName << "\n";
                            continue;
                        }
                        std::string photonHistName = histPrefix + photonTrigger;
                        TH1* photonHist = (TH1*)pDir->Get(photonHistName.c_str());
                        if (!photonHist) {
                            std::cerr << "[WARN] Hist => " << photonHistName
                                      << " not found => " << rootFileName << "\n";
                            continue;
                        }

                        TH1* ratioHist = (TH1*)photonHist->Clone(("ratio_" + photonTrigger).c_str());
                        ratioHist->SetDirectory(0);
                        ratioHist->Divide(photonHist, minbiasHist, 1.0,1.0,"B");

                        int color = kBlack;
                        auto itColor = triggerColorMap.find(photonTrigger);
                        if (itColor != triggerColorMap.end()) color = itColor->second;

                        ratioHist->SetMarkerStyle(20);
                        ratioHist->SetMarkerColor(color);
                        ratioHist->SetLineColor(color);
                        ratioHist->SetTitle(("Turn-On for " + combinationName).c_str());
                        ratioHist->GetXaxis()->SetTitle(xAxisTitleTurnOn.c_str());
                        ratioHist->GetYaxis()->SetTitle("Ratio to MBD");
                        ratioHist->GetYaxis()->SetRangeUser(0, 2.0);

                        if (firstDrawTurnOn) {
                            ratioHist->Draw("E1");
                            firstDrawTurnOn = false;
                        } else {
                            ratioHist->Draw("E1 SAME");
                        }

                        TLine* line = new TLine(ratioHist->GetXaxis()->GetXmin(),1,
                                                ratioHist->GetXaxis()->GetXmax(),1);
                        line->SetLineStyle(1);
                        line->SetLineColor(kBlack);
                        line->Draw("SAME");

                        // Fit if available
                        std::string displayPhoton = photonTrigger;
                        auto it_nm = triggerNameMap.find(photonTrigger);
                        if (it_nm != triggerNameMap.end()) {
                            displayPhoton = it_nm->second;
                        }
                        std::ostringstream legendEntry;
                        legendEntry << displayPhoton;

                        // Search fit parameters
                        std::pair<std::string, std::string> key = std::make_pair(combinationName, photonTrigger);
                        auto it_fitParams = TriggerConfig::triggerFitParameters.find(key);
                        if (it_fitParams == TriggerConfig::triggerFitParameters.end()) {
                            std::string comboNoFirmware = Utils::stripFirmwareTag(combinationName);
                            key = std::make_pair(comboNoFirmware, photonTrigger);
                            it_fitParams = TriggerConfig::triggerFitParameters.find(key);
                        }
                        if (it_fitParams == TriggerConfig::triggerFitParameters.end()) {
                            key = std::make_pair("", photonTrigger);
                            it_fitParams = TriggerConfig::triggerFitParameters.find(key);
                        }

                        if (enableFits && it_fitParams != TriggerConfig::triggerFitParameters.end()) {
                            DataStructures::FitParameters params = it_fitParams->second;
                            TF1* fitFunc = nullptr;
                            if (fitFunctionType == "sigmoid") {
                                fitFunc = Utils::sigmoidFit(("fit_" + photonTrigger).c_str(),
                                                            0.0, 20.0,
                                                            params.amplitudeEstimate,
                                                            params.slopeEstimate,
                                                            params.xOffsetEstimate,
                                                            params.amplitudeMin,
                                                            params.amplitudeMax,
                                                            params.slopeMin,
                                                            params.slopeMax,
                                                            params.xOffsetMin,
                                                            params.xOffsetMax);
                            }
                            else if (fitFunctionType == "erf") {
                                fitFunc = Utils::erfFit(("fit_" + photonTrigger).c_str(),
                                                        0.0, 20.0,
                                                        params.amplitudeEstimate,
                                                        params.xOffsetEstimate,
                                                        params.sigmaEstimate,
                                                        params.amplitudeMin,
                                                        params.amplitudeMax,
                                                        params.xOffsetMin,
                                                        params.xOffsetMax,
                                                        params.sigmaMin,
                                                        params.sigmaMax);
                            }
                            else {
                                // default to sigmoid
                                fitFunc = Utils::sigmoidFit(("fit_" + photonTrigger).c_str(),
                                                            0.0, 20.0,
                                                            params.amplitudeEstimate,
                                                            params.slopeEstimate,
                                                            params.xOffsetEstimate,
                                                            params.amplitudeMin,
                                                            params.amplitudeMax,
                                                            params.slopeMin,
                                                            params.slopeMax,
                                                            params.xOffsetMin,
                                                            params.xOffsetMax);
                            }

                            fitFunc->SetLineColor(color);
                            ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
                            ratioHist->Fit(fitFunc, "R");
                            fitFunc->Draw("SAME");

                            double A      = fitFunc->GetParameter(0);
                            double A_error= fitFunc->GetParError(0);

                            double x99      = 0;
                            double x99_error= 0;
                            if (fitFunctionType == "sigmoid") {
                                double k      = fitFunc->GetParameter(1);
                                double x0     = fitFunc->GetParameter(2);
                                double k_error= fitFunc->GetParError(1);
                                double x0_error= fitFunc->GetParError(2);

                                x99 = x0 + (std::log(99) / k);
                                x99_error = std::sqrt((x0_error*x0_error) +
                                                      (std::pow((std::log(99)/(k*k)),2)*k_error*k_error));
                                std::cout << "[FIT] Sigmoid for " << photonTrigger
                                          << ": A=" << A << " ± " << A_error
                                          << ", k=" << k << " ± " << k_error
                                          << ", x0="<< x0 << " ± " << x0_error
                                          << ", x99="<< x99 << " ± " << x99_error << " GeV\n";
                            }
                            else if (fitFunctionType == "erf") {
                                double x0 = fitFunc->GetParameter(1);
                                double sigma = fitFunc->GetParameter(2);
                                double x0_error = fitFunc->GetParError(1);
                                double sigma_error = fitFunc->GetParError(2);
                                double erfInvVal = TMath::ErfInverse(0.98);
                                x99 = x0 + sqrt(2)*sigma*erfInvVal;
                                x99_error = std::sqrt(x0_error*x0_error +
                                    (sqrt(2)*erfInvVal*sigma_error)*
                                    (sqrt(2)*erfInvVal*sigma_error));
                                std::cout << "[FIT] Erf for " << photonTrigger
                                          << ": A=" << A << " ± " << A_error
                                          << ", x0=" << x0 << " ± " << x0_error
                                          << ", sigma=" << sigma << " ± " << sigma_error
                                          << ", x99=" << x99 << " ± " << x99_error << " GeV\n";
                            }

                            triggerEfficiencyPoints[photonTrigger] = x99;
                            combinationToTriggerEfficiencyPoints[combinationName][photonTrigger] = x99;

                            legendEntry << ", 99% eff = "
                                        << std::fixed << std::setprecision(2)
                                        << x99 << " GeV";

                            if (x99 > ratioHist->GetXaxis()->GetXmin() &&
                                x99 < ratioHist->GetXaxis()->GetXmax()) {
                                TLine* verticalLine = new TLine(x99, 0, x99, 1);
                                verticalLine->SetLineStyle(2);
                                verticalLine->SetLineColor(color);
                                verticalLine->SetLineWidth(5);
                                verticalLine->Draw("SAME");
                            }
                        }
                        else {
                            std::cerr << "[INFO] No fit parameters found => "
                                      << photonTrigger
                                      << " in combination => " << combinationName << "\n";
                        }

                        legPhotonTurnOn->AddEntry(ratioHist, legendEntry.str().c_str(), "p");
                        std::cout << "[DEBUG] Legend Entry Added => " << legendEntry.str() << "\n";
                    }

                    legPhotonTurnOn->Draw();

                    TLegend* legendEfficiencyLine = new TLegend(0.18, 0.62, 0.38, 0.72);
                    legendEfficiencyLine->SetTextSize(0.03);
                    legendEfficiencyLine->SetBorderSize(0);
                    legendEfficiencyLine->SetFillStyle(0);
                    TLine* dummyLine = new TLine(0,0,0,0);
                    dummyLine->SetLineStyle(2);
                    dummyLine->SetLineColor(kGray+1);
                    dummyLine->SetLineWidth(2);
                    legendEfficiencyLine->AddEntry(dummyLine, "99% Efficiency Point", "l");
                    legendEfficiencyLine->Draw();

                    if (!firmwareStatus.empty()) {
                        cPhotonTurnOn->cd();
                        TLatex firmwareStatusText;
                        firmwareStatusText.SetNDC();
                        firmwareStatusText.SetTextAlign(22);
                        firmwareStatusText.SetTextSize(0.03);
                        firmwareStatusText.SetTextColor(kBlack);
                        firmwareStatusText.DrawLatex(0.5, 0.96, firmwareStatus.c_str());
                    }

                    cPhotonTurnOn->Modified();
                    cPhotonTurnOn->Update();

                    std::string outputTurnOnFileName = plotDirectory + "/" + histPrefixForFile + "_TurnOn.png";
                    cPhotonTurnOn->SaveAs(outputTurnOnFileName.c_str());
                    std::cout << "[INFO] Saved turn-on plot => " << outputTurnOnFileName << "\n";

                    delete cPhotonTurnOn;
                }
            }
        }

        // ----------------------------------------------------------------------
        // 9) Jet triggers => always use h_jet_energy_ prefix
        // ----------------------------------------------------------------------
        {
            // Identify jet triggers in this combination
            std::vector<std::string> jetTriggersInCombination;
            for (const auto& trig : triggers)
            {
                if (trig == "Jet_8_GeV_plus_MBD_NS_geq_1" ||
                    trig == "Jet_10_GeV_plus_MBD_NS_geq_1" ||
                    trig == "Jet_12_GeV_plus_MBD_NS_geq_1")
                {
                    jetTriggersInCombination.push_back(trig);
                }
            }

            // If any jet triggers are present, we do TWO steps:
            //  (A) Overlay the raw jet-energy distributions
            //  (B) Make the ratio-based turn-on plot vs. MBD
            if (!jetTriggersInCombination.empty())
            {
                //----------------------------------------------------------------------
                // (A) JET RAW OVERLAY (similar to photon overlay)
                //----------------------------------------------------------------------
                {
                    TCanvas* cJetOverlay = new TCanvas("cJetOverlay", "Jet Energy Overlay", 800, 600);
                    cJetOverlay->SetLogy();

                    TLegend* legJetOverlay = new TLegend(0.5, 0.58, 0.8, 0.88);
                    legJetOverlay->SetTextSize(0.028);

                    // Also include MBD for the raw overlay
                    std::vector<std::string> jetsPlusMBD = jetTriggersInCombination;
                    jetsPlusMBD.push_back("MBD_NandS_geq_1");

                    // We'll store hist clones to find global maxima
                    std::vector<TH1*> clonedHists;
                    double globalMaxBinContent = 0.0;
                    double globalMaxXValue     = 0.0;

                    // 1) First pass: read & clone => track global maxima for y and x
                    for (const auto& possibleTrig : jetsPlusMBD)
                    {
                        TDirectory* jDir = inputFile->GetDirectory(possibleTrig.c_str());
                        if (!jDir)
                        {
                            std::cerr << "[WARN][JetOverlay] Directory '" << possibleTrig
                                      << "' not found => " << rootFileName << "\n";
                            continue;
                        }

                        // Build histogram name, always h_jet_energy_ for jets/MBD
                        std::string jetHName = jetHistPrefix + possibleTrig;
                        TH1* jetHist = (TH1*)jDir->Get(jetHName.c_str());
                        if (!jetHist)
                        {
                            std::cerr << "[WARN][JetOverlay] Hist => " << jetHName
                                      << " not found => " << rootFileName << "\n";
                            continue;
                        }

                        // Clone so we can manipulate
                        TH1* histClone = (TH1*)jetHist->Clone();
                        histClone->SetDirectory(nullptr);

                        // find local maxima for Y + furthest non-empty bin for X
                        double localMax = histClone->GetMaximum();
                        if (localMax > globalMaxBinContent) {
                            globalMaxBinContent = localMax;
                        }
                        int lastNonEmptyBin = histClone->FindLastBinAbove(0.0);
                        double maxXval      = histClone->GetXaxis()->GetBinUpEdge(lastNonEmptyBin);
                        if (maxXval > globalMaxXValue) {
                            globalMaxXValue = maxXval;
                        }

                        clonedHists.push_back(histClone);
                    }

                    // Decide on final X range: at least 50 or largest bin edge
                    double finalXmax = (globalMaxXValue < 50.0) ? 50.0 : (globalMaxXValue + 2.0);

                    // Decide on final Y range for log scale => ~30% above global max
                    double finalYmax = globalMaxBinContent * 1.3;

                    // 2) Second pass: draw them
                    bool firstJetDraw = true;
                    for (auto* histClone : clonedHists)
                    {
                        // Helper to see if 'name' ends with 'suffix'
                        auto endsWith = [&](const std::string& fullName, const std::string& suffix) {
                            return (fullName.size() >= suffix.size()) &&
                                   (0 == fullName.compare(fullName.size() - suffix.size(), suffix.size(), suffix));
                        };

                        // Figure out which actual trigger name from the histogram name:
                        std::string hName = histClone->GetName();  // e.g. "h_jet_energy_Jet_8_GeV_plus_MBD_NS_geq_1"

                        // We'll guess a fallback trigger is MBD:
                        std::string matchedTrigger = "MBD_NandS_geq_1";
                        // But check if it ends with each possible jet name:
                        if (endsWith(hName, "Jet_8_GeV_plus_MBD_NS_geq_1"))   matchedTrigger = "Jet_8_GeV_plus_MBD_NS_geq_1";
                        else if (endsWith(hName, "Jet_10_GeV_plus_MBD_NS_geq_1")) matchedTrigger = "Jet_10_GeV_plus_MBD_NS_geq_1";
                        else if (endsWith(hName, "Jet_12_GeV_plus_MBD_NS_geq_1")) matchedTrigger = "Jet_12_GeV_plus_MBD_NS_geq_1";

                        // Assign color based on matchedTrigger
                        int color = kBlack;
                        auto itCol = triggerColorMap.find(matchedTrigger);
                        if (itCol != triggerColorMap.end()) color = itCol->second;

                        // Set drawing style
                        histClone->SetLineColor(color);
                        histClone->SetLineWidth(2);
                        histClone->SetTitle(("Jet Overlay for " + combinationName).c_str());
                        histClone->GetXaxis()->SetTitle("Jet Energy [GeV]");
                        histClone->GetYaxis()->SetTitle("Prescaled Counts");

                        // Force the X and Y ranges
                        histClone->GetXaxis()->SetRangeUser(0.0, finalXmax);
                        histClone->SetMaximum(finalYmax);

                        // Draw
                        if (firstJetDraw) {
                            histClone->Draw("HIST");
                            firstJetDraw = false;
                        }
                        else {
                            histClone->Draw("HIST SAME");
                        }

                        // Legend label
                        std::string displayJet = matchedTrigger; // fallback
                        auto it_nm = triggerNameMap.find(matchedTrigger);
                        if (it_nm != triggerNameMap.end()) {
                            displayJet = it_nm->second; // nicer label
                        }
                        legJetOverlay->AddEntry(histClone, displayJet.c_str(), "l");
                    }

                    legJetOverlay->Draw();

                    // Possibly run numbers
                    auto it_runJetOverlay = combinationToValidRuns.find(combinationName);
                    if (it_runJetOverlay != combinationToValidRuns.end()) {
                        drawRunNumbersOnCanvas(it_runJetOverlay->second, firmwareStatus);
                    }

                    cJetOverlay->Modified();
                    cJetOverlay->Update();

                    // Save
                    std::string outJetOverlay = plotDirectory + "/" + jetHistPrefix + "_JetOverlay.png";
                    cJetOverlay->SaveAs(outJetOverlay.c_str());
                    std::cout << "[INFO] Saved Jet overlay => " << outJetOverlay << "\n";

                    // Clean up
                    for (auto* h : clonedHists) delete h;
                    delete cJetOverlay;
                }

                //----------------------------------------------------------------------
                // (B) Jet turn-on ratio vs. MBD
                //----------------------------------------------------------------------
                TDirectory* jetMinbiasDir = inputFile->GetDirectory("MBD_NandS_geq_1");
                if (!jetMinbiasDir)
                {
                    std::cerr << "[WARN][Jet] MBD_NandS_geq_1 dir not found => " << rootFileName
                              << "\n=> skip Jet turn-on\n";
                }
                else
                {
                    // minbias hist => "h_jet_energy_MBD_NandS_geq_1"
                    std::string mbJetHistName = jetHistPrefix + "MBD_NandS_geq_1";
                    TH1* mbJetHist = (TH1*)jetMinbiasDir->Get(mbJetHistName.c_str());
                    if (!mbJetHist)
                    {
                        std::cerr << "[WARN][Jet] Hist => " << mbJetHistName
                                  << " not found => " << rootFileName
                                  << "\n=> skip Jet turn-on\n";
                    }
                    else
                    {
                        // Ratio-based "Jet Turn-On"
                        TCanvas* cJetTurnOn = new TCanvas("cJetTurnOn", "Jet Turn-On", 800, 600);
                        TLegend* legJetTurnOn = new TLegend(0.18, 0.72, 0.45, 0.9);
                        legJetTurnOn->SetTextSize(0.028);
                        bool firstJet = true;

                        for (const auto& jetTrig : jetTriggersInCombination)
                        {
                            TDirectory* jDir = inputFile->GetDirectory(jetTrig.c_str());
                            if (!jDir)
                            {
                                std::cerr << "[WARN][Jet] Dir => " << jetTrig
                                          << " missing => " << rootFileName << "\n";
                                continue;
                            }

                            std::string jetHName = jetHistPrefix + jetTrig;
                            TH1* jetHist = (TH1*)jDir->Get(jetHName.c_str());
                            if (!jetHist)
                            {
                                std::cerr << "[WARN][Jet] Hist => " << jetHName
                                          << " not found => " << rootFileName << "\n";
                                continue;
                            }

                            // Clone + ratio vs. MBD
                            TH1* ratioJet = (TH1*)jetHist->Clone(("ratioJet_" + jetTrig).c_str());
                            ratioJet->SetDirectory(0);
                            ratioJet->Divide(jetHist, mbJetHist, 1.0, 1.0, "B");
                            ratioJet->SetTitle(("Jet Turn-On => " + combinationName).c_str());
                            ratioJet->GetXaxis()->SetTitle("Jet Energy [GeV]");
                            ratioJet->GetYaxis()->SetTitle("Ratio to MBD");
                            ratioJet->GetYaxis()->SetRangeUser(0, 2.0);

                            int color = kGray+3;
                            auto itCol = triggerColorMap.find(jetTrig);
                            if (itCol != triggerColorMap.end()) color = itCol->second;

                            ratioJet->SetMarkerStyle(20);
                            ratioJet->SetMarkerColor(color);
                            ratioJet->SetLineColor(color);

                            if (firstJet)
                            {
                                ratioJet->Draw("E1");
                                firstJet = false;
                            }
                            else
                            {
                                ratioJet->Draw("E1 SAME");
                            }

                            // dashed line @ y=1
                            TLine* ln = new TLine(ratioJet->GetXaxis()->GetXmin(), 1,
                                                  ratioJet->GetXaxis()->GetXmax(), 1);
                            ln->SetLineStyle(2);
                            ln->SetLineColor(kBlack);
                            ln->Draw("SAME");

                            // Legend
                            std::string displayJet = jetTrig;
                            auto it_nm = triggerNameMap.find(jetTrig);
                            if (it_nm != triggerNameMap.end())
                            {
                                displayJet = it_nm->second;
                            }
                            legJetTurnOn->AddEntry(ratioJet, displayJet.c_str(), "p");
                        }

                        legJetTurnOn->Draw();

                        // Possibly run numbers
                        auto it_runJet = combinationToValidRuns.find(combinationName);
                        if (it_runJet != combinationToValidRuns.end())
                        {
                            drawRunNumbersOnCanvas(it_runJet->second, firmwareStatus);
                        }

                        cJetTurnOn->Modified();
                        cJetTurnOn->Update();

                        // Save
                        std::string outJetTurnOn = plotDirectory + "/" + jetHistPrefix + "_JetTurnOn.png";
                        cJetTurnOn->SaveAs(outJetTurnOn.c_str());
                        delete cJetTurnOn;
                    }
                } // end else of MBD directory check
            } // end if (!jetTriggersInCombination.empty())
        }
        // ----------------------------------------------------------------------
        // 9) If fitOnly => skip rest
        // ----------------------------------------------------------------------
        if (fitOnly) {
            std::cout << "[DEBUG] fitOnly == true => skipping rest.\n";
            inputFile->Close();
            delete inputFile;

            auto it_run = combinationToValidRuns.find(combinationName);
            if (it_run != combinationToValidRuns.end()) {
                const std::vector<int>& validRuns = it_run->second;
                PlotRunByRunHistograms(outputDirectory, plotDirectory, triggers,
                                       validRuns, triggerColorMap, triggerNameMap,
                                       firmwareStatus);
            }
            else {
                std::cout << "[DEBUG] No valid runs found for => " << combinationName << "\n";
            }
            continue;
        }

        // ----------------------------------------------------------------------
        // 10) Invariant Mass
        // ----------------------------------------------------------------------
        ProcessInvariantMassHistograms(inputFile, plotDirectory, triggers,
                                       triggerColorMap, combinationName);

        // ----------------------------------------------------------------------
        // 11) Meson Mass vs Pt
        // ----------------------------------------------------------------------
        ProcessMesonMassVsPt(plotDirectory, combinationName, triggers,
                             triggerEfficiencyPoints, DataStructures::pT_bins);

        // ----------------------------------------------------------------------
        // 12) Isolation Energy
        // ----------------------------------------------------------------------
        std::cout << "[DEBUG] About to call ProcessIsolationEnergyHistogramsWithCuts => '"
                  << combinationName << "'...\n";

        try {
            ProcessIsolationEnergyHistogramsWithCuts(inputFile, plotDirectory,
                                                     triggers, combinationName);
        }
        catch (...) {
            std::cerr << "[FATAL] Exception in ProcessIsolationEnergyHistogramsWithCuts for => "
                      << combinationName << "\n";
        }

        // Cleanup
        inputFile->Close();
        delete inputFile;

        // ----------------------------------------------------------------------
        // 13) Run-by-run overlays
        // ----------------------------------------------------------------------
        auto it_run = combinationToValidRuns.find(combinationName);
        if (it_run != combinationToValidRuns.end()) {
            const std::vector<int>& validRuns = it_run->second;
            PlotRunByRunHistograms(outputDirectory, plotDirectory, triggers,
                                   validRuns, triggerColorMap, triggerNameMap,
                                   firmwareStatus);
        }
        else {
            std::cout << "[DEBUG] No valid runs found for => " << combinationName << "\n";
        }
    } // end for each combinedRootFile

    // --------------------------------------------------------------------------
    // 14) Final isolation data steps
    // --------------------------------------------------------------------------
    std::string csvOutputPath = "/Users/patsfan753/Desktop/isolation_data.csv";
    WriteIsolationDataToCSV(
      csvOutputPath,
      dataMap_inMassWindow,
      dataMap_outsideMassWindow
    );
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

    std::cout << "[DEBUG] PlotCombinedHistograms() completed successfully.\n";
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

void ProcessAllShowerShapeQA(
    const std::string& outputDirectory,
    const std::vector<std::string>& combinedRootFiles,
    const std::map<std::string, std::vector<int>>& combinationToValidRuns )
{
    //--------------------------------------------------------------------------
    // Helper: check if 'str' ends with 'suffix'
    //--------------------------------------------------------------------------
    auto endsWith = [](const std::string& str, const std::string& suffix) {
        if (suffix.size() > str.size()) return false;
        return std::equal(suffix.rbegin(), suffix.rend(), str.rbegin());
    };

    //--------------------------------------------------------------------------
    // Helper: extract combination name from file name
    // e.g. "MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Combined.root"
    // => "MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1"
    //--------------------------------------------------------------------------
    auto extractCombinationName = [&](const std::string& filename) -> std::string {
        const std::string suffix = "_Combined.root";
        if (endsWith(filename, suffix)) {
            return filename.substr(0, filename.size() - suffix.size());
        }
        return filename; // fallback if doesn't match
    };

    //--------------------------------------------------------------------------
    // List of known "shower-shape" base names we want to overlay
    //--------------------------------------------------------------------------
    std::vector<std::string> showerShapeBases = {
        "E3by7_over_E7by7",
        "w72",
        "E1x1_over_ClusterE",
        "E1x1_over_E3x3",
        "E3x2_over_E3x5",
        "E1by7_over_E7by7",
        "E3x3_over_ClusterE",
        "E3x3_over_E3x7",
        "weta",
        "wphi"
    };

    //--------------------------------------------------------------------------
    // For debugging
    //--------------------------------------------------------------------------
    std::cout << "[DEBUG] Entering ProcessAllShowerShapeQA() ...\n";
    std::cout << "        # combinedRootFiles = " << combinedRootFiles.size() << "\n\n";

    //--------------------------------------------------------------------------
    // We'll define a base directory for the "Plots"
    //--------------------------------------------------------------------------
    std::string basePlotDir = "/Users/patsfan753/Desktop/DirectPhotonAna/Plots";

    //--------------------------------------------------------------------------
    // Loop over each combined root file
    //--------------------------------------------------------------------------
    for (const auto& fileName : combinedRootFiles)
    {
        // 1) Build the full path + extract combinationName
        std::string rootFilePath    = outputDirectory + "/" + fileName;
        std::string combinationName = extractCombinationName(fileName);

        // 2) Make sure we have an output folder for this combination
        std::string comboPlotDir = basePlotDir + "/" + combinationName;
        gSystem->mkdir(comboPlotDir.c_str(), true);

        std::string showerQA_Dir = comboPlotDir + "/showerShapeQA";
        gSystem->mkdir(showerQA_Dir.c_str(), true);

        // 3) Open the combined file
        std::cout << "[INFO] Processing combination: " << combinationName
                  << "\n       => " << rootFilePath << std::endl;

        TFile* inputFile = TFile::Open(rootFilePath.c_str(), "READ");
        if (!inputFile || inputFile->IsZombie()) {
            std::cerr << "[ERROR] Could NOT open file => " << rootFilePath << "\n";
            if (inputFile) { delete inputFile; }
            continue;
        }

        // 4) Get top-level directories => triggers
        TList* dirList = inputFile->GetListOfKeys();
        if (!dirList) {
            std::cerr << "[WARN] No top-level directories found in => "
                      << fileName << "\n";
            inputFile->Close();
            delete inputFile;
            continue;
        }

        // Optionally we can see if we have runNumbers:
        // auto itRun = combinationToValidRuns.find(combinationName);
        // if (itRun != combinationToValidRuns.end()) {
        //     const std::vector<int>& runNumbers = itRun->second;
        // }

        // 5) Loop over each TKey => each is presumably a directory for a trigger
        TIter nextDir(dirList);
        TKey* dirKey;
        while ((dirKey = (TKey*)nextDir())) {
            // read object
            TObject* dirObj = dirKey->ReadObj();
            if (!dirObj) {
                std::cerr << "[WARN] dirKey->ReadObj() returned null for key="
                          << dirKey->GetName() << "\n";
                continue;
            }

            // Check if it is indeed a TDirectory
            if (!dirObj->InheritsFrom(TDirectory::Class())) {
                // If it's not a TDirectory, we don't want it for shower-shape QA
                std::cout << "[DEBUG] Object '" << dirKey->GetName()
                          << "' is not a TDirectory => skipping.\n";
                delete dirObj;
                continue;
            }

            // If it is a directory:
            TDirectory* trigDir = dynamic_cast<TDirectory*>(dirObj);
            if (!trigDir) {
                std::cerr << "[ERROR] Could not cast object to TDirectory => "
                          << dirKey->GetName() << "\n";
                delete dirObj;
                continue;
            }

            // We'll keep 'trigDir' around. DO NOT DELETE 'dirObj' here,
            // or you'll lose 'trigDir'.
            std::string triggerName = trigDir->GetName();
            std::cout << "[DEBUG] Found trigger directory => " << triggerName << "\n";

            // Make an output subfolder for each trigger
            std::string trigPlotDir = showerQA_Dir + "/" + triggerName;
            gSystem->mkdir(trigPlotDir.c_str(), true);

            // 6) For each "shower-shape base" => look up the two hist names
            for (const auto& base : showerShapeBases) {

                std::string noCutName   = base + "_NoShowerShapeCuts_"   + triggerName;
                std::string withCutName = base + "_withShowerShapeCuts_" + triggerName;

                TH1* hNoCut   = dynamic_cast<TH1*>(trigDir->Get(noCutName.c_str()));
                TH1* hWithCut = dynamic_cast<TH1*>(trigDir->Get(withCutName.c_str()));

                if (!hNoCut && !hWithCut) {
                    // no relevant hists => skip
                    // std::cout << "[DEBUG] For base='" << base << "', no hist found => skipping.\n";
                    continue;
                }

                if (hNoCut) {
                    std::cout << "[DEBUG] Found => " << noCutName
                              << " (#entries=" << hNoCut->GetEntries() << ")\n";
                } else {
                    std::cerr << "[WARN] Missing => " << noCutName << "\n";
                }

                if (hWithCut) {
                    std::cout << "[DEBUG] Found => " << withCutName
                              << " (#entries=" << hWithCut->GetEntries() << ")\n";
                } else {
                    std::cerr << "[WARN] Missing => " << withCutName << "\n";
                }

                // 7) Make a canvas + overlay
                TCanvas c("c", ("Overlay: " + base).c_str(), 800, 600);
                TLegend leg(0.58, 0.65, 0.88, 0.82);
                leg.SetBorderSize(0);
                leg.SetTextSize(0.04);

                bool drawnSomething = false;

                if (hNoCut) {
                    hNoCut->SetLineColor(kBlack);
                    hNoCut->SetLineWidth(2);
                    hNoCut->SetStats(false);
                    hNoCut->SetTitle(base.c_str());

                    hNoCut->Draw("HIST");
                    leg.AddEntry(hNoCut, "NoShowerShapeCuts", "l");
                    drawnSomething = true;
                }

                if (hWithCut) {
                    hWithCut->SetLineColor(kRed);
                    hWithCut->SetLineWidth(2);
                    hWithCut->SetStats(false);

                    if (!drawnSomething) {
                        hWithCut->SetTitle(base.c_str());
                        hWithCut->Draw("HIST");
                        leg.AddEntry(hWithCut, "withShowerShapeCuts", "l");
                        drawnSomething = true;
                    } else {
                        hWithCut->Draw("HIST SAME");
                        leg.AddEntry(hWithCut, "withShowerShapeCuts", "l");
                    }
                }

                if (!drawnSomething) {
                    // theoretically won't happen because we continue if both are null
                    continue;
                }

                leg.Draw();
                c.Modified();
                c.Update();

                // 8) Save
                std::string outName = trigPlotDir + "/" + base + "_Overlay.png";
                c.SaveAs(outName.c_str());
                std::cout << "[INFO] => Saved overlay => " << outName << "\n";
            }

            // We do *not* delete 'trigDir' or 'dirObj' because we want
            // them to remain valid within the scope of this file read.

        } // end while (dirKey)

        // 9) Close the file
        inputFile->Close();
        delete inputFile;
    }

    std::cout << "\n[INFO] Finished ProcessAllShowerShapeQA.\n\n";
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
    
//    overrideTriggerStatus[46697]["Photon_3_GeV_plus_MBD_NS_geq_1"] = "OFF";
  
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
    
    ProcessAllShowerShapeQA(
        outputDirectory,
        combinedRootFiles,
        combinationToValidRuns  // If you do not need the run-number map, you can omit or pass an empty map
    );
}
