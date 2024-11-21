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




// Function to analyze combinations from the CSV file
std::map<std::set<std::string>, std::vector<int>> AnalyzeWhatTriggerGroupsAvailable(const std::string& csvFilePath) {
    const std::vector<std::string>& allTriggers = TriggerConfig::allTriggers;
    const std::vector<std::string>& photonTriggers = TriggerConfig::photonTriggers;

    
    // Map from trigger name to column index in the CSV
    std::map<std::string, int> triggerToIndex;
    
    // Map from run number to set of triggers that are 'ON'
    std::map<int, std::set<std::string>> runToActiveTriggers;
    
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
                if (status == "ON") {
                    activeTriggers.insert(trigger);
                }
            }
        }
        // Store active triggers for this run
        runToActiveTriggers[runNumber] = activeTriggers;
    }
    
    file.close();
    
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
    std::map<std::set<std::string>, std::vector<int>> combinationToRuns;
    
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
                    combinationToRuns[combination].push_back(runNumber);
                }
            }
        }
    }
    // Now, identify combinations with identical run lists
    // Map from sorted run lists to sets of trigger combinations
    std::map<std::vector<int>, std::vector<std::set<std::string>>> runListToCombinations;

    for (const auto& entry : combinationToRuns) {
        const std::set<std::string>& combination = entry.first;
        std::vector<int> runList = entry.second;
        std::sort(runList.begin(), runList.end()); // Ensure run list is sorted for comparison

        runListToCombinations[runList].push_back(combination);
    }

    // For each run list, find the largest combination(s) and remove subsets
    std::map<std::set<std::string>, std::vector<int>> filteredCombinationToRuns;

    for (const auto& entry : runListToCombinations) {
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
                filteredCombinationToRuns[combo] = runList;
            }
        }
    }

    // Now 'filteredCombinationToRuns' contains only the largest combinations per unique run list
    // You can return this map instead of 'combinationToRuns'
    return filteredCombinationToRuns;
}



void PrintSortedCombinations(const std::map<std::set<std::string>, std::vector<int>>& combinationToRuns) {
    // Sort combinations based on size and triggers
    std::vector<std::pair<std::set<std::string>, std::vector<int>>> sortedCombinations(
        combinationToRuns.begin(), combinationToRuns.end());

    std::sort(sortedCombinations.begin(), sortedCombinations.end(),
              [](const auto& a, const auto& b) {
                  // Compare based on the number of triggers active (excluding 'MBD_NandS_geq_1')
                  int countA = a.first.size();
                  int countB = b.first.size();
                  if (countA != countB) {
                      return countA < countB;
                  } else {
                      // If same number of triggers, sort alphabetically
                      return a.first < b.first;
                  }
              });

    // Print the sorted combinations
    for (const auto& kv : sortedCombinations) {
        const std::set<std::string>& triggers = kv.first;
        const std::vector<int>& runs = kv.second;

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


void ProcessAndMergeRootFiles(
    const std::map<std::set<std::string>, std::vector<int>>& combinationToRuns,
    const std::string& outputDirectory, std::map<std::string, std::vector<int>>& combinationToValidRuns) {
    
    std::cout << "Starting ProcessAndMergeRootFiles" << std::endl;
    // Disable automatic addition of histograms to directories
    TH1::AddDirectory(false);
    
    // Iterate over each trigger combination
    for (const auto& kv : combinationToRuns) {
        const std::set<std::string>& triggers = kv.first;
        const std::vector<int>& runs = kv.second;
        
        // Create a string to represent the combination for naming
        std::string combinationName;
        for (const auto& trigger : triggers) {
            combinationName += trigger + "_";
        }
        // Remove the trailing underscore
        if (!combinationName.empty()) {
            combinationName.pop_back();
        }
        
        std::cout << "\nProcessing combination: " << combinationName << std::endl;
        std::cout << "Runs in this combination: ";
        for (const auto& runNumber : runs) {
            std::cout << runNumber << " ";
        }
        std::cout << std::endl;
        
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
            continue;
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
        
        // If no valid runs were found, skip this combination
        if (validRuns.empty()) {
            std::cout << "No valid runs found for combination: " << combinationName << ". Skipping this combination." << std::endl;
            continue;
        }
        
        // If no histograms were merged, skip writing the output ROOT file
        if (mergedHistograms.empty()) {
            std::cout << "No histograms were merged for combination: " << combinationName << ". Skipping writing output ROOT file." << std::endl;
            continue;
        }
        
        // Write all merged histograms to the final ROOT file
        std::cout << "Writing merged histograms to final ROOT file: " << finalRootFilePath << std::endl;
        TFile finalFile(finalRootFilePath.c_str(), "RECREATE");
        if (finalFile.IsZombie() || !finalFile.IsOpen()) {
            std::cerr << "Failed to create final ROOT file: " << finalRootFilePath << std::endl;
            mergedHistograms.clear();
            continue;
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
    
    // Re-enable automatic addition of histograms to directories
    TH1::AddDirectory(true);
}


std::vector<std::string> ExtractTriggersFromFilename(const std::string& filename, const std::vector<std::string>& allTriggers) {
    // Remove '_Combined.root' suffix
    std::string baseName = filename;
    std::string suffix = "_Combined.root";
    if (Utils::EndsWith(baseName, suffix)) {
        baseName = baseName.substr(0, baseName.length() - suffix.length());
    }

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

    // Pi0 estimates
    double sigmaPi0Estimate = 0.021;  // Estimate for pi0 sigma
    double meanPi0Estimate = 0.135;   // Estimate for pi0 mean
    double amplitudePi0Estimate = hPi0Mass->GetBinContent(hPi0Mass->GetXaxis()->FindBin(meanPi0Estimate));

    // Eta estimates
    double sigmaEtaEstimate = 0.04;   // Eta sigma estimate
    double meanEtaEstimate = 0.62;    // Eta mean estimate
    double amplitudeEtaEstimate = 0.1 * amplitudePi0Estimate; // Eta amplitude much smaller than pi0

    // Define the totalFit function as two Gaussians (pi0 and eta) plus a fourth-order polynomial (pol4)
    totalFit = new TF1("totalFit", "gaus(0) + gaus(3) + pol4(6)", fitStart, fitEnd);
    totalFit->SetLineColor(kRed);

    // Set initial parameters for pi0 and eta Gaussian components
    totalFit->SetParameters(amplitudePi0Estimate, meanPi0Estimate, sigmaPi0Estimate, amplitudeEtaEstimate, meanEtaEstimate, sigmaEtaEstimate);

    // Set limits on eta mean and sigma to help fit convergence
    totalFit->SetParLimits(4, 0.55, 0.65);   // Eta mean constrained between 550 and 650 MeV
    totalFit->SetParLimits(5, 0.03, 0.05);   // Eta sigma constrained between 30 and 50 MeV

    // Perform the fit
    TFitResultPtr fitResult = hPi0Mass->Fit("totalFit", "SR+");

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
    for (int i = 6; i < 11; i++) {  // Parameters 6 to 10 correspond to the pol4 background
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
        ptRangeLabel << "pT: " << Utils::formatToThreeSigFigs(data.cuts.pTMin) << " - " << Utils::formatToThreeSigFigs(data.cuts.pTMax) << " GeV";
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
    labelText.DrawLatex(0.45, 0.9, "Active Trigger Group:");
    valueText.DrawLatex(0.65, 0.9, triggerGroupName.c_str());

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
        valueText.DrawLatex(0.65, 0.65, ptRangeLabel.str().c_str());
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

    // Call printHistogramData to write the CSV file
    printHistogramData(histogramDataVector, plotDirectory, combinationName);
}

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


void generateMesonPlotVsPt(
    const std::vector<double>& pTCenters,
    const std::vector<double>& meanValues,
    const std::vector<double>& meanErrors,
    const std::string& yAxisLabel,
    const std::string& outputFilePath,
    const std::string& triggerName,
    const std::string& cutCombination,
    int markerStyle,
    int markerColor,
    double clusECore,
    double chi,
    double asymmetry,
    double yMin = std::numeric_limits<double>::quiet_NaN(),
    double yMax = std::numeric_limits<double>::quiet_NaN(),
    double xMin = std::numeric_limits<double>::quiet_NaN(),
    double xMax = std::numeric_limits<double>::quiet_NaN()) {
    if (pTCenters.empty()) {
        return; // Nothing to plot if pTCenters is empty
    }
    
    TGraphErrors* graph = new TGraphErrors(pTCenters.size());
    for (size_t i = 0; i < pTCenters.size(); ++i) {
        graph->SetPoint(i, pTCenters[i], meanValues[i]);
        graph->SetPointError(i, 0, meanErrors[i]);
    }
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(1);
    graph->SetLineWidth(2);
    graph->SetMarkerColor(markerColor);
    graph->SetLineColor(markerColor);

    // *** Move xMin and xMax calculation here ***
    // Automatically calculate x-axis range if xMin or xMax is NaN
    if (std::isnan(xMin) || std::isnan(xMax)) {
        xMin = 2;  // Default starting point for x-axis
        xMax = *std::max_element(pTCenters.begin(), pTCenters.end()) * 1.1;  // Adding 10% margin
    }

    // Automatically calculate y-axis range if yMin or yMax is not provided
    if (std::isnan(yMin) || std::isnan(yMax)) {
        yMin = std::numeric_limits<double>::max();
        yMax = std::numeric_limits<double>::lowest();
        for (size_t i = 0; i < meanValues.size(); ++i) {
            double xVal = pTCenters[i];
            if (xVal >= xMin && xVal <= xMax) {
                double yVal = meanValues[i];
                double yErr = meanErrors[i];
                yMin = std::min(yMin, yVal - yErr);
                yMax = std::max(yMax, yVal + yErr);
            }
        }
        // Add a 5% margin around yMin and yMax
        double yMargin = 0.05 * (yMax - yMin);
        yMin -= yMargin;
        yMax += yMargin;
    }

    // Create canvas
    TCanvas canvas;

    // Draw graph
    TH1F* hFrame = canvas.DrawFrame(xMin, yMin, xMax, yMax, (";Leading Cluster p_{T} [GeV];" + yAxisLabel).c_str());
    hFrame->GetXaxis()->SetNdivisions(5);
    hFrame->GetXaxis()->SetLimits(xMin, xMax);
    hFrame->GetXaxis()->CenterLabels(false);

    graph->Draw("P SAME");

    // Add trigger and cut information on the top right of the plot
    TLatex labelText;
    labelText.SetNDC();
    labelText.SetTextSize(0.042);

    TLatex valueText;
    valueText.SetNDC();
    valueText.SetTextSize(0.042);

    // Use the map to translate the trigger name if available
    std::string displayTriggerName = triggerName;
    if (TriggerCombinationNames::triggerCombinationNameMap.find(triggerName) != TriggerCombinationNames::triggerCombinationNameMap.end()) {
        displayTriggerName = TriggerCombinationNames::triggerCombinationNameMap.at(triggerName);
    }

    labelText.DrawLatex(0.2, 0.87, "#font[62]{Trigger:}");
    valueText.DrawLatex(0.3, 0.87, displayTriggerName.c_str());
    
    
    labelText.DrawLatex(0.2, 0.8, "#font[62]{ECore #geq}");
    std::ostringstream eCoreWithUnit;
    eCoreWithUnit << clusECore << "   GeV";
    valueText.DrawLatex(0.42, 0.8, eCoreWithUnit.str().c_str());

    labelText.DrawLatex(0.2, 0.73, "#font[62]{#chi^{2} <}");
    std::ostringstream chiStr;
    chiStr << chi;
    valueText.DrawLatex(0.42, 0.73, chiStr.str().c_str());

    labelText.DrawLatex(0.2, 0.66, "#font[62]{Asymmetry <}");
    std::ostringstream asymmetryStr;
    asymmetryStr << asymmetry;
    valueText.DrawLatex(0.42, 0.66, asymmetryStr.str().c_str());

    // Ensure the directory exists
    std::string outputDirPath = outputFilePath.substr(0, outputFilePath.find_last_of("/"));
    gSystem->mkdir(outputDirPath.c_str(), true);

    // Save plot
    canvas.SaveAs(outputFilePath.c_str());
    std::cout << "Saved plot: " << outputFilePath << std::endl;

    // Clean up
    delete graph;
}


void ProcessMesonMassVsPt(const std::string& plotDirectory,
                          const std::string& combinationName,
                          const std::vector<std::string>& triggers,
                          const std::map<std::string, double>& triggerEfficiencyPoints) {
    // Construct the CSV file path
    std::string csvFilePath = plotDirectory + "/InvariantMassInformation_" + combinationName + ".csv";

    // Open the CSV file
    std::ifstream csvFile(csvFilePath);
    if (!csvFile.is_open()) {
        std::cerr << "Error: Could not open CSV file " << csvFilePath << std::endl;
        return;
    }

    // Read the CSV data
    std::vector<DataStructures::HistogramData> dataVector;

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

        dataVector.push_back(data);
    }

    csvFile.close();

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

    // For each cut combination, process data
    for (const auto& cutPair : dataByCuts) {
        const std::string& cutCombination = cutPair.first;
        const std::vector<DataStructures::HistogramData>& dataList = cutPair.second;

        // Collect pT centers and mean values for the main plot
        std::vector<double> pTCentersPi0;
        std::vector<double> meanPi0Values;
        std::vector<double> meanPi0Errors;

        std::vector<double> pTCentersEta;
        std::vector<double> meanEtaValues;
        std::vector<double> meanEtaErrors;

        double clusECore = dataList[0].cuts.clusECore;
        double chi = dataList[0].cuts.chi;
        double asymmetry = dataList[0].cuts.asymmetry;

        // Organize data per pT bin and trigger
        std::map<std::pair<double, double>, std::map<std::string, DataStructures::HistogramData>> dataPerPtBin;

        for (const auto& data : dataList) {
            if (data.cuts.pTMin == -1 || data.cuts.pTMax == -1) {
                continue; // Skip data without pT bins
            }

            std::pair<double, double> pTBin = std::make_pair(data.cuts.pTMin, data.cuts.pTMax);
            dataPerPtBin[pTBin][data.cuts.triggerName] = data;
        }

        // Maps to store data for the overlay plots
        std::set<std::string> triggersInData; // To keep track of triggers present in data

        // For π⁰
        std::map<std::string, std::vector<double>> triggerToPtCentersPi0;
        std::map<std::string, std::vector<double>> triggerToMeanPi0Values;
        std::map<std::string, std::vector<double>> triggerToMeanPi0Errors;

        // For η
        std::map<std::string, std::vector<double>> triggerToPtCentersEta;
        std::map<std::string, std::vector<double>> triggerToMeanEtaValues;
        std::map<std::string, std::vector<double>> triggerToMeanEtaErrors;

        // Now, for each pT bin, collect data from all triggers
        for (const auto& ptBinData : dataPerPtBin) {
            double pTMin = ptBinData.first.first;
            double pTMax = ptBinData.first.second;
            double pTCenter = (pTMin + pTMax) / 2.0;

            const auto& triggerDataMap = ptBinData.second;

            // For the overlay plots, collect data from all triggers
            for (const auto& triggerDataPair : triggerDataMap) {
                const std::string& triggerName = triggerDataPair.first;
                const DataStructures::HistogramData& data = triggerDataPair.second;

                // Validate meanPi0 before adding
                bool validPi0 = data.meanPi0 > 0.0 && data.meanPi0 < 0.4;
                bool validEta = data.meanEta > 0.3 && data.meanEta < 1.0;

                triggersInData.insert(triggerName);

                if (validPi0) {
                    triggerToPtCentersPi0[triggerName].push_back(pTCenter);
                    triggerToMeanPi0Values[triggerName].push_back(data.meanPi0);
                    triggerToMeanPi0Errors[triggerName].push_back(data.meanPi0Error);
                }

                if (validEta) {
                    triggerToPtCentersEta[triggerName].push_back(pTCenter);
                    triggerToMeanEtaValues[triggerName].push_back(data.meanEta);
                    triggerToMeanEtaErrors[triggerName].push_back(data.meanEtaError);
                }
            }

            // Decide which trigger to use for the main plot
            std::string triggerToUse = "MBD_NandS_geq_1";

            if (!triggerEfficiencyPoints.empty()) {
                for (const auto& photonTrigger : triggers) {
                    if (photonTrigger != "MBD_NandS_geq_1") {
                        auto it = triggerEfficiencyPoints.find(photonTrigger);
                        if (it != triggerEfficiencyPoints.end()) {
                            double x99 = it->second;
                            if (pTCenter >= x99) {
                                triggerToUse = photonTrigger;
                                break; // Use the first photon trigger that becomes efficient
                            }
                        }
                    }
                }
            }

            // For the main plot, check if data from the selected trigger is available
            auto dataIt = triggerDataMap.find(triggerToUse);
            if (dataIt != triggerDataMap.end()) {
                const DataStructures::HistogramData& selectedData = dataIt->second;

                // Validate meanPi0 and meanEta before adding
                bool validPi0 = selectedData.meanPi0 > 0.0 && selectedData.meanPi0 < 0.4;
                bool validEta = selectedData.meanEta > 0.3 && selectedData.meanEta < 1.0;

                if (!validPi0) {
                    std::cerr << "Warning: Invalid meanPi0 (" << selectedData.meanPi0
                              << ") for pT bin " << pTMin << " - " << pTMax
                              << ". Skipping this point for pi0." << std::endl;
                }

                if (!validEta) {
                    std::cerr << "Warning: Invalid meanEta (" << selectedData.meanEta
                              << ") for pT bin " << pTMin << " - " << pTMax
                              << ". Skipping this point for eta." << std::endl;
                }

                if (validPi0) {
                    pTCentersPi0.push_back(pTCenter);
                    meanPi0Values.push_back(selectedData.meanPi0);
                    meanPi0Errors.push_back(selectedData.meanPi0Error);
                    std::cout << "Added valid Pi0: pT = " << pTCenter
                              << ", meanPi0 = " << selectedData.meanPi0 << std::endl;
                }

                if (validEta) {
                    pTCentersEta.push_back(pTCenter);
                    meanEtaValues.push_back(selectedData.meanEta);
                    meanEtaErrors.push_back(selectedData.meanEtaError);
                    std::cout << "Added valid Eta: pT = " << pTCenter
                              << ", meanEta = " << selectedData.meanEta << std::endl;
                }
            } else {
                std::cerr << "Warning: Data not found for pT bin " << pTMin << " - " << pTMax
                          << " and trigger " << triggerToUse << std::endl;
            }
        }

        // Generate π⁰ mass vs pT plot if we have valid data
        if (!pTCentersPi0.empty()) {
            std::string outputFilePathPi0 = plotDirectory + "/" + cutCombination + "/meanPi0_vs_pT.png";
            generateMesonPlotVsPt(pTCentersPi0, meanPi0Values, meanPi0Errors,
                                  "Mean #pi^{0} Mass [GeV]", outputFilePathPi0, combinationName,
                                  cutCombination, 20, kRed, clusECore, chi, asymmetry, 0.14, 0.17, 2.0, 6.0);
        } else {
            std::cout << "No valid π⁰ data to plot for cut combination " << cutCombination << std::endl;
        }

        // Generate η mass vs pT plot if we have valid data
        if (!pTCentersEta.empty()) {
            std::string outputFilePathEta = plotDirectory + "/" + cutCombination + "/meanEta_vs_pT.png";
            generateMesonPlotVsPt(pTCentersEta, meanEtaValues, meanEtaErrors,
                                  "Mean #eta Mass [GeV]", outputFilePathEta, combinationName,
                                  cutCombination, 21, kBlue, clusECore, chi, asymmetry, 0.55, 0.72, 2.0, 6.0);
        } else {
            std::cout << "No valid η data to plot for cut combination " << cutCombination << std::endl;
        }

        // Generate overlay plot for π⁰ mass vs pT
        if (!triggerToPtCentersPi0.empty()) {
            // Create canvas
            TCanvas canvas;

            // Create multigraph
            TMultiGraph* mg = new TMultiGraph();

            // Create legend
            TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
            legend->SetTextSize(0.03);

            // Offset per trigger to avoid overlapping points
            double offsetStep = 0.1;
            int triggerIndex = 0;
            int numTriggers = triggersInData.size();

            for (const auto& triggerName : triggersInData) {
                // Get the data
                const std::vector<double>& pTCenters = triggerToPtCentersPi0[triggerName];
                const std::vector<double>& meanValues = triggerToMeanPi0Values[triggerName];
                const std::vector<double>& meanErrors = triggerToMeanPi0Errors[triggerName];

                // Create a new TGraphErrors
                TGraphErrors* graph = new TGraphErrors(pTCenters.size());
                for (size_t i = 0; i < pTCenters.size(); ++i) {
                    // Apply offset to pT center
                    double pTOffset = -offsetStep * (numTriggers - 1) / 2.0 + offsetStep * triggerIndex;
                    double pT = pTCenters[i] + pTOffset;
                    graph->SetPoint(i, pT, meanValues[i]);
                    graph->SetPointError(i, 0, meanErrors[i]);
                }

                // Set marker style and color
                int markerStyle = 20; // Use the same marker style for all triggers
                int markerColor = kBlack;
                // Get color from triggerColorMap if available
                auto it_color = TriggerConfig::triggerColorMap.find(triggerName);
                if (it_color != TriggerConfig::triggerColorMap.end()) {
                    markerColor = it_color->second;
                }
                graph->SetMarkerStyle(markerStyle);
                graph->SetMarkerSize(1);
                graph->SetLineWidth(2);
                graph->SetMarkerColor(markerColor);
                graph->SetLineColor(markerColor);

                // Add to multigraph
                mg->Add(graph, "P");

                // Add to legend
                std::string displayTriggerName = triggerName;
                if (TriggerConfig::triggerNameMap.find(triggerName) != TriggerConfig::triggerNameMap.end()) {
                    displayTriggerName = TriggerConfig::triggerNameMap.at(triggerName);
                }
                legend->AddEntry(graph, displayTriggerName.c_str(), "p");

                ++triggerIndex;
            }

            // Draw multigraph
            mg->Draw("A");
            mg->SetTitle(("Mean #pi^{0} Mass vs p_{T} (" + cutCombination + ")").c_str());
            mg->GetXaxis()->SetTitle("Leading Cluster p_{T} [GeV]");
            mg->GetYaxis()->SetTitle("Mean #pi^{0} Mass [GeV]");
            mg->GetXaxis()->SetLimits(2.0, 6.0);
            mg->GetXaxis()->SetNdivisions(505);
            mg->GetYaxis()->SetRangeUser(0.14, 0.19);

            // Draw legend
            legend->Draw();

            // Add cut combination information to the canvas
            TLatex labelText;
            labelText.SetNDC();
            labelText.SetTextSize(0.035);

            TLatex valueText;
            valueText.SetNDC();
            valueText.SetTextSize(0.042);

            labelText.DrawLatex(0.2, 0.87, "#font[62]{ECore #geq}");
            std::ostringstream eCoreWithUnit;
            eCoreWithUnit << clusECore << "   GeV";
            valueText.DrawLatex(0.35, 0.87, eCoreWithUnit.str().c_str());

            labelText.DrawLatex(0.2, 0.82, "#font[62]{#chi^{2} <}");
            std::ostringstream chiStr;
            chiStr << chi;
            valueText.DrawLatex(0.35, 0.82, chiStr.str().c_str());

            labelText.DrawLatex(0.2, 0.77, "#font[62]{Asymmetry <}");
            std::ostringstream asymmetryStr;
            asymmetryStr << asymmetry;
            valueText.DrawLatex(0.38, 0.77, asymmetryStr.str().c_str());

            // Ensure the directory exists
            std::string outputDirPath = plotDirectory + "/" + cutCombination;
            gSystem->mkdir(outputDirPath.c_str(), true);

            // Save plot
            std::string outputFilePath = outputDirPath + "/meanPi0_vs_pT_Overlay.png";
            canvas.SaveAs(outputFilePath.c_str());
            std::cout << "Saved overlay plot: " << outputFilePath << std::endl;

            // Clean up
            delete mg;
            delete legend;
        }

        // Generate overlay plot for η mass vs pT
        if (!triggerToPtCentersEta.empty()) {
            // Create canvas
            TCanvas canvas;

            // Create multigraph
            TMultiGraph* mgEta = new TMultiGraph();

            // Create legend
            TLegend* legendEta = new TLegend(0.6, 0.68, 0.9, 0.88);
            legendEta->SetTextSize(0.035);

            // Offset per trigger to avoid overlapping points
            double offsetStep = 0.1;
            int triggerIndex = 0;
            int numTriggers = triggersInData.size();

            for (const auto& triggerName : triggersInData) {
                // Get the data
                const std::vector<double>& pTCenters = triggerToPtCentersEta[triggerName];
                const std::vector<double>& meanValues = triggerToMeanEtaValues[triggerName];
                const std::vector<double>& meanErrors = triggerToMeanEtaErrors[triggerName];

                if (pTCenters.empty()) {
                    continue; // Skip if no data for this trigger
                }

                // Create a new TGraphErrors
                TGraphErrors* graph = new TGraphErrors(pTCenters.size());
                for (size_t i = 0; i < pTCenters.size(); ++i) {
                    // Apply offset to pT center
                    double pTOffset = -offsetStep * (numTriggers - 1) / 2.0 + offsetStep * triggerIndex;
                    double pT = pTCenters[i] + pTOffset;
                    graph->SetPoint(i, pT, meanValues[i]);
                    graph->SetPointError(i, 0, meanErrors[i]);
                }

                // Set marker style and color
                int markerStyle = 20; // Use the same marker style for all triggers
                int markerColor = kBlack;
                // Get color from triggerColorMap if available
                auto it_color = TriggerConfig::triggerColorMap.find(triggerName);
                if (it_color != TriggerConfig::triggerColorMap.end()) {
                    markerColor = it_color->second;
                }
                graph->SetMarkerStyle(markerStyle);
                graph->SetMarkerSize(1);
                graph->SetLineWidth(2);
                graph->SetMarkerColor(markerColor);
                graph->SetLineColor(markerColor);

                // Add to multigraph
                mgEta->Add(graph, "P");

                // Add to legend
                std::string displayTriggerName = triggerName;
                if (TriggerConfig::triggerNameMap.find(triggerName) != TriggerConfig::triggerNameMap.end()) {
                    displayTriggerName = TriggerConfig::triggerNameMap.at(triggerName);
                }
                legendEta->AddEntry(graph, displayTriggerName.c_str(), "p");

                ++triggerIndex;
            }

            // Draw multigraph
            mgEta->Draw("A");
            mgEta->SetTitle(("Mean #eta Mass vs p_{T} (" + cutCombination + ")").c_str());
            mgEta->GetXaxis()->SetTitle("Leading Cluster p_{T} [GeV]");
            mgEta->GetYaxis()->SetTitle("Mean #eta Mass [GeV]");
            mgEta->GetXaxis()->SetLimits(2.0, 6.0);
            mgEta->GetXaxis()->SetNdivisions(505);
            mgEta->GetYaxis()->SetRangeUser(0.45, 0.75);

            // Draw legend
            legendEta->Draw();

            // Add cut combination information to the canvas
            TLatex labelText;
            labelText.SetNDC();
            labelText.SetTextSize(0.035);

            TLatex valueText;
            valueText.SetNDC();
            valueText.SetTextSize(0.035);

            labelText.DrawLatex(0.2, 0.9, "#font[62]{ECore #geq}");
            std::ostringstream eCoreWithUnit;
            eCoreWithUnit << clusECore << "   GeV";
            valueText.DrawLatex(0.4, 0.9, eCoreWithUnit.str().c_str());

            labelText.DrawLatex(0.2, 0.85, "#font[62]{#chi^{2} <}");
            std::ostringstream chiStr;
            chiStr << chi;
            valueText.DrawLatex(0.4, 0.85, chiStr.str().c_str());

            labelText.DrawLatex(0.2, 0.8, "#font[62]{Asymmetry <}");
            std::ostringstream asymmetryStr;
            asymmetryStr << asymmetry;
            valueText.DrawLatex(0.4, 0.8, asymmetryStr.str().c_str());

            // Ensure the directory exists
            std::string outputDirPath = plotDirectory + "/" + cutCombination;
            gSystem->mkdir(outputDirPath.c_str(), true);

            // Save plot
            std::string outputFilePath = outputDirPath + "/meanEta_vs_pT_Overlay.png";
            canvas.SaveAs(outputFilePath.c_str());
            std::cout << "Saved η overlay plot: " << outputFilePath << std::endl;

            // Clean up
            delete mgEta;
            delete legendEta;
        }
    }
}


void PlotRunByRunHistograms(
    const std::string& outputDirectory,
    const std::string& plotDirectory,
    const std::vector<std::string>& triggers,
    const std::vector<int>& runNumbers,
    const std::map<std::string, int>& triggerColorMap,
    const std::map<std::string, std::string>& triggerNameMap) {

    // Create the run-by-run overlays directory
    std::string runByRunDir = plotDirectory + "/runByRun8by8overlays";
    gSystem->mkdir(runByRunDir.c_str(), true);

    // Loop over each run number
    for (int runNumber : runNumbers) {
        std::string runFileName = std::to_string(runNumber) + "_HistOutput.root";
        std::string runFilePath = outputDirectory + "/" + runFileName;

        // Open the ROOT file for the run
        TFile* runFile = TFile::Open(runFilePath.c_str(), "READ");
        if (!runFile || runFile->IsZombie()) {
            std::cerr << "Error: Could not open run file " << runFilePath << std::endl;
            continue;
        }

        // Create a canvas for the overlay plot
        TCanvas* canvas = new TCanvas("canvas", "Run-by-Run Overlay Plot", 800, 600);
        TLegend* legend = new TLegend(0.55, 0.4, 0.8, 0.8);
        legend->SetTextSize(0.035);
        canvas->SetLogy();

        bool firstDraw = true;

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

            int color = kBlack; // Default color
            auto it = TriggerConfig::triggerColorMap.find(trigger);
            if (it != TriggerConfig::triggerColorMap.end()) {
                color = it->second;
            }

            histClone->SetLineColor(color);
            histClone->SetLineWidth(2);

            std::string canvasTitle = "Run " + std::to_string(runNumber) + " Overlay";
            histClone->SetTitle(canvasTitle.c_str());
            histClone->GetXaxis()->SetTitle("Maximum 8x8 Energy Sum (EMCal) [GeV]");
            histClone->GetYaxis()->SetTitle("Events");

            if (firstDraw) {
                histClone->Draw("HIST");
                firstDraw = false;
            } else {
                histClone->Draw("HIST SAME");
            }

            // Use the map to translate the trigger name if available
            std::string displayTriggerName = trigger;
            if (triggerNameMap.find(trigger) != triggerNameMap.end()) {
                displayTriggerName = triggerNameMap.at(trigger);
            }

            // Add histogram to legend with display name
            legend->AddEntry(histClone, displayTriggerName.c_str(), "l");
        }

        legend->Draw();
        
        TLatex runNumberText;
        runNumberText.SetNDC(); // Use normalized device coordinates (0-1)
        runNumberText.SetTextAlign(31); // Align right and top
        runNumberText.SetTextSize(0.04); // Adjust text size as needed
        runNumberText.SetTextColor(kBlack); // Set text color

        // Construct the run number string
        std::ostringstream runNumberStr;
        runNumberStr << "Run Number: " << runNumber;

        // Draw the run number at the top right corner
        runNumberText.DrawLatex(0.9, 0.85, runNumberStr.str().c_str());

        // Update the canvas
        canvas->Modified();
        canvas->Update();

        // Save the canvas
        std::ostringstream outputFileName;
        outputFileName << runByRunDir << "/Run_" << runNumber << "_h8by8TowerEnergySum_Overlay.png";
        canvas->SaveAs(outputFileName.str().c_str());
        std::cout << "Saved run-by-run overlay plot to " << outputFileName.str() << std::endl;

        // Clean up
        delete canvas;
        runFile->Close();
        delete runFile;
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
        std::string chiLabel = "Chi2/NDF < " + Utils::formatToThreeSigFigs(chi);
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
    float isolatedYield_in;
    float isolatedYieldError_in;
    float isolatedYield_out;
    float isolatedYieldError_out;
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
    const std::map<std::string, double>& triggerEfficiencyPoints,
    const std::vector<std::pair<float, float>>& exclusionRanges,
    std::map<SpectraGroupKey, std::map<float, CombinedSpectraData>>& combinedSpectraDataMap
) {
    std::cout << "\033[1;34m[INFO]\033[0m Starting SortAndCombineSpectraData function.\n";

    // Step 1: Map each trigger to its priority based on TriggerConfig::allTriggers
    std::map<std::string, int> triggerPriorityMap;
    const std::vector<std::string>& allTriggers = TriggerConfig::allTriggers;
    for (size_t i = 0; i < allTriggers.size(); ++i) {
        triggerPriorityMap[allTriggers[i]] = static_cast<int>(i);
    }
    std::cout << "\033[1;34m[INFO]\033[0m Mapped triggers to their priorities based on TriggerConfig::allTriggers.\n";

    // Organize data into groups
    std::map<SpectraGroupKey, std::map<float, std::vector<std::tuple<std::string, float, DataStructures::IsolationData, DataStructures::IsolationData>>>> groupedData;

    std::cout << "\033[1;34m[INFO]\033[0m Organizing data into groups.\n";

    for (const auto& entry : dataMap_inMassWindow) {
        const auto& key = entry.first;
        const auto& isoData_in = entry.second;
        std::string triggerGroupName = std::get<0>(key);
        std::string triggerName = std::get<1>(key);
        float eCore = std::get<2>(key);
        float chi = std::get<3>(key);
        float asymmetry = std::get<4>(key);
        float ptMin = std::get<5>(key);
        float ptMax = std::get<6>(key);
        float pTCenter = (ptMin + ptMax) / 2.0;
        float isoMin = std::get<7>(key);
        float isoMax = std::get<8>(key);
        std::string massWindowLabel = std::get<9>(key);

        // Check for exclusion ranges
        std::pair<float, float> isoEtRange = {isoMin, isoMax};
        if (std::find(exclusionRanges.begin(), exclusionRanges.end(), isoEtRange) != exclusionRanges.end()) {
            std::cout << "\033[1;33m[WARNING]\033[0m Excluding isoEtRange [" << isoMin << ", " << isoMax << "] due to exclusionRanges.\n";
            continue; // Exclude this isoEtRange
        }

        SpectraGroupKey spectraGroupKey = {triggerGroupName, eCore, chi, asymmetry};

        // Find corresponding outsideMassWindow data
        auto outsideKey = key;
        std::get<9>(outsideKey) = "outsideMassWindow";

        auto it_outside = dataMap_outsideMassWindow.find(outsideKey);
        if (it_outside == dataMap_outsideMassWindow.end()) {
            std::cout << "\033[1;31m[ERROR]\033[0m No outsideMassWindow data found for key. Skipping this entry.\n";
            continue;
        }
        const auto& isoData_out = it_outside->second;

        // Add data to groupedData
        groupedData[spectraGroupKey][pTCenter].emplace_back(triggerName, pTCenter, isoData_in, isoData_out);
    }

    std::cout << "\033[1;34m[INFO]\033[0m Data organized into " << groupedData.size() << " groups.\n";

    // Now, for each group, sort triggers and combine data
    for (const auto& [spectraGroupKey, ptDataMap] : groupedData) {
        std::cout << "\033[1;34m[INFO]\033[0m Processing group: TriggerGroupName=" << spectraGroupKey.triggerGroupName
                  << ", ECore=" << spectraGroupKey.eCore << ", Chi=" << spectraGroupKey.chi
                  << ", Asymmetry=" << spectraGroupKey.asymmetry << "\n";

        // Get list of triggers in this group
        std::set<std::string> triggerSet;
        for (const auto& [pTCenter, dataVec] : ptDataMap) {
            for (const auto& dataEntry : dataVec) {
                const std::string& triggerName = std::get<0>(dataEntry);
                triggerSet.insert(triggerName);
            }
        }
        // Convert to vector and sort triggers
        std::vector<std::string> triggerList(triggerSet.begin(), triggerSet.end());

        std::sort(triggerList.begin(), triggerList.end(),
            [&](const std::string& a, const std::string& b) -> bool {
                bool aHasEff = triggerEfficiencyPoints.find(a) != triggerEfficiencyPoints.end();
                bool bHasEff = triggerEfficiencyPoints.find(b) != triggerEfficiencyPoints.end();

                if (aHasEff && bHasEff) {
                    double effA = triggerEfficiencyPoints.at(a);
                    double effB = triggerEfficiencyPoints.at(b);
                    if (effA != effB)
                        return effA > effB; // Higher efficiency threshold first
                    else
                        return triggerPriorityMap[a] < triggerPriorityMap[b];
                } else if (aHasEff) {
                    return true;  // a has efficiency, b does not
                } else if (bHasEff) {
                    return false; // b has efficiency, a does not
                } else {
                    return triggerPriorityMap[a] < triggerPriorityMap[b];
                }
            });

        std::cout << "\033[1;34m[INFO]\033[0m Sorted triggers in group: ";
        for (const auto& triggerName : triggerList) {
            std::cout << triggerName << " ";
        }
        std::cout << "\n";

        // For each pT bin, select appropriate trigger data
        for (const auto& [pTCenter, dataVec] : ptDataMap) {
            bool triggerAssigned = false;
            CombinedSpectraData combinedData;

            std::cout << "\033[1;34m[INFO]\033[0m Processing pT bin centered at " << pTCenter << " GeV.\n";

            for (const auto& triggerName : triggerList) {
                // Check efficiency threshold
                double efficiencyThreshold = 0.0;
                auto effIt = triggerEfficiencyPoints.find(triggerName);
                if (effIt != triggerEfficiencyPoints.end()) {
                    efficiencyThreshold = effIt->second;
                }
                if (pTCenter >= efficiencyThreshold) {
                    // Find data for this trigger
                    auto dataIt = std::find_if(dataVec.begin(), dataVec.end(),
                        [&](const auto& dataEntry) {
                            return std::get<0>(dataEntry) == triggerName;
                        });
                    if (dataIt != dataVec.end()) {
                        combinedData.pTCenter = pTCenter;
                        combinedData.isolatedYield_in = std::get<2>(*dataIt).isolatedYield;
                        combinedData.isolatedYieldError_in = std::get<2>(*dataIt).isolatedYieldError;
                        combinedData.isolatedYield_out = std::get<3>(*dataIt).isolatedYield;
                        combinedData.isolatedYieldError_out = std::get<3>(*dataIt).isolatedYieldError;
                        triggerAssigned = true;
                        std::cout << "\033[1;34m[INFO]\033[0m Assigned trigger '" << triggerName << "' for pT bin " << pTCenter << " GeV.\n";
                        break;
                    } else {
                        std::cout << "\033[1;33m[WARNING]\033[0m No data found for trigger '" << triggerName << "' at pT bin " << pTCenter << " GeV.\n";
                    }
                } else {
                    std::cout << "\033[1;33m[WARNING]\033[0m pT bin " << pTCenter << " GeV is below efficiency threshold " << efficiencyThreshold << " GeV for trigger '" << triggerName << "'.\n";
                }
            }
            if (triggerAssigned) {
                combinedSpectraDataMap[spectraGroupKey][pTCenter] = combinedData;
            } else {
                std::cout << "\033[1;31m[ERROR]\033[0m No suitable trigger assigned for pT bin " << pTCenter << " GeV in group '" << spectraGroupKey.triggerGroupName << "'.\n";
            }
        }
    }
    std::cout << "\033[1;34m[INFO]\033[0m Finished SortAndCombineSpectraData function.\n";
}


void GenerateCombinedSpectraPlots(
    const std::map<SpectraGroupKey, std::map<float, CombinedSpectraData>>& combinedSpectraDataMap,
    const std::string& basePlotDirectory,
    const std::map<std::string, std::string>& triggerCombinationNameMap
) {
    std::cout << "\033[1;34m[INFO]\033[0m Starting GenerateCombinedSpectraPlots function.\n";

    for (const auto& [spectraGroupKey, ptDataMap] : combinedSpectraDataMap) {
        const std::string& triggerGroupName = spectraGroupKey.triggerGroupName;
        float eCore = spectraGroupKey.eCore;
        float chi = spectraGroupKey.chi;
        float asym = spectraGroupKey.asymmetry;

        // Map to human-readable names
        std::string readableTriggerGroupName = Utils::getTriggerCombinationName(
            triggerGroupName, triggerCombinationNameMap);

        std::cout << "\033[1;34m[INFO]\033[0m Processing plot for Trigger Group: " << readableTriggerGroupName
                  << ", ECore > " << eCore << " GeV, Chi2/NDF < " << chi << ", Asymmetry < " << asym << ".\n";

        // Define output directory (same as per-trigger spectra plots)
        std::ostringstream dirStream;
        dirStream << basePlotDirectory << "/" << triggerGroupName
                  << "/E" << Utils::formatToThreeSigFigs(eCore)
                  << "_Chi" << Utils::formatToThreeSigFigs(chi)
                  << "_Asym" << Utils::formatToThreeSigFigs(asym)
                  << "/Spectra/Overlay";
        std::string dirPath = dirStream.str();
        if (gSystem->mkdir(dirPath.c_str(), true) != 0) {
            std::cout << "\033[1;31m[ERROR]\033[0m Failed to create directory: " << dirPath << "\n";
            continue;
        } else {
            std::cout << "\033[1;34m[INFO]\033[0m Created directory: " << dirPath << "\n";
        }

        // Create a TCanvas
        TCanvas* canvas = new TCanvas("canvas", "Combined Isolated Photon Spectra", 800, 600);
        canvas->SetLogy();

        // Prepare data vectors
        std::vector<double> pTValues, isolatedYields_in, isolatedYieldsError_in;
        std::vector<double> isolatedYields_out, isolatedYieldsError_out;

        for (const auto& [pTCenter, combinedData] : ptDataMap) {
            pTValues.push_back(combinedData.pTCenter);
            isolatedYields_in.push_back(combinedData.isolatedYield_in);
            isolatedYieldsError_in.push_back(combinedData.isolatedYieldError_in);
            isolatedYields_out.push_back(combinedData.isolatedYield_out);
            isolatedYieldsError_out.push_back(combinedData.isolatedYieldError_out);
        }

        if (pTValues.empty()) {
            std::cout << "\033[1;31m[ERROR]\033[0m No data points available for plotting in group '" << readableTriggerGroupName << "'. Skipping plot.\n";
            delete canvas;
            continue;
        }

        // Create TGraphErrors
        TGraphErrors* graphIn = new TGraphErrors(pTValues.size(),
                                                 pTValues.data(),
                                                 isolatedYields_in.data(),
                                                 nullptr,
                                                 isolatedYieldsError_in.data());
        graphIn->SetMarkerStyle(20);
        graphIn->SetMarkerColor(kBlue);
        graphIn->SetLineColor(kBlue);
        graphIn->SetLineWidth(2);

        TGraphErrors* graphOut = new TGraphErrors(pTValues.size(),
                                                  pTValues.data(),
                                                  isolatedYields_out.data(),
                                                  nullptr,
                                                  isolatedYieldsError_out.data());
        graphOut->SetMarkerStyle(21);
        graphOut->SetMarkerColor(kRed);
        graphOut->SetLineColor(kRed);
        graphOut->SetLineWidth(2);

        // Create a TMultiGraph
        TMultiGraph* multiGraph = new TMultiGraph();
        multiGraph->Add(graphIn, "P");
        multiGraph->Add(graphOut, "P");

        // Set titles and axis labels
        std::string plotTitle = "Combined Isolated Photon Spectra for Trigger Group: " + readableTriggerGroupName;
        multiGraph->SetTitle(plotTitle.c_str());
        multiGraph->GetXaxis()->SetTitle("Cluster p_{T} [GeV]");
        multiGraph->GetYaxis()->SetTitle("Yield");

        // Set Y-axis range
        double globalMaxY_in = *std::max_element(isolatedYields_in.begin(), isolatedYields_in.end());
        double globalMaxY_out = *std::max_element(isolatedYields_out.begin(), isolatedYields_out.end());
        double globalMaxY = std::max(globalMaxY_in, globalMaxY_out);

        if (globalMaxY <= 0.0) {
            std::cout << "\033[1;31m[ERROR]\033[0m Maximum Y value is non-positive for group '" << readableTriggerGroupName << "'. Skipping plot.\n";
            delete graphIn;
            delete graphOut;
            delete multiGraph;
            delete canvas;
            continue;
        }

        multiGraph->SetMinimum(1e-1);
        multiGraph->SetMaximum(globalMaxY * 1.5);

        // Draw the multigraph
        multiGraph->Draw("A");

        // Create and draw legend
        TLegend* legend = new TLegend(0.38, 0.78, 0.7, 0.9);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.024);
        legend->AddEntry(graphIn, "Isolated Photons from Meson Decay Yield", "p");
        legend->AddEntry(graphOut, "Prompt Photon Candidate Yield", "p");
        legend->Draw();

        // Add labels using TLatex
        TLatex labelText;
        labelText.SetNDC();
        labelText.SetTextSize(0.025);
        labelText.SetTextColor(kBlack);

        double xStart = 0.55;
        double yStartLabel = 0.72;
        double yStepLabel = 0.05;

        // Prepare label strings
        std::string triggerGroupLabel = "Trigger Group: " + readableTriggerGroupName;
        std::string eCoreLabel = "ECore > " + Utils::formatToThreeSigFigs(eCore) + " GeV";
        std::string chiLabel = "Chi2/NDF < " + Utils::formatToThreeSigFigs(chi);
        std::string asymLabel = "Asymmetry < " + Utils::formatToThreeSigFigs(asym);

        // Draw labels
        labelText.DrawLatex(xStart, yStartLabel, triggerGroupLabel.c_str());
        labelText.DrawLatex(xStart, yStartLabel - yStepLabel, eCoreLabel.c_str());
        labelText.DrawLatex(xStart, yStartLabel - 2 * yStepLabel, chiLabel.c_str());
        labelText.DrawLatex(xStart, yStartLabel - 3 * yStepLabel, asymLabel.c_str());

        // Update canvas and save
        canvas->Modified();
        canvas->Update();

        std::string outputFilePath = dirPath + "/CombinedOverlaySpectra.png";
        canvas->SaveAs(outputFilePath.c_str());
        std::cout << "\033[1;34m[INFO]\033[0m Saved combined overlay spectra plot to " << outputFilePath << std::endl;

        // Clean up
        delete legend;
        delete graphIn;
        delete graphOut;
        delete multiGraph;
        delete canvas;
    }

    std::cout << "\033[1;34m[INFO]\033[0m Finished GenerateCombinedSpectraPlots function.\n";
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
        TLegend legend(0.18, 0.75, 0.48, 0.9); // Adjust as needed
        legend.SetBorderSize(0);
        legend.SetTextSize(0.024);

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
        labelText.SetTextSize(0.024);       // Adjusted text size
        labelText.SetTextColor(kBlack);    // Ensured text color is black
        double xStart = 0.7; // Starting x-coordinate (left side)
        double yStartLabel = 0.9; // Starting y-coordinate
        double yStepLabel = 0.035;  // Vertical spacing between lines

        // Prepare label strings
        std::string triggerGroupLabel = "Trigger Group: " + readableTriggerGroupName;
        std::string triggerNameLabel = "Trigger: " + readableTriggerName;
        std::string eCoreLabel = "ECore > " + Utils::formatToThreeSigFigs(eCore) + " GeV";
        std::string chiLabel = "Chi2/NDF < " + Utils::formatToThreeSigFigs(chi);
        std::string asymLabel = "Asymmetry < " + Utils::formatToThreeSigFigs(asym);
        std::string massWindowLabelStr = "Mass Window: " + massWindowLabel;
        std::string coneRadiusLabel = "#Delta R_{cone} < 0.3";

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
    const std::map<std::string, double>& triggerEfficiencyPoints,
    std::map<std::string, std::vector<std::string>>& sortedTriggersByGroupName,
    std::map<std::string, std::map<std::pair<float, float>,
    std::vector<DataStructures::IsolationDataWithPt>>>& combinedTriggerDataMap) {
    
    // Step 1: Map each trigger to its priority based on TriggerConfig::allTriggers
    std::map<std::string, int> triggerPriorityMap;
    const std::vector<std::string>& allTriggers = TriggerConfig::allTriggers;
    for (size_t i = 0; i < allTriggers.size(); ++i) {
        triggerPriorityMap[allTriggers[i]] = static_cast<int>(i);
    }

    // Step 2: Populate sortedTriggersByGroupName without duplication
    std::cout << "\033[34m[INFO]\033[0m Populating sortedTriggersByGroupName...\n";
    for (const auto& [groupKey, isoEtMap] : groupedData) {
        const std::string& triggerGroupName = std::get<0>(groupKey);
        const std::string& triggerName = std::get<1>(groupKey);

        // Initialize the group if it doesn't exist
        if (sortedTriggersByGroupName.find(triggerGroupName) == sortedTriggersByGroupName.end()) {
            sortedTriggersByGroupName[triggerGroupName] = {};
        }

        // Add the trigger to the group if it's not already present
        if (std::find(sortedTriggersByGroupName[triggerGroupName].begin(),
                      sortedTriggersByGroupName[triggerGroupName].end(),
                      triggerName) == sortedTriggersByGroupName[triggerGroupName].end()) {
            sortedTriggersByGroupName[triggerGroupName].push_back(triggerName);
            std::cout << "\033[34m[INFO]\033[0m Added trigger '" << triggerName
                      << "' to group '" << triggerGroupName << "'\n";
        }
    }

    // Step 3: Sort triggers within each group based on their efficiency thresholds in descending order
    std::cout << "\033[34m[INFO]\033[0m Sorting triggers within each group based on efficiency thresholds...\n";
    for (auto& [triggerGroupName, triggerList] : sortedTriggersByGroupName) {
        std::sort(triggerList.begin(), triggerList.end(),
                  [&](const std::string& a, const std::string& b) -> bool {
                      bool aHasEff = triggerEfficiencyPoints.find(a) != triggerEfficiencyPoints.end();
                      bool bHasEff = triggerEfficiencyPoints.find(b) != triggerEfficiencyPoints.end();
                      
                      if (aHasEff && bHasEff) {
                          double effA = triggerEfficiencyPoints.at(a);
                          double effB = triggerEfficiencyPoints.at(b);
                          if (effA != effB)
                              return effA > effB; // Higher efficiency threshold first
                          else
                              return triggerPriorityMap[a] < triggerPriorityMap[b];
                      } else if (aHasEff) {
                          return true;  // a has efficiency, b does not
                      } else if (bHasEff) {
                          return false; // b has efficiency, a does not
                      } else {
                          return triggerPriorityMap[a] < triggerPriorityMap[b];
                      }
                  });
        
        // Debugging output
        std::cout << "\033[1;32mTrigger Group Name\033[0m: " << triggerGroupName << "\n";
        std::cout << "\033[1;32mSorted Trigger List\033[0m: ";
        for (const auto& trigger : triggerList) {
            if (triggerEfficiencyPoints.find(trigger) != triggerEfficiencyPoints.end()) {
                std::cout << trigger << " (Eff Threshold: " << triggerEfficiencyPoints.at(trigger) << "), ";
            } else {
                std::cout << trigger << " (Eff Threshold: N/A), ";
            }
        }
        std::cout << "\n";
    }



    // Step 4: Combine triggers for each group and isoEtRange
    std::cout << "\033[34m[INFO]\033[0m Combining triggers within each group...\n";

    for (const auto& [triggerGroupName, sortedTriggerList] : sortedTriggersByGroupName) {
        std::cout << "\033[33m[PROCESSING]\033[0m Combining triggers for group: "
                  << "\033[1m" << triggerGroupName << "\033[0m\n";

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
            std::cerr << "\033[31m[ERROR]\033[0m No isoEt ranges found for group '"
                      << triggerGroupName << "'\n";
            continue;
        }

        // Iterate over each isoEtRange
        for (const auto& isoEtRange : allIsoEtRanges) {
            std::cout << "\033[34m[INFO]\033[0m Processing isoEtRange: [" << isoEtRange.first
                      << ", " << isoEtRange.second << "]\n";

            // Collect all unique pT bins across triggers for this isoEtRange
            std::set<std::pair<float, float>> allPtBins;
            for (const auto& triggerName : sortedTriggerList) {
                // Find the corresponding groupKey
                bool foundGroupKey = false;
                GroupKey currentGroupKey;
                for (const auto& [gk, isoEtMap] : groupedData) {
                    if (std::get<0>(gk) == triggerGroupName && std::get<1>(gk) == triggerName) {
                        currentGroupKey = gk;
                        foundGroupKey = true;
                        break;
                    }
                }
                if (!foundGroupKey) {
                    std::cerr << "\033[31m[ERROR]\033[0m Group key not found for trigger '"
                              << triggerName << "'\n";
                    continue;
                }

                // Get isoDataList for this isoEtRange
                auto isoIt = groupedData.at(currentGroupKey).find(isoEtRange);
                if (isoIt != groupedData.at(currentGroupKey).end()) {
                    for (const auto& isoData : isoIt->second) {
                        allPtBins.emplace(std::make_pair(isoData.ptMin, isoData.ptMax));
                    }
                }
            }

            if (allPtBins.empty()) {
                std::cerr << "\033[31m[ERROR]\033[0m No pT bins found for isoEt range: ["
                          << isoEtRange.first << ", " << isoEtRange.second << "]\n";
                continue;
            }

            std::cout << "\033[34m[INFO]\033[0m Found " << allPtBins.size() << " pT bins\n";

            // Iterate over each pT bin and select appropriate trigger data
            std::vector<DataStructures::IsolationDataWithPt> selectedDataPoints;

            // Iterate over each pT bin and select appropriate trigger data
            for (const auto& ptBin : allPtBins) {
                std::cout << "\033[32m[DEBUG]\033[0m Processing pT bin: [" << ptBin.first
                          << ", " << ptBin.second << "]\n";
                double pTCenter = (ptBin.first + ptBin.second) / 2.0;
                bool triggerAssigned = false;
                DataStructures::IsolationDataWithPt selectedIsoData;

                std::cout << "\033[32m[DEBUG]\033[0m pTCenter: " << pTCenter << " GeV\n";

                // Iterate triggers in sorted order (now descending efficiency thresholds)
                for (const auto& triggerName : sortedTriggerList) {
                    // Declare efficiencyThreshold at the beginning of the loop
                    double efficiencyThreshold = 0.0; // Initialize to 0.0 or another appropriate default

                    // When checking for efficiency thresholds
                    auto effIt = triggerEfficiencyPoints.find(triggerName);
                    if (effIt == triggerEfficiencyPoints.end()) {
                        std::cout << "\033[31m[WARNING]\033[0m No efficiency threshold for trigger '"
                                  << triggerName << "'. Assigning default value 0.0.\n";
                        // efficiencyThreshold is already initialized to 0.0
                    } else {
                        efficiencyThreshold = effIt->second;
                    }

                    std::cout << "\033[32m[DEBUG]\033[0m Evaluating trigger '" << triggerName
                              << "' with efficiency threshold " << efficiencyThreshold
                              << " GeV against pTCenter " << pTCenter << " GeV\n";


                    if (pTCenter >= efficiencyThreshold) {
                        // Find the isoData for this trigger and pT bin
                        bool foundGroupKey = false;
                        GroupKey currentGroupKey;
                        for (const auto& [gk, isoEtMap] : groupedData) {
                            if (std::get<0>(gk) == triggerGroupName && std::get<1>(gk) == triggerName) {
                                currentGroupKey = gk;
                                foundGroupKey = true;
                                break;
                            }
                        }
                        if (!foundGroupKey) {
                            std::cerr << "\033[31m[ERROR]\033[0m Group key not found for trigger '"
                                      << triggerName << "'\n";
                            continue;
                        }

                        // Get the isoDataList for this isoEtRange
                        auto isoIt = groupedData.at(currentGroupKey).find(isoEtRange);
                        if (isoIt != groupedData.at(currentGroupKey).end()) {
                            // Find the isoData that contains pTCenter
                            auto dataIt = std::find_if(isoIt->second.begin(), isoIt->second.end(),
                                [&](const DataStructures::IsolationDataWithPt& id) {
                                    return id.ptMin <= pTCenter && pTCenter < id.ptMax;
                                });
                            if (dataIt != isoIt->second.end()) {
                                selectedIsoData = *dataIt;
                                triggerAssigned = true;
                                std::cout << "\033[32m[DEBUG]\033[0m Assigned to trigger '"
                                          << triggerName << "'\n";
                                break; // Trigger selected for this pT bin
                            } else {
                                std::cout << "\033[31m[DEBUG]\033[0m Data not found for trigger '"
                                          << triggerName << "' with pT bin ["
                                          << ptBin.first << ", " << ptBin.second << "]\n";
                            }
                        }
                    } else {
                        std::cout << "\033[33m[DEBUG]\033[0m pTCenter " << pTCenter
                                  << " GeV is below threshold for trigger '" << triggerName << "'\n";
                    }
                }
                // If no trigger met the efficiency threshold, use "Minbias" if available
                if (!triggerAssigned) {
                    const std::string minbiasTrigger = "MBD_NandS_geq_1";
                    if (std::find(sortedTriggerList.begin(), sortedTriggerList.end(), minbiasTrigger) != sortedTriggerList.end()) {
                        // Find the isoData for "Minbias"
                        bool foundGroupKey = false;
                        GroupKey currentGroupKey;
                        for (const auto& [gk, isoEtMap] : groupedData) {
                            if (std::get<0>(gk) == triggerGroupName && std::get<1>(gk) == minbiasTrigger) {
                                currentGroupKey = gk;
                                foundGroupKey = true;
                                break;
                            }
                        }
                        if (foundGroupKey) {
                            // Get the isoDataList for this isoEtRange
                            auto isoIt = groupedData.at(currentGroupKey).find(isoEtRange);
                            if (isoIt != groupedData.at(currentGroupKey).end()) {
                                // Find the isoData with this pT bin
                                auto dataIt = std::find_if(isoIt->second.begin(), isoIt->second.end(),
                                    [&](const DataStructures::IsolationDataWithPt& id) {
                                        return std::abs(id.ptMin - ptBin.first) < 1e-6 && std::abs(id.ptMax - ptBin.second) < 1e-6;
                                    });
                                if (dataIt != isoIt->second.end()) {
                                    selectedIsoData = *dataIt;
                                    triggerAssigned = true;
                                    std::cout << "\033[32m[DEBUG]\033[0m Assigned to Minbias trigger '"
                                              << minbiasTrigger << "'\n";
                                }
                            }
                        }
                    }
                }

                // If still not assigned, skip this pT bin
                if (triggerAssigned) {
                    selectedDataPoints.push_back(selectedIsoData);
                } else {
                    std::cout << "\033[31m[WARNING]\033[0m No suitable trigger found for pT bin ["
                              << ptBin.first << ", " << ptBin.second << "] in isoEtRange ["
                              << isoEtRange.first << ", " << isoEtRange.second << "]\n";
                }
            }

            // Assign selectedDataPoints to combinedTriggerDataMap
            combinedTriggerDataMap[triggerGroupName][isoEtRange] = selectedDataPoints;

            // Debugging output for combined data
            std::cout << "\033[36m[DEBUG]\033[0m Combined data points for group '"
                      << triggerGroupName << "', isoEtRange [" << isoEtRange.first << ", "
                      << isoEtRange.second << "]: " << selectedDataPoints.size() << " points\n";
        }
    }

    std::cout << "\033[34m[INFO]\033[0m Trigger sorting and combining completed.\n";
}

void GenerateCombinedRatioPlot(
    const std::map<std::string, std::map<std::pair<float, float>, std::vector<DataStructures::IsolationDataWithPt>>>& combinedTriggerDataMap,
    const std::map<GroupKey, std::map<std::pair<float, float>, std::vector<DataStructures::IsolationDataWithPt>>>& groupedData,
    const std::string& basePlotDirectory) {
    // -----------------------------
    // ** Combined Plot Generation **
    // -----------------------------
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
            // Assuming all dataPoints have the same groupKey attributes
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

            // Map triggerGroupName and triggerName to human-readable names
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

            // Create a TMultiGraph for the combined data
            TMultiGraph* combinedMultiGraph = new TMultiGraph();

            // Create a legend
            TLegend combinedLegend(0.18, 0.75, 0.48, 0.9);
            combinedLegend.SetBorderSize(0);
            combinedLegend.SetTextSize(0.024);

            // Group data points by triggerName for coloring
            std::map<std::string, std::vector<DataStructures::IsolationDataWithPt>> dataByTrigger;
            for (const auto& isoData : dataPoints) {
                dataByTrigger[isoData.triggerName].push_back(isoData);
            }

            // Iterate over each trigger's data points
            for (const auto& [triggerName, triggerDataPoints] : dataByTrigger) {
                std::vector<double> weightedPts;
                std::vector<double> ratios;
                std::vector<double> errors;

                for (const auto& isoData : triggerDataPoints) {
                    weightedPts.push_back(isoData.weightedPt);
                    ratios.push_back(isoData.ratio);
                    errors.push_back(isoData.error);
                }
                if (weightedPts.empty()) {
                    std::cerr << "\033[31m[WARNING]\033[0m No valid data points for trigger: " << triggerName
                              << ". Skipping.\n";
                    continue;
                }

                // Create a TGraphErrors for this trigger
                TGraphErrors* graph = new TGraphErrors(weightedPts.size(),
                                                      weightedPts.data(),
                                                      ratios.data(),
                                                      nullptr,
                                                      errors.data());
                graph->SetMarkerStyle(20);
                // Assign color based on triggerName
                auto colorIt = TriggerConfig::triggerColorMap.find(triggerName);
                if (colorIt != TriggerConfig::triggerColorMap.end()) {
                    graph->SetMarkerColor(colorIt->second);
                    graph->SetLineColor(colorIt->second);
                } else {
                    graph->SetMarkerColor(kBlack);
                    graph->SetLineColor(kBlack);
                }
                graph->SetLineWidth(2);

                // Add to the multigraph
                combinedMultiGraph->Add(graph, "P");

                // Add entry to legend
                std::string readableTriggerName = triggerName;
                auto triggerNameIt = TriggerConfig::triggerNameMap.find(triggerName);
                if (triggerNameIt != TriggerConfig::triggerNameMap.end()) {
                    readableTriggerName = triggerNameIt->second;
                }
                combinedLegend.AddEntry(graph, readableTriggerName.c_str(), "p");
                
                std::cout << "\033[32m[DEBUG]\033[0m Added data for trigger: " << triggerName
                          << " (" << readableTriggerName << ") with " << weightedPts.size() << " points.\n";
            }

            // Set titles
            std::string plotTitle = "Isolation Ratio vs pT (" + massWindowLabel + ")";
            combinedMultiGraph->SetTitle(plotTitle.c_str());
            combinedMultiGraph->GetXaxis()->SetTitle("Weighted p_{T} of Leading Cluster [GeV]");
            const std::string yAxisTitle = (massWindowLabel == "inMassWindow") ?
                "#frac{Isolated Photons from #pi^{0}/#eta Decays}{All Photons from #pi^{0}/#eta Decays}" :
                "#frac{Isolated Prompt Photons}{All Prompt Photons}";
            
            combinedMultiGraph->GetYaxis()->SetTitle(yAxisTitle.c_str());
            combinedMultiGraph->GetYaxis()->SetRangeUser(0, 2.0);
            combinedMultiGraph->GetXaxis()->SetLimits(2.0, 25.0);

            // Draw multigraph
            combinedMultiGraph->Draw("A P");

            // Draw a dashed line at y = 1
            TLine* combinedLine = new TLine(2, 1, 25, 1);
            combinedLine->SetLineStyle(2); // Dashed line
            combinedLine->Draw();

            // Draw the legend
            combinedLegend.Draw();

            // Add labels using TLatex in the top-left corner
            TLatex labelText;
            labelText.SetNDC();
            labelText.SetTextSize(0.026);       // Adjust text size as needed
            labelText.SetTextColor(kBlack);    // Ensure text color is black
            
            double xStart = 0.55; // Starting x-coordinate (left side)
            double yStartLabel = 0.9; // Starting y-coordinate
            double yStepLabel = 0.04;  // Vertical spacing between lines
            
            // Prepare label strings
            std::string triggerGroupLabel = "Trigger Group: " + readableTriggerGroupName;
            std::string eCoreLabel = "ECore > " + Utils::formatToThreeSigFigs(eCore) + " GeV";
            std::string chiLabel = "Chi2/NDF < " + Utils::formatToThreeSigFigs(chi);
            std::string asymLabel = "Asymmetry < " + Utils::formatToThreeSigFigs(asym);
            std::string massWindowLabelStr = "Mass Window: " + massWindowLabel;
            std::string coneRadiusLabel = "#Delta R_{cone} < 0.3";
            // Draw labels
            labelText.DrawLatex(xStart, yStartLabel, triggerGroupLabel.c_str());
            labelText.DrawLatex(xStart, yStartLabel - yStepLabel, eCoreLabel.c_str());
            labelText.DrawLatex(xStart, yStartLabel - 2 * yStepLabel, chiLabel.c_str());
            labelText.DrawLatex(xStart, yStartLabel - 3 * yStepLabel, asymLabel.c_str());
            labelText.DrawLatex(xStart, yStartLabel - 4 * yStepLabel, massWindowLabelStr.c_str());
            labelText.DrawLatex(xStart, yStartLabel - 5 * yStepLabel, coneRadiusLabel.c_str());

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
            delete combinedMultiGraph;
            delete combinedLine;
            // The graphs are managed by TMultiGraph and ROOT's memory management
        }
    }
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
    const std::vector<std::pair<float, float>>& exclusionRanges = {},
    const std::map<std::string, double>& triggerEfficiencyPoints = {},
    bool drawRefA = false,
    bool drawRefB = false) {
    std::map<SpectraGroupKey, std::map<float, CombinedSpectraData>> combinedSpectraDataMap;

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
    
    SortAndCombineTriggers(groupedData, triggerEfficiencyPoints, sortedTriggersByGroupName, combinedTriggerDataMap);

    GenerateCombinedRatioPlot(combinedTriggerDataMap, groupedData, basePlotDirectory);
    
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

    // Sort and combine spectra data
    SortAndCombineSpectraData(
        dataMap_inMassWindow,
        dataMap_outsideMassWindow,
        triggerEfficiencyPoints,
        exclusionRanges,
        combinedSpectraDataMap
    );

    // Call GenerateCombinedSpectraPlots
    GenerateCombinedSpectraPlots(
        combinedSpectraDataMap,
        basePlotDirectory,
        TriggerCombinationNames::triggerCombinationNameMap
    );
    std::cout << "\033[33m[INFO]\033[0m Finished processing isolation data." << std::endl;
}



void PlotCombinedHistograms(
    const std::string& outputDirectory,
    const std::vector<std::string>& combinedRootFiles,
    const std::map<std::string, std::vector<int>>& combinationToValidRuns) {
    
    // List of all triggers
    const std::vector<std::string>& allTriggers = TriggerConfig::allTriggers;
    const std::vector<std::string>& photonTriggers = TriggerConfig::photonTriggers;
    const auto& triggerColorMap = TriggerConfig::triggerColorMap;
    const auto& triggerNameMap = TriggerConfig::triggerNameMap;
    std::map<std::string, double> triggerEfficiencyPoints;
    
    // Base directory for plots
    std::string basePlotDirectory = "/Users/patsfan753/Desktop/DirectPhotonAna/Plots";

    // Create the base plot directory if it doesn't exist
    gSystem->mkdir(basePlotDirectory.c_str(), true);

    // Function to draw run numbers on canvas
    auto drawRunNumbersOnCanvas = [](const std::vector<int>& runNumbers) {
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

        // Define a separate font size for the header
        double headerTextSize = 0.045;

        // Draw the header with a different font size
        runNumbersLatex.SetTextSize(headerTextSize);
        std::ostringstream headerText;
        headerText << "Run Numbers (" << runNumbers.size() << " runs)";
        runNumbersLatex.DrawLatex(0.5, 0.9, headerText.str().c_str());

        // Skip plotting run numbers for specific size
        if (runNumbers.size() == 697 || runNumbers.size() == 700 || runNumbers.size() == 585 || runNumbers.size() == 728 || runNumbers.size() == 104) {
            std::cout << "[INFO] Skipping run number plotting for runNumbers.size() = 692." << std::endl;
            return;
        }

        // Reset the text size for the run numbers
        runNumbersLatex.SetTextSize(textSize);

        // Calculate the number of rows based on columns
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

        std::cout << "Processing file: " << rootFileName << std::endl;
        std::cout << "Triggers found: ";
        for (const auto& trigger : triggers) {
            std::cout << trigger << " ";
        }
        std::cout << std::endl;

        // Create a combination name for the folder
        std::string combinationName;
        for (const auto& trigger : triggers) {
            combinationName += trigger + "_";
        }
        // Remove the trailing underscore
        if (!combinationName.empty()) {
            combinationName.pop_back();
        }

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

        // -------------------- Overlay Plot --------------------
        // Create a canvas
        TCanvas* canvas = new TCanvas("canvas", "Overlay Plot", 800, 600);
        TLegend* legend = new TLegend(0.65, 0.53, 0.85, 0.88);
        legend->SetTextSize(0.03);
        canvas->SetLogy();

        bool firstDraw = true;

        // Loop over triggers and plot histograms
        for (const auto& trigger : triggers) {
            // Get the trigger directory
            TDirectory* triggerDir = inputFile->GetDirectory(trigger.c_str());
            if (!triggerDir) {
                std::cerr << "Trigger directory '" << trigger << "' not found in file " << rootFileName << std::endl;
                continue;
            }
            
            // Construct histogram name
            std::string histName = "h8by8TowerEnergySum_" + trigger;

            // Get histogram from file
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
            histClone->GetXaxis()->SetTitle("Maximum 8x8 Energy Sum (EMCal) [GeV]");
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

        // Display run numbers on the plot
        auto it = combinationToValidRuns.find(combinationName);
        if (it != combinationToValidRuns.end()) {
            const std::vector<int>& validRuns = it->second;
            // Draw run numbers on canvas using the drawRunNumbersOnCanvas function
            drawRunNumbersOnCanvas(validRuns);
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
                    TLegend* legendTurnOn = new TLegend(0.18, 0.75, 0.45, 0.9);
                    legendTurnOn->SetTextSize(0.025);
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
                        
                        // Perform the fit and append the parameters to the legend entry
                        std::pair<std::string, std::string> key = std::make_pair(combinationName, photonTrigger);
                        auto it_fitParams = TriggerConfig::triggerFitParameters.find(key);
                        if (it_fitParams == TriggerConfig::triggerFitParameters.end()) {
                            // Try with empty combinationName
                            key = std::make_pair("", photonTrigger);
                            it_fitParams = TriggerConfig::triggerFitParameters.find(key);
                        }

                        if (enableFits && it_fitParams != TriggerConfig::triggerFitParameters.end()) {
                            DataStructures::FitParameters params = it_fitParams->second;
                            TF1* fitFunc = Utils::sigmoidFit(("fit_" + photonTrigger).c_str(), 0.0, 20.0,
                                                      params.amplitudeEstimate, params.slopeEstimate, params.xOffsetEstimate,
                                                      params.amplitudeMin, params.amplitudeMax,
                                                      params.slopeMin, params.slopeMax,
                                                      params.xOffsetMin, params.xOffsetMax);
                            fitFunc->SetLineColor(color);
                            ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
                            ratioHist->Fit(fitFunc, "R");
                            fitFunc->Draw("SAME");
                            
                            // Retrieve the fit parameters regardless of convergence
                            double A = fitFunc->GetParameter(0);
                            double k = fitFunc->GetParameter(1);
                            double x0 = fitFunc->GetParameter(2);
                            
                            double A_error = fitFunc->GetParError(0); // Error on Amplitude
                            double k_error = fitFunc->GetParError(1); // Error on Slope
                            double x0_error = fitFunc->GetParError(2); // Error on XOffset
                            
                            // Calculate the 99% efficiency point and its error
                            double x99 = x0 + (std::log(99) / k);
                            
                            // Propagate error for x99 using partial derivatives
                            double x99_error = std::sqrt(
                                                         (x0_error * x0_error) +
                                                         (std::pow((std::log(99) / (k * k)), 2) * k_error * k_error)
                                                         );
                            
                            // Store the 99% efficiency point
                            triggerEfficiencyPoints[photonTrigger] = x99;
                            
                            // Append fit parameters to the legend entry with errors
                            legendEntry << ", Amp = " << std::fixed << std::setprecision(2) << A
                            << ", x(y = 0.99) = " << std::fixed << std::setprecision(2) << x99 << " GeV";
                            
                            // Debugging output to confirm correct parameter retrieval
                            std::cout << "Fit for " << photonTrigger
                            << ": A = " << A << " ± " << A_error
                            << ", k = " << k << " ± " << k_error
                            << ", x0 = " << x0 << " ± " << x0_error
                            << ", x99 = " << x99 << " ± " << x99_error << " GeV" << std::endl;
                            
                            
                            // Draw a dashed vertical line at x99
                            if (x99 > ratioHist->GetXaxis()->GetXmin() && x99 < ratioHist->GetXaxis()->GetXmax()) {
                                TLine* verticalLine = new TLine(x99, 0, x99, 1); // Draw line up to y = 1
                                verticalLine->SetLineStyle(2); // Dashed line
                                verticalLine->SetLineColor(color); // Use the color of the current trigger
                                verticalLine->SetLineWidth(2); // Set the line width (3 is thicker)
                                verticalLine->Draw("SAME");
                            }
                        }
                        
                        // Add the updated legend entry with the fit parameters
                        legendTurnOn->AddEntry(ratioHist, legendEntry.str().c_str(), "p");
                        std::cout << "Legend Entry Added: " << legendEntry.str() << std::endl;  // Debugging statement
                    }
                    
                    // Draw the legend
                    legendTurnOn->Draw();
                    
                    // Add a separate legend for the 99% efficiency line
                    TLegend* legendEfficiencyLine = new TLegend(0.18, 0.62, 0.38, 0.72); // Adjust position as needed
                    legendEfficiencyLine->SetTextSize(0.03);
                    legendEfficiencyLine->SetBorderSize(0);
                    legendEfficiencyLine->SetFillStyle(0);
                    
                    // Create a dummy line for the 99% efficiency point
                    TLine* dummyLine = new TLine(0, 0, 0, 0); // The coordinates do not matter since it is a dummy
                    dummyLine->SetLineStyle(2); // Dashed line
                    dummyLine->SetLineColor(kGray + 1); // Gray color
                    dummyLine->SetLineWidth(2);
                    
                    // Add the dummy line to the legend
                    legendEfficiencyLine->AddEntry(dummyLine, "99% Efficiency Point", "l"); // "l" for line
                    legendEfficiencyLine->Draw();
                    
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
        ProcessInvariantMassHistograms(inputFile, plotDirectory, triggers, triggerColorMap, combinationName);

        ProcessMesonMassVsPt(plotDirectory, combinationName, triggers, triggerEfficiencyPoints);
        
        ProcessIsolationEnergyHistogramsWithCuts(inputFile, plotDirectory, triggers, combinationName);

        
        inputFile->Close();
        delete inputFile;
        
        // Now, generate run-by-run overlays
        auto it_run = combinationToValidRuns.find(combinationName);
        if (it_run != combinationToValidRuns.end()) {
            const std::vector<int>& validRuns = it_run->second;
            PlotRunByRunHistograms(outputDirectory, plotDirectory, triggers, validRuns, triggerColorMap, triggerNameMap);
        } else {
            std::cout << "No valid runs found for combination: " << combinationName << std::endl;
        }
    }
    // After processing all histograms, write the CSV file
    std::string csvOutputPath = "/Users/patsfan753/Desktop/isolation_data.csv";
    WriteIsolationDataToCSV(csvOutputPath);
    readDataFromCSV(csvOutputPath, dataMap_inMassWindow, dataMap_outsideMassWindow);
    
    std::vector<std::pair<float, float>> exclusionRanges = {
    };

    ProcessIsolationData(dataMap_inMassWindow, basePlotDirectory, exclusionRanges, triggerEfficiencyPoints, false, true);
    ProcessIsolationData(dataMap_outsideMassWindow, basePlotDirectory, exclusionRanges, triggerEfficiencyPoints, true, false);
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


void AnalyzeTriggerGroupings() {
    gROOT->LoadMacro("sPhenixStyle.C");
    SetsPhenixStyle();

    std::string csvFilePath = "/Users/patsfan753/Desktop/DirectPhotonAna/triggerAnalysisCombined.csv";

    std::string outputDirectory = "/Users/patsfan753/Desktop/DirectPhotonAna/output";

    // Get the map of trigger combinations to run numbers
    std::map<std::set<std::string>, std::vector<int>> combinationToRuns = AnalyzeWhatTriggerGroupsAvailable(csvFilePath);

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

        // Create a string to represent the combination for naming
        std::string combinationName;
        for (const auto& trigger : triggers) {
            combinationName += trigger + "_";
        }
        // Remove the trailing underscore
        if (!combinationName.empty()) {
            combinationName.pop_back();
        }

        // Define the final output ROOT file path
        std::string finalRootFilePath = outputDirectory + "/" + combinationName + "_Combined.root";
        // Define the valid runs text file path
        std::string validRunsFilePath = outputDirectory + "/" + combinationName + "_ValidRuns.txt";

        // Check if both files exist
        bool rootFileExists = !gSystem->AccessPathName(finalRootFilePath.c_str());
        bool validRunsFileExists = !gSystem->AccessPathName(validRunsFilePath.c_str());

        if (!rootFileExists || !validRunsFileExists) {
            // At least one file does not exist
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

                // Output the size of the list of run numbers
                std::cout << "Combination: " << combinationName << ", Number of valid runs: " << validRuns.size() << std::endl;
            } else {
                std::cerr << "Failed to open valid runs file: " << validRunsFilePath << std::endl;
                allOutputFilesExist = false;
                break;
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
            combinedRootFiles.push_back(fileName);
        }
    }
    gSystem->FreeDirectory(dirp);

    // Now plot the combined histograms
    PlotCombinedHistograms(outputDirectory, combinedRootFiles, combinationToValidRuns);
    
    ProcessAllIsolationEnergies(outputDirectory, combinedRootFiles, TriggerCombinationNames::triggerCombinationNameMap);
}
