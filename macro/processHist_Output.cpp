#include <iostream>
#include <vector>
#include <string>
#include <TSystem.h>
#include <dirent.h>
#include <sys/stat.h>
#include <algorithm>

// Define ANSI color codes
#define ANSI_RESET  "\x1b[0m"
#define ANSI_RED    "\x1b[31m"
#define ANSI_GREEN  "\x1b[32m"
#define ANSI_YELLOW "\x1b[33m"
#define ANSI_CYAN   "\x1b[36m"


void get_scaledowns(int runnumber, int scaledowns[]) {
    // Connect to the PostgreSQL server
    TSQLServer *db = TSQLServer::Connect("pgsql://sphnxdaqdbreplica:5432/daq", "phnxro", "");
    
    // Check if the connection was successful
    if (db) {
        printf("Server info: %s\n", db->ServerInfo()); // Print server information
    } else {
        printf("Failed to connect to the database\n"); // Print error message if connection fails
        return; // Exit the function if the connection fails
    }

    // Variables to store query results
    TSQLRow *row;
    TSQLResult *res;
    char sql[1000]; // Buffer to hold SQL query strings

    // Loop over each of the 64 scaledown values
    for (int is = 0; is < 64; is++) {
        // Format the SQL query to retrieve the scaledown value for the given run number and bit position
        sprintf(sql, "SELECT scaledown%02d FROM gl1_scaledown WHERE runnumber = %d;", is, runnumber);
        printf("Executing query: %s\n", sql); // Print the SQL query for debugging

        // Execute the query
        res = db->Query(sql);

        // Check the number of rows and fields returned by the query
        int nrows = res->GetRowCount();
        int nfields = res->GetFieldCount();

        // Iterate over the rows of the result set
        for (int i = 0; i < nrows; i++) {
            row = res->Next(); // Get the next row
            // Iterate over the fields in the current row
            for (int j = 0; j < nfields; j++) {
                scaledowns[is] = stoi(row->GetField(j)); // Convert the field value to an integer and store it in the scaledowns array
            }
            delete row; // Clean up the row object
        }
        delete res; // Clean up the result set object
    }

    delete db; // Close the database connection
}



void scale_histogram(TH1* hist, int scaledown, int triggerIndex) {
    if (!hist) {
        std::cerr << ANSI_RED << "Error: Histogram is null for trigger index " << triggerIndex << ANSI_RESET << std::endl;
        return;
    }

    // Apply the scaling if the scaledown factor is valid
    if (scaledown > 0) {
        hist->Scale(scaledown + 1.0);  // Scale by (scaledown + 1)
        std::cout << ANSI_YELLOW << "Scaled histogram \"" << hist->GetName()
                  << "\" for trigger index " << triggerIndex
                  << " by factor " << scaledown + 1
                  << " (scale-down factor: " << scaledown << ")" << ANSI_RESET << std::endl;
    } else if (scaledown == 0) {
        std::cout << ANSI_CYAN << "No scaling applied to histogram \""
                  << hist->GetName() << "\" for trigger index "
                  << triggerIndex << " (scale-down factor is 0)" << ANSI_RESET << std::endl;
    } else {
        std::cout << ANSI_RED << "Skipping histogram \"" << hist->GetName()
                  << "\" for trigger index " << triggerIndex
                  << " because the trigger was off (scale-down factor = -1)" << ANSI_RESET << std::endl;
    }
}

void processQAHistograms(const std::string& outputFileName, int runNumber) {
    // Retrieve scale-down factors for the current run
    int scaledowns[64] = {0};
    get_scaledowns(runNumber, scaledowns);
    
    // Open the merged ROOT file
    TFile* file = TFile::Open(outputFileName.c_str(), "UPDATE");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open merged ROOT file: " << outputFileName << std::endl;
        return;
    }
    
    // Navigate to the QA directory
    TDirectory* qaDir = (TDirectory*)file->Get("QA");
    if (!qaDir) {
        std::cerr << "Error: QA directory not found in the ROOT file: " << outputFileName << std::endl;
        file->Close();
        return;
    }

    for (int triggerIndex = 0; triggerIndex < 64; ++triggerIndex) {
        std::cout << "Processing QA Trigger Directory: Trigger" << triggerIndex << std::endl;
        std::string triggerDirName = "Trigger" + std::to_string(triggerIndex);
        TDirectory* triggerDir = (TDirectory*)qaDir->Get(triggerDirName.c_str());
        if (!triggerDir) {
            std::cout << "Trigger directory not found: " << triggerDirName << std::endl;
            continue;
        }
        
        triggerDir->cd();
        
        std::vector<std::string> histNames = {
            "h8by8TowerEnergySum_" + std::to_string(triggerIndex),
            "h_jet_energy_" + std::to_string(triggerIndex),
            "h_hcal_energy_" + std::to_string(triggerIndex),
            "h_jet_emcal_energy_" + std::to_string(triggerIndex),
            "h_jet_hcalin_energy_" + std::to_string(triggerIndex),
            "h_jet_hcalout_energy_" + std::to_string(triggerIndex),
            "hCluster_maxECore_" + std::to_string(triggerIndex)
        };

        // Scale each histogram in the list
        for (const auto& histName : histNames) {
            TH1* hist = (TH1*)triggerDir->Get(histName.c_str());
            if (hist) {
                scale_histogram(hist, scaledowns[triggerIndex], triggerIndex);
                hist->Write("", TObject::kOverwrite);
            } else {
                std::cout << "Histogram not found: " << histName << std::endl;
            }
        }
        qaDir->cd();
    }
    
    // PhotonAnalysis Directory
    TDirectory* photonDir = (TDirectory*)file->Get("PhotonAnalysis");
    if (photonDir) {
        for (int triggerIndex = 0; triggerIndex < 64; ++triggerIndex) {
            std::cout << "Processing PhotonAnalysis Trigger Directory: Trigger" << triggerIndex << std::endl;
            std::string triggerDirName = "Trigger" + std::to_string(triggerIndex);
            TDirectory* triggerDir = (TDirectory*)photonDir->Get(triggerDirName.c_str());
            if (!triggerDir) {
                std::cout << "Trigger directory not found: " << triggerDirName << std::endl;
                continue;
            }
            triggerDir->cd();

            // List of histogram prefixes to scale
            std::vector<std::string> histogramPrefixes = {
                "isolatedPhotonCount_E",
                "allPhotonCount_E",
                "ptPhoton_E"
            };

            // Get the list of keys and make a copy of the histogram names
            TList* keyList = triggerDir->GetListOfKeys();
            int totalHistograms = keyList->GetEntries();
            std::cout << "Total histograms in Trigger" << triggerIndex << ": " << totalHistograms << std::endl;

            std::vector<std::string> histNames;
            TIter next(keyList);
            TKey* key;
            while ((key = (TKey*)next())) {
                histNames.push_back(key->GetName());
            }

            // Counters for statistics
            int histogramsProcessed = 0;
            int histogramsSkippedZeroEntries = 0;
            int histogramsSkippedTriggerOff = 0;
            int histogramsSkippedNotMatching = 0;
            int histogramsSkippedNotTH1 = 0;

            for (const auto& histName : histNames) {
                TObject* obj = triggerDir->Get(histName.c_str());
                if (!obj) continue;
                std::string className = obj->ClassName();

                // Check if the object is a TH1 histogram
                if (className.find("TH1") != std::string::npos) {
                    // Check if histogram name starts with any of the prefixes and ends with '_' + triggerIndex
                    bool matchesPattern = false;
                    for (const auto& prefix : histogramPrefixes) {
                        std::string suffix = "_" + std::to_string(triggerIndex);
                        if (histName.find(prefix) == 0 && histName.size() >= suffix.size() &&
                            histName.compare(histName.size() - suffix.size(), suffix.size(), suffix) == 0) {
                            matchesPattern = true;
                            break; // Found a matching prefix
                        }
                    }
                    if (matchesPattern) {
                        // Histogram matches the pattern
                        TH1* hist = (TH1*)obj;
                        if (hist) {
                            if (hist->GetEntries() == 0) {
                                histogramsSkippedZeroEntries++;
                                continue;
                            }
                            // Check the scaledown factor
                            if (scaledowns[triggerIndex] == -1) {
                                histogramsSkippedTriggerOff++;
                                continue;
                            }
                            scale_histogram(hist, scaledowns[triggerIndex], triggerIndex);
                            hist->Write("", TObject::kOverwrite);
                            histogramsProcessed++;
                        }
                    } else {
                        histogramsSkippedNotMatching++;
                    }
                } else {
                    histogramsSkippedNotTH1++;
                }

                // Optional: Print progress every 100 histograms
                if (histogramsProcessed % 100 == 0 && histogramsProcessed > 0) {
                    std::cout << "Processed " << histogramsProcessed << " histograms in Trigger" << triggerIndex << std::endl;
                }
            }

            // Report statistics for this trigger
            std::cout << "Trigger " << triggerIndex << " processing summary:" << std::endl;
            std::cout << "  Total histograms: " << totalHistograms << std::endl;
            std::cout << "  Histograms processed (scaled): " << histogramsProcessed << std::endl;
            std::cout << "  Histograms skipped (zero entries): " << histogramsSkippedZeroEntries << std::endl;
            std::cout << "  Histograms skipped (trigger off): " << histogramsSkippedTriggerOff << std::endl;
            std::cout << "  Histograms skipped (not matching pattern): " << histogramsSkippedNotMatching << std::endl;
            std::cout << "  Histograms skipped (not TH1): " << histogramsSkippedNotTH1 << std::endl;

            photonDir->cd();
        }
    }
    file->Close();
    std::cout << "Successfully processed and scaled histograms for run " << runNumber << std::endl;
}

void mergeRunFiles(const std::string& runNumber, const std::string& baseDir, const std::string& outputDir) {
    std::string runPath = baseDir + runNumber + "/";
    std::string outputFileName = outputDir + runNumber + "_HistOutput.root";

    std::cout << "Processing run number: " << runNumber << std::endl;

    // Open the directory containing the ROOT files
    DIR* dir = opendir(runPath.c_str());
    if (!dir) {
        std::cerr << "Error: Cannot open directory " << runPath << std::endl;
        return;
    }

    std::vector<std::string> rootFiles;
    struct dirent* entry;
    while ((entry = readdir(dir))) {
        std::string fileName = entry->d_name;
        if (fileName.find(".root") != std::string::npos) {
            std::string filePath = runPath + fileName;
            rootFiles.push_back(filePath);
        }
    }
    closedir(dir);

    if (rootFiles.empty()) {
        std::cerr << "No ROOT files found in: " << runPath << std::endl;
        return;
    }

    // Construct the hadd command to merge all segment files into one output file
    std::string haddCommand = "hadd -f -T " + outputFileName;
    for (const auto& file : rootFiles) {
        haddCommand += " " + file;
    }

    // Execute the hadd command to merge the files
    std::cout << "Merging files for run: " << runNumber << std::endl;
    int haddResult = gSystem->Exec(haddCommand.c_str());
    if (haddResult != 0) {
        std::cerr << "Error: hadd failed with exit code " << haddResult << " for run " << runNumber << std::endl;
    } else {
        std::cout << "Successfully merged files into " << outputFileName << std::endl;

        // After merging, process and scale histograms in the QA directory
        processQAHistograms(outputFileName, std::stoi(runNumber));
    }
}

void mergeAllRuns(std::vector<std::string>& runNumbers, const std::string& outputDir) {
    if (runNumbers.empty()) {
        std::cerr << "Error: No run numbers provided to merge." << std::endl;
        return;
    }

    // Ensure run numbers are sorted lexicographically
    std::sort(runNumbers.begin(), runNumbers.end());

    std::string runNumberLowest = runNumbers.front();
    std::string runNumberHighest = runNumbers.back();

    std::string finalOutputFile = outputDir + "Final_Merged_Hists_runnumber" + runNumberLowest + "_runnumber" + runNumberHighest + ".root";

    std::string haddCommand = "hadd -f -T " + finalOutputFile;
    for (const auto& runNumber : runNumbers) {
        std::string runOutputFile = outputDir + runNumber + "_HistOutput.root";
        haddCommand += " " + runOutputFile;
    }

    std::cout << "Merging run number ROOT files into final output file: " << finalOutputFile << std::endl;
    int haddResult = gSystem->Exec(haddCommand.c_str());
    if (haddResult != 0) {
        std::cerr << "Error: hadd failed with exit code " << haddResult << " during final merging." << std::endl;
    } else {
        std::cout << "Successfully merged run number files into " << finalOutputFile << std::endl;
    }
}

// Main function to process all runs by merging their segment files
void processHist_Output() {
    std::string baseDir = "/sphenix/tg/tg01/bulk/jbennett/DirectPhotons/output/";
    std::string outputDir = "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/output/";
    std::vector<std::string> runNumbers;

    DIR* dir = opendir(baseDir.c_str());
    if (!dir) {
        std::cerr << "Error: Cannot open base directory " << baseDir << std::endl;
        return;
    }

    struct dirent* entry;
    struct stat s;
    while ((entry = readdir(dir))) {
        std::string folderName = entry->d_name;
        std::string folderPath = baseDir + folderName;

        if (stat(folderPath.c_str(), &s) == 0 && S_ISDIR(s.st_mode) && folderName != "." && folderName != "..") {
            runNumbers.push_back(folderName);
        }
    }
    closedir(dir);

    if (runNumbers.empty()) {
        std::cerr << "No run directories found in base directory: " << baseDir << std::endl;
        return;
    }

    // Process each run number found
    for (const auto& runNumber : runNumbers) {
        mergeRunFiles(runNumber, baseDir, outputDir);
    }
    mergeAllRuns(runNumbers, outputDir);
}
