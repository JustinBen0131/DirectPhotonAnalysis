#include <iostream>
#include <vector>
#include <string>
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <dirent.h>  // For directory reading
#include <sys/stat.h> // For checking if path is a directory

#define ANSI_BOLD "\033[1m"
#define ANSI_RED "\033[1;31m"
#define ANSI_GREEN "\033[1;32m"
#define ANSI_CYAN "\033[1;36m"
#define ANSI_YELLOW "\033[1;33m"
#define ANSI_RESET "\033[0m"

// Struct to hold parsed cut values
struct CutValues {
    float clusECore = -1.0;
    float chi = -1.0;
    float asymmetry = -1.0;
    float pT_min = -1.0;
    float pT_max = -1.0;
    int triggerIndex = -1;
};

// Function to parse histogram names and extract cut values and pT range
CutValues parseHistName(const std::string& histName) {
    CutValues cuts;

    // Regex to extract the relevant values: E, Chi, Asym, pT range, and trigger index
    std::regex re("E([0-9]+(?:point[0-9]*)?)_Chi([0-9]+(?:point[0-9]*)?)_Asym([0-9]+(?:point[0-9]*)?)_pT_([0-9]+(?:point[0-9]*)?)to([0-9]+(?:point[0-9]*)?)_(\\d+)");
    std::smatch match;

    // Lambda to convert strings with 'point' to float values
    auto convert = [](const std::string& input) -> float {
        std::string temp = input;
        size_t pointPos = temp.find("point");
        if (pointPos != std::string::npos) {
            temp.replace(pointPos, 5, ".");
        }
        try {
            return std::stof(temp);
        } catch (const std::exception&) {
            return -1.0f;
        }
    };

    // If the regex matches the histogram name
    if (std::regex_search(histName, match, re) && match.size() == 7) {
        cuts.clusECore = convert(match[1].str());
        cuts.chi = convert(match[2].str());
        cuts.asymmetry = convert(match[3].str());
        cuts.pT_min = convert(match[4].str());
        cuts.pT_max = convert(match[5].str());
        cuts.triggerIndex = std::stoi(match[6].str());

        std::cout << "Parsed histogram: " << histName << std::endl;
        std::cout << "  ECore: " << cuts.clusECore << ", Chi: " << cuts.chi
                  << ", Asymmetry: " << cuts.asymmetry << ", pT_min: " << cuts.pT_min
                  << ", pT_max: " << cuts.pT_max << ", Trigger Index: " << cuts.triggerIndex << std::endl;
    } else {
        std::cerr << "Error: Failed to parse histogram name: " << histName << std::endl;
    }

    return cuts;
}

// Function to reorganize histograms in the final output file
void reorganizeHistograms(TFile* outputFile) {
    // Go to PhotonAnalysis directory
    TDirectory* photonDir = outputFile->GetDirectory("PhotonAnalysis");
    if (!photonDir) {
        std::cerr << "PhotonAnalysis directory not found!" << std::endl;
        return;
    }

    // Iterate over Trigger directories in PhotonAnalysis
    TIter next_trigger(photonDir->GetListOfKeys());
    TKey* triggerKey;
    while ((triggerKey = (TKey*)next_trigger())) {
        TDirectory* triggerDir = dynamic_cast<TDirectory*>(triggerKey->ReadObj());
        if (!triggerDir || std::string(triggerDir->GetName()).find("Trigger") != 0) continue;

        std::string triggerName = triggerDir->GetName();
        std::cout << "Processing " << triggerName << std::endl;

        // Iterate over histograms in the Trigger directory
        TIter next_hist(triggerDir->GetListOfKeys());
        TKey* histKey;
        while ((histKey = (TKey*)next_hist())) {
            TH1* hist = dynamic_cast<TH1*>(histKey->ReadObj());
            if (!hist) continue;

            std::string histName = hist->GetName();

            // Parse the histogram name to extract the cut values and pT range
            CutValues cuts = parseHistName(histName);

            if (cuts.clusECore < 0 || cuts.chi < 0 || cuts.asymmetry < 0 || cuts.pT_min < 0 || cuts.pT_max < 0 || cuts.triggerIndex < 0) {
                std::cerr << "Skipping histogram due to parsing error: " << histName << std::endl;
                continue;
            }

            // Create a directory structure for the cuts
            std::stringstream cutDirName;
            cutDirName << "PhotonAnalysis/" << triggerName << "/E" << cuts.clusECore
                       << "_Chi" << cuts.chi << "_Asym" << cuts.asymmetry;
            outputFile->cd();
            if (!gDirectory->GetDirectory(cutDirName.str().c_str())) {
                outputFile->mkdir(cutDirName.str().c_str());
            }

            // Create a directory structure for the pT range within the cuts
            std::stringstream pTDirName;
            pTDirName << cutDirName.str() << "/pT_" << cuts.pT_min << "to" << cuts.pT_max;
            outputFile->cd();
            if (!gDirectory->GetDirectory(pTDirName.str().c_str())) {
                outputFile->mkdir(pTDirName.str().c_str());
            }

            // Move the histogram into the correct directory
            outputFile->cd(pTDirName.str().c_str());
            hist->Write(histName.c_str(), TObject::kOverwrite);

            std::cout << "Moved histogram: " << histName << " to " << pTDirName.str() << std::endl;

            delete hist;  // Clean up memory
        }
    }

    outputFile->Write();  // Save changes to the file
    std::cout << "Reorganization complete." << std::endl;
}

// Function to get scale-down factors from the database
void get_scaledowns(int runnumber, int scaledowns[]) {
    std::cout << ANSI_CYAN << "Connecting to the database to retrieve scale-down factors for run number: " << runnumber << ANSI_RESET << std::endl;
    TSQLServer *db = TSQLServer::Connect("pgsql://sphnxdaqdbreplica:5432/daq", "phnxro", "");
    if (!db) {
        std::cerr << ANSI_RED << "Error: Could not connect to the database" << ANSI_RESET << std::endl;
        return;
    }

    TSQLResult *res = nullptr;
    TSQLRow *row = nullptr;
    char sql[1000];

    for (int is = 0; is < 64; is++) {
        sprintf(sql, "select scaledown%02d from gl1_scaledown where runnumber = %d;", is, runnumber);
        res = db->Query(sql);
        if (!res) {
            std::cerr << ANSI_RED << "Error: Failed to execute query for scale-down factor at index " << is << ": " << sql << ANSI_RESET << std::endl;
            scaledowns[is] = 0; // Default to 0 if query fails
            continue;
        }

        if (res->GetRowCount() == 0) {
            std::cout << ANSI_YELLOW << "No scale-down factor found for trigger bit " << is << ". Defaulting to 0." << ANSI_RESET << std::endl;
            scaledowns[is] = 0;
        } else {
            row = res->Next();
            if (row) {
                try {
                    scaledowns[is] = std::stoi(row->GetField(0));
                    std::cout << ANSI_GREEN << "Successfully retrieved scale-down factor for trigger bit " << is << ": " << scaledowns[is] << ANSI_RESET << std::endl;
                } catch (const std::exception &e) {
                    std::cerr << ANSI_RED << "Error: Failed to convert scale-down factor to integer for trigger bit " << is << ": " << e.what() << ANSI_RESET << std::endl;
                    scaledowns[is] = 0;
                }
                delete row;
            }
        }
        delete res;
    }
    delete db;
    std::cout << ANSI_CYAN << "Completed retrieval of scale-down factors for run number: " << runnumber << ANSI_RESET << std::endl;
}

// Function to scale histograms based on the scale-down factors
void scale_histogram(TH1* hist, int scaledown, int triggerIndex) {
    if (scaledown > 0) {
        hist->Scale(scaledown + 1);
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

// Function to scale histograms inside QA directories
void scale_histograms_in_directory(TDirectory* dir, int scaledowns[]) {
    TIter next_key(dir->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next_key())) {
        std::string dirName = key->GetName();
        TDirectory* subDir = dynamic_cast<TDirectory*>(key->ReadObj());
        if (!subDir) continue;

        // Check for sub-directories named "TriggerX"
        if (dirName.find("Trigger") == 0) {
            std::cout << ANSI_GREEN << "Accessing directory: " << dirName << ANSI_RESET << std::endl;
            int triggerIndex = std::stoi(dirName.substr(7)); // Extract trigger index from "TriggerX"

            TIter next_hist(subDir->GetListOfKeys());
            TKey* histKey;
            while ((histKey = (TKey*)next_hist())) {
                std::string histName = histKey->GetName();

                // Check if the histogram matches the specified patterns
                if (histName.find("h8by8TowerEnergySum_") == 0 ||
                    histName.find("h_hcal_energy_") == 0 ||
                    histName.find("h_jet_energy_") == 0 ||
                    histName.find("hCluster_maxECore_") == 0) {
                    
                    TH1* hist = dynamic_cast<TH1*>(histKey->ReadObj());
                    if (hist && triggerIndex >= 0 && triggerIndex < 64) {
                        std::cout << ANSI_YELLOW << "Found histogram \"" << histName
                                  << "\" in trigger directory: " << dirName << ANSI_RESET << std::endl;
                        scale_histogram(hist, scaledowns[triggerIndex], triggerIndex);
                        hist->Write("", TObject::kOverwrite);
                        delete hist;
                    }
                }
            }
            delete subDir;
        }
    }
}

// Function to process ROOT files in the output directory and merge histograms
void processRunByRun_histOutput() {
    std::string outputDir = "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/output/";
    std::vector<std::string> rootFiles;
    DIR* dir = opendir(outputDir.c_str());
    if (!dir) {
        std::cerr << ANSI_RED << "Error: Cannot open directory " << outputDir << ANSI_RESET << std::endl;
        return;
    }

    struct dirent* entry;
    while ((entry = readdir(dir))) {
        std::string fileName = entry->d_name;
        if (fileName.find("_HistOutput.root") != std::string::npos) {
            rootFiles.push_back(outputDir + fileName);
            std::cout << ANSI_GREEN << "Found ROOT file: " << fileName << ANSI_RESET << std::endl;
        }
    }
    closedir(dir);

    if (rootFiles.empty()) {
        std::cerr << ANSI_RED << "No ROOT files found in: " << outputDir << ANSI_RESET << std::endl;
        return;
    }

    for (const auto& rootFile : rootFiles) {
        std::cout << ANSI_CYAN << "Processing file: " << rootFile << ANSI_RESET << std::endl;

        // Extract run number from the file name
        size_t start = rootFile.find_last_of("/") + 1;
        size_t end = rootFile.find("_HistOutput.root");
        int runNumber = 0;
        try {
            runNumber = std::stoi(rootFile.substr(start, end - start));
            std::cout << ANSI_GREEN << "Extracted run number: " << runNumber << ANSI_RESET << std::endl;
        } catch (const std::exception &e) {
            std::cerr << ANSI_RED << "Error: Failed to extract run number from file name: " << e.what() << ANSI_RESET << std::endl;
            continue;
        }

        int scaledowns[64] = {0};
        get_scaledowns(runNumber, scaledowns);

        // Open the ROOT file in update mode to scale histograms
        TFile* file = TFile::Open(rootFile.c_str(), "UPDATE");
        if (!file || file->IsZombie()) {
            std::cerr << ANSI_RED << "Error: Failed to open file " << rootFile << ANSI_RESET << std::endl;
            delete file;
            continue;
        }

        // Access QA directory
        TDirectory* qaDir = file->GetDirectory("QA");
        if (qaDir) {
            std::cout << ANSI_CYAN << "Accessing QA directory in " << rootFile << ANSI_RESET << std::endl;
            scale_histograms_in_directory(qaDir, scaledowns);
        } else {
            std::cerr << ANSI_RED << "QA directory not found in " << rootFile << ANSI_RESET << std::endl;
        }

        file->Close();
        delete file;
    }

    // Now merge all the processed files
    std::string outputFileName = outputDir + "Final_Merged_HistOutput.root";
    std::string haddCommand = "hadd -f " + outputFileName;
    for (const auto& rootFile : rootFiles) {
        haddCommand += " " + rootFile;
    }
    std::cout << ANSI_GREEN << "Merging files into: " << outputFileName << ANSI_RESET << std::endl;
    int haddResult = gSystem->Exec(haddCommand.c_str());
    if (haddResult != 0) {
        std::cerr << ANSI_RED << "Error: hadd failed for merging ROOT files" << ANSI_RESET << std::endl;
    } else {
        std::cout << ANSI_BOLD << "Successfully merged ROOT files into " << outputFileName << ANSI_RESET << std::endl;

        // Now reorganize the merged ROOT file
        TFile* finalFile = TFile::Open(outputFileName.c_str(), "UPDATE");
        if (finalFile && !finalFile->IsZombie()) {
            reorganizeHistograms(finalFile);
            finalFile->Close();
            delete finalFile;
        } else {
            std::cerr << ANSI_RED << "Error: Could not open the final merged ROOT file for reorganization" << ANSI_RESET << std::endl;
        }
    }
}

