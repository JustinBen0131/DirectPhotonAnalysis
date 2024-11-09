#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <fstream>
#include <algorithm>
#include <sys/stat.h>

// ROOT headers
#include <TSystem.h>

// Function to load run numbers from a text file into an unordered_set for fast lookups
std::unordered_set<std::string> loadRunNumbersFromFile(const std::string& filename) {
    std::unordered_set<std::string> validRunNumbers;
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        validRunNumbers.insert(line);  // Add each line (run number) to the set
    }
    file.close();
    return validRunNumbers;
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

// Main function to merge per-run ROOT files into final output file
void mergeRunFilesToFinalOutput() {
    std::string outputDir = "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/output/";
    std::string inputFile = "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/runListTriggerLUTv1.txt";  // Path to the input text file
    std::unordered_set<std::string> validRunNumbers = loadRunNumbersFromFile(inputFile);
    std::vector<std::string> runNumbers;

    // Add all valid run numbers to runNumbers vector if per-run ROOT file exists
    for (const auto& runNumber : validRunNumbers) {
        std::string runOutputFile = outputDir + runNumber + "_HistOutput.root";
        struct stat buffer;
        if (stat(runOutputFile.c_str(), &buffer) == 0) {
            runNumbers.push_back(runNumber);
        } else {
            std::cout << "Per-run ROOT file not found for run " << runNumber << ". Skipping." << std::endl;
        }
    }

    if (runNumbers.empty()) {
        std::cerr << "No per-run ROOT files found in output directory: " << outputDir << std::endl;
        return;
    }

    std::sort(runNumbers.begin(), runNumbers.end());

    mergeAllRuns(runNumbers, outputDir);
}

