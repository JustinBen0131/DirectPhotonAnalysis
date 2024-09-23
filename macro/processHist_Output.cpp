#include <iostream>
#include <vector>
#include <string>
#include <TSystem.h>
#include <TFile.h>
#include <dirent.h>  // For directory reading
#include <sys/stat.h> // For checking if path is a directory

#define ANSI_BOLD "\033[1m"
#define ANSI_RED "\033[1;31m"
#define ANSI_GREEN "\033[1;32m"
#define ANSI_CYAN "\033[1;36m"
#define ANSI_YELLOW "\033[1;33m"
#define ANSI_RESET "\033[0m"

// Function to merge segment ROOT files for a given run number using hadd
void mergeRunFiles(const std::string& runNumber, const std::string& baseDir, const std::string& outputDir) {
    // Construct paths
    std::string runPath = baseDir + runNumber + "/";
    std::string outputFileName = outputDir + runNumber + "_HistOutput.root";

    // Open the directory containing the ROOT files
    DIR* dir = opendir(runPath.c_str());
    if (!dir) {
        std::cerr << ANSI_RED << "Error: Cannot open directory " << runPath << ANSI_RESET << std::endl;
        return;
    }

    std::vector<std::string> rootFiles;
    struct dirent* entry;
    while ((entry = readdir(dir))) {
        std::string fileName = entry->d_name;
        if (fileName.find(".root") != std::string::npos) {
            rootFiles.push_back(runPath + fileName);
        }
    }
    closedir(dir);

    if (rootFiles.empty()) {
        std::cerr << ANSI_RED << "No ROOT files found in: " << runPath << ANSI_RESET << std::endl;
        return;
    }

    // Construct the hadd command to merge all segment files into one output file
    std::string haddCommand = "hadd -f " + outputFileName;
    for (const auto& file : rootFiles) {
        haddCommand += " " + file;
    }

    // Execute the hadd command
    std::cout << ANSI_GREEN << "Merging files for run: " << runNumber << ANSI_RESET << std::endl;
    int haddResult = gSystem->Exec(haddCommand.c_str());
    if (haddResult != 0) {
        std::cerr << ANSI_RED << "Error: hadd failed for run " << runNumber << ANSI_RESET << std::endl;
        return;
    }

    // Verify the merged output to ensure no duplicates
    TFile* outputFile = TFile::Open(outputFileName.c_str(), "READ");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << ANSI_RED << "Error: Failed to open output file " << outputFileName << ANSI_RESET << std::endl;
        return;
    }

    std::cout << ANSI_BOLD << "Successfully merged run number: " << runNumber << " into " << outputFileName << ANSI_RESET << std::endl;
    outputFile->Close();
    delete outputFile;
}

// Main function to process all runs by merging their segment files
void processHist_Output() {
    std::string baseDir = "/sphenix/tg/tg01/bulk/jbennett/DirectPhotons/output/";
    std::string outputDir = "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/output/";
    std::vector<std::string> runNumbers;

    std::cout << ANSI_BOLD << "Scanning base directory: " << baseDir << ANSI_RESET << std::endl;
    
    // Open the base directory to find all run directories
    DIR* dir = opendir(baseDir.c_str());
    if (!dir) {
        std::cerr << ANSI_RED << "Error: Cannot open base directory " << baseDir << ANSI_RESET << std::endl;
        return;
    }

    struct dirent* entry;
    struct stat s;
    while ((entry = readdir(dir))) {
        std::string folderName = entry->d_name;
        std::string folderPath = baseDir + folderName;

        // Check if the entry is a directory
        if (stat(folderPath.c_str(), &s) == 0 && S_ISDIR(s.st_mode) && folderName != "." && folderName != "..") {
            runNumbers.push_back(folderName);
        }
    }
    closedir(dir);

    // Process each run number found
    for (const auto& runNumber : runNumbers) {
        mergeRunFiles(runNumber, baseDir, outputDir);
    }
}
