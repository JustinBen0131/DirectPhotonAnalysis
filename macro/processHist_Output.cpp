#include <iostream>
#include <vector>
#include <string>
#include <TSystem.h>
#include <dirent.h>  // For directory reading
#include <sys/stat.h> // For checking if path is a directory

// Function to merge ROOT files for a given run number using hadd
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
    std::string haddCommand = "hadd -f " + outputFileName;
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
}
