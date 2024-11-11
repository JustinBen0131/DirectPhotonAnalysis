#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TKey.h>
#include <TColor.h>
#include <sys/stat.h>
#include <TRandom3.h>
#include <TLatex.h>
#include <regex>
#include <TASImage.h>
#include <TImage.h>

bool enableFits = false; // Set to true if you want to enable the fits

// Photon triggers of interest
std::vector<int> photonTriggers = {25, 26, 27};

// Numerator triggers for h8by8TowerEnergySum turn-on curves
std::vector<int> h8by8NumeratorTriggers = photonTriggers;

// Denominator trigger
int denominatorTrigger = 10;

std::map<int, std::map<int, bool>> triggerStatusMap;
// Mapping trigger indices to names based on the provided trigger list
std::map<int, std::string> triggerNameMap = {
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

// Helper function to get the trigger name based on index
std::string getTriggerName(int triggerIndex) {
    if (triggerNameMap.find(triggerIndex) != triggerNameMap.end()) {
        return triggerNameMap[triggerIndex];
    }
    return "Unknown Trigger";
}

std::map<int, int> triggerColorMap = {
    {10, kRed},
    {25, kBlue},
    {26, kGreen + 2},
    {27, kMagenta},
    {29, kBlue},
    {30, kGreen + 2},
    {31, kMagenta}
};



// Helper function to normalize trigger names for comparison
std::string normalizeTriggerName(const std::string& triggerName) {
    std::string normalized = triggerName;

    // Replace special characters and standardize spaces/underscores
    std::replace(normalized.begin(), normalized.end(), '_', ' ');
    normalized.erase(std::remove(normalized.begin(), normalized.end(), '>'), normalized.end());
    normalized.erase(std::remove(normalized.begin(), normalized.end(), '='), normalized.end());

    // Convert to lowercase for case-insensitive matching
    std::transform(normalized.begin(), normalized.end(), normalized.begin(), ::tolower);
    return normalized;
}

// Updated function to get the trigger index from its normalized name
int getTriggerIndexFromName(const std::string& triggerName) {
    std::string normalizedHeader = normalizeTriggerName(triggerName);

    for (const auto& pair : triggerNameMap) {
        std::string normalizedMapName = normalizeTriggerName(pair.second);
        if (normalizedMapName == normalizedHeader) {
            return pair.first;
        }
    }
    return -1; // Return -1 if not found
}

std::string trim(const std::string& str) {
    size_t start = str.find_first_not_of(" \t\n\r\f\v");
    size_t end = str.find_last_not_of(" \t\n\r\f\v");
    return (start == std::string::npos) ? "" : str.substr(start, end - start + 1);
}

std::map<int, std::map<int, bool>> readTriggerStatusFromCSV(const std::string& csvFilePath) {
    std::map<int, std::map<int, bool>> triggerStatusMap;
    std::ifstream file(csvFilePath);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open CSV file at path: " << csvFilePath << std::endl;
        return triggerStatusMap;
    }

    std::string line;

    // Read the header line to get trigger columns
    if (!std::getline(file, line)) {
        std::cerr << "Error: Failed to read the header line from CSV." << std::endl;
        return triggerStatusMap;
    }

    std::vector<std::string> headers;
    std::istringstream headerStream(line);
    std::string header;

    // Extract headers
    std::cout << "Parsing headers from CSV..." << std::endl;
    while (std::getline(headerStream, header, ',')) {
        headers.push_back(trim(header));
        std::cout << "Header found: " << headers.back() << std::endl;
    }

    // Map the header names to trigger indices
    std::vector<int> triggerIndices(headers.size(), -1);
    for (size_t i = 1; i < headers.size(); ++i) { // Skip the first column (runNumber)
        triggerIndices[i] = getTriggerIndexFromName(headers[i]);
        if (triggerIndices[i] == -1) {
            std::cerr << "Warning: Unrecognized trigger name in header: " << headers[i] << std::endl;
        } else {
            std::cout << "Mapped header '" << headers[i] << "' to trigger index " << triggerIndices[i] << std::endl;
        }
    }

    // Read the data lines
    int lineNumber = 1;
    while (std::getline(file, line)) {
        lineNumber++;
        std::istringstream lineStream(line);
        std::string cell;

        // Read and trim the run number
        std::getline(lineStream, cell, ',');
        int runNumber = -1;
        try {
            runNumber = std::stoi(trim(cell));
            std::cout << "Processing run number: " << runNumber << std::endl;
        } catch (const std::invalid_argument&) {
            std::cerr << "Error: Invalid run number '" << cell << "' on line " << lineNumber << std::endl;
            continue;
        }

        // Read the trigger status for each trigger bit
        std::map<int, bool> triggerMap;
        for (size_t i = 1; i < headers.size(); ++i) {
            if (!std::getline(lineStream, cell, ',')) {
                std::cerr << "Warning: Missing cell data on line " << lineNumber << ", column " << i << std::endl;
                break;
            }
            cell = trim(cell); // Trim whitespace

            // Only process if we found a valid trigger index
            if (triggerIndices[i] != -1) {
                triggerMap[triggerIndices[i]] = (cell == "ON");
                std::cout << "Trigger status for '" << headers[i] << "' (index " << triggerIndices[i] << "): " << cell << std::endl;
            }
        }
        triggerStatusMap[runNumber] = triggerMap;
    }

    std::cout << "Finished reading CSV. Total runs processed: " << triggerStatusMap.size() << std::endl;

    // Debug: Print out all run numbers in the map
    for (const auto& entry : triggerStatusMap) {
        std::cout << "Run number in map: " << entry.first << std::endl;
    }

    return triggerStatusMap;
}

// List of run numbers to ignore
std::vector<int> ignoredRunNumbers = {44616, 45034, 45035, 45048, 45061, 45052, 45100, 45103, 45105, 45106, 45153, 45155, 45178, 45181, 45196, 45199, 45203, 45246, 45248, 45288, 45290, 45292, 45315, 45318, 45443, 45485, 45489, 45645, 45856};


std::map<std::string, std::vector<int>> groupRunsByTriggerStatus(const std::map<int, std::map<int, bool>>& triggerStatusMap) {
    std::map<std::string, std::vector<int>> groupedRuns;

    // Define a set for efficient checking of ignored run numbers
    std::set<int> ignoredRunSet(ignoredRunNumbers.begin(), ignoredRunNumbers.end());

    for (const auto& entry : triggerStatusMap) {
        int runNumber = entry.first;
        const auto& triggerMap = entry.second;

        // Check if the run is in the ignored list
        if (ignoredRunSet.count(runNumber) > 0) {
            std::cout << "Skipping run: " << runNumber << " (ignored)" << std::endl;
            continue; // Skip this run
        }

        // Ensure the denominator trigger is ON for all runs
        if (triggerMap.count(denominatorTrigger) == 0 || !triggerMap.at(denominatorTrigger)) {
            continue;
        }

        // Generate a human-readable key for grouping based on photon triggers
        std::string groupKey = "photonCurves_";
        for (size_t i = 0; i < photonTriggers.size(); ++i) {
            int triggerBit = photonTriggers[i];
            bool isTriggerOn = triggerMap.count(triggerBit) ? triggerMap.at(triggerBit) : false;
            groupKey += (isTriggerOn ? std::to_string(triggerBit) + "on_" : std::to_string(triggerBit) + "off_");
        }
        
        // Remove trailing underscore
        if (!groupKey.empty() && groupKey.back() == '_') {
            groupKey.pop_back();
        }

        // Add the run number to the corresponding group
        groupedRuns[groupKey].push_back(runNumber);

        // Additional grouping for each photon trigger individually
        for (const auto& triggerBit : photonTriggers) {
            if (triggerMap.count(triggerBit) && triggerMap.at(triggerBit)) {
                std::string individualGroupKey = "Trigger" + std::to_string(denominatorTrigger) + "_" + std::to_string(triggerBit) + "on";
                groupedRuns[individualGroupKey].push_back(runNumber);
            }
        }
    }

    return groupedRuns;
}


// Function to parse the groupKey and get trigger statuses
std::map<int, bool> parseGroupKey(const std::string& groupKey) {
    std::map<int, bool> triggerStatus;
    // The groupKey can be in formats like "photonCurves_25on_26on_27on" or "Trigger10_25on"

    std::vector<std::string> tokens;
    std::stringstream ss(groupKey);
    std::string token;
    while (std::getline(ss, token, '_')) {
        tokens.push_back(token);
    }

    for (const auto& trigToken : tokens) {
        // The token is like '25on', '25off', 'Trigger10', or '27on'
        std::regex re("(Trigger)?(\\d+)(on|off)?");
        std::smatch match;
        if (std::regex_match(trigToken, match, re)) {
            int triggerIndex = std::stoi(match[2]);
            bool isOn = false;
            if (match[3].matched) {
                isOn = (match[3] == "on");
            } else if (match[1].matched) {
                // If 'Trigger' is present without 'on' or 'off', assume it's 'on'
                isOn = true;
            }
            triggerStatus[triggerIndex] = isOn;
        }
    }
    return triggerStatus;
}



void drawRunNumbersOnCanvas(const std::vector<int>& runNumbers, const std::string& groupKey) {
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
    if (runNumbers.size() == 1) {
        numColumns = 1;
        textSize = 0.04;
        xStart = 0.82; xEnd = 0.92;
        yStart = 0.93; yEnd = 0.83;
        xSpacingFactor = 1.0;
        ySpacingFactor = 1.0;
    } else if (runNumbers.size() == 14) {
        numColumns = 2;
        textSize = 0.028;
        xStart = 0.7; xEnd = 0.9;
        yStart = 0.92; yEnd = 0.52;
        xSpacingFactor = 0.9;
        ySpacingFactor = 1.2;
    } else if (runNumbers.size() == 14) {
        numColumns = 2;
        textSize = 0.028;
        xStart = 0.7; xEnd = 0.9;
        yStart = 0.92; yEnd = 0.52;
        xSpacingFactor = 0.9;
        ySpacingFactor = 1.2;
    } else if (runNumbers.size() == 2) {
        numColumns = 2;
        textSize = 0.035;
        xStart = 0.7; xEnd = 0.9;
        yStart = 0.92; yEnd = 0.82;
        xSpacingFactor = 0.9;
        ySpacingFactor = 1.2;
    } else if (runNumbers.size() == 9) {
        numColumns = 3;
        textSize = 0.038;
        xStart = 0.55; xEnd = 0.9;
        yStart = 0.92; yEnd = 0.62;
        xSpacingFactor = 1.1;
        ySpacingFactor = 1.1;
    } else if (runNumbers.size() == 103) {
        numColumns = 10;
        textSize = 0.018;
        xStart = 0.47; xEnd = 0.93;
        yStart = 0.93; yEnd = 0.42;
        xSpacingFactor = 1.0;
        ySpacingFactor = 1.0;
    } else if (runNumbers.size() == 21) {
        numColumns = 3;
        textSize = 0.025;
        xStart = 0.52; xEnd = 0.92;
        yStart = 0.93; yEnd = 0.42;
        xSpacingFactor = 1.0;
        ySpacingFactor = 1.0;
    } else if (runNumbers.size() == 76) {
        numColumns = 8;
        textSize = 0.02;
        xStart = 0.48; xEnd = 0.91;
        yStart = 0.93; yEnd = 0.45;
        xSpacingFactor = 1.0;
        ySpacingFactor = 1.0;
    } else if (runNumbers.size() == 104) {
        numColumns = 8;
        textSize = 0.02;
        xStart = 0.48; xEnd = 0.91;
        yStart = 0.93; yEnd = 0.45;
        xSpacingFactor = 1.0;
        ySpacingFactor = 1.0;
    } else if (runNumbers.size() == 112) {
        numColumns = 8;
        textSize = 0.02;
        xStart = 0.48; xEnd = 0.91;
        yStart = 0.93; yEnd = 0.45;
        xSpacingFactor = 1.0;
        ySpacingFactor = 1.0;
    } else if (runNumbers.size() == 75) {
        numColumns = 5;
        textSize = 0.02;
        xStart = 0.48; xEnd = 0.91;
        yStart = 0.93; yEnd = 0.45;
        xSpacingFactor = 1.0;
        ySpacingFactor = 1.0;
    } else if (runNumbers.size() == 84) {
        numColumns = 5;
        textSize = 0.02;
        xStart = 0.48; xEnd = 0.91;
        yStart = 0.93; yEnd = 0.45;
        xSpacingFactor = 1.0;
        ySpacingFactor = 1.0;
    } else if (runNumbers.size() == 87) {
        numColumns = 6;
        textSize = 0.02;
        xStart = 0.48; xEnd = 0.91;
        yStart = 0.93; yEnd = 0.45;
        xSpacingFactor = 1.0;
        ySpacingFactor = 1.0;
    } else if (runNumbers.size() == 115) {
        numColumns = 8;
        textSize = 0.02;
        xStart = 0.48; xEnd = 0.91;
        yStart = 0.93; yEnd = 0.45;
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

    // Set the text size
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
    double ySpacing = ySpacingFactor * (yStart - yEnd) / (numRows + 1); // +1 for the "Runs:" label


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
}

TF1* sigmoidFit(const std::string& name, double xmin, double xmax,
                double amplitude, double slope, double xOffset,
                double amplitudeMin, double amplitudeMax,
                double slopeMin, double slopeMax,
                double xOffsetMin, double xOffsetMax) {
    // Define a sigmoid function for fitting
    TF1* fitFunc = new TF1(name.c_str(), "[0]/(1+exp(-[1]*(x-[2])))", xmin, xmax);
    fitFunc->SetParNames("Amplitude", "Slope", "XOffset");

    // Set the manual parameters
    fitFunc->SetParameter(0, amplitude);  // Amplitude
    fitFunc->SetParameter(1, slope);      // Slope
    fitFunc->SetParameter(2, xOffset);    // XOffset

    // Set the manual limits for the parameters
    fitFunc->SetParLimits(0, amplitudeMin, amplitudeMax);  // Amplitude limits
    fitFunc->SetParLimits(1, slopeMin, slopeMax);          // Slope limits
    fitFunc->SetParLimits(2, xOffsetMin, xOffsetMax);      // XOffset limits

    return fitFunc;
}

void plotGroupedTriggerTurnOnCurves(const std::string& inputDir, const std::string& outputDir, const std::string& groupKey, const std::vector<int>& runNumbers, const std::map<int, bool>& groupTriggerStatus) {
    // Define the histogram base name
    std::string histNameBase = "h8by8TowerEnergySum_";

    std::vector<int> numeratorTriggers;

    // Determine numerator triggers based on the groupKey or groupTriggerStatus dynamically
    for (int trig : photonTriggers) {
        // Check if the groupKey matches any specific trigger being ON
        std::string keyCheck = "Trigger10_" + std::to_string(trig) + "on";
        if (groupKey == keyCheck) {
            numeratorTriggers.push_back(trig);
        }
    }

    // If the groupKey does not match a specific trigger, use groupTriggerStatus to check which triggers are ON
    if (numeratorTriggers.empty()) {
        for (int trig : photonTriggers) {
            if (groupTriggerStatus.count(trig) && groupTriggerStatus.at(trig)) {
                numeratorTriggers.push_back(trig);
            }
        }
    }

    // If no numerator triggers are 'ON', skip plotting
    if (numeratorTriggers.empty()) {
        std::cout << "No numerator triggers are ON in group " << groupKey << ", skipping turn-on curve plotting." << std::endl;
        return;
    }
    
    // Sum the denominator histogram over the runs
    TH1* denominatorHist = nullptr;

    for (const int runNumber : runNumbers) {
        // Construct the input file path
        std::string inputFilePath = inputDir + "/" + std::to_string(runNumber) + "_HistOutput.root";

        TFile* inputFile = TFile::Open(inputFilePath.c_str(), "READ");

        if (!inputFile || inputFile->IsZombie()) {
            std::cerr << "Error: Could not open file " << inputFilePath << std::endl;
            continue;
        }

        std::string dirName = "QA/Trigger" + std::to_string(denominatorTrigger);
        TDirectory* dir = (TDirectory*)inputFile->Get(dirName.c_str());

        if (!dir) {
            std::cerr << "Warning: Directory " << dirName << " not found in file " << inputFilePath << std::endl;
            inputFile->Close();
            continue;
        }

        std::string histName = histNameBase + std::to_string(denominatorTrigger);
        TH1* hist = (TH1*)dir->Get(histName.c_str());

        if (hist) {
            if (!denominatorHist) {
                denominatorHist = (TH1*)hist->Clone("combined_denominator");
                denominatorHist->SetDirectory(0); // Detach from file
            } else {
                denominatorHist->Add(hist);
            }
        }

        inputFile->Close();
    }

    if (!denominatorHist || denominatorHist->GetEntries() == 0) {
        std::cerr << "Denominator histogram is empty for group " << groupKey << ", cannot plot turn-on curves." << std::endl;
        return;
    }

    // Now sum the numerator histograms
    std::map<int, TH1*> numeratorHists;

    for (int numeratorTrigger : numeratorTriggers) {
        TH1* combinedNumeratorHist = nullptr;

        for (const int runNumber : runNumbers) {
            // Construct the input file path
            std::string inputFilePath = inputDir + "/" + std::to_string(runNumber) + "_HistOutput.root";

            TFile* inputFile = TFile::Open(inputFilePath.c_str(), "READ");

            if (!inputFile || inputFile->IsZombie()) {
                std::cerr << "Error: Could not open file " << inputFilePath << std::endl;
                continue;
            }

            std::string dirName = "QA/Trigger" + std::to_string(numeratorTrigger);
            TDirectory* dir = (TDirectory*)inputFile->Get(dirName.c_str());

            if (!dir) {
                std::cerr << "Warning: Directory " << dirName << " not found in file " << inputFilePath << std::endl;
                inputFile->Close();
                continue;
            }

            std::string histName = histNameBase + std::to_string(numeratorTrigger);
            TH1* hist = (TH1*)dir->Get(histName.c_str());

            if (hist) {
                if (!combinedNumeratorHist) {
                    combinedNumeratorHist = (TH1*)hist->Clone(("combined_numerator_" + std::to_string(numeratorTrigger)).c_str());
                    combinedNumeratorHist->SetDirectory(0); // Detach from file
                } else {
                    combinedNumeratorHist->Add(hist);
                }
            }

            inputFile->Close();
        }

        if (combinedNumeratorHist && combinedNumeratorHist->GetEntries() > 0) {
            numeratorHists[numeratorTrigger] = combinedNumeratorHist;
        } else {
            std::cerr << "Combined numerator histogram is empty for trigger " << numeratorTrigger << " in group " << groupKey << std::endl;
        }
    }

    // Now we can compute the ratios and plot
    TCanvas* canvas = new TCanvas("canvas_turnon", "Trigger Turn-On", 800, 600);
    TLegend* legend = new TLegend(0.2, 0.7, 0.4, 0.9);
    legend->SetTextSize(0.03);

    std::vector<int> colors = {kBlue, kGreen + 2, kMagenta, kRed};

    int colorIndex = 0;
    bool firstDraw = true;

    for (const auto& pair : numeratorHists) {
        int numeratorTrigger = pair.first;
        TH1* numeratorHist = pair.second;

        // Compute the ratio
        TH1* ratioHist = (TH1*)numeratorHist->Clone(("ratio_" + std::to_string(numeratorTrigger)).c_str());
        ratioHist->Divide(denominatorHist);

        // Set styles
        int color = kBlack; // Default color
        if (triggerColorMap.find(numeratorTrigger) != triggerColorMap.end()) {
            color = triggerColorMap[numeratorTrigger];
        }
        ratioHist->SetLineColor(color);
        ratioHist->SetMarkerColor(color);
        ratioHist->SetLineWidth(2);
        ratioHist->SetTitle((histNameBase + " Turn-On Curve").c_str());
        ratioHist->GetYaxis()->SetTitle("TriggerX/Trigger10 (MBD N&S>=1)");
        ratioHist->GetXaxis()->SetTitle("Energy [GeV]");
        ratioHist->GetYaxis()->SetRangeUser(0, 1.5);  // Adjust y-axis range for better view

        // Draw the ratio histogram
        if (firstDraw) {
            ratioHist->Draw("E");
            firstDraw = false;
        } else {
            ratioHist->Draw("E SAME");
        }

        // Add entry to legend
        std::string triggerName = getTriggerName(numeratorTrigger);
        legend->AddEntry(ratioHist, triggerName.c_str(), "l");

        colorIndex++;
    }

    // Draw the dashed line at y = 1
    TLine* line = new TLine(denominatorHist->GetXaxis()->GetXmin(), 1, denominatorHist->GetXaxis()->GetXmax(), 1);
    line->SetLineStyle(2);  // Dashed line
    line->SetLineColor(kBlack);
    line->Draw("SAME");

    // Draw the legend and save the canvas
    legend->Draw();

    // Save the canvas with a filename indicating the group
    std::string outputFileName = outputDir + "h8by8TowerEnergySum_TurnOn.png";
    canvas->SaveAs(outputFileName.c_str());

    delete canvas;
}

// Function to sort run numbers numerically
void sortRunNumbers(std::vector<std::string>& runNumbers) {
    // Convert run numbers to integers
    std::vector<int> runNumbersInt;
    for (const auto& runStr : runNumbers) {
        try {
            runNumbersInt.push_back(std::stoi(runStr));
        } catch (const std::exception& e) {
            std::cerr << "Error converting run number " << runStr << " to integer: " << e.what() << std::endl;
        }
    }

    // Sort the integers
    std::sort(runNumbersInt.begin(), runNumbersInt.end());

    // Convert back to strings
    runNumbers.clear();
    for (const auto& runInt : runNumbersInt) {
        runNumbers.push_back(std::to_string(runInt));
    }
}

void plotGroupTabulatedHistograms(const std::string& inputDir, const std::string& outputBaseDir, const std::vector<int>& runNumbers, const std::string& groupKey) {
    // First, sort the runNumbers vector to have run numbers in sequential order
    std::vector<int> sortedRunNumbers = runNumbers;
    std::sort(sortedRunNumbers.begin(), sortedRunNumbers.end());

    // Convert run numbers to strings
    std::vector<std::string> runNumberStrings;
    for (const auto& runNumber : sortedRunNumbers) {
        runNumberStrings.push_back(std::to_string(runNumber));
    }

    const int runsPerImage = 32; // 32 runs per image
    const int columns = 8;
    const int rows = 4;
    const int totalImages = (runNumberStrings.size() + runsPerImage - 1) / runsPerImage;

    for (int imageIndex = 0; imageIndex < totalImages; ++imageIndex) {
        // Create a canvas with adjusted dimensions
        TCanvas* canvas = new TCanvas(("canvasTabulated_" + groupKey + "_" + std::to_string(imageIndex)).c_str(), "Tabulated Plot", 1600, 800);
        canvas->Divide(columns, rows);

        std::vector<TImage*> images; // Store images to keep them in memory

        for (int i = 0; i < runsPerImage; ++i) {
            int runIndex = imageIndex * runsPerImage + i;
            if (runIndex >= runNumberStrings.size()) break;

            std::string runNumber = runNumberStrings[runIndex];

            // Construct the path to the image file
            std::string imagePath = outputBaseDir + runNumber + "/h8by8TowerEnergySum.png";

            // Check if the image file exists
            if (gSystem->AccessPathName(imagePath.c_str())) {
                std::cerr << "Warning: Image file " << imagePath << " does not exist." << std::endl;
                continue;
            }

            canvas->cd(i + 1);
            if (!gPad) {
                std::cerr << "Error: Could not access pad " << i + 1 << " on canvas" << std::endl;
                continue;
            }
            gPad->SetFillColor(0);
            gPad->Clear();

            // Load the image
            TImage* img = TImage::Open(imagePath.c_str());
            if (!img) {
                std::cerr << "Error: Could not open image " << imagePath << std::endl;
                continue;
            }

            // Adjust pad margins
            gPad->SetLeftMargin(0.0);
            gPad->SetRightMargin(0.0);
            gPad->SetTopMargin(0.0);
            gPad->SetBottomMargin(0.0);

            // Draw the image
            img->Draw("X");

            // Add run number label
            TLatex latex;
            latex.SetTextSize(0.04);
            latex.SetNDC();
            latex.DrawLatex(0.05, 0.9, ("Run: " + runNumber).c_str());

            gPad->Modified();
            gPad->Update();

            // Store the image to keep it in memory
            images.push_back(img);
        }

        // Save the tabulated plot as a PNG file
        std::string outputFileName = outputBaseDir + groupKey + "/h8by8TowerEnergySum_Tabulated_" + std::to_string(imageIndex + 1) + ".png";
        canvas->SaveAs(outputFileName.c_str());
        delete canvas;

        // Delete the images after saving the canvas
        for (auto img : images) {
            delete img;
        }
        images.clear();
    }
}

void plotGroupedHistograms(const std::string& inputDir, const std::string& baseOutputDir, const std::map<std::string, std::vector<int>>& groupedRuns) {
    // Define the histogram to be used
    std::string histNameBase = "h8by8TowerEnergySum_";
    
    std::vector<int> defaultTriggerIndices = photonTriggers;
    defaultTriggerIndices.insert(defaultTriggerIndices.begin(), denominatorTrigger); // Include denominator trigger


    for (const auto& group : groupedRuns) {
        const std::string& groupKey = group.first;
        const std::vector<int>& runNumbers = group.second;

        std::cout << "Processing group: " << groupKey << " with " << runNumbers.size() << " runs." << std::endl;

        // Create a unique folder for the current group
        std::string outputDir = baseOutputDir + groupKey + "/";
        gSystem->mkdir(outputDir.c_str(), true);

        std::vector<int> triggerIndices;
         // Dynamically determine triggers to include based on the groupKey
         if (groupKey.find("Trigger10_") != std::string::npos) {
             for (int trig : photonTriggers) {
                 std::string keyCheck = "Trigger10_" + std::to_string(trig) + "on";
                 if (groupKey == keyCheck) {
                     triggerIndices.push_back(denominatorTrigger);
                     triggerIndices.push_back(trig);
                 }
             }
         } else {
             // Use default triggers including the denominator trigger
             triggerIndices = photonTriggers;
             triggerIndices.insert(triggerIndices.begin(), denominatorTrigger);
         }
        

        TCanvas* canvas = new TCanvas("canvas", "Grouped Overlay Plot", 800, 600);
        TLegend* legend = new TLegend(0.2, 0.2, 0.4, 0.35);
        legend->SetTextSize(0.035);
        canvas->SetLogy();

        std::vector<int> colors = {kRed, kBlue, kGreen + 2, kMagenta};
        bool firstDraw = true;

        // Loop over the trigger indices and sum histograms
        for (size_t j = 0; j < triggerIndices.size(); ++j) {
            int triggerIndex = triggerIndices[j];

            TH1* combinedHist = nullptr;

            // Loop through all runs in the current group
            for (const int runNumber : runNumbers) {
                // Construct the input file path
                std::string inputFilePath = inputDir + "/" + std::to_string(runNumber) + "_HistOutput.root";

                TFile* inputFile = TFile::Open(inputFilePath.c_str(), "READ");

                if (!inputFile || inputFile->IsZombie()) {
                    std::cerr << "Error: Could not open file " << inputFilePath << std::endl;
                    continue;
                }

                std::string dirName = "QA/Trigger" + std::to_string(triggerIndex);
                TDirectory* dir = (TDirectory*)inputFile->Get(dirName.c_str());

                if (!dir) {
                    std::cerr << "Warning: Directory " << dirName << " not found in file " << inputFilePath << std::endl;
                    inputFile->Close();
                    continue;
                }

                std::string histName = histNameBase + std::to_string(triggerIndex);
                TH1* hist = (TH1*)dir->Get(histName.c_str());

                if (hist) {
                    if (!combinedHist) {
                        combinedHist = (TH1*)hist->Clone(("combined_" + std::to_string(triggerIndex)).c_str());
                        combinedHist->SetDirectory(0); // Detach from file
                    } else {
                        combinedHist->Add(hist);
                    }
                }

                inputFile->Close();
            }

            // Draw the combined histogram if it has content
            if (combinedHist && combinedHist->GetEntries() > 0) {
                combinedHist->SetLineWidth(2);
                int color = kBlack; // Default color
                if (triggerColorMap.find(triggerIndex) != triggerColorMap.end()) {
                    color = triggerColorMap[triggerIndex];
                }
                combinedHist->SetLineColor(color);
                
                
                combinedHist->SetTitle(("Grouped Overlay for " + histNameBase).c_str());
                combinedHist->GetXaxis()->SetTitle("Energy [GeV]");
                combinedHist->GetYaxis()->SetTitle("Events");

                if (firstDraw) {
                    combinedHist->Draw("HIST");
                    firstDraw = false;
                } else {
                    combinedHist->Draw("HIST SAME");
                }

                // Only add triggers that are actually 'ON' to the legend
                legend->AddEntry(combinedHist, getTriggerName(triggerIndex).c_str(), "l");
            }
        }

        legend->Draw();

        // After drawing histograms and legend, draw run numbers
        std::cout << "Drawing run numbers on canvas for group: " << groupKey << std::endl;
        drawRunNumbersOnCanvas(runNumbers, groupKey);

        // Update the canvas to ensure everything is drawn
        canvas->Modified();
        canvas->Update();

        // Save the canvas with a filename indicating the trigger combination
        std::string outputFileName = outputDir + "h8by8TowerEnergySum_Group.png";
        canvas->SaveAs(outputFileName.c_str());

        delete canvas;

        // Parse the groupKey to get the trigger statuses
        std::map<int, bool> groupTriggerStatus = parseGroupKey(groupKey);

        // Call the function to plot the turn-on curves
        plotGroupedTriggerTurnOnCurves(inputDir, outputDir, groupKey, runNumbers, groupTriggerStatus);
        
        plotGroupTabulatedHistograms(inputDir, baseOutputDir, runNumbers, groupKey);
    }
}




void plotTriggerTurnOnCurves(TFile* inputFile, const std::string& histNameBase, const std::vector<int>& numeratorTriggers, int denominatorTrigger, const std::string& outputFilePath, const std::string& runNumber) {
    TCanvas* canvas = new TCanvas("canvas", "Trigger Turn-On", 800, 600);
    TLegend* legend = new TLegend(0.6, 0.7, 0.8, 0.9);
    legend->SetHeader("Turn-On Curves", "C");
    legend->SetTextSize(0.03);

    std::string denominatorHistName = histNameBase + std::to_string(denominatorTrigger);
    std::string denominatorDirName = "QA/Trigger" + std::to_string(denominatorTrigger);
    TDirectory* denominatorDir = (TDirectory*)inputFile->Get(denominatorDirName.c_str());
    TH1* denominatorHist = (TH1*)denominatorDir->Get(denominatorHistName.c_str());

    if (!denominatorHist) {
        std::cerr << "Error: Denominator histogram not found for trigger " << denominatorTrigger << std::endl;
        return;
    }

    for (size_t i = 0; i < numeratorTriggers.size(); ++i) {
        int numeratorTrigger = numeratorTriggers[i];
        std::string numeratorHistName = histNameBase + std::to_string(numeratorTrigger);
        std::string numeratorDirName = "QA/Trigger" + std::to_string(numeratorTrigger);
        TDirectory* numeratorDir = (TDirectory*)inputFile->Get(numeratorDirName.c_str());
        TH1* numeratorHist = (TH1*)numeratorDir->Get(numeratorHistName.c_str());

        if (!numeratorHist) {
            std::cerr << "Error: Numerator histogram not found for trigger " << numeratorTrigger << std::endl;
            continue;
        }

        // Clone the numerator and compute the ratio
        TH1* ratioHist = (TH1*)numeratorHist->Clone(("ratio_" + std::to_string(numeratorTrigger)).c_str());
        ratioHist->Divide(denominatorHist);

        int color = kBlack; // Default color
        if (triggerColorMap.find(numeratorTrigger) != triggerColorMap.end()) {
            color = triggerColorMap[numeratorTrigger];
        }
        ratioHist->SetLineColor(color);
        ratioHist->SetMarkerColor(color);
        ratioHist->SetLineWidth(2);
        ratioHist->SetTitle((histNameBase + " Turn-On Curve").c_str());
        ratioHist->GetYaxis()->SetTitle("TriggerX/Trigger10 (MBD N&S>=1)");
        ratioHist->GetXaxis()->SetTitle("Energy [GeV]");
        ratioHist->GetYaxis()->SetRangeUser(0, 2);  // Adjust y-axis range for better view

        // Draw the ratio histogram
        if (i == 0) {
            ratioHist->Draw("E");
        } else {
            ratioHist->Draw("E SAME");
        }

        // Add entry to legend
        std::string triggerName = getTriggerName(numeratorTrigger);
        legend->AddEntry(ratioHist, triggerName.c_str(), "l");

        // Manually set sigmoid parameters and limits for each trigger
        double amplitude, slope, xOffset;
        double amplitudeMin, amplitudeMax, slopeMin, slopeMax, xOffsetMin, xOffsetMax;
        
        if (numeratorTrigger == 25) {
            amplitude = 1.1; slope = 0.58; xOffset = 6;
            amplitudeMin = 1.0; amplitudeMax = 1.3;
            slopeMin = 0.55; slopeMax = 0.6;
            xOffsetMin = 5.5; xOffsetMax = 6.5;
            
        } else if (numeratorTrigger == 26) {
            amplitude = 1.2; slope = 0.5; xOffset = 6.5;
            amplitudeMin = 1.0; amplitudeMax = 1.2;
            slopeMin = 0.4; slopeMax = 0.6;
            xOffsetMin = 6; xOffsetMax = 7;
            
        } else if (numeratorTrigger == 27) {
            amplitude = 1.2; slope = 0.52; xOffset = 7.2;
            amplitudeMin = 1.0; amplitudeMax = 1.3;
            slopeMin = 0.4; slopeMax = 0.6;
            xOffsetMin = 6.5; xOffsetMax = 7.5;
        }
        
        if (numeratorTrigger == 21) {
            amplitude = 1.0; slope = 0.85; xOffset = 12.5;
            amplitudeMin = 0.95; amplitudeMax = 1.05;
            slopeMin = 0.8; slopeMax = 0.9;
            xOffsetMin = 12; xOffsetMax = 13;
        } else if (numeratorTrigger == 22) {
            amplitude = 1.2; slope = 0.82; xOffset = 12.2;
            amplitudeMin = 0.95; amplitudeMax = 1.05;
            slopeMin = 0.8; slopeMax = 0.85;
            xOffsetMin = 12; xOffsetMax = 12.5;
        } else if (numeratorTrigger == 23) {
            amplitude = 1.0; slope = 0.78; xOffset = 12.7;
            amplitudeMin = 0.95; amplitudeMax = 1.05;
            slopeMin = 0.75; slopeMax = 0.85;
            xOffsetMin = 12.2; xOffsetMax = 13.2;
        }
        if (enableFits) {
            // Fit the ratio with a sigmoid function using the manual parameters and limits
            double xmin = ratioHist->GetXaxis()->GetXmin();
            double xmax = ratioHist->GetXaxis()->GetXmax();
            TF1* fitFunc = sigmoidFit("sigmoidFit_" + std::to_string(numeratorTrigger), xmin, xmax, amplitude, slope, xOffset,
                                      amplitudeMin, amplitudeMax, slopeMin, slopeMax, xOffsetMin, xOffsetMax);
            int color = kBlack; // Default color
            if (triggerColorMap.find(numeratorTrigger) != triggerColorMap.end()) {
                color = triggerColorMap[numeratorTrigger];
            }
            fitFunc->SetLineColor(color);  // Match the color of the fit to the histogram
            ratioHist->Fit(fitFunc, "R");

            // Draw the run number on the top-right corner of the canvas
            TLatex latex;
            latex.SetTextSize(0.04);
            latex.SetNDC();
            latex.SetTextAlign(33); // Align top-right
            latex.DrawLatex(0.9, 0.95, ("Run: " + runNumber).c_str());

            
            // Draw the fit
            fitFunc->Draw("SAME");
        }
    }

    // Draw the dashed line at y = 1
    TLine* line = new TLine(denominatorHist->GetXaxis()->GetXmin(), 1, denominatorHist->GetXaxis()->GetXmax(), 1);
    line->SetLineStyle(2);  // Dashed line
    line->SetLineColor(kBlack);
    line->Draw("SAME");

    TLatex latex;
    latex.SetTextSize(0.045);
    latex.SetNDC(); // Normalized Device Coordinates
    latex.SetTextAlign(33); // Align top-right
    latex.DrawLatex(0.9, 0.25, ("Run: " + runNumber).c_str());
    
    // Draw the legend and save the canvas
    legend->Draw();
    canvas->SaveAs(outputFilePath.c_str());
    delete canvas;
}


bool isNumeric(const std::string& str) {
    return !str.empty() && std::all_of(str.begin(), str.end(), ::isdigit);
}


void plotOverlayHistograms(const std::string& inputFilePath, const std::string& outputDir, const std::string& runNumber, const std::map<int, std::map<int, bool>>& triggerStatusMap) {
    
    // Check if the runNumber is numeric before attempting conversion
    int runNum = -1;
    if (isNumeric(runNumber)) {
        try {
            runNum = std::stoi(runNumber);
        } catch (const std::invalid_argument&) {
            std::cerr << "Error: Invalid run number '" << runNumber << "' encountered during stoi conversion." << std::endl;
            return;
        }
    } else {
        std::cout << "Debug: Non-numeric run number '" << runNumber << "', skipping integer conversion." << std::endl;
    }

    
    // Set the input file path
    TFile* inputFile = TFile::Open(inputFilePath.c_str(), "READ");

    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open the file " << inputFilePath << std::endl;
        return;
    }
    
    // Only attempt to use the trigger status map if runNum is valid
    std::map<int, bool> triggerMap;
    if (runNum != -1) {
        if (triggerStatusMap.find(runNum) != triggerStatusMap.end()) {
            triggerMap = triggerStatusMap.at(runNum);
        } else {
            std::cerr << "Run number " << runNum << " not found in the trigger status map." << std::endl;
        }
    }
    // Set up histogram names and titles
    std::vector<std::string> histNames = {
        "h8by8TowerEnergySum_",
        "h_hcal_energy_",
        "h_jet_energy_",
        "hCluster_maxECore_"
    };

    std::vector<std::string> titles = {
        "Max 8x8 Tower Energy Sum; Max 8x8 Tower Energy Sum [GeV]; Events",
        "Max HCal Tower Energy Sums; Energy [GeV]; Events",
        "Maximum 0.8x0.8 Energy Sum (EMCAL + HCAL); Maximum 0.8x0.8 Energy Sum (EMCAL + HCAL) [GeV]; Events",
        "Max ECore Cluster; Max ECore [GeV]; Events"
    };

    gSystem->mkdir(outputDir.c_str(), true);

    std::vector<std::string> outputFiles = {
        outputDir + "h8by8TowerEnergySum.png",
        outputDir + "h_hcal_energy.png",
        outputDir + "h_jet_energy.png",
        outputDir + "hCluster_maxECore.png"
    };

    // Set the style for the plots
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);

    // Loop over the histogram types
    for (size_t i = 0; i < histNames.size(); ++i) {
        TCanvas* canvas = new TCanvas(("canvas" + std::to_string(i)).c_str(), "Overlay Plot", 800, 600);
        canvas->SetLogy();
        TLegend* legend = new TLegend(0.5, 0.65, 0.7, 0.85);
        legend->SetHeader("Triggers", "C");
        legend->SetTextSize(0.04);

        bool firstDraw = true;

        // Determine the correct set of trigger indices based on the histogram name
        std::vector<int> triggerIndices;
        if (histNames[i] == "h8by8TowerEnergySum_") {
            // Use the photonTriggers vector
            triggerIndices = photonTriggers;
            triggerIndices.insert(triggerIndices.begin(), denominatorTrigger); // Include denominator trigger
        } else if (histNames[i] == "hCluster_maxECore_") {
            triggerIndices = {10, 29, 30, 31};
        } else if (histNames[i] == "h_hcal_energy_") {
            triggerIndices = {10, 21, 22, 23};
        } else if (histNames[i] == "h_jet_energy_") {
            triggerIndices = {10, 21, 22, 23};
        }
        // Loop over the trigger indices and plot the histograms
        for (size_t j = 0; j < triggerIndices.size(); ++j) {
            int triggerIndex = triggerIndices[j];

            // Skip if the trigger is 'OFF'
            if (triggerMap.find(triggerIndex) != triggerMap.end() && !triggerMap[triggerIndex]) {
                continue;
            }
            
            // Get the trigger name based on the index
            std::string triggerName = getTriggerName(triggerIndex);

            // Navigate to the correct directory in the ROOT file
            std::string dirName = "QA/Trigger" + std::to_string(triggerIndex);
            TDirectory* dir = (TDirectory*)inputFile->Get(dirName.c_str());
            if (!dir) {
                std::cerr << "Warning: Directory " << dirName << " not found in the file." << std::endl;
                continue;
            }

            std::string histName = histNames[i] + std::to_string(triggerIndex);
            TH1* hist = (TH1*)dir->Get(histName.c_str());

            if (!hist) {
                std::cerr << "Warning: Histogram " << histName << " not found in the directory " << dirName << "." << std::endl;
                continue;
            }

            hist->SetLineWidth(2);
            // Use the triggerColorMap to set the color
            int color = kBlack; // Default color
            if (triggerColorMap.find(triggerIndex) != triggerColorMap.end()) {
                color = triggerColorMap[triggerIndex];
            }
            hist->SetLineColor(color);

            hist->SetTitle(titles[i].c_str());

            // Draw histograms with different styles for hCluster_maxECore_
            if (histNames[i] == "hCluster_maxECore_") {
                hist->SetMarkerStyle(20); // Use points instead of lines
                int color = kBlack; // Default color
                if (triggerColorMap.find(triggerIndex) != triggerColorMap.end()) {
                    color = triggerColorMap[triggerIndex];
                }
                hist->SetMarkerColor(color);
                if (firstDraw) {
                    hist->Draw("HIST"); // Draw points
                    firstDraw = false;
                } else {
                    hist->Draw("HIST SAME");
                }
            } else {
                if (firstDraw) {
                    hist->Draw("HIST");
                    firstDraw = false;
                } else {
                    hist->Draw("HIST SAME");
                }
            }

            // Add the actual trigger name to the legend
            legend->AddEntry(hist, triggerName.c_str(), "l");
        }

        legend->Draw();
        TLatex latex;
        latex.SetTextSize(0.065);
        latex.SetNDC();
        latex.SetTextAlign(33); // Align top-right
        latex.DrawLatex(0.8, 0.9, ("Run: " + runNumber).c_str());
        
        canvas->SaveAs(outputFiles[i].c_str());
        std::cout << "Saved plot: " << outputFiles[i] << std::endl;
        delete canvas;
    }

    inputFile->Close();
    delete inputFile;
    
    // After closing inputFile, call the new function for the specified histograms
     TFile* inputFileNew = TFile::Open(inputFilePath.c_str(), "READ");
     if (!inputFileNew || inputFileNew->IsZombie()) {
         std::cerr << "Error: Could not reopen the file " << inputFilePath << std::endl;
         return;
     }

    plotTriggerTurnOnCurves(inputFileNew, "h8by8TowerEnergySum_", h8by8NumeratorTriggers, denominatorTrigger, outputDir + "h8by8TowerEnergySum_TurnOn.png", runNumber);

     // Plot ratios for h_jet_energy (triggers 17/10, 18/10, 19/10)
     std::vector<int> hJetNumeratorTriggers = {21, 22, 23};
     plotTriggerTurnOnCurves(inputFileNew, "h_jet_energy_", hJetNumeratorTriggers, 10, outputDir + "h_jet_energy_TurnOn.png", runNumber);

     inputFileNew->Close();
     delete inputFileNew;
}

// Function to process a single run file
void processRunFile(const std::string& inputFilePath, const std::string& outputDir, const std::string& runNumber, const std::map<int, std::map<int, bool>>& triggerStatusMap) {
    gSystem->mkdir(outputDir.c_str(), true); // Create output directory if it doesn't exist
    plotOverlayHistograms(inputFilePath, outputDir, runNumber, triggerStatusMap);
}


void plotTabulatedHistograms(const std::string& inputDir, const std::string& outputDir, std::vector<std::string>& runNumbers) {
    // First, sort the runNumbers vector to have run numbers in sequential order
    sortRunNumbers(runNumbers);

    const int runsPerImage = 34; // 34 runs per image
    const int columns = 8;
    const int rows = 5;
    const int totalPads = columns * rows; // 40 pads
    const int totalImages = (runNumbers.size() + runsPerImage - 1) / runsPerImage;

    std::cout << "Debug: Number of run numbers to process: " << runNumbers.size() << std::endl;
    std::cout << "Debug: Total images to generate: " << totalImages << std::endl;

    for (int imageIndex = 0; imageIndex < totalImages; ++imageIndex) {
        // Adjust canvas dimensions to be wider than tall
        TCanvas* canvas = new TCanvas(("canvasTabulated" + std::to_string(imageIndex)).c_str(), "Tabulated Plot", 1600, 1000);
        canvas->Divide(columns, rows);

        std::cout << "Debug: Creating canvas for image " << imageIndex + 1 << "/" << totalImages << std::endl;

        std::vector<TImage*> images; // Store images to keep them in memory

        for (int i = 0; i < runsPerImage; ++i) {
            int runIndex = imageIndex * runsPerImage + i;
            if (runIndex >= runNumbers.size()) break;

            std::string runNumber = runNumbers[runIndex];

            // Construct the path to the image file
            std::string imagePath = outputDir + runNumber + "/h8by8TowerEnergySum.png";
            std::cout << "Debug: Attempting to load image: " << imagePath << std::endl;

            // Check if the image file exists
            if (gSystem->AccessPathName(imagePath.c_str())) {
                std::cerr << "Warning: Image file " << imagePath << " does not exist." << std::endl;
                continue;
            }

            canvas->cd(i + 1);
            if (!gPad) {
                std::cerr << "Error: Could not access pad " << i + 1 << " on canvas" << std::endl;
                continue;
            }
            gPad->SetFillColor(0);
            gPad->Clear();

            // Load the image
            TImage* img = TImage::Open(imagePath.c_str());
            if (!img) {
                std::cerr << "Error: Could not open image " << imagePath << std::endl;
                continue;
            }

            // Adjust pad margins
            gPad->SetLeftMargin(0.0);
            gPad->SetRightMargin(0.0);
            gPad->SetTopMargin(0.0);
            gPad->SetBottomMargin(0.0);

            // Draw the image
            img->Draw("X");

            gPad->Modified();
            gPad->Update();

            // Store the image to keep it in memory
            images.push_back(img);
        }

        // Save the tabulated plot as a PNG file
        std::string outputFileName = outputDir + "h8by8TowerEnergySum_Tabulated_" + std::to_string(imageIndex + 1) + ".png";
        std::cout << "Debug: Saving tabulated plot as " << outputFileName << std::endl;
        canvas->SaveAs(outputFileName.c_str());
        delete canvas;

        // Delete the images after saving the canvas
        for (auto img : images) {
            delete img;
        }
        images.clear();
    }
}

// Main function to run both processes
void AnalyzeTriggerEfficencies() {
    gROOT->LoadMacro("sPhenixStyle.C");
    SetsPhenixStyle();

    std::string baseDir = "/Users/patsfan753/Desktop/DirectPhotonAna/";
    std::string combinedInputFilePath = baseDir + "runListTriggerLutv1_combinedHistOutput.root";
    std::string combinedOutputDir = baseDir + "Plots/";
    std::string csvFilePath = baseDir + "runListTriggerLutv1/triggerAnalysisLUTv1.csv";
    
    auto triggerStatusMap = readTriggerStatusFromCSV(csvFilePath);
    

    // Ensure the combined input file exists
    TFile* combinedFile = TFile::Open(combinedInputFilePath.c_str(), "READ");
    if (!combinedFile || combinedFile->IsZombie()) {
        std::cerr << "Error: Could not open the combined input file: " << combinedInputFilePath << std::endl;
        return;
    }
    combinedFile->Close();
    delete combinedFile;

    // Create the output directory if it doesn't exist
    gSystem->mkdir(combinedOutputDir.c_str(), true);

    // Process the combined file for plotting
    std::cout << "Processing combined input file: " << combinedInputFilePath << std::endl;
    std::cout << "Saving plots to directory: " << combinedOutputDir << std::endl;
    plotOverlayHistograms(combinedInputFilePath, combinedOutputDir, "Combined", triggerStatusMap);
    

    // Now process individual run files
    std::string inputDir = baseDir + "runListTriggerLutv1/output/";
    std::string outputBaseDir = baseDir + "runListTriggerLutv1/Plots/";

    // Use TSystemDirectory to list files
    TSystemDirectory dir("outputDir", inputDir.c_str());
    TList *files = dir.GetListOfFiles();
    std::vector<std::string> runNumbers;
    
    if (files) {
        TSystemFile *file;
        TString fname;
        TIter next(files);
        while ((file=(TSystemFile*)next())) {
            fname = file->GetName();
            if (!file->IsDirectory()) {
                // Check if fname matches pattern *_HistOutput.root
                if (fname.EndsWith("_HistOutput.root")) {
                    // Extract the run number from the filename
                    std::string filename = fname.Data();
                    size_t pos = filename.find("_HistOutput.root");
                    std::string runNumber = filename.substr(0, pos);
                    
                    runNumbers.push_back(runNumber);
                    // Create output directory
                    std::string outputDir = outputBaseDir + runNumber + "/";
                    // Full input file path
                    std::string inputFilePath = inputDir + filename;

                    processRunFile(inputFilePath, outputDir, runNumber, triggerStatusMap);

                }
            }
        }
        plotTabulatedHistograms(inputDir, outputBaseDir, runNumbers);
    }
    auto groupedRuns = groupRunsByTriggerStatus(triggerStatusMap);

    // Plot histograms for grouped runs
    plotGroupedHistograms(inputDir, outputBaseDir, groupedRuns);
}
