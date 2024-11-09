#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

#include <iostream>
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


struct CutValues {
    float clusECore = 0;
    float asymmetry = 0;
    float chi = 0;
    int triggerIndex = 0;
    float pTMin = -1;  // Default to -1 indicating no pT bin
    float pTMax = -1;  // Default to -1 indicating no pT bin
    float isoMin = -1; // Default to -1 indicating no isoEt range
    float isoMax = -1; // Default to -1 indicating no isoEt range
    bool invalid = false; // Flag to indicate parsing failure
};



// Mapping trigger indices to names based on the provided trigger list
std::map<int, std::string> triggerNameMap = {
    {0, "Clock"},
    {1, "ZDC South"},
    {2, "ZDC North"},
    {3, "ZDC Coincidence"},
    {4, "HCAL Singles"},
    {5, "HCAL Coincidence"},
    {8, "MBD S >= 1"},
    {9, "MBD N >= 1"},
    {10, "MBD N&S >= 1"},
    {11, "MBD N&S >= 2"},
    {12, "MBD N&S >= 1, vtx < 10 cm"},
    {13, "MBD N&S >= 1, vtx < 30 cm"},
    {14, "MBD N&S >= 1, vtx < 60 cm"},
    {15, "HCAL Singles + MBD NS >= 1"},
    {16, "Jet 6 GeV + MBD NS >= 1"},
    {17, "Jet 8 GeV + MBD NS >= 1"},
    {18, "Jet 10 GeV + MBD NS >= 1"},
    {19, "Jet 12 GeV + MBD NS >= 1"},
    {20, "Jet 6 GeV"},
    {21, "Jet 8 GeV"},
    {22, "Jet 10 GeV"},
    {23, "Jet 12 GeV"},
    {24, "Photon 1 GeV + MBD NS >= 1"},
    {25, "Photon 2 GeV + MBD NS >= 1"},
    {26, "Photon 3 GeV + MBD NS >= 1"},
    {27, "Photon 4 GeV + MBD NS >= 1"},
    {28, "Photon 1 GeV"},
    {29, "Photon 2 GeV"},
    {30, "Photon 3 GeV"},
    {31, "Photon 4 GeV"}
};

// Helper function to get the trigger name based on index
std::string getTriggerName(int triggerIndex) {
    if (triggerNameMap.find(triggerIndex) != triggerNameMap.end()) {
        return triggerNameMap[triggerIndex];
    }
    return "Unknown Trigger";
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


void plotTriggerTurnOnCurves(TFile* inputFile, const std::string& histNameBase, const std::vector<int>& numeratorTriggers, int denominatorTrigger, const std::string& outputFilePath, const std::string& runNumber) {
    TCanvas* canvas = new TCanvas("canvas", "Trigger Turn-On", 800, 600);
    TLegend* legend = new TLegend(0.2, 0.7, 0.4, 0.9);
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

    // Loop over the numerator triggers and calculate the ratio
    std::vector<int> colors = {kRed, kBlue, kGreen + 2, kMagenta};
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

        // Set styles
        ratioHist->SetLineColor(colors[i % colors.size()]);
        ratioHist->SetMarkerColor(colors[i % colors.size()]);
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
        
        if (numeratorTrigger == 29) {
            amplitude = 1.1; slope = 0.58; xOffset = 6;
            amplitudeMin = 1.0; amplitudeMax = 1.3;
            slopeMin = 0.55; slopeMax = 0.6;
            xOffsetMin = 5.5; xOffsetMax = 6.5;
            
        } else if (numeratorTrigger == 30) {
            amplitude = 1.2; slope = 0.5; xOffset = 6.5;
            amplitudeMin = 1.0; amplitudeMax = 1.2;
            slopeMin = 0.4; slopeMax = 0.6;
            xOffsetMin = 6; xOffsetMax = 7;
            
        } else if (numeratorTrigger == 31) {
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

        // Fit the ratio with a sigmoid function using the manual parameters and limits
        double xmin = ratioHist->GetXaxis()->GetXmin();
        double xmax = ratioHist->GetXaxis()->GetXmax();
        TF1* fitFunc = sigmoidFit("sigmoidFit_" + std::to_string(numeratorTrigger), xmin, xmax, amplitude, slope, xOffset,
                                  amplitudeMin, amplitudeMax, slopeMin, slopeMax, xOffsetMin, xOffsetMax);
        fitFunc->SetLineColor(colors[i % colors.size()]);  // Match the color of the fit to the histogram
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

    // Draw the dashed line at y = 1
    TLine* line = new TLine(denominatorHist->GetXaxis()->GetXmin(), 1, denominatorHist->GetXaxis()->GetXmax(), 1);
    line->SetLineStyle(2);  // Dashed line
    line->SetLineColor(kBlack);
    line->Draw("SAME");

    // Draw the legend and save the canvas
    legend->Draw();
    canvas->SaveAs(outputFilePath.c_str());
    delete canvas;
}


void plotOverlayHistograms(const std::string& inputFilePath, const std::string& outputDir, const std::string& runNumber) {
    // Set the input file path
    TFile* inputFile = TFile::Open(inputFilePath.c_str(), "READ");

    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open the file " << inputFilePath << std::endl;
        return;
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

    // Define colors for each trigger index (ensure colors are distinct and easily visible)
    std::vector<int> colors = {kRed, kBlue, kGreen + 2, kMagenta, kCyan + 1};

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
            triggerIndices = {10, 29, 30, 31};
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
            hist->SetLineColor(colors[j % colors.size()]); // Set unique color for each trigger
            hist->SetTitle(titles[i].c_str());

            // Draw histograms with different styles for hCluster_maxECore_
            if (histNames[i] == "hCluster_maxECore_") {
                hist->SetMarkerStyle(20); // Use points instead of lines
                hist->SetMarkerColor(colors[j % colors.size()]);
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
        latex.SetTextSize(0.04);
        latex.SetNDC();
        latex.SetTextAlign(33); // Align top-right
        latex.DrawLatex(0.9, 0.95, ("Run: " + runNumber).c_str());
        
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

     // Plot ratios for h8by8TowerEnergySum (triggers 25/10, 26/10, 27/10)
     std::vector<int> h8by8NumeratorTriggers = {29, 30, 31};
    plotTriggerTurnOnCurves(inputFileNew, "h8by8TowerEnergySum_", h8by8NumeratorTriggers, 10, outputDir + "h8by8TowerEnergySum_TurnOn.png", runNumber);

     // Plot ratios for h_jet_energy (triggers 17/10, 18/10, 19/10)
     std::vector<int> hJetNumeratorTriggers = {21, 22, 23};
     plotTriggerTurnOnCurves(inputFileNew, "h_jet_energy_", hJetNumeratorTriggers, 10, outputDir + "h_jet_energy_TurnOn.png", runNumber);

     inputFileNew->Close();
     delete inputFileNew;
}

// Function to process a single run file
void processRunFile(const std::string& inputFilePath, const std::string& outputDir, const std::string& runNumber) {
    gSystem->mkdir(outputDir.c_str(), true); // Create output directory if it doesn't exist
    plotOverlayHistograms(inputFilePath, outputDir, runNumber);
}



// Main function to run both processes
void AnalyzeTriggerEfficencies() {
    gROOT->LoadMacro("sPhenixStyle.C");
    SetsPhenixStyle();

    std::string baseDir = "/Users/patsfan753/Desktop/DirectPhotonAna/";
    std::string combinedInputFilePath = baseDir + "runListTriggerLutv1_combinedHistOutput.root";
    std::string combinedOutputDir = baseDir + "Plots/";

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
    plotOverlayHistograms(combinedInputFilePath, combinedOutputDir, "Combined");
    

    // Now process individual run files
    std::string inputDir = baseDir + "runListTriggerLutv1/output/";
    std::string outputBaseDir = baseDir + "runListTriggerLutv1/Plots/";

    // Use TSystemDirectory to list files
    TSystemDirectory dir("outputDir", inputDir.c_str());
    TList *files = dir.GetListOfFiles();
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
                    // Create output directory
                    std::string outputDir = outputBaseDir + runNumber + "/";
                    // Full input file path
                    std::string inputFilePath = inputDir + filename;
                    // Process the run file
                    processRunFile(inputFilePath, outputDir, runNumber);
                }
            }
        }
    }
}
