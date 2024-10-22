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

bool saveHistogramsFlag = false;  // Set to true to save histograms, false to skip

std::string inputDir = "/Users/patsfan753/Desktop/DirectPhotonAna/output/";
std::string outputDir = "/Users/patsfan753/Desktop/DirectPhotonAna/Plots/";
std::string inputFilePath = inputDir + "Final_Merged_Hists_runnumber44686_runnumber44707.root";


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
    {24, "Photon 2 GeV + MBD NS >= 1"},
    {25, "Photon 3 GeV + MBD NS >= 1"},
    {26, "Photon 4 GeV + MBD NS >= 1"},
    {27, "Photon 5 GeV + MBD NS >= 1"},
    {28, "Photon 2 GeV"},
    {29, "Photon 3 GeV"},
    {30, "Photon 4 GeV"},
    {31, "Photon 5 GeV"}
};


// Helper function to get the trigger name based on index
std::string getTriggerName(int triggerIndex) {
    if (triggerNameMap.find(triggerIndex) != triggerNameMap.end()) {
        return triggerNameMap[triggerIndex];
    }
    return "Unknown Trigger";
}


// Recursive function to save histograms in the ROOT file as PNGs in the same directory structure
void saveHistograms(TDirectory* dir, const std::string& currentPath) {
    // Get the list of keys in the current directory
    TIter next(dir->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) {
        TObject* obj = key->ReadObj();
        if (obj->InheritsFrom(TDirectory::Class())) {
            // If the object is a directory, create a corresponding folder and recursively process it
            std::string subDirName = currentPath + "/" + obj->GetName();
            gSystem->mkdir(subDirName.c_str(), true);
            saveHistograms((TDirectory*)obj, subDirName);
        } else if (obj->InheritsFrom(TH1::Class())) {
            // If the object is a histogram, save it as a PNG
            TH1* hist = (TH1*)obj;
            std::string outputPath = currentPath + "/" + hist->GetName() + ".png";

            // Create a canvas for drawing
            TCanvas canvas;

            // Check if the histogram name matches the isoEt patterns to add labels
            std::string histName = hist->GetName();
            if (histName.find("h2_isoEtE_") == 0 || histName.find("h1_isoEt_") == 0) {
                // Extract trigger index and Delta R value from the histogram name
                size_t underscore1 = histName.find('_');
                size_t underscore2 = histName.find('_', underscore1 + 1);
                size_t underscore3 = histName.find('_', underscore2 + 1);
                std::string triggerIndexStr = histName.substr(underscore2 + 1, underscore3 - underscore2 - 1);
                std::cout << "Extracted trigger index string: " << triggerIndexStr << std::endl;
                // Updated code inside saveHistograms() to handle special cases like "TowersInCone"
                if (histName.find("TowersInCone") != std::string::npos) {
                    std::cerr << "Warning: Histogram name contains 'TowersInCone', skipping trigger index extraction." << std::endl;
                    hist->Draw();
                    canvas.SaveAs(outputPath.c_str());
                    std::cout << "Saved: " << outputPath << std::endl;
                    delete obj;
                    continue; // Skip further processing for "TowersInCone" histograms
                }

                int triggerIndex = std::stoi(histName.substr(underscore2 + 1, underscore3 - underscore2 - 1));
                
                std::string dR_cut_str = histName.substr(underscore3 + 1);

                // Get the corresponding trigger name
                std::string triggerName = getTriggerName(triggerIndex);

                // Draw the histogram
                hist->Draw("COLZ");

                // Add Delta R cut and trigger name on the canvas
                TLatex text;
                text.SetNDC();
                text.SetTextSize(0.03);
                text.DrawLatex(0.15, 0.85, ("Trigger: " + triggerName).c_str());
                text.DrawLatex(0.15, 0.80, ("#Delta R < " + dR_cut_str).c_str());
            } else {
                // For other histograms, draw normally
                hist->Draw();
            }

            // Save the canvas as a PNG file
            canvas.SaveAs(outputPath.c_str());
            std::cout << "Saved: " << outputPath << std::endl;
        }
        delete obj;
    }
}

// Function to save all histograms in a ROOT file as PNGs in the same directory structure
void saveHistogramsAsPNGs(const std::string& inputFilePath, const std::string& outputBaseDir) {
    TFile* inputFile = TFile::Open(inputFilePath.c_str(), "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open the file " << inputFilePath << std::endl;
        return;
    }

    // Extract the run number from the file name
    std::string runNumber = inputFilePath.substr(inputFilePath.find_last_of("/") + 1, inputFilePath.find("_HistOutput.root") - inputFilePath.find_last_of("/") - 1);

    // Create the base output directory for this run number
    std::string runOutputDir = outputBaseDir + runNumber;
    gSystem->mkdir(runOutputDir.c_str(), true);

    // Start the recursive histogram saving from the ROOT file's top directory
    saveHistograms(inputFile, runOutputDir);

    inputFile->Close();
    delete inputFile;
}

// Function to process each run file in the output directory
void processRunFiles() {
    // Open the input directory
    void* dir = gSystem->OpenDirectory(inputDir.c_str());
    if (!dir) {
        std::cerr << "Error: Cannot open directory " << inputDir << std::endl;
        return;
    }

    const char* entry;
    while ((entry = gSystem->GetDirEntry(dir))) {
        std::string fileName = entry;
        if (fileName.find("_HistOutput.root") != std::string::npos) {
            std::string filePath = inputDir + fileName;
            saveHistogramsAsPNGs(filePath, outputDir);
        }
    }
    gSystem->FreeDirectory(dir);
}

int getDistinctColor(int index) {
    // Define a fixed set of visually distinct colors
    static std::vector<int> colors = {
        kRed, kBlue, kGreen + 2, kMagenta, kCyan + 1, kOrange + 7, kViolet + 1,
        kPink + 9, kSpring + 4, kTeal, kAzure + 7, kGray + 2, kBlack
    };

    // Cycle through the colors based on the index to avoid running out of colors
    return colors[index % colors.size()];
}

void OverlayCaloQA() {
    // Define trigger index for this specific analysis
    int triggerIndex = 10;
    
    // Get the corresponding trigger name
    std::string triggerName = getTriggerName(triggerIndex);

    // Set output file paths for the overlay plots, now including trigger name in file names
    std::string outputFileChi2 = outputDir + "Overlay_hClusterChi2_Trigger" + std::to_string(triggerIndex) + "_" + triggerName + ".png";
    std::string outputFileEEMCal = outputDir + "Overlay_hTotalCaloEEMCal_Trigger" + std::to_string(triggerIndex) + "_" + triggerName + ".png";
    std::string outputFileClusterPt = outputDir + "Overlay_hClusterPt_Trigger" + std::to_string(triggerIndex) + "_" + triggerName + ".png";

    // Create canvases for the three overlays
    TCanvas* canvasChi2 = new TCanvas("OverlayCaloQA_Chi2", "Overlay Calo QA Chi2", 800, 600);
    TCanvas* canvasEEMCal = new TCanvas("OverlayCaloQA_EEMCal", "Overlay Calo QA EEMCal", 800, 600);
    TCanvas* canvasClusterPt = new TCanvas("OverlayCaloQA_ClusterPt", "Overlay Calo QA Cluster Pt", 800, 600);
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);

    // Set logarithmic y-axis scale
    canvasChi2->SetLogy();
    canvasEEMCal->SetLogy();
    canvasClusterPt->SetLogy();

    // Create legends with properly set positions and sizes
    TLegend* legendChi2 = new TLegend(0.55, 0.55, 0.9, 0.9);
    TLegend* legendEEMCal = new TLegend(0.55, 0.55, 0.9, 0.9);
    TLegend* legendClusterPt = new TLegend(0.55, 0.55, 0.9, 0.9);
    legendChi2->SetNColumns(2);
    legendEEMCal->SetNColumns(2);
    legendClusterPt->SetNColumns(2);
    legendChi2->SetTextSize(0.025);
    legendEEMCal->SetTextSize(0.025);
    legendClusterPt->SetTextSize(0.025);
    legendChi2->SetTextFont(42);
    legendEEMCal->SetTextFont(42);
    legendClusterPt->SetTextFont(42);

    // Vectors to keep track of dummy lines and run numbers
    std::vector<TLine*> dummyLinesChi2;
    std::vector<TLine*> dummyLinesEEMCal;
    std::vector<TLine*> dummyLinesClusterPt;
    std::vector<std::pair<std::string, TLine*>> runEntriesChi2;
    std::vector<std::pair<std::string, TLine*>> runEntriesEEMCal;
    std::vector<std::pair<std::string, TLine*>> runEntriesClusterPt;

    // Open the input directory
    void* dir = gSystem->OpenDirectory(inputDir.c_str());
    if (!dir) {
        std::cerr << "Error: Cannot open directory " << inputDir << std::endl;
        return;
    }

    const char* entry;
    bool firstDrawChi2 = true;
    bool firstDrawEEMCal = true;
    bool firstDrawClusterPt = true;
    int colorIndex = 0; // Track the current color index

    // Iterate through each file in the directory
    while ((entry = gSystem->GetDirEntry(dir))) {
        std::string fileName = entry;

        // Skip the Final_Merged_HistOutput.root file
        if (fileName == "Final_Merged_HistOutput.root") {
            std::cout << "Skipping file: " << fileName << std::endl;
            continue;
        }

        if (fileName.find("_HistOutput.root") != std::string::npos) {
            std::string filePath = inputDir + fileName;
            std::cout << "Processing file: " << filePath << std::endl;

            // Open the ROOT file
            TFile* inputFile = TFile::Open(filePath.c_str(), "READ");
            if (!inputFile || inputFile->IsZombie()) {
                std::cerr << "Error: Could not open the file " << filePath << std::endl;
                continue;
            }

            // Construct the directory and histogram names for trigger index 10
            std::string dirName = "QA/Trigger" + std::to_string(triggerIndex);
            TDirectory* qaDir = (TDirectory*)inputFile->Get(dirName.c_str());
            if (!qaDir) {
                std::cerr << "Warning: Directory " << dirName << " not found in the file " << fileName << "." << std::endl;
                inputFile->Close();
                continue;
            }

            // Extract run number from the file name
            std::string runNumber = fileName.substr(0, fileName.find("_HistOutput.root"));

            // Overlay hClusterChi2_10
            TH1F* histChi2 = (TH1F*)qaDir->Get("hClusterChi2_10");
            if (histChi2 && histChi2->InheritsFrom(TH1::Class())) {
                histChi2->SetLineWidth(2);
                int color = getDistinctColor(colorIndex++);
                histChi2->SetLineColor(color);
                std::string titleChi2 = "Overlay of Cluster Chi2 Distributions for Trigger: " + triggerName;
                histChi2->SetTitle(titleChi2.c_str());
                histChi2->GetXaxis()->SetTitle("#chi^{2}");
                histChi2->GetYaxis()->SetTitle("Counts");

                // Draw histograms on the Chi2 canvas
                histChi2->SetDirectory(0); // Detach the histogram from the file
                canvasChi2->cd();
                if (firstDrawChi2) {
                    histChi2->Draw("HIST");
                    firstDrawChi2 = false;
                } else {
                    histChi2->Draw("HIST SAME");
                }

                // Create a persistent dummy line for the Chi2 legend
                TLine* dummyLine = new TLine();
                dummyLine->SetLineColor(color);
                dummyLine->SetLineWidth(2);
                dummyLinesChi2.push_back(dummyLine);
                runEntriesChi2.emplace_back(runNumber, dummyLine);
            }

            // Overlay hTotalCaloEEMCal_10
            TH1F* histEEMCal = (TH1F*)qaDir->Get("hTotalCaloEEMCal_10");
            if (histEEMCal && histEEMCal->InheritsFrom(TH1::Class())) {
                histEEMCal->SetLineWidth(2);
                int color = getDistinctColor(colorIndex++);
                histEEMCal->SetLineColor(color);
                std::string titleEEMCal = "Overlay of Total Calo EEMCal Distributions for Trigger: " + triggerName;
                histEEMCal->SetTitle(titleEEMCal.c_str());
                histEEMCal->GetXaxis()->SetTitle("Energy");
                histEEMCal->GetYaxis()->SetTitle("Counts");

                // Draw histograms on the EEMCal canvas
                histEEMCal->SetDirectory(0); // Detach the histogram from the file
                canvasEEMCal->cd();
                if (firstDrawEEMCal) {
                    histEEMCal->Draw("HIST");
                    firstDrawEEMCal = false;
                } else {
                    histEEMCal->Draw("HIST SAME");
                }

                // Create a persistent dummy line for the EEMCal legend
                TLine* dummyLine = new TLine();
                dummyLine->SetLineColor(color);
                dummyLine->SetLineWidth(2);
                dummyLinesEEMCal.push_back(dummyLine);
                runEntriesEEMCal.emplace_back(runNumber, dummyLine);
            }

            // Overlay hClusterPt_10
            TH1F* histClusterPt = (TH1F*)qaDir->Get("hClusterPt_10");
            if (histClusterPt && histClusterPt->InheritsFrom(TH1::Class())) {
                histClusterPt->SetLineWidth(2);
                int color = getDistinctColor(colorIndex++);
                histClusterPt->SetLineColor(color);
                std::string titleClusterPt = "Overlay of Cluster Pt Distributions for Trigger: " + triggerName;
                histClusterPt->SetTitle(titleClusterPt.c_str());
                histClusterPt->GetXaxis()->SetTitle("Pt [GeV]");
                histClusterPt->GetYaxis()->SetTitle("Counts");
                histClusterPt->GetXaxis()->SetRangeUser(0, 30); // Set the x-axis range to 0 - 30

                // Draw histograms on the Cluster Pt canvas
                histClusterPt->SetDirectory(0); // Detach the histogram from the file
                canvasClusterPt->cd();
                if (firstDrawClusterPt) {
                    histClusterPt->Draw("HIST");
                    firstDrawClusterPt = false;
                } else {
                    histClusterPt->Draw("HIST SAME");
                }

                // Create a persistent dummy line for the Cluster Pt legend
                TLine* dummyLine = new TLine();
                dummyLine->SetLineColor(color);
                dummyLine->SetLineWidth(2);
                dummyLinesClusterPt.push_back(dummyLine);
                runEntriesClusterPt.emplace_back(runNumber, dummyLine);
            }

            // Clean up the input file
            inputFile->Close();
        }
    }
    gSystem->FreeDirectory(dir);

    // Sort run entries by run number in increasing numerical order
    std::sort(runEntriesChi2.begin(), runEntriesChi2.end());
    std::sort(runEntriesEEMCal.begin(), runEntriesEEMCal.end());
    std::sort(runEntriesClusterPt.begin(), runEntriesClusterPt.end());

    // Add sorted entries to the legends
    for (const auto& entry : runEntriesChi2) {
        legendChi2->AddEntry(entry.second, ("Run: " + entry.first).c_str(), "l");
    }
    for (const auto& entry : runEntriesEEMCal) {
        legendEEMCal->AddEntry(entry.second, ("Run: " + entry.first).c_str(), "l");
    }
    for (const auto& entry : runEntriesClusterPt) {
        legendClusterPt->AddEntry(entry.second, ("Run: " + entry.first).c_str(), "l");
    }

    // Draw the legends and save the plots
    if (legendChi2->GetNRows() > 0) {
        canvasChi2->cd();
        legendChi2->Draw();
        canvasChi2->Update();
        canvasChi2->SaveAs(outputFileChi2.c_str());
        std::cout << "Saved overlay plot: " << outputFileChi2 << std::endl;
    } else {
        std::cerr << "Warning: Legend for Chi2 has no entries to draw." << std::endl;
    }

    if (legendEEMCal->GetNRows() > 0) {
        canvasEEMCal->cd();
        legendEEMCal->Draw();
        canvasEEMCal->Update();
        canvasEEMCal->SaveAs(outputFileEEMCal.c_str());
        std::cout << "Saved overlay plot: " << outputFileEEMCal << std::endl;
    } else {
        std::cerr << "Warning: Legend for EEMCal has no entries to draw." << std::endl;
    }

    if (legendClusterPt->GetNRows() > 0) {
        canvasClusterPt->cd();
        legendClusterPt->Draw();
        canvasClusterPt->Update();
        canvasClusterPt->SaveAs(outputFileClusterPt.c_str());
        std::cout << "Saved overlay plot: " << outputFileClusterPt << std::endl;
    } else {
        std::cerr << "Warning: Legend for Cluster Pt has no entries to draw." << std::endl;
    }

    // Clean up
    delete legendChi2;
    delete canvasChi2;
    delete legendEEMCal;
    delete canvasEEMCal;
    delete legendClusterPt;
    delete canvasClusterPt;

    // Delete all dummy lines
    for (auto line : dummyLinesChi2) {
        delete line;
    }
    for (auto line : dummyLinesEEMCal) {
        delete line;
    }
    for (auto line : dummyLinesClusterPt) {
        delete line;
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


void plotTriggerTurnOnCurves(TFile* inputFile, const std::string& histNameBase, const std::vector<int>& numeratorTriggers, int denominatorTrigger, const std::string& outputFilePath) {
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
        
        if (numeratorTrigger == 25) {
            amplitude = 1.2; slope = 1.0; xOffset = 2.5;
            amplitudeMin = 1.0; amplitudeMax = 1.3;
            slopeMin = 0.8; slopeMax = 1.2;
            xOffsetMin = 2.0; xOffsetMax = 3.0;
        } else if (numeratorTrigger == 26) {
            amplitude = 1.2; slope = 0.85; xOffset = 4.2;
            amplitudeMin = 1.0; amplitudeMax = 1.2;
            slopeMin = 0.8; slopeMax = 1.0;
            xOffsetMin = 4.0; xOffsetMax = 5.0;
        } else if (numeratorTrigger == 27) {
            amplitude = 1.1; slope = 0.75; xOffset = 5.5;
            amplitudeMin = 0.9; amplitudeMax = 1.2;
            slopeMin = 0.5; slopeMax = 0.8;
            xOffsetMin = 5.0; xOffsetMax = 6.0;
        }
        
        if (numeratorTrigger == 17) {
            amplitude = 1.2; slope = 1.1; xOffset = 10;
            amplitudeMin = 1.0; amplitudeMax = 1.3;
            slopeMin = 0.8; slopeMax = 1.2;
            xOffsetMin = 9.0; xOffsetMax = 11.0;
        } else if (numeratorTrigger == 18) {
            amplitude = 1.2; slope = 1.1; xOffset = 10;
            amplitudeMin = 1.0; amplitudeMax = 1.3;
            slopeMin = 0.8; slopeMax = 1.2;
            xOffsetMin = 9.0; xOffsetMax = 11.0;
        } else if (numeratorTrigger == 19) {
            amplitude = 1.2; slope = 1.1; xOffset = 10;
            amplitudeMin = 1.0; amplitudeMax = 1.3;
            slopeMin = 0.8; slopeMax = 1.2;
            xOffsetMin = 9.0; xOffsetMax = 11.0;
        }

        // Fit the ratio with a sigmoid function using the manual parameters and limits
        double xmin = ratioHist->GetXaxis()->GetXmin();
        double xmax = ratioHist->GetXaxis()->GetXmax();
        TF1* fitFunc = sigmoidFit("sigmoidFit_" + std::to_string(numeratorTrigger), xmin, xmax, amplitude, slope, xOffset,
                                  amplitudeMin, amplitudeMax, slopeMin, slopeMax, xOffsetMin, xOffsetMax);
        fitFunc->SetLineColor(colors[i % colors.size()]);  // Match the color of the fit to the histogram
        ratioHist->Fit(fitFunc, "R");

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


void plotOverlayHistograms() {
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

    // Set output file paths
    std::string outputDir = "/Users/patsfan753/Desktop/DirectPhotonAna/Plots/";
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
            triggerIndices = {10, 25, 26, 27};
        } else if (histNames[i] == "hCluster_maxECore_") {
            triggerIndices = {10, 25, 26, 27};
        } else if (histNames[i] == "h_hcal_energy_") {
            triggerIndices = {10, 17, 18, 19};
        } else if (histNames[i] == "h_jet_energy_") {
            triggerIndices = {10, 17, 18, 19};
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
     std::vector<int> h8by8NumeratorTriggers = {25, 26, 27};
     plotTriggerTurnOnCurves(inputFileNew, "h8by8TowerEnergySum_", h8by8NumeratorTriggers, 10, outputDir + "h8by8TowerEnergySum_TurnOn.png");

     // Plot ratios for h_jet_energy (triggers 17/10, 18/10, 19/10)
     std::vector<int> hJetNumeratorTriggers = {17, 18, 19};
     plotTriggerTurnOnCurves(inputFileNew, "h_jet_energy_", hJetNumeratorTriggers, 10, outputDir + "h_jet_energy_TurnOn.png");

     inputFileNew->Close();
     delete inputFileNew;
}



// Structure to hold the parsed cut values
struct CutValues {
    float clusECore = 0;
    float asymmetry = 0;
    float chi = 0;
    int triggerIndex = 0;
    float pTMin = -1;  // Default to -1 indicating no pT bin
    float pTMax = -1;  // Default to -1 indicating no pT bin
};

// Function to parse histogram names and extract cut values, trigger index, and optional pT range
CutValues parseHistName(const std::string& histName) {
    CutValues cuts;

    // Regex pattern that handles both cases: with and without pT bins
    std::regex re("invMass(?:_noPtBins)?_E([0-9]+(?:point[0-9]*)?)_Chi([0-9]+(?:point[0-9]*)?)_Asym([0-9]+(?:point[0-9]*)?)(?:_pT_([0-9]+(?:point[0-9]*)?)to([0-9]+(?:point[0-9]*)?))?_(\\d+)");
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

            cuts.triggerIndex = std::stoi(match[6].str());

            // Diagnostic prints
            std::cout << "Parsed histogram: " << histName << std::endl;
            std::cout << "  clusECore: " << cuts.clusECore << ", Chi: " << cuts.chi << ", Asymmetry: " << cuts.asymmetry
                      << ", pTMin: " << cuts.pTMin << ", pTMax: " << cuts.pTMax << ", Trigger Index: " << cuts.triggerIndex << std::endl;
        }
    } else {
        std::cerr << "Error: Failed to parse histogram name: " << histName << std::endl;
    }

    return cuts;
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


// Function to format numbers to three significant figures without scientific notation
std::string formatToThreeSigFigs(double value) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(3) << value; // Use fixed notation with three decimal places
    return stream.str();
}

void saveAnnotatedInvariantMassHistograms(const std::string& inputFilePath) {
    TFile* inputFile = TFile::Open(inputFilePath.c_str(), "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open the file " << inputFilePath << std::endl;
        return;
    }

    // Get the 'PhotonAnalysis' directory
    TDirectory* invMassDir = (TDirectory*)inputFile->Get("PhotonAnalysis");
    if (!invMassDir) {
        std::cerr << "Error: 'PhotonAnalysis' directory not found." << std::endl;
        inputFile->Close();
        return;
    }

    // Iterate through trigger directories
    TIter nextTriggerDir(invMassDir->GetListOfKeys());
    TKey* keyTriggerDir;
    while ((keyTriggerDir = (TKey*)nextTriggerDir())) {
        TObject* objTriggerDir = keyTriggerDir->ReadObj();
        if (objTriggerDir->InheritsFrom(TDirectory::Class())) {
            TDirectory* triggerDir = (TDirectory*)objTriggerDir;
            std::string triggerName = triggerDir->GetName();
            std::string outputTriggerDir = outputDir + "InvMass/" + triggerName;
            gSystem->mkdir(outputTriggerDir.c_str(), true);

            // Iterate through histograms in the trigger subdirectory
            TIter nextHist(triggerDir->GetListOfKeys());
            TKey* keyHist;
            while ((keyHist = (TKey*)nextHist())) {
                TObject* objHist = keyHist->ReadObj();
                if (objHist->InheritsFrom(TH1::Class())) {
                    TH1* hist = (TH1*)objHist;
                    std::string histName = hist->GetName();

                    // Check if the histogram has any entries
                    if (hist->GetEntries() == 0) {
                        std::cerr << "Skipping empty histogram: " << histName << std::endl;
                        delete objHist;
                        continue;  // Skip this histogram
                    }

                    // Filter histograms: only process those with 'invMass' in the name
                    if (histName.find("invMass_E") != 0 && histName.find("invMass_noPtBins_E") != 0) {
                        std::cerr << "Skipping histogram: " << histName << " (does not match invMass naming pattern)" << std::endl;
                        delete objHist;
                        continue;
                    }

                    // Parse the histogram name to get cut values and trigger index
                    CutValues cuts = parseHistName(histName);

                    // Construct the correct output directory based on parsed triggerIndex
                    std::string outputTriggerDir = outputDir + "InvMass/Trigger" + std::to_string(cuts.triggerIndex);
                    gSystem->mkdir(outputTriggerDir.c_str(), true);  // Create trigger-specific directory

                    // Construct the directory for the specific cut
                    std::ostringstream cutDirStream;
                    cutDirStream << outputTriggerDir << "/E" << formatToThreeSigFigs(cuts.clusECore)
                                 << "_Chi" << formatToThreeSigFigs(cuts.chi)
                                 << "_Asym" << formatToThreeSigFigs(cuts.asymmetry);
                    std::string cutDirPath = cutDirStream.str();
                    gSystem->mkdir(cutDirPath.c_str(), true);  // Create cut directory

                    // If there's a pT range, create the pT folder inside the cut-specific directory
                    std::string outputDirPath = cutDirPath;
                    if (cuts.pTMin != -1 && cuts.pTMax != -1) {
                        std::ostringstream ptDirStream;
                        ptDirStream << cutDirPath << "/pT_" << formatToThreeSigFigs(cuts.pTMin)
                                    << "_to_" << formatToThreeSigFigs(cuts.pTMax);
                        outputDirPath = ptDirStream.str();
                        gSystem->mkdir(outputDirPath.c_str(), true);  // Create pT directory
                    }

                    // Construct the output file path for the histogram PNG
                    std::string outputFilePath = outputDirPath + "/" + histName + ".png";


                    // Set axis labels
                    hist->GetXaxis()->SetTitle("M_{#gamma#gamma} [GeV]");
                    hist->GetYaxis()->SetTitle("Counts");

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
                    double sigmaPi0 = totalFit->GetParameter(2);
                    double meanEta = totalFit->GetParameter(4);
                    double sigmaEta = totalFit->GetParameter(5);

                    double massRatio = meanEta / meanPi0;
                    
                    
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
                        ptRangeLabel << "pT: " << formatToThreeSigFigs(cuts.pTMin) << " - " << formatToThreeSigFigs(cuts.pTMax) << " GeV";
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
                    std::string triggerNameLabel = getTriggerName(cuts.triggerIndex);
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

                    // Format in three columns: Trigger, pi0, and eta

                    // First column: Trigger information
                    labelText.DrawLatex(0.5, 0.90, "Trigger:");
                    valueText.DrawLatex(0.62, 0.90, triggerNameLabel.c_str());

                    labelText.DrawLatex(0.5, 0.85, "clusECore:");
                    valueText.DrawLatex(0.62, 0.85, formatToThreeSigFigs(cuts.clusECore).c_str());

                    labelText.DrawLatex(0.5, 0.80, "Chi2:");
                    valueText.DrawLatex(0.62, 0.80, formatToThreeSigFigs(cuts.chi).c_str());

                    labelText.DrawLatex(0.5, 0.75, "Asymmetry:");
                    valueText.DrawLatex(0.62, 0.75, formatToThreeSigFigs(cuts.asymmetry).c_str());

                    // If pT range is available, add it to the legend
                    if (!ptRangeLabel.str().empty()) {
                        labelText.DrawLatex(0.5, 0.70, "pT Range:");
                        valueText.DrawLatex(0.62, 0.70, ptRangeLabel.str().c_str());
                    }

                    // Second column: Pi0 information
                    labelText.DrawLatex(0.6, 0.55, "#pi^{0}:");
                    labelText.DrawLatex(0.8, 0.55, "#eta:");

                    labelText.DrawLatex(0.50, 0.5, "Mass:");
                    valueText.DrawLatex(0.60, 0.5, formatToThreeSigFigs(meanPi0).c_str());
                    valueText.DrawLatex(0.80, 0.5, formatToThreeSigFigs(meanEta).c_str());

                    labelText.DrawLatex(0.50, 0.45, "Sigma:");
                    valueText.DrawLatex(0.60, 0.45, formatToThreeSigFigs(sigmaPi0).c_str());
                    valueText.DrawLatex(0.80, 0.45, formatToThreeSigFigs(sigmaEta).c_str());

                    labelText.DrawLatex(0.50, 0.4, "S/B Ratio:");
                    valueText.DrawLatex(0.60, 0.4, formatToThreeSigFigs(signalToBackgroundPi0Ratio).c_str());
                    valueText.DrawLatex(0.80, 0.4, formatToThreeSigFigs(signalToBackgroundEtaRatio).c_str());

                    // Third row: Mass ratio
                    labelText.DrawLatex(0.5, 0.35, "Mass Ratio (#eta/#pi^{0}):");
                    valueText.DrawLatex(0.68, 0.35, formatToThreeSigFigs(massRatio).c_str());

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
                delete objHist;
            }
        }
        delete objTriggerDir;
    }

    inputFile->Close();
    delete inputFile;
}

void isolationEnergies(const std::string& inputFilePath, const std::string& outputDir) {
    // Open the input file
    TFile* inputFile = TFile::Open(inputFilePath.c_str(), "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open the file " << inputFilePath << std::endl;
        return;
    }

    // Get the 'QA' directory
    TDirectory* qaDir = (TDirectory*)inputFile->Get("QA");
    if (!qaDir) {
        std::cerr << "Error: 'QA' directory not found in file " << inputFilePath << std::endl;
        inputFile->Close();
        return;
    }

    std::cout << "QA directory found. Processing trigger subdirectories...\n";

    // Create the base output directory 'IsolationEnergies'
    std::string isolationEnergiesDir = outputDir + "/IsolationEnergies";
    gSystem->mkdir(isolationEnergiesDir.c_str(), true);

    // Iterate through trigger subdirectories in the 'QA' directory
    TIter nextTriggerDir(qaDir->GetListOfKeys());
    TKey* keyTriggerDir;
    while ((keyTriggerDir = (TKey*)nextTriggerDir())) {
        TObject* objTriggerDir = keyTriggerDir->ReadObj();

        if (objTriggerDir->InheritsFrom(TDirectory::Class())) {
            TDirectory* triggerDir = (TDirectory*)objTriggerDir;
            std::string triggerDirName = triggerDir->GetName();

            // Extract the trigger index from the directory name (e.g., "Trigger3")
            if (triggerDirName.find("Trigger") != 0) {
                std::cout << "Skipping non-trigger directory: " << triggerDirName << std::endl;
                delete objTriggerDir;
                continue;
            }

            int triggerIndex = std::stoi(triggerDirName.substr(7));  // Extract the index after "Trigger"
            std::string triggerName = getTriggerName(triggerIndex);

            std::cout << "Processing directory: " << triggerDirName << " for trigger " << triggerIndex << " (" << triggerName << ")\n";

            // Create the output directory for this trigger index
            std::string triggerOutputDir = isolationEnergiesDir + "/" + triggerDirName;
            gSystem->mkdir(triggerOutputDir.c_str(), true);

            // Iterate through histograms in the trigger directory
            TIter nextHist(triggerDir->GetListOfKeys());
            TKey* keyHist;
            while ((keyHist = (TKey*)nextHist())) {
                TObject* objHist = keyHist->ReadObj();

                // Only process histograms with specific naming patterns
                std::string histName = objHist->GetName();
                std::string prefix2D = "h2_cluster_iso_Ecore_";
                std::string prefix1D = "h1_isoEt_";

                if (histName.find(prefix2D) == 0) {
                        // It's a 2D histogram
                        std::string outputFilePath = triggerOutputDir + "/" + histName + ".png";

                        TH2* h2 = (TH2*)objHist;
                        
                        // Set the title for 2D histograms
                        h2->SetTitle(("Isolation Energy vs Ecore for " + triggerName).c_str());

                        // Disable the statistics box for 2D histograms
                        gStyle->SetOptStat(0);

                        // Set Z-axis range as requested
                        h2->GetZaxis()->SetRangeUser(0, 10000);

                        // Set Z-axis title
                        h2->GetZaxis()->SetTitle("Counts");

                        // Create the canvas and set the right margin to fit the color palette
                        TCanvas canvas;
                        canvas.SetRightMargin(0.18);
                        canvas.SetLogz();
                    
                    
                        // Draw the 2D histogram with color palette (Z-axis)
                        h2->Draw("COLZ");

                        TLatex latex;
                        latex.SetNDC();  // Use normalized coordinates (0-1 range)
                        latex.SetTextSize(0.04);  // Set text size
                        latex.DrawLatex(0.5, 0.8, ("Trigger: " + getTriggerName(triggerIndex)).c_str()); // Concatenate string


                        // Save the canvas as a PNG file
                        canvas.SaveAs(outputFilePath.c_str());
                        std::cout << "Saved 2D histogram: " << outputFilePath << std::endl;

                    } else if (histName.find(prefix1D) == 0) {
                        // It's a 1D histogram
                        std::string outputFilePath = triggerOutputDir + "/" + histName + ".png";

                        TH1* h1 = (TH1*)objHist;
                        // Set the title for 1D histograms
                        h1->SetTitle(("Isolation Energy for " + triggerName).c_str());
                        h1->GetXaxis()->SetRangeUser(-10, 10);
                        // Enable the statistics box for 1D histograms
                        gStyle->SetOptStat(1);
                        TCanvas canvas;

                        h1->Draw();
                        // Add TLatex for the trigger name (concatenate string manually)
                        TLatex latex;
                        latex.SetNDC();  // Use normalized coordinates (0-1 range)
                        latex.SetTextSize(0.04);  // Set text size
                        latex.DrawLatex(0.65, 0.8, ("Trigger: " + getTriggerName(triggerIndex)).c_str());  // Concatenate string


                        
                        canvas.SaveAs(outputFilePath.c_str());
                        std::cout << "Saved 1D histogram: " << outputFilePath << std::endl;
                    } else {
                        std::cout << "Skipping histogram: " << histName << " (does not match pattern)" << std::endl;
                    }

                    // Clean up the read object (if necessary)
                    delete objHist;
                }
            } else {
                std::cout << "Skipping non-directory object: " << objTriggerDir->GetName() << std::endl;
            }

            // Clean up the trigger directory object
            delete objTriggerDir;
        }

    // Close the input file
    inputFile->Close();
    delete inputFile;
}
void savePhotonAnalysisHistograms(const std::string& inputFilePath, const std::string& outputDir) {
    TFile* inputFile = TFile::Open(inputFilePath.c_str(), "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open the file " << inputFilePath << std::endl;
        return;
    }

    // Get the 'PhotonAnalysis' directory
    TDirectory* photonAnalysisDir = (TDirectory*)inputFile->Get("PhotonAnalysis");
    if (!photonAnalysisDir) {
        std::cerr << "Error: 'PhotonAnalysis' directory not found." << std::endl;
        inputFile->Close();
        return;
    }

    // Create the base output directory 'IsolationEnergies'
    std::string baseOutputDir = outputDir + "/IsolationEnergies";
    gSystem->mkdir(baseOutputDir.c_str(), true);

    // Iterate through trigger directories in the 'PhotonAnalysis' directory
    TIter nextTriggerDir(photonAnalysisDir->GetListOfKeys());
    TKey* keyTriggerDir;
    while ((keyTriggerDir = (TKey*)nextTriggerDir())) {
        TObject* objTriggerDir = keyTriggerDir->ReadObj();
        if (objTriggerDir->InheritsFrom(TDirectory::Class())) {
            TDirectory* triggerDir = (TDirectory*)objTriggerDir;
            std::string triggerName = triggerDir->GetName();
            int triggerIndex = std::stoi(triggerName.substr(7)); // Extract trigger index from "TriggerX"
            std::string outputTriggerDir = baseOutputDir + "/" + triggerName;
            gSystem->mkdir(outputTriggerDir.c_str(), true);

            // Iterate through histograms in the trigger subdirectory
            TIter nextHist(triggerDir->GetListOfKeys());
            TKey* keyHist;
            while ((keyHist = (TKey*)nextHist())) {
                TObject* objHist = keyHist->ReadObj();
                if (objHist->InheritsFrom(TH1::Class())) {
                    TH1* hist = (TH1*)objHist;
                    std::string histName = hist->GetName();

                    // Filter histograms: Process only those with isolatedPhotonCount, allPhotonCount, or ptPhoton in the name
                    if (histName.find("isolatedPhotonCount_E") != 0 &&
                        histName.find("allPhotonCount_E") != 0 &&
                        histName.find("ptPhoton_E") != 0) {
                        std::cerr << "Skipping histogram: " << histName << " (does not match expected naming pattern)" << std::endl;
                        delete objHist;
                        continue;
                    }

                    // Parse the histogram name to get cut values, pT range, and trigger index
                    CutValues cuts = parseHistName(histName);

                    // Construct the correct output directory based on parsed triggerIndex
                    std::string outputCutDir = outputTriggerDir + "/E" + formatToThreeSigFigs(cuts.clusECore)
                                              + "_Chi" + formatToThreeSigFigs(cuts.chi)
                                              + "_Asym" + formatToThreeSigFigs(cuts.asymmetry);
                    gSystem->mkdir(outputCutDir.c_str(), true);  // Create cut-specific directory

                    // Construct the pT directory inside the cut directory
                    std::string outputPtDir = outputCutDir + "/pT_" + formatToThreeSigFigs(cuts.pTMin)
                                             + "_to_" + formatToThreeSigFigs(cuts.pTMax);
                    gSystem->mkdir(outputPtDir.c_str(), true);  // Create pT-specific directory

                    // Construct the output file path for the histogram PNG
                    std::string outputFilePath = outputPtDir + "/" + histName + ".png";

                    // Create a canvas and draw the histogram
                    TCanvas canvas;
                    canvas.SetRightMargin(0.15);  // Leave space for Z-axis color bar if needed

                    // For 1D histograms, just draw and save the canvas
                    hist->Draw();

                    // Add TLatex to annotate the trigger index
                    TLatex latex;
                    latex.SetNDC();  // Use normalized coordinates (0-1 range)
                    latex.SetTextSize(0.04);  // Set text size
                    latex.DrawLatex(0.7, 0.85, ("Trigger: " + getTriggerName(cuts.triggerIndex)).c_str());  // Annotate with trigger name

                    // Save the canvas as a PNG file
                    canvas.SaveAs(outputFilePath.c_str());
                    std::cout << "Saved histogram: " << outputFilePath << std::endl;

                    // Clean up the histogram object
                    delete objHist;
                }
            }

            // Clean up the trigger directory object
            delete objTriggerDir;
        } else {
            std::cout << "Skipping non-directory object: " << objTriggerDir->GetName() << std::endl;
        }
    }

    // Close the input file
    inputFile->Close();
    delete inputFile;
}




// Main function to run both processes
void plotOutput() {
    gROOT->LoadMacro("sPhenixStyle.C");
    SetsPhenixStyle();
    
    if (saveHistogramsFlag) {
        // Run the function to process individual run files and save histograms as PNGs
        processRunFiles();
    }

    // Run the function to overlay histograms from the merged file
    plotOverlayHistograms();
    OverlayCaloQA();
//    saveAnnotatedInvariantMassHistograms(inputFilePath);
    isolationEnergies(inputFilePath, outputDir);
    savePhotonAnalysisHistograms(inputFilePath, outputDir);
}
