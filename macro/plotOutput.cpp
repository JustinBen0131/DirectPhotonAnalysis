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

std::string inputDir = "/Users/patsfan753/Desktop/DirectPhotonAna/";
std::string outputDir = "/Users/patsfan753/Desktop/DirectPhotonAna/Plots/";
std::string inputFilePath = inputDir + "Final_Merged_Hists_runnumber46623_runnumber47230.root";
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

std::string formatToThreeSigFigs(double value) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(3) << value; // Use fixed notation with three decimal places
    return stream.str();
}

// Function to parse photon histogram names and extract cut values, pT range, isoEt range, and trigger index
CutValues parsePhotonHistName(const std::string& histName) {
    CutValues cuts;

    // Updated regex pattern to match negative numbers and ensure proper grouping
    std::regex re("(isolatedPhotonCount|allPhotonCount|ptPhoton)_E(-?[0-9]+(?:point[0-9]*)?)_Chi(-?[0-9]+(?:point[0-9]*)?)_Asym(-?[0-9]+(?:point[0-9]*)?)"
                  "(?:_isoEt_(-?[0-9]+(?:point[0-9]*)?)to(-?[0-9]+(?:point[0-9]*)?))?_pT_(-?[0-9]+(?:point[0-9]*)?)to(-?[0-9]+(?:point[0-9]*)?)_(\\d+)");
    std::smatch match;

    // Lambda function to convert strings with 'point' and handle negative numbers
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
        if (match.size() >= 10) {
            cuts.clusECore = convert(match[2].str());
            cuts.chi = convert(match[3].str());
            cuts.asymmetry = convert(match[4].str());

            // If isoEt range is present, assign isoMin and isoMax, otherwise keep them as -1
            if (match[5].matched && match[6].matched) {
                cuts.isoMin = convert(match[5].str());
                cuts.isoMax = convert(match[6].str());
            }

            // pT bin range is always present
            cuts.pTMin = convert(match[7].str());
            cuts.pTMax = convert(match[8].str());

            // Trigger index
            cuts.triggerIndex = std::stoi(match[9].str());
        }
    } else {
        std::cerr << "Error: Failed to parse histogram name: " << histName << std::endl;
        cuts.invalid = true;
    }

    return cuts;
}


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
    std::string outputFileVtxZ = outputDir + "Overlay_hVtxZ_Trigger" + std::to_string(triggerIndex) + "_" + triggerName + ".png";
    std::string outputFileET = outputDir + "Overlay_hET_Trigger" + std::to_string(triggerIndex) + "_" + triggerName + ".png";


    // Create canvases for the three overlays
    TCanvas* canvasChi2 = new TCanvas("OverlayCaloQA_Chi2", "Overlay Calo QA Chi2", 800, 600);
    TCanvas* canvasEEMCal = new TCanvas("OverlayCaloQA_EEMCal", "Overlay Calo QA EEMCal", 800, 600);
    TCanvas* canvasClusterPt = new TCanvas("OverlayCaloQA_ClusterPt", "Overlay Calo QA Cluster Pt", 800, 600);
    TCanvas* canvasVtxZ = new TCanvas("OverlayCaloQA_VtxZ", "Overlay Calo QA VtxZ", 800, 600);
    TCanvas* canvasET = new TCanvas("OverlayCaloQA_ET", "Overlay Calo QA ET", 800, 600);
     
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);

    // Set logarithmic y-axis scale
    canvasChi2->SetLogy();
    canvasEEMCal->SetLogy();
    canvasClusterPt->SetLogy();
    canvasET->SetLogy();

    // Create legends with properly set positions and sizes
    TLegend* legendChi2 = new TLegend(0.55, 0.55, 0.9, 0.9);
    TLegend* legendEEMCal = new TLegend(0.55, 0.55, 0.9, 0.9);
    TLegend* legendClusterPt = new TLegend(0.55, 0.55, 0.9, 0.9);
    TLegend* legendVtxZ = new TLegend(0.18, 0.55, 0.5, 0.9);
    TLegend* legendET = new TLegend(0.55, 0.55, 0.9, 0.9);
    

    legendChi2->SetNColumns(2);
    legendEEMCal->SetNColumns(2);
    legendClusterPt->SetNColumns(2);
    legendVtxZ->SetNColumns(2);
    legendET->SetNColumns(2);

    legendChi2->SetTextSize(0.025);
    legendEEMCal->SetTextSize(0.025);
    legendClusterPt->SetTextSize(0.025);
    legendVtxZ->SetTextSize(0.025);
    legendET->SetTextSize(0.025);

    legendChi2->SetTextFont(42);
    legendEEMCal->SetTextFont(42);
    legendClusterPt->SetTextFont(42);
    legendVtxZ->SetTextFont(42);
    legendET->SetTextFont(42);

    // Vectors to keep track of dummy lines and run numbers
    std::vector<TLine*> dummyLinesChi2, dummyLinesEEMCal, dummyLinesClusterPt, dummyLinesVtxZ, dummyLinesET;
    std::vector<std::pair<std::string, TLine*>> runEntriesChi2, runEntriesEEMCal, runEntriesClusterPt, runEntriesVtxZ, runEntriesET;

    // Open the input directory
    void* dir = gSystem->OpenDirectory(inputDir.c_str());
    if (!dir) {
        std::cerr << "Error: Cannot open directory " << inputDir << std::endl;
        return;
    }

    const char* entry;
    bool firstDrawChi2 = true, firstDrawEEMCal = true, firstDrawClusterPt = true, firstDrawVtxZ = true, firstDrawET = true;
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
                histChi2->GetXaxis()->SetTitle("#chi^{2}");
                histChi2->GetYaxis()->SetTitle("Counts");

                // Draw histograms on the Chi2 canvas
                histChi2->SetDirectory(0);
                canvasChi2->cd();
                if (firstDrawChi2) {
                    histChi2->Draw("HIST");
                    firstDrawChi2 = false;
                } else {
                    histChi2->Draw("HIST SAME");
                }

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
                histEEMCal->GetXaxis()->SetTitle("Energy");
                histEEMCal->GetYaxis()->SetTitle("Counts");

                // Draw histograms on the EEMCal canvas
                histEEMCal->SetDirectory(0);
                canvasEEMCal->cd();
                if (firstDrawEEMCal) {
                    histEEMCal->Draw("HIST");
                    firstDrawEEMCal = false;
                } else {
                    histEEMCal->Draw("HIST SAME");
                }

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
                histClusterPt->GetXaxis()->SetTitle("Pt [GeV]");
                histClusterPt->GetYaxis()->SetTitle("Counts");
                histClusterPt->GetXaxis()->SetRangeUser(0, 30);

                // Draw histograms on the Cluster Pt canvas
                histClusterPt->SetDirectory(0);
                canvasClusterPt->cd();
                if (firstDrawClusterPt) {
                    histClusterPt->Draw("HIST");
                    firstDrawClusterPt = false;
                } else {
                    histClusterPt->Draw("HIST SAME");
                }

                TLine* dummyLine = new TLine();
                dummyLine->SetLineColor(color);
                dummyLine->SetLineWidth(2);
                dummyLinesClusterPt.push_back(dummyLine);
                runEntriesClusterPt.emplace_back(runNumber, dummyLine);
            }

            // Overlay h_ET_10
            TH1F* histET = (TH1F*)qaDir->Get("h_ET_10");
            if (histET && histET->InheritsFrom(TH1::Class())) {
                histET->SetLineWidth(2);
                int color = getDistinctColor(colorIndex++);
                histET->SetLineColor(color);
                histET->GetXaxis()->SetTitle("Energy [GeV]");
                histET->GetYaxis()->SetTitle("Counts");

                // Draw histograms on the ET canvas
                histET->SetDirectory(0);
                canvasET->cd();
                if (firstDrawET) {
                    histET->Draw("HIST");
                    firstDrawET = false;
                } else {
                    histET->Draw("HIST SAME");
                }

                TLine* dummyLine = new TLine();
                dummyLine->SetLineColor(color);
                dummyLine->SetLineWidth(2);
                dummyLinesET.push_back(dummyLine);
                runEntriesET.emplace_back(runNumber, dummyLine);
            }

            // Overlay hVtxZ_10
             TH1F* histVtxZ = (TH1F*)qaDir->Get("hVtxZ_10");
             if (histVtxZ && histVtxZ->InheritsFrom(TH1::Class())) {
                 histVtxZ->SetLineWidth(2);
                 int color = getDistinctColor(colorIndex++);
                 histVtxZ->SetLineColor(color);
                 histVtxZ->GetXaxis()->SetTitle("z [cm]");
                 histVtxZ->GetYaxis()->SetTitle("Counts");
                 
                 // Adjust y-axis to ensure all plots fit
                 histVtxZ->SetMaximum(histVtxZ->GetMaximum() * 1.75); // Increase y-axis range by 20%

                 // Draw histograms on the VtxZ canvas
                 histVtxZ->SetDirectory(0);
                 canvasVtxZ->cd();
                 if (firstDrawVtxZ) {
                     histVtxZ->Draw("HIST");
                     firstDrawVtxZ = false;
                 } else {
                     histVtxZ->Draw("HIST SAME");
                 }

                 TLine* dummyLine = new TLine();
                 dummyLine->SetLineColor(color);
                 dummyLine->SetLineWidth(2);
                 dummyLinesVtxZ.push_back(dummyLine);
                 runEntriesVtxZ.emplace_back(runNumber, dummyLine);
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
    std::sort(runEntriesVtxZ.begin(), runEntriesVtxZ.end());
    std::sort(runEntriesET.begin(), runEntriesET.end());

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
    for (const auto& entry : runEntriesVtxZ) {
        legendVtxZ->AddEntry(entry.second, ("Run: " + entry.first).c_str(), "l");
    }
    for (const auto& entry : runEntriesET) {
        legendET->AddEntry(entry.second, ("Run: " + entry.first).c_str(), "l");
    }

    // Draw the legends and save the plots for each canvas
    if (legendChi2->GetNRows() > 0) {
        canvasChi2->cd();
        legendChi2->Draw();
        canvasChi2->Update();
        canvasChi2->SaveAs(outputFileChi2.c_str());
        std::cout << "Saved overlay plot: " << outputFileChi2 << std::endl;
    }

    if (legendEEMCal->GetNRows() > 0) {
        canvasEEMCal->cd();
        legendEEMCal->Draw();
        canvasEEMCal->Update();
        canvasEEMCal->SaveAs(outputFileEEMCal.c_str());
        std::cout << "Saved overlay plot: " << outputFileEEMCal << std::endl;
    }

    if (legendClusterPt->GetNRows() > 0) {
        canvasClusterPt->cd();
        legendClusterPt->Draw();
        canvasClusterPt->Update();
        canvasClusterPt->SaveAs(outputFileClusterPt.c_str());
        std::cout << "Saved overlay plot: " << outputFileClusterPt << std::endl;
    }

    if (legendVtxZ->GetNRows() > 0) {
        canvasVtxZ->cd();
        legendVtxZ->Draw();
        canvasVtxZ->Update();
        canvasVtxZ->SaveAs(outputFileVtxZ.c_str());
        std::cout << "Saved overlay plot: " << outputFileVtxZ << std::endl;
    }

    if (legendET->GetNRows() > 0) {
        canvasET->cd();
        legendET->Draw();
        canvasET->Update();
        canvasET->SaveAs(outputFileET.c_str());
        std::cout << "Saved overlay plot: " << outputFileET << std::endl;
    }

    // Clean up
    delete legendChi2;
    delete canvasChi2;
    delete legendEEMCal;
    delete canvasEEMCal;
    delete legendClusterPt;
    delete canvasClusterPt;
    delete legendVtxZ;
    delete canvasVtxZ;
    delete legendET;
    delete canvasET;

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
    for (auto line : dummyLinesVtxZ) {
        delete line;
    }
    for (auto line : dummyLinesET) {
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
   // OverlayCaloQA();
}
