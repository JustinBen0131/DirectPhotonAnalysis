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

std::string inputDir = "/Users/patsfan753/Desktop/DirectPhotonAna/";
std::string outputDir = "/Users/patsfan753/Desktop/DirectPhotonAna/Plots/";
std::string inputFilePath = inputDir + "Final_Merged_Hists_runnumber46623_runnumber47230.root";

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


std::string formatToThreeSigFigs(double value) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(3) << value; // Use fixed notation with three decimal places
    return stream.str();
}

struct CutValues {
    float clusECore = -1;   // Default to -1 indicating no specific value
    float chi = -1;         // Default to -1 indicating no specific value
    float asymmetry = -1;   // Default to -1 indicating no specific value
    float isoMin = -1;      // Default to -1 indicating no isoEt range
    float isoMax = -1;      // Default to -1 indicating no isoEt range
    float pTMin = -1;       // Default to -1 indicating no pT range
    float pTMax = -1;       // Default to -1 indicating no pT range
    int triggerIndex = -1;  // Default to -1 indicating no trigger index
    std::string massWindowLabel;
};

CutValues parseIsolationQAHistName(const std::string& histName) {
    CutValues cuts;
    
    std::regex re("(?:_E([0-9]+(?:point[0-9]*)?)_Chi([0-9]+(?:point[0-9]*)?)_Asym([0-9]+(?:point[0-9]*)?)"
                  "(?:_(inMassWindow|outsideMassWindow))?"
                  "(?:_isoEt_(-?[0-9]+(?:point[0-9]*)?)to(-?[0-9]+(?:point[0-9]*)?))?)?_pT_([0-9]+(?:point[0-9]*)?)to([0-9]+(?:point[0-9]*)?)_(\\d+)|_(\\d+)");

    std::smatch match;

    // Lambda function to convert strings with 'point' to float values
    auto convert = [](const std::string& input) -> float {
        std::string temp = input;
        size_t pointPos = temp.find("point");
        if (pointPos != std::string::npos) {
            temp.replace(pointPos, 5, ".");
        }
        return std::stof(temp);
    };

    // Check if the regex matches the histogram name
    if (std::regex_search(histName, match, re)) {
        if (match.size() >= 10) {
            // Handle detailed histograms with cuts, mass window, pT range, and triggerIndex
            if (match[1].matched) {
                cuts.clusECore = convert(match[1].str());
                cuts.chi = convert(match[2].str());
                cuts.asymmetry = convert(match[3].str());

                if (match[4].matched) {
                    cuts.massWindowLabel = match[4].str();  // Capture massWindowLabel if present
                }

                if (match[5].matched && match[6].matched) {
                    cuts.isoMin = convert(match[5].str());
                    cuts.isoMax = convert(match[6].str());
                }

                cuts.pTMin = convert(match[7].str());
                cuts.pTMax = convert(match[8].str());
                cuts.triggerIndex = std::stoi(match[9].str());
            }
            // Handle simple histograms with just pT range and triggerIndex
            else if (match[7].matched && match[8].matched) {
                cuts.pTMin = convert(match[7].str());
                cuts.pTMax = convert(match[8].str());
                cuts.triggerIndex = std::stoi(match[9].str());
            }
            // Handle histograms with just triggerIndex
            else if (match[10].matched) {
                cuts.triggerIndex = std::stoi(match[10].str());
            }

            // Diagnostic prints
            std::cout << "Parsed histogram: " << histName << std::endl;
            std::cout << "  clusECore: " << cuts.clusECore << ", Chi: " << cuts.chi
                      << ", Asymmetry: " << cuts.asymmetry << ", isoMin: " << cuts.isoMin
                      << ", isoMax: " << cuts.isoMax << ", pTMin: " << cuts.pTMin
                      << ", pTMax: " << cuts.pTMax << ", Trigger Index: " << cuts.triggerIndex
                      << ", Mass Window Label: " << cuts.massWindowLabel << std::endl;
        }
    } else {
        std::cerr << "Error: Failed to parse histogram name: " << histName << std::endl;
    }

    return cuts;
}


// Helper function to get the trigger name based on index
std::string getTriggerName(int triggerIndex) {
    if (triggerNameMap.find(triggerIndex) != triggerNameMap.end()) {
        return triggerNameMap[triggerIndex];
    }
    return "Unknown Trigger";
}



void saveHistogram(TObject* objHist, const std::string& outputFilePath, const std::string& triggerName, const CutValues& cuts) {
    TCanvas canvas;
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);

    if (objHist->InheritsFrom(TH2::Class())) {
        TH2* h2 = (TH2*)objHist;
        h2->SetTitle(("Isolation Energy vs Ecore for " + triggerName).c_str());
        gStyle->SetOptStat(0);
        h2->GetZaxis()->SetRangeUser(0, 10000);
        h2->GetZaxis()->SetTitle("Counts");
        canvas.SetRightMargin(0.18);
        canvas.SetLogz();
        h2->Draw("COLZ");

        // Output Trigger Name
        latex.DrawLatex(0.2, 0.8, ("Trigger: " + triggerName).c_str());

        // Output pT range if valid
        if (cuts.pTMin >= 0 && cuts.pTMax >= 0) {
            latex.DrawLatex(0.2, 0.75, ("pT: " + formatToThreeSigFigs(cuts.pTMin) + " to " + formatToThreeSigFigs(cuts.pTMax)).c_str());
        }

        // Output cut combination if valid
        if (cuts.clusECore >= 0 && cuts.chi >= 0 && cuts.asymmetry >= 0) {
            latex.DrawLatex(0.2, 0.7, ("Cuts: ECore=" + formatToThreeSigFigs(cuts.clusECore) + ", Chi=" + formatToThreeSigFigs(cuts.chi) +
                                      ", Asym=" + formatToThreeSigFigs(cuts.asymmetry)).c_str());
        }

    } else if (objHist->InheritsFrom(TH1::Class())) {
        TH1* h1 = (TH1*)objHist;
        h1->SetTitle(("Isolation Energy for " + triggerName).c_str());
        h1->GetXaxis()->SetRangeUser(-10, 10);
        gStyle->SetOptStat(1);
        h1->Draw();

        // Output Trigger Name
        latex.DrawLatex(0.18, 0.8, ("Trigger: " + triggerName).c_str());

        // Output pT range if valid
        if (cuts.pTMin >= 0 && cuts.pTMax >= 0) {
            latex.DrawLatex(0.18, 0.75, ("pT: " + formatToThreeSigFigs(cuts.pTMin) + " to " + formatToThreeSigFigs(cuts.pTMax)).c_str());
        }

        // Output cut combination if valid
        if (cuts.clusECore >= 0 && cuts.chi >= 0 && cuts.asymmetry >= 0) {
            latex.DrawLatex(0.18, 0.7, ("Cuts: ECore=" + formatToThreeSigFigs(cuts.clusECore) + ", Chi=" + formatToThreeSigFigs(cuts.chi) +
                                      ", Asym=" + formatToThreeSigFigs(cuts.asymmetry)).c_str());
        }
    }

    canvas.SaveAs(outputFilePath.c_str());
    std::cout << "Saved histogram: " << outputFilePath << std::endl;
}

struct IsolatedPhotonLog {
    int triggerIndex;
    float clusECore;
    float chi;
    float asymmetry;
    float pTMin;
    float pTMax;
    float isoMin;
    float isoMax;
    int isolatedEntries;
    std::string massWindowLabel;
};

struct TotalPhotonLog {
    int triggerIndex;
    float clusECore;
    float chi;
    float asymmetry;
    float pTMin;
    float pTMax;
    int totalEntries;
    std::string massWindowLabel;
};

struct PtWeightingLog {
    int triggerIndex;
    float clusECore;
    float chi;
    float asymmetry;
    float pTMin;
    float pTMax;
    double weightedAveragePt;
    std::string massWindowLabel;
};

std::map<std::tuple<int, float, float, float, float, float, float, float, std::string>, IsolatedPhotonLog> isolatedPhotonMap;
std::map<std::tuple<int, float, float, float, float, float, std::string>, TotalPhotonLog> totalPhotonMap;
std::map<std::tuple<int, float, float, float, float, float, std::string>, PtWeightingLog> pTweightingMap;


void processDirectory(TDirectory* baseDir, const std::string& isolationEnergiesDir) {
    std::string dirName = baseDir->GetName();
    std::cout << "Processing directory: " << dirName << "\n";

    // Iterate through trigger subdirectories in the base directory
    TIter nextTriggerDir(baseDir->GetListOfKeys());
    TKey* keyTriggerDir;
    while ((keyTriggerDir = (TKey*)nextTriggerDir())) {
        TObject* objTriggerDir = keyTriggerDir->ReadObj();
        // Check if the object is a directory
        if (objTriggerDir->InheritsFrom(TDirectory::Class())) {
            TDirectory* triggerDir = (TDirectory*)objTriggerDir;
            std::string triggerDirName = triggerDir->GetName();

            // Extract the trigger index from the directory name (e.g., "Trigger3")
            if (triggerDirName.find("Trigger") != 0) {
                std::cout << "Skipping non-trigger directory: " << triggerDirName << std::endl;
                delete objTriggerDir;
                continue;
            }
            int triggerIndex = std::stoi(triggerDirName.substr(7));
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

                // Get the histogram name
                std::string histName = objHist->GetName();

                // Only process histograms that follow the naming scheme
                if (histName.find("h2_cluster_iso_Ecore_") == 0 || histName.find("h1_isoEt_") == 0 || histName.find("isolatedPhotonCount_E") == 0 || histName.find("allPhotonCount_E") == 0 || histName.find("ptPhoton_E") == 0) {
                    std::cout << "Processing histogram: " << histName << std::endl;

                    // Parse the histogram name for the cuts, pT range, and isoEt range (for isolatedPhotonCount_E)
                    CutValues cuts = parseIsolationQAHistName(histName);
                    // Print debug information about the parsed cut values
                    std::cout << "Parsed Cut Values: "
                              << "clusECore: " << cuts.clusECore
                              << ", chi: " << cuts.chi
                              << ", asymmetry: " << cuts.asymmetry
                              << ", pTMin: " << cuts.pTMin
                              << ", pTMax: " << cuts.pTMax
                              << ", isoMin: " << cuts.isoMin
                              << ", isoMax: " << cuts.isoMax << "\n";


                    auto baseKey = std::make_tuple(triggerIndex, cuts.clusECore, cuts.chi, cuts.asymmetry, cuts.pTMin, cuts.pTMax, cuts.massWindowLabel);
                    auto isoKey = std::make_tuple(triggerIndex, cuts.clusECore, cuts.chi, cuts.asymmetry, cuts.pTMin, cuts.pTMax, cuts.isoMin, cuts.isoMax, cuts.massWindowLabel);

                    if (histName.find("isolatedPhotonCount_E") == 0) {
                        if (objHist->InheritsFrom(TH1::Class())) {
                            IsolatedPhotonLog& log = isolatedPhotonMap[isoKey];
                            log.triggerIndex = triggerIndex;
                            log.clusECore = cuts.clusECore;
                            log.chi = cuts.chi;
                            log.asymmetry = cuts.asymmetry;
                            log.pTMin = cuts.pTMin;
                            log.pTMax = cuts.pTMax;
                            log.isoMin = cuts.isoMin;
                            log.isoMax = cuts.isoMax;
                            log.massWindowLabel = cuts.massWindowLabel;
                            log.isolatedEntries = ((TH1*)objHist)->GetEntries();

                            std::cout << "  Isolated Photon Count: " << log.isolatedEntries << "\n";
                            std::cout << "  Mass Window Label: " << log.massWindowLabel << "\n";
                        }
                    }

                    if (histName.find("allPhotonCount_E") == 0) {
                        if (objHist->InheritsFrom(TH1::Class())) {
                            TotalPhotonLog& log = totalPhotonMap[baseKey];
                            log.triggerIndex = triggerIndex;
                            log.clusECore = cuts.clusECore;
                            log.chi = cuts.chi;
                            log.asymmetry = cuts.asymmetry;
                            log.pTMin = cuts.pTMin;
                            log.pTMax = cuts.pTMax;
                            log.massWindowLabel = cuts.massWindowLabel;
                            log.totalEntries = ((TH1*)objHist)->GetEntries();

                            std::cout << "  Total Photon Count: " << log.totalEntries << "\n";
                            std::cout << "  Mass Window Label: " << log.massWindowLabel << "\n";
                        }
                    }
                    if (histName.find("ptPhoton_E") == 0) {
                          double weightedSumPt = 0;
                          double totalPhotonCounts = 0;
                          int nBins = ((TH1*)objHist)->GetNbinsX();

                          for (int bin = 1; bin <= nBins; ++bin) {
                              double photonCount = ((TH1*)objHist)->GetBinContent(bin);
                              double pt = ((TH1*)objHist)->GetBinCenter(bin);
                              weightedSumPt += pt * photonCount;
                              totalPhotonCounts += photonCount;
                          }

                          PtWeightingLog& log = pTweightingMap[baseKey];
                          log.triggerIndex = triggerIndex;
                          log.clusECore = cuts.clusECore;
                          log.chi = cuts.chi;
                          log.asymmetry = cuts.asymmetry;
                          log.pTMin = cuts.pTMin;
                          log.pTMax = cuts.pTMax;
                          log.massWindowLabel = cuts.massWindowLabel;
                          log.weightedAveragePt = (totalPhotonCounts > 0) ? weightedSumPt / totalPhotonCounts : 0;

                          std::cout << "  Weighted Average pT: " << log.weightedAveragePt << "\n";
                          std::cout << "  Mass Window Label: " << log.massWindowLabel << "\n";
                    }
                    // Check if we have a valid cut variation
                    if (cuts.clusECore >= 0 && cuts.chi >= 0 && cuts.asymmetry >= 0) {
                        // Folder structure for the cut variation
                        std::string cutDir = triggerOutputDir + "/E_" + formatToThreeSigFigs(cuts.clusECore) +
                                             "_Chi" + formatToThreeSigFigs(cuts.chi) +
                                             "_Asym" + formatToThreeSigFigs(cuts.asymmetry);
                        gSystem->mkdir(cutDir.c_str(), true);

                        // Check if we have a valid pT range
                        if (cuts.pTMin >= 0 && cuts.pTMax >= 0) {
                            // Subdirectory for the pT range under the cut variation folder
                            std::string pTBinDir = cutDir + "/pT_" + formatToThreeSigFigs(cuts.pTMin) + "to" + formatToThreeSigFigs(cuts.pTMax);
                            gSystem->mkdir(pTBinDir.c_str(), true);

                            // Handle isolatedPhotonCount_E histograms with isoEt-specific subfolder
                            if (histName.find("isolatedPhotonCount_E") == 0 && cuts.isoMin != -1 && cuts.isoMax != -1) {
                                std::string isoDir = pTBinDir + "/isoEt_" + formatToThreeSigFigs(cuts.isoMin) + "to" + formatToThreeSigFigs(cuts.isoMax);
                                gSystem->mkdir(isoDir.c_str(), true);

                                // Save the isolated photon histogram in the isoEt subfolder
                                std::string outputFilePath = isoDir + "/" + histName + ".png";
                                saveHistogram(objHist, outputFilePath, triggerName, cuts);

                            } else if (histName.find("allPhotonCount_E") == 0 || histName.find("ptPhoton_E") == 0) {
                                // Handle allPhotonCount_E and ptPhoton_E histograms (no isoEt)
                                std::string outputFilePath = pTBinDir + "/" + histName + ".png";
                                saveHistogram(objHist, outputFilePath, triggerName, cuts);

                            } else {
                                // Save other histograms (non-isoEt) in the pT range folder
                                std::string outputFilePath = pTBinDir + "/" + histName + ".png";
                                saveHistogram(objHist, outputFilePath, triggerName, cuts);
                            }
                        } else {
                            // Save the cut variation histogram directly if no pT range is provided
                            std::string outputFilePath = cutDir + "/" + histName + ".png";
                            saveHistogram(objHist, outputFilePath, triggerName, cuts);
                        }
                    } else {
                        std::string outputFilePath = triggerOutputDir + "/" + histName + ".png";
                        saveHistogram(objHist, outputFilePath, triggerName, cuts);
                    }
                } else {
                    std::cout << "Skipping histogram: " << histName << " (does not match required naming pattern)" << std::endl;
                }
                delete objHist;  // Clean up after processing
            }
        } else {
            std::cout << "Skipping non-directory object: " << objTriggerDir->GetName() << std::endl;
        }

        delete objTriggerDir;
    }
}
void outputHistogramLogs(const std::string& outputFilePath) {
    std::cout << "\n================= Isolated Photon Map =================\n";
    std::cout << "TriggerIndex | clusECore | Chi | Asymmetry | pTMin | pTMax | isoMin | isoMax | IsolatedEntries | MassWindowLabel\n";
    std::cout << "----------------------------------------------------------------------------------------------\n";
    for (const auto& entry : isolatedPhotonMap) {
        auto key = entry.first;
        const IsolatedPhotonLog& log = entry.second;
        std::cout << std::get<0>(key) << " | "
                  << std::get<1>(key) << " | "
                  << std::get<2>(key) << " | "
                  << std::get<3>(key) << " | "
                  << std::get<4>(key) << " | "
                  << std::get<5>(key) << " | "
                  << std::get<6>(key) << " | "
                  << std::get<7>(key) << " | "
                  << log.isolatedEntries << " | "
                  << log.massWindowLabel << "\n";
    }

    std::cout << "\n================= Total Photon Map =================\n";
    std::cout << "TriggerIndex | clusECore | Chi | Asymmetry | pTMin | pTMax | TotalEntries | MassWindowLabel\n";
    std::cout << "------------------------------------------------------------------------------------\n";
    for (const auto& entry : totalPhotonMap) {
        auto key = entry.first;
        const TotalPhotonLog& log = entry.second;
        std::cout << std::get<0>(key) << " | "
                  << std::get<1>(key) << " | "
                  << std::get<2>(key) << " | "
                  << std::get<3>(key) << " | "
                  << std::get<4>(key) << " | "
                  << std::get<5>(key) << " | "
                  << log.totalEntries << " | "
                  << log.massWindowLabel << "\n";
    }

    std::cout << "\n================= Weighted pT Map =================\n";
    std::cout << "TriggerIndex | clusECore | Chi | Asymmetry | pTMin | pTMax | WeightedAveragePt | MassWindowLabel\n";
    std::cout << "------------------------------------------------------------------------------------\n";
    for (const auto& entry : pTweightingMap) {
        auto key = entry.first;
        const PtWeightingLog& log = entry.second;
        std::cout << std::get<0>(key) << " | "
                  << std::get<1>(key) << " | "
                  << std::get<2>(key) << " | "
                  << std::get<3>(key) << " | "
                  << std::get<4>(key) << " | "
                  << std::get<5>(key) << " | "
                  << log.weightedAveragePt << " | "
                  << log.massWindowLabel << "\n";
    }

    std::ofstream outFile(outputFilePath);

    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open the file " << outputFilePath << " for writing.\n";
        return;
    }

    // Write CSV header
    outFile << "Trigger,ECore,Chi,Asymmetry,pT Min,pT Max,isoMin,isoMax,Isolated Counts,Total Counts,Isolated/Total,Statistical Error,Weighted pT,MassWindowLabel\n";

    // Iterate through isolatedPhotonMap and correlate with totalPhotonMap and pTweightingMap
    for (const auto& isoEntry : isolatedPhotonMap) {
        auto isoKey = isoEntry.first;
        const IsolatedPhotonLog& isoLog = isoEntry.second;

        // Extract the common key for totalPhotonMap and pTweightingMap (including massWindowLabel)
        auto commonKey = std::make_tuple(
            std::get<0>(isoKey),  // triggerIndex
            std::get<1>(isoKey),  // clusECore
            std::get<2>(isoKey),  // chi
            std::get<3>(isoKey),  // asymmetry
            std::get<4>(isoKey),  // pTMin
            std::get<5>(isoKey),  // pTMax
            isoLog.massWindowLabel // massWindowLabel
        );

        // Find corresponding entries in totalPhotonMap and pTweightingMap
        auto totalEntry = totalPhotonMap.find(commonKey);
        auto pTweightingEntry = pTweightingMap.find(commonKey);

        // Ensure corresponding entries exist in both maps
        if (totalEntry != totalPhotonMap.end() && pTweightingEntry != pTweightingMap.end()) {
            const TotalPhotonLog& totalLog = totalEntry->second;
            const PtWeightingLog& pTLog = pTweightingEntry->second;

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

            // Write CSV row
            outFile << isoLog.triggerIndex << ","
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
                    << isoLog.massWindowLabel << "\n";
        }
    }

    outFile.close();
    std::cout << "CSV file successfully written to " << outputFilePath << "\n";
}



void isolationEnergies(const std::string& inputFilePath, const std::string& outputDir) {
    // Open the input file
    TFile* inputFile = TFile::Open(inputFilePath.c_str(), "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open the file " << inputFilePath << std::endl;
        return;
    }

    // Get both 'QA' and 'PhotonAnalysis' directories
    TDirectory* qaDir = (TDirectory*)inputFile->Get("QA");
    TDirectory* photonDir = (TDirectory*)inputFile->Get("PhotonAnalysis");

    if (!qaDir && !photonDir) {
        std::cerr << "Error: Neither 'QA' nor 'PhotonAnalysis' directories found in file " << inputFilePath << std::endl;
        inputFile->Close();
        return;
    }

    // Create the base output directory 'IsolationEnergies'
    std::string isolationEnergiesDir = outputDir + "/IsolationEnergies";
    gSystem->mkdir(isolationEnergiesDir.c_str(), true);

    // Process both directories if they exist
    if (qaDir) {
        processDirectory(qaDir, isolationEnergiesDir);
    }
    if (photonDir) {
        processDirectory(photonDir, isolationEnergiesDir);
    }

    inputFile->Close();
    delete inputFile;
}
// Function to format floating point values to three significant figures
std::string reFormatToSave(float value) {
    std::ostringstream oss;
    oss.precision(3);
    oss << std::fixed << value;
    return oss.str();
}
// Function to read CSV data into separate maps for inMassWindow and outsideMassWindow
void readDataFromCSV(const std::string& filename,
                     std::map<std::string, std::vector<std::tuple<float, float, float, float, float, float, float>>>& dataMap_inMassWindow,
                     std::map<std::string, std::vector<std::tuple<float, float, float, float, float, float, float>>>& dataMap_outsideMassWindow) {

    std::ifstream file(filename);
    std::string line;
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Skip the header line
    std::getline(file, line);
    std::cout << "Skipping header: " << line << std::endl;

    // Read CSV data
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token, massWindowLabel;
        int trigger;
        float eCore, chi, asym, ptMin, ptMax, isoMin, isoMax, ratio = 0, error = 0, weightedPt;

        // Parse values from CSV
        std::getline(ss, token, ',');
        trigger = std::stoi(token);

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

        // Skip the "Isolated Counts" and "Total Counts" columns
        std::getline(ss, token, ',');
        std::getline(ss, token, ',');

        // Parse the ratio (Isolated/Total)
        std::getline(ss, token, ',');
        ratio = std::stof(token);

        // Parse the statistical error
        std::getline(ss, token, ',');
        error = std::stof(token);

        // Parse the weighted pT value
        std::getline(ss, token, ',');
        weightedPt = std::stof(token);

        // Parse the MassWindowLabel
        std::getline(ss, massWindowLabel, ',');

        // Key for grouping: Trigger + Cut combination (ECore, Chi, Asymmetry)
        std::string key = "Trigger" + std::to_string(trigger) + "_E_" + reFormatToSave(eCore) + "_Chi" + reFormatToSave(chi) +
                          "_Asym" + reFormatToSave(asym);

        // Add data to the appropriate map based on MassWindowLabel
        if (massWindowLabel == "inMassWindow") {
            dataMap_inMassWindow[key].emplace_back(ptMin, ptMax, isoMin, isoMax, ratio, error, weightedPt);
        } else if (massWindowLabel == "outsideMassWindow") {
            dataMap_outsideMassWindow[key].emplace_back(ptMin, ptMax, isoMin, isoMax, ratio, error, weightedPt);
        }
    }
    file.close();
}


// Function to plot Ratio vs pT with fixed isoEt ranges, with an option to exclude certain ranges
void plotRatioVsPt(const std::map<std::string, std::vector<std::tuple<float, float, float, float, float, float, float>>>& dataMap,
                   const std::string& outputDir, const std::string& label,
                   const std::vector<std::pair<float, float>>& exclusionRanges = {},
                   bool drawRefA = false, bool drawRefB = false) {
    // Define isoEt ranges and assign colors
    std::vector<std::pair<float, float>> isoEtRanges = {
        {-5, 0},
        {0, 2},
        {2, 5},
        {5, 10},
        {-10, 0},
        {0, 10}
    };

//    std::vector<std::pair<float, float>> isoEtRanges = {
//        {-5, 0},
//        {0, 5},
//        {-2, 0},
//        {0, 2},
//        {-2, -5},
//        {2, 5},
//        {-5, -10},
//        {5, 10},
//        {-10, 0},
//        {0, 10}
//    };
    
    
    std::vector<int> colors = {kRed + 4, kGreen, kBlue, kMagenta, kRed, kOrange};

    const std::string plotTitle = "Ratio of Isolated Photons from Meson Decays to All Clusters from Meson Decays Compared to PHENIX Data";

    // Create a single dummy legend
    TLegend* dummyLegend = new TLegend(0.18, 0.84, 0.38, 0.89); // Enlarged legend size
    dummyLegend->SetBorderSize(0);  // Remove legend border
    dummyLegend->SetTextSize(0.03); // Set text size for the legend
     
    for (size_t i = 0; i < isoEtRanges.size(); ++i) {
        std::pair<float, float> currentRange = isoEtRanges[i];
        if (std::find(exclusionRanges.begin(), exclusionRanges.end(), currentRange) != exclusionRanges.end()) {
            continue;
        }

        // Format values to two decimal places
        std::ostringstream formattedRange;
        formattedRange << "#font[62]{sPHENIX, Isolation Criteria:}"
                       << std::fixed << std::setprecision(2)
                       << isoEtRanges[i].first << " #leq E_{T, iso} < " << isoEtRanges[i].second << " GeV, #Delta R_{cone} < 0.3";


        std::string legendEntry = formattedRange.str();
        TGraph* dummyGraph = new TGraph();
        dummyGraph->SetMarkerStyle(20);
        dummyGraph->SetMarkerColor(colors[i]);
        dummyLegend->AddEntry(dummyGraph, legendEntry.c_str(), "p");
    }

    
    std::vector<double> referencePTGamma = {3.36, 4.39, 5.41, 6.42, 7.43, 8.44, 9.80, 11.83, 14.48};
    std::vector<double> referenceRatio = {0.594, 0.664, 0.626, 0.658, 0.900, 0.715, 0.872, 0.907, 0.802};
    std::vector<double> referenceStatError = {0.014, 0.028, 0.043, 0.061, 0.113, 0.130, 0.120, 0.190, 0.290};

    // Define second reference data (from the second table - Reference Two)
    std::vector<double> referenceTwoPTGamma = {3.34, 4.38, 5.40, 6.41, 7.42, 8.43, 9.78, 11.81, 14.41};
    std::vector<double> referenceTwoRatio = {0.477, 0.455, 0.448, 0.430, 0.338, 0.351, 0.400, 0.286, 0.371};
    std::vector<double> referenceTwoStatError = {0.0020, 0.0060, 0.012, 0.021, 0.032, 0.053, 0.070, 0.130, 0.180};

    TLegend* refLegend = new TLegend(0.18, 0.77, 0.38, 0.81);  // Adjust position for reference points
    refLegend->SetBorderSize(0);  // Remove legend border
    refLegend->SetTextSize(0.03); // Set text size for the legend

    // Add dummy markers and labels for reference data
    if (drawRefA) {
        TGraph* refGraphOneDummy = new TGraph();  // Dummy graph for Reference 1
        refGraphOneDummy->SetMarkerStyle(22);
        refGraphOneDummy->SetMarkerColor(kBlack);
        refLegend->AddEntry(refGraphOneDummy, "2003 pp Dataset, PHENIX Measurement", "p");
    }

    if (drawRefB) {
        TGraph* refGraphTwoDummy = new TGraph();  // Dummy graph for Reference 2
        refGraphTwoDummy->SetMarkerStyle(20);
        refGraphTwoDummy->SetMarkerColor(kBlue);
        refLegend->AddEntry(refGraphTwoDummy, "#font[62]{PHENIX (2003 pp RHIC Run):} Ratio = #frac{Isolated Photons from #pi^{0} Decays}{All Photons from #pi^{0} Decays}", "p");
    }

    // Iterate over the map and create a plot for each unique combination
    for (const auto& entry : dataMap) {
        const std::string& cutKey = entry.first;
        const auto& points = entry.second;

        std::cout << "Creating plot for cut combination: " << cutKey << std::endl;

        // Extract trigger and cut information from cutKey
        size_t pos = cutKey.find("_E_");
        std::string triggerStr = cutKey.substr(0, pos);
        int triggerIndex = std::stoi(triggerStr.substr(triggerStr.find("Trigger") + 7));
        std::string cutFolder = cutKey.substr(pos + 1);

        // Format trigger name and cut values
        std::string triggerName = getTriggerName(triggerIndex);
        std::string eCore = cutFolder.substr(cutFolder.find("E_") + 2, 5);
        std::string chi = cutFolder.substr(cutFolder.find("Chi") + 3, 5);
        std::string asymmetry = cutFolder.substr(cutFolder.find("Asym") + 4, 5);
        
        
        std::string savePath = "/Users/patsfan753/Desktop/DirectPhotonAna/Plots/IsolationEnergies/" + triggerStr + "/" + cutFolder + "_" + label;

        std::cout << "Saving to folder: " << savePath << std::endl;

        // Create the canvas and graphs
        TCanvas* canvas = new TCanvas(("canvas_" + label).c_str(), ("Ratio vs pT " + label).c_str(), 800, 600);
        TMultiGraph* multiGraph = new TMultiGraph();
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
 
        for (size_t i = 0; i < isoEtRanges.size(); ++i) {
            // Create a new TGraphErrors for each isoEt range
            TGraphErrors* graphWithError = new TGraphErrors();
            int color = colors[i];
            std::pair<float, float> currentRange = isoEtRanges[i];

            // Add points to the graph for the current isoEt range
            for (const auto& point : points) {
                float ptMin = std::get<0>(point);
                float ptMax = std::get<1>(point);
                float isoMin = std::get<2>(point);
                float isoMax = std::get<3>(point);
                float ratio = std::get<4>(point);
                float error = std::get<5>(point);
                float weightedPt = std::get<6>(point);

                // Check if the current isoEt range is in the exclusion list
                std::pair<float, float> isoRange = {isoMin, isoMax};
                if (std::find(exclusionRanges.begin(), exclusionRanges.end(), isoRange) != exclusionRanges.end()) {
                    std::cout << "Excluding isoEt range: isoMin = " << isoMin << ", isoMax = " << isoMax << std::endl;
                    continue;  // Skip the points belonging to the excluded ranges
                }

                // Only add points that belong to the current isoEt range
                if (isoMin == currentRange.first && isoMax == currentRange.second) {
                    int index = graphWithError->GetN();
                    graphWithError->SetPoint(index, weightedPt, ratio);
                    graphWithError->SetPointError(index, 0, error);
                    graphWithError->SetMarkerStyle(20);
                    graphWithError->SetMarkerColor(color); // Set the correct color for this range
                    graphWithError->SetLineColor(color);
                }
            }

            // Add the graph for the current isoEt range to the multiGraph
            multiGraph->Add(graphWithError);
        }

        //multiGraph->SetTitle(plotTitle.c_str());
        multiGraph->GetXaxis()->SetRangeUser(0, 15);  // X-axis range
        multiGraph->GetXaxis()->SetTitleSize(0.045);
        multiGraph->GetXaxis()->SetTitle("Weighted Average p_{T} of all Photons from #pi^{0}/#eta Decay");
        multiGraph->GetYaxis()->SetRangeUser(0, 1.5);  // Y-axis range
        multiGraph->GetYaxis()->SetTitleSize(0.045);
        multiGraph->GetYaxis()->SetTitle("#frac{Isolated Photons from #pi^{0}/#eta Decay}{All Photons from #pi^{0}/#eta Decay}");
        multiGraph->Draw("AP");
        
        // Draw a dashed line at y = 1
        TLine* line = new TLine(2, 1, 15, 1);
        line->SetLineStyle(2); // Dashed line
        line->Draw();
        
        if (drawRefA) {
            // Create reference TGraphErrors for the first set of reference points
            TGraphErrors* refGraphOne = new TGraphErrors(referencePTGamma.size());
            for (size_t i = 0; i < referencePTGamma.size(); ++i) {
                refGraphOne->SetPoint(i, referencePTGamma[i], referenceRatio[i]);
                refGraphOne->SetPointError(i, 0, referenceStatError[i]);
            }
            refGraphOne->SetMarkerStyle(22);
            refGraphOne->SetMarkerColor(kBlack);
            refGraphOne->SetLineColor(kBlack);
            refGraphOne->Draw("P SAME");
        }

        if (drawRefB) {
            // Create reference TGraphErrors for the second set of reference points
            TGraphErrors* refGraphTwo = new TGraphErrors(referenceTwoPTGamma.size());
            for (size_t i = 0; i < referenceTwoPTGamma.size(); ++i) {
                refGraphTwo->SetPoint(i, referenceTwoPTGamma[i], referenceTwoRatio[i]);
                refGraphTwo->SetPointError(i, 0, referenceTwoStatError[i]);
            }
            refGraphTwo->SetMarkerStyle(20);
            refGraphTwo->SetMarkerColor(kBlue);
            refGraphTwo->SetLineColor(kBlue);
            refGraphTwo->Draw("P SAME");
        }

        // Draw reference legend if either refA or refB is drawn
        if (drawRefA || drawRefB) {
            refLegend->Draw();
        }

        // Add trigger and cut information on the top right of the plot
        TLatex labelText;
        labelText.SetNDC();
        labelText.SetTextSize(0.03);

        TLatex valueText;
        valueText.SetNDC();
        valueText.SetTextSize(0.03);

        labelText.DrawLatex(0.21, 0.65, "#font[62]{Trigger:}");
        valueText.DrawLatex(0.31, 0.65, triggerName.c_str());

        labelText.DrawLatex(0.21, 0.6, "#font[62]{ECore #geq}");
        std::string eCoreWithUnit = eCore + " GeV";
        valueText.DrawLatex(0.31, 0.6, eCoreWithUnit.c_str());

        labelText.DrawLatex(0.21, 0.55, "#font[62]{#chi^{2} <}");
        valueText.DrawLatex(0.26, 0.55, chi.c_str());

        labelText.DrawLatex(0.21, 0.5, "#font[62]{Asymmetry <}");
        valueText.DrawLatex(0.34, 0.5, asymmetry.c_str());
        
        // Add the dummy legend and latex for labels
        dummyLegend->Draw();  // Paste the fixed dummy legend

        // Save the plot to the appropriate directory
        gSystem->mkdir(savePath.c_str(), true);
        canvas->SaveAs((savePath + "/RatioVsPt_" + label + ".png").c_str());

        // Cleanup
        delete multiGraph;
        delete canvas;
    }

    // Cleanup the dummy legend
    delete dummyLegend;
}

// Function to overlay Ratio vs pT plots for both inMassWindow and outsideMassWindow with exclusion range
void overlayRatioInOutMassWindow(const std::map<std::string, std::vector<std::tuple<float, float, float, float, float, float, float>>>& dataMap_inMassWindow,
                                 const std::map<std::string, std::vector<std::tuple<float, float, float, float, float, float, float>>>& dataMap_outsideMassWindow,
                                 const std::string& outputDir,
                                 const std::vector<std::pair<float, float>>& exclusionRanges = {},
                                 bool drawRefA = false, bool drawRefB = false) {
    // Define isoEt ranges and assign colors
    std::vector<std::pair<float, float>> isoEtRanges = {
        {-5, 0},
        {0, 2},
        {2, 5},
        {5, 10},
        {-10, 0},
        {0, 10}
    };

    // Colors for inMassWindow and outsideMassWindow datasets
    std::vector<int> colors_inMassWindow = {kRed + 4, kGreen, kBlue, kMagenta, kRed, kOrange};
    std::vector<int> colors_outsideMassWindow = {kRed - 3, kGreen - 6, kBlue - 7, kMagenta - 6, kRed - 2, kOrange - 5};

    // Create a single legend for isoEt ranges, adding labels for inMassWindow and outsideMassWindow
    TLegend* legend = new TLegend(0.16, 0.72, 0.45, 0.88);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.03);

    for (size_t i = 0; i < isoEtRanges.size(); ++i) {
        std::pair<float, float> currentRange = isoEtRanges[i];

        // Skip adding the legend entry if the range is excluded
        if (std::find(exclusionRanges.begin(), exclusionRanges.end(), currentRange) != exclusionRanges.end()) {
            continue;
        }

        std::string legendEntry_inMassWindow = "Clusters in Meson Mass Window (#pi^{0}, #eta): " + std::to_string(isoEtRanges[i].first) + " #leq E_{T, iso} < " + std::to_string(isoEtRanges[i].second);
        TGraph* dummyGraph_inMassWindow = new TGraph();
        dummyGraph_inMassWindow->SetMarkerStyle(20);
        dummyGraph_inMassWindow->SetMarkerColor(colors_inMassWindow[i]);
        legend->AddEntry(dummyGraph_inMassWindow, legendEntry_inMassWindow.c_str(), "p");

        std::string legendEntry_outsideMassWindow = "Clusters outside Meson Mass Window (#pi^{0}, #eta): " + std::to_string(isoEtRanges[i].first) + " #leq E_{T, iso} < " + std::to_string(isoEtRanges[i].second);
        TGraph* dummyGraph_outsideMassWindow = new TGraph();
        dummyGraph_outsideMassWindow->SetMarkerStyle(24);
        dummyGraph_outsideMassWindow->SetMarkerColor(colors_outsideMassWindow[i]);
        legend->AddEntry(dummyGraph_outsideMassWindow, legendEntry_outsideMassWindow.c_str(), "p");
    }

    // Reference data setup
    std::vector<double> referencePTGamma = {3.36, 4.39, 5.41, 6.42, 7.43, 8.44, 9.80, 11.83, 14.48};
    std::vector<double> referenceRatio = {0.594, 0.664, 0.626, 0.658, 0.900, 0.715, 0.872, 0.907, 0.802};
    std::vector<double> referenceStatError = {0.014, 0.028, 0.043, 0.061, 0.113, 0.130, 0.120, 0.190, 0.290};

    std::vector<double> referenceTwoPTGamma = {3.34, 4.38, 5.40, 6.41, 7.42, 8.43, 9.78, 11.81, 14.41};
    std::vector<double> referenceTwoRatio = {0.477, 0.455, 0.448, 0.430, 0.338, 0.351, 0.400, 0.286, 0.371};
    std::vector<double> referenceTwoStatError = {0.0020, 0.0060, 0.012, 0.021, 0.032, 0.053, 0.070, 0.130, 0.180};

    TLegend* refLegend = new TLegend(0.14, 0.69, 0.45, 0.7);
    refLegend->SetBorderSize(0);
    refLegend->SetTextSize(0.027);

    if (drawRefA) {
        TGraph* refGraphOneDummy = new TGraph();
        refGraphOneDummy->SetMarkerStyle(22);
        refGraphOneDummy->SetMarkerColor(kBlack);
        refLegend->AddEntry(refGraphOneDummy, "2003 pp Dataset, PHENIX Measurement", "p");
    }

    if (drawRefB) {
        TGraph* refGraphTwoDummy = new TGraph();
        refGraphTwoDummy->SetMarkerStyle(21);
        refGraphTwoDummy->SetMarkerColor(kGray + 2);
        refLegend->AddEntry(refGraphTwoDummy, "2003 pp Dataset, PHENIX Measurement", "p");
    }

    for (const auto& entry : dataMap_inMassWindow) {
        const std::string& cutKey = entry.first;
        const auto& points_inMassWindow = entry.second;

        // Check if the corresponding key exists in the outsideMassWindow data
        auto outsideIter = dataMap_outsideMassWindow.find(cutKey);
        if (outsideIter == dataMap_outsideMassWindow.end()) continue;
        const auto& points_outsideMassWindow = outsideIter->second;
        
        // Extract trigger and cut details from cutKey
        size_t pos = cutKey.find("_E_");
        std::string triggerStr = cutKey.substr(0, pos);
        std::string cutFolder = cutKey.substr(pos + 1);


        std::string savePath = outputDir + "IsolationEnergies/" + triggerStr + "/" + cutFolder;
        gSystem->mkdir(savePath.c_str(), true);


        TCanvas* canvas = new TCanvas(("canvas_" + cutKey).c_str(), ("Overlay Ratio vs pT for " + cutKey).c_str(), 800, 600);
        TMultiGraph* multiGraph = new TMultiGraph();

        // Plot points for both inMassWindow and outsideMassWindow
        for (size_t i = 0; i < isoEtRanges.size(); ++i) {
            int color_inMassWindow = colors_inMassWindow[i];
            int color_outsideMassWindow = colors_outsideMassWindow[i];
            std::pair<float, float> currentRange = isoEtRanges[i];

            // Graphs for inMassWindow and outsideMassWindow
            TGraphErrors* graph_inMassWindow = new TGraphErrors();
            TGraphErrors* graph_outsideMassWindow = new TGraphErrors();

            for (const auto& point : points_inMassWindow) {
                float isoMin = std::get<2>(point);
                float isoMax = std::get<3>(point);
                float ratio = std::get<4>(point);
                float error = std::get<5>(point);
                float weightedPt = std::get<6>(point);

                if (std::find(exclusionRanges.begin(), exclusionRanges.end(), std::make_pair(isoMin, isoMax)) != exclusionRanges.end()) {
                    continue;
                }

                if (isoMin == currentRange.first && isoMax == currentRange.second) {
                    int index = graph_inMassWindow->GetN();
                    graph_inMassWindow->SetPoint(index, weightedPt, ratio);
                    graph_inMassWindow->SetPointError(index, 0, error);
                    graph_inMassWindow->SetMarkerStyle(20);
                    graph_inMassWindow->SetMarkerColor(color_inMassWindow);
                }
            }

            for (const auto& point : points_outsideMassWindow) {
                float isoMin = std::get<2>(point);
                float isoMax = std::get<3>(point);
                float ratio = std::get<4>(point);
                float error = std::get<5>(point);
                float weightedPt = std::get<6>(point);

                if (std::find(exclusionRanges.begin(), exclusionRanges.end(), std::make_pair(isoMin, isoMax)) != exclusionRanges.end()) {
                    continue;
                }

                if (isoMin == currentRange.first && isoMax == currentRange.second) {
                    int index = graph_outsideMassWindow->GetN();
                    graph_outsideMassWindow->SetPoint(index, weightedPt, ratio);
                    graph_outsideMassWindow->SetPointError(index, 0, error);
                    graph_outsideMassWindow->SetMarkerStyle(24);
                    graph_outsideMassWindow->SetMarkerColor(color_outsideMassWindow);
                }
            }

            multiGraph->Add(graph_inMassWindow);
            multiGraph->Add(graph_outsideMassWindow);
        }

        // Draw the multiGraph
        multiGraph->SetTitle(("Overlay Ratio vs pT for " + cutKey).c_str());
        multiGraph->GetXaxis()->SetRangeUser(0, 15);
        multiGraph->GetXaxis()->SetTitle("Weighted p_{T} of all Photons from Meson Decay");
        multiGraph->GetYaxis()->SetRangeUser(0, 1.5);
        multiGraph->GetYaxis()->SetTitle("Ratio (Isolated/Total)");
        multiGraph->Draw("AP");

        // Draw a dashed line at y = 1
        TLine* line = new TLine(2, 1, 15, 1);
        line->SetLineStyle(2);
        line->Draw();

        if (drawRefA) {
            TGraphErrors* refGraphOne = new TGraphErrors(referencePTGamma.size());
            for (size_t i = 0; i < referencePTGamma.size(); ++i) {
                refGraphOne->SetPoint(i, referencePTGamma[i], referenceRatio[i]);
                refGraphOne->SetPointError(i, 0, referenceStatError[i]);
            }
            refGraphOne->SetMarkerStyle(22);
            refGraphOne->SetMarkerColor(kBlack);
            refGraphOne->SetLineColor(kBlack);
            refGraphOne->Draw("P SAME");
        }

        if (drawRefB) {
            TGraphErrors* refGraphTwo = new TGraphErrors(referenceTwoPTGamma.size());
            for (size_t i = 0; i < referenceTwoPTGamma.size(); ++i) {
                refGraphTwo->SetPoint(i, referenceTwoPTGamma[i], referenceTwoRatio[i]);
                refGraphTwo->SetPointError(i, 0, referenceTwoStatError[i]);
            }
            refGraphTwo->SetMarkerStyle(21);
            refGraphTwo->SetMarkerColor(kGray + 2);
            refGraphTwo->SetLineColor(kGray + 2);
            refGraphTwo->Draw("P SAME");
        }

        if (drawRefA || drawRefB) {
            refLegend->Draw();
        }

        legend->Draw();

        gSystem->mkdir(savePath.c_str(), true);
        canvas->SaveAs((savePath + "/Overlay_RatioVsPt_InsideOutsideMassWindow.png").c_str());

        delete multiGraph;
        delete canvas;
    }

    delete legend;
}



void analyzeClusterIso() {
    gROOT->LoadMacro("sPhenixStyle.C");
    SetsPhenixStyle();
    std::string csvFilePath = "/Users/patsfan753/Desktop/DirectPhotonInformation.csv";
//    isolationEnergies(inputFilePath, outputDir);
//

//    outputHistogramLogs(csvFilePath);
    
    std::vector<std::pair<float, float>> exclusionRanges = {{-5, 0}, {0, 2}, {2, 5}, {5, 10}, {0, 10}};
    // Maps to store data for inMassWindow and outsideMassWindow
    std::map<std::string, std::vector<std::tuple<float, float, float, float, float, float, float>>> dataMap_inMassWindow;
    std::map<std::string, std::vector<std::tuple<float, float, float, float, float, float, float>>> dataMap_outsideMassWindow;

    // Read the data from CSV into maps
    readDataFromCSV(csvFilePath, dataMap_inMassWindow, dataMap_outsideMassWindow);

    // Plot for inMassWindow data
    plotRatioVsPt(dataMap_inMassWindow, outputDir, "inMassWindow", exclusionRanges, false, true);

    // Plot for outsideMassWindow data
    plotRatioVsPt(dataMap_outsideMassWindow, outputDir, "outsideMassWindow", exclusionRanges, false, true);
    
    overlayRatioInOutMassWindow(dataMap_inMassWindow, dataMap_outsideMassWindow, outputDir, exclusionRanges, false, false);
}

