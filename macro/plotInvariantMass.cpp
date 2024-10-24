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

struct HistogramData {
    CutValues cuts;  // Holds triggerIndex, clusECore, chi, asymmetry, pTMin, pTMax
    std::string histName;  // Name of the histogram

    // Fitted parameters and their errors
    double meanPi0;
    double meanPi0Error;
    double sigmaPi0;
    double sigmaPi0Error;
    double meanEta;
    double meanEtaError;
    double sigmaEta;
    double sigmaEtaError;

    // Mass ratio and its error
    double massRatio;
    double massRatioError;

    // Signal and background yields for pi0
    double signalPi0Yield;
    double signalPi0Error;
    double backgroundPi0Yield;
    double backgroundPi0Error;
    double signalToBackgroundPi0Ratio;
    double signalToBackgroundPi0Error;

    // Signal and background yields for eta
    double signalEtaYield;
    double signalEtaError;
    double backgroundEtaYield;
    double backgroundEtaError;
    double signalToBackgroundEtaRatio;
    double signalToBackgroundEtaError;
};

std::vector<HistogramData> saveAnnotatedInvariantMassHistograms(const std::string& inputFilePath) {
    std::vector<HistogramData> histogramDataVector;
    
    TFile* inputFile = TFile::Open(inputFilePath.c_str(), "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open the file " << inputFilePath << std::endl;
        return histogramDataVector;
    }

    // Get the 'PhotonAnalysis' directory
    TDirectory* invMassDir = (TDirectory*)inputFile->Get("PhotonAnalysis");
    if (!invMassDir) {
        std::cerr << "Error: 'PhotonAnalysis' directory not found." << std::endl;
        inputFile->Close();
        return histogramDataVector;
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
                    double meanPi0error = totalFit->GetParError(1);
                    double sigmaPi0 = totalFit->GetParameter(2);
                    double sigmaPi0error = totalFit->GetParError(2);
                    double meanEta = totalFit->GetParameter(4);
                    double meanEtaError = totalFit->GetParError(4);
                    double sigmaEta = totalFit->GetParameter(5);
                    double sigmaEtaError = totalFit->GetParError(5);

                    double massRatio = meanEta / meanPi0;
                    
                    // Calculate the relative errors
                    double relativeErrorMeanEta = meanEtaError / meanEta;
                    double relativeErrorMeanPi0 = meanPi0error / meanPi0;

                    // Calculate the relative error of massRatio
                    double relativeErrorMassRatio = sqrt(
                        pow(relativeErrorMeanEta, 2) +
                        pow(relativeErrorMeanPi0, 2)
                    );

                    // Compute the absolute error of massRatio
                    double massRatioError = massRatio * relativeErrorMassRatio;
                    
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
                    
                    // Store the data in the HistogramData structure
                    HistogramData data;
                    data.cuts = cuts;
                    data.histName = histName;
                    data.meanPi0 = meanPi0;
                    data.meanPi0Error = meanPi0error;
                    data.sigmaPi0 = sigmaPi0;
                    data.sigmaPi0Error = sigmaPi0error;
                    data.meanEta = meanEta;
                    data.meanEtaError = meanEtaError;
                    data.sigmaEta = sigmaEta;
                    data.sigmaEtaError = sigmaEtaError;
                    data.massRatio = massRatio;
                    data.massRatioError = massRatioError;
                    data.signalPi0Yield = signalPi0Yield;
                    data.signalPi0Error = signalPi0Error;
                    data.backgroundPi0Yield = backgroundPi0Yield;
                    data.backgroundPi0Error = backgroundPi0Error;
                    data.signalToBackgroundPi0Ratio = signalToBackgroundPi0Ratio;
                    data.signalToBackgroundPi0Error = signalToBackgroundPi0Error;
                    data.signalEtaYield = signalEtaYield;
                    data.signalEtaError = signalEtaError;
                    data.backgroundEtaYield = backgroundEtaYield;
                    data.backgroundEtaError = backgroundEtaError;
                    data.signalToBackgroundEtaRatio = signalToBackgroundEtaRatio;
                    data.signalToBackgroundEtaError = signalToBackgroundEtaError;

                    // Add the data to the vector
                    histogramDataVector.push_back(data);
                    
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
    
    return histogramDataVector; // Return the collected data
}

void printHistogramData(const std::vector<HistogramData>& histogramDataVector) {
    // Organize the data by trigger index, then cut combinations, then pT bins
    // For that, best use nested maps
    // map<TriggerIndex, map<CutCombinationString, vector<HistogramData>>>
    std::map<int, std::map<std::string, std::vector<HistogramData>>> dataMap;
    
    // First, organize the data
    for (const auto& data : histogramDataVector) {
        int triggerIndex = data.cuts.triggerIndex;
        // Create a string representing the cut combination
        std::ostringstream cutCombinationStream;
        cutCombinationStream << "E" << data.cuts.clusECore
                             << "_Chi" << data.cuts.chi
                             << "_Asym" << data.cuts.asymmetry;
        std::string cutCombination = cutCombinationStream.str();
        
        dataMap[triggerIndex][cutCombination].push_back(data);
    }
    
    // Now, print the data
    for (const auto& triggerPair : dataMap) {
        int triggerIndex = triggerPair.first;
        std::string triggerName = getTriggerName(triggerIndex);
        std::cout << "Trigger Index: " << triggerIndex << " (" << triggerName << ")\n";
        
        for (const auto& cutPair : triggerPair.second) {
            std::string cutCombination = cutPair.first;
            std::cout << "  Cut Combination: " << cutCombination << "\n";
            
            // Now print the data for each pT bin
            // We'll create a table header
            std::cout << "    pTMin - pTMax | meanPi0 ± error | sigmaPi0 ± error | meanEta ± error | sigmaEta ± error | massRatio ± error\n";
            std::cout << "    ----------------------------------------------------------------------------------------------------------\n";
            for (const auto& data : cutPair.second) {
                // If pTMin is -1, it means no pT bin, so we can indicate that
                std::string pTBin;
                if (data.cuts.pTMin == -1) {
                    pTBin = "No pT bin";
                } else {
                    pTBin = formatToThreeSigFigs(data.cuts.pTMin) + " - " + formatToThreeSigFigs(data.cuts.pTMax);
                }
                std::cout << "    " << std::setw(14) << pTBin << " | "
                          << formatToThreeSigFigs(data.meanPi0) << " ± " << formatToThreeSigFigs(data.meanPi0Error) << " | "
                          << formatToThreeSigFigs(data.sigmaPi0) << " ± " << formatToThreeSigFigs(data.sigmaPi0Error) << " | "
                          << formatToThreeSigFigs(data.meanEta) << " ± " << formatToThreeSigFigs(data.meanEtaError) << " | "
                          << formatToThreeSigFigs(data.sigmaEta) << " ± " << formatToThreeSigFigs(data.sigmaEtaError) << " | "
                          << formatToThreeSigFigs(data.massRatio) << " ± " << formatToThreeSigFigs(data.massRatioError) << "\n";
            }
            std::cout << "\n"; // Add an empty line between cut combinations
        }
        std::cout << "\n"; // Add an empty line between triggers
    }
}

void plotMesonMeanVsPt(const std::vector<HistogramData>& histogramDataVector) {
    // Organize data by trigger index and cut combination
    std::map<int, std::map<std::string, std::vector<HistogramData>>> dataMap;

    std::map<int, std::vector<std::string>> cutVariations;
    std::map<int, std::vector<double>> meanPi0sPerTrigger;
    std::map<int, std::vector<double>> meanPi0ErrorsPerTrigger;
    std::map<int, std::vector<double>> meanEtasPerTrigger;
    std::map<int, std::vector<double>> meanEtaErrorsPerTrigger;
    
    
    for (const auto& data : histogramDataVector) {
        int triggerIndex = data.cuts.triggerIndex;
        // Create a string representing the cut combination with formatted values
        std::ostringstream cutCombinationStream;
        cutCombinationStream << "E" << formatToThreeSigFigs(data.cuts.clusECore)
                             << "_Chi" << formatToThreeSigFigs(data.cuts.chi)
                             << "_Asym" << formatToThreeSigFigs(data.cuts.asymmetry);
        std::string cutCombination = cutCombinationStream.str();

        dataMap[triggerIndex][cutCombination].push_back(data);
        
        if (data.cuts.pTMin == -1) {
            // Initialize per-trigger data if not already done
            if (std::find(cutVariations[triggerIndex].begin(), cutVariations[triggerIndex].end(), cutCombination) == cutVariations[triggerIndex].end()) {
                cutVariations[triggerIndex].push_back(cutCombination);

                // Initialize with NaN or invalid value
                meanPi0sPerTrigger[triggerIndex].push_back(std::numeric_limits<double>::quiet_NaN());
                meanPi0ErrorsPerTrigger[triggerIndex].push_back(0);
                meanEtasPerTrigger[triggerIndex].push_back(std::numeric_limits<double>::quiet_NaN());
                meanEtaErrorsPerTrigger[triggerIndex].push_back(0);
            }
            size_t index = cutVariations[triggerIndex].size() - 1;

            // Check for invalid fits
            bool invalidPi0Fit = (data.meanPi0Error > 1.0) || (data.meanPi0Error == 0);
            bool invalidEtaFit = (data.meanEtaError > 1.0) || (data.meanEtaError == 0);

            if (!invalidPi0Fit) {
                meanPi0sPerTrigger[triggerIndex][index] = data.meanPi0;
                meanPi0ErrorsPerTrigger[triggerIndex][index] = data.meanPi0Error;
            }

            if (!invalidEtaFit) {
                meanEtasPerTrigger[triggerIndex][index] = data.meanEta;
                meanEtaErrorsPerTrigger[triggerIndex][index] = data.meanEtaError;
            }
        }
    }

    // Iterate over triggers and cut combinations
    for (const auto& triggerPair : dataMap) {
        int triggerIndex = triggerPair.first;
        std::string triggerName = getTriggerName(triggerIndex);

        for (const auto& cutPair : triggerPair.second) {
            std::string cutCombination = cutPair.first;
            const auto& dataList = cutPair.second;

            // Prepare vectors for plotting π⁰
            std::vector<double> pTCentersPi0;
            std::vector<double> meanPi0s;
            std::vector<double> meanPi0Errors;

            // Prepare vectors for plotting η
            std::vector<double> pTCentersEta;
            std::vector<double> meanEtas;
            std::vector<double> meanEtaErrors;

            // Handle pT bins and mean values
            for (const auto& data : dataList) {
                if (data.cuts.pTMin == -1) {
                    continue;  // Skip entries without pT bins
                }

                double pTCenter = (data.cuts.pTMin + data.cuts.pTMax) / 2.0;

                // Criteria for invalid fits (adjust as needed)
                bool invalidPi0Fit = (data.meanPi0Error > 1.0) || (data.meanPi0Error == 0);
                bool invalidEtaFit = (data.meanEtaError > 1.0) || (data.meanEtaError == 0);

                if (!invalidPi0Fit) {
                    pTCentersPi0.push_back(pTCenter);
                    meanPi0s.push_back(data.meanPi0);
                    meanPi0Errors.push_back(data.meanPi0Error);
                }

                if (!invalidEtaFit) {
                    pTCentersEta.push_back(pTCenter);
                    meanEtas.push_back(data.meanEta);
                    meanEtaErrors.push_back(data.meanEtaError);
                }
            }

            // Generate Mean π⁰ Mass vs. pT Plot
            if (!pTCentersPi0.empty()) {
                // Create TGraphErrors for π⁰
                TGraphErrors* graphPi0 = new TGraphErrors(pTCentersPi0.size());
                for (size_t i = 0; i < pTCentersPi0.size(); ++i) {
                    graphPi0->SetPoint(i, pTCentersPi0[i], meanPi0s[i]);
                    graphPi0->SetPointError(i, 0, meanPi0Errors[i]);
                }
                graphPi0->SetMarkerStyle(21);
                graphPi0->SetMarkerSize(1);
                graphPi0->SetLineWidth(2);
                graphPi0->SetMarkerColor(kBlue);
                graphPi0->SetLineColor(kBlue);

                // Determine y-axis range automatically for π⁰, including error bars
                double minY = std::numeric_limits<double>::max();
                double maxY = std::numeric_limits<double>::lowest();
                for (size_t i = 0; i < meanPi0s.size(); ++i) {
                    double yVal = meanPi0s[i];
                    double yErr = meanPi0Errors[i];
                    minY = std::min(minY, yVal - yErr);
                    maxY = std::max(maxY, yVal + yErr);
                }
                double yMargin = 0.05 * (maxY - minY);  // Add 5% margin
                minY -= yMargin;
                maxY += yMargin;

                // Create canvas
                TCanvas canvas;
                canvas.SetGrid();

                // Draw graph for π⁰
                TH1F* hFrame = canvas.DrawFrame(0, minY, 20, maxY, ";p_{T} [GeV];Mean #pi^{0} Mass [GeV]");
                graphPi0->Draw("P SAME");

                // Add labels
                TLatex latex;
                latex.SetNDC();
                latex.SetTextSize(0.04);
                latex.DrawLatex(0.6, 0.85, ("Trigger: " + triggerName).c_str());
                latex.DrawLatex(0.6, 0.80, ("Cut: " + cutCombination).c_str());

                // Save plot to the corresponding folder
                std::string outputDirPath = outputDir + "InvMass/Trigger" + std::to_string(triggerIndex) + "/" + cutCombination + "/";
                std::string outputFilePath = outputDirPath + "MeanPi0_vs_pT.png";

                // Ensure the directory exists
                gSystem->mkdir(outputDirPath.c_str(), true);

                canvas.SaveAs(outputFilePath.c_str());
                std::cout << "Saved plot: " << outputFilePath << std::endl;

                // Clean up
                delete graphPi0;
            }

            // Generate Mean η Mass vs. pT Plot
            if (!pTCentersEta.empty()) {
                // Create TGraphErrors for η
                TGraphErrors* graphEta = new TGraphErrors(pTCentersEta.size());
                for (size_t i = 0; i < pTCentersEta.size(); ++i) {
                    graphEta->SetPoint(i, pTCentersEta[i], meanEtas[i]);
                    graphEta->SetPointError(i, 0, meanEtaErrors[i]);
                }
                graphEta->SetMarkerStyle(20);
                graphEta->SetMarkerSize(1);
                graphEta->SetLineWidth(2);
                graphEta->SetMarkerColor(kRed);
                graphEta->SetLineColor(kRed);

                // Determine y-axis range automatically for η, including error bars
                double minY = std::numeric_limits<double>::max();
                double maxY = std::numeric_limits<double>::lowest();
                for (size_t i = 0; i < meanEtas.size(); ++i) {
                    double yVal = meanEtas[i];
                    double yErr = meanEtaErrors[i];
                    minY = std::min(minY, yVal - yErr);
                    maxY = std::max(maxY, yVal + yErr);
                }
                double yMargin = 0.05 * (maxY - minY);  // Add 5% margin
                minY -= yMargin;
                maxY += yMargin;

                // Create canvas
                TCanvas canvas;
                canvas.SetGrid();

                // Draw graph for η
                TH1F* hFrame = canvas.DrawFrame(0, minY, 20, maxY, ";p_{T} [GeV];Mean #eta Mass [GeV]");
                graphEta->Draw("P SAME");

                // Add labels
                TLatex latex;
                latex.SetNDC();
                latex.SetTextSize(0.04);
                latex.DrawLatex(0.6, 0.85, ("Trigger: " + triggerName).c_str());
                latex.DrawLatex(0.6, 0.80, ("Cut: " + cutCombination).c_str());

                // Save plot to the corresponding folder
                std::string outputDirPath = outputDir + "InvMass/Trigger" + std::to_string(triggerIndex) + "/" + cutCombination + "/";
                std::string outputFilePath = outputDirPath + "MeanEta_vs_pT.png";

                // Ensure the directory exists
                gSystem->mkdir(outputDirPath.c_str(), true);

                canvas.SaveAs(outputFilePath.c_str());
                std::cout << "Saved plot: " << outputFilePath << std::endl;

                // Clean up
                delete graphEta;
            }
        }
    }
    
    for (const auto& triggerPair : dataMap) {
        int triggerIndex = triggerPair.first;
        std::string triggerName = getTriggerName(triggerIndex);

        // Retrieve collected data for this trigger
        const auto& cuts = cutVariations[triggerIndex];
        const auto& pi0Means = meanPi0sPerTrigger[triggerIndex];
        const auto& pi0Errors = meanPi0ErrorsPerTrigger[triggerIndex];
        const auto& etaMeans = meanEtasPerTrigger[triggerIndex];
        const auto& etaErrors = meanEtaErrorsPerTrigger[triggerIndex];

        // Generate Mean π⁰ Mass vs. Cut Variation Plot
        {
            std::vector<double> xValues;
            std::vector<double> xErrors;
            std::vector<double> yValues;
            std::vector<double> yErrors;
            std::vector<std::string> xLabels;

            for (size_t i = 0; i < pi0Means.size(); ++i) {
                if (std::isnan(pi0Means[i])) continue; // Skip invalid fits
                xValues.push_back(i + 1);
                xErrors.push_back(0);
                yValues.push_back(pi0Means[i]);
                yErrors.push_back(pi0Errors[i]);
                xLabels.push_back(cuts[i]);
            }

            if (!xValues.empty()) {
                TGraphErrors* graphPi0 = new TGraphErrors(xValues.size(), &xValues[0], &yValues[0], &xErrors[0], &yErrors[0]);
                graphPi0->SetMarkerStyle(21);
                graphPi0->SetMarkerSize(1);
                graphPi0->SetLineWidth(2);

                // Determine y-axis range including error bars
                double minY = std::numeric_limits<double>::max();
                double maxY = std::numeric_limits<double>::lowest();
                for (size_t i = 0; i < yValues.size(); ++i) {
                    double yVal = yValues[i];
                    double yErr = yErrors[i];
                    minY = std::min(minY, yVal - yErr);
                    maxY = std::max(maxY, yVal + yErr);
                }
                double yMargin = 0.1 * (maxY - minY);
                minY -= yMargin;
                maxY += yMargin;

                // Create canvas
                TCanvas canvas;
                canvas.SetGrid();

                // Draw graph
                double xMin = 0;
                double xMax = xValues.size() + 1;
                TH1F* hFrame = canvas.DrawFrame(xMin, minY, xMax, maxY, ";Cut Variation;Mean #pi^{0} Mass [GeV]");
                graphPi0->Draw("P SAME");

                // Set x-axis labels
                TAxis* axis = hFrame->GetXaxis();
                axis->SetNdivisions(xValues.size(), false);
                axis->SetLabelSize(0.03);
                axis->SetTickLength(0);
                axis->SetLabelOffset(999); // Hide default labels

                // Create custom labels
                for (size_t i = 0; i < xValues.size(); ++i) {
                    double x = xValues[i];
                    TLatex label;
                    label.SetTextAlign(33); // Center align
                    label.SetTextSize(0.02);
                    label.SetTextAngle(90); // Vertical text
                    label.DrawLatex(x, minY - (maxY - minY)*0.02, xLabels[i].c_str());
                }

                // Add labels
                TLatex latex;
                latex.SetNDC();
                latex.SetTextSize(0.04);
                latex.DrawLatex(0.6, 0.85, ("Trigger: " + triggerName).c_str());

                // Save plot
                std::string outputDirPath = outputDir + "InvMass/Trigger" + std::to_string(triggerIndex) + "/";
                std::string outputFilePath = outputDirPath + "MeanPi0_vs_CutVariation.png";
                gSystem->mkdir(outputDirPath.c_str(), true);
                canvas.SaveAs(outputFilePath.c_str());
                std::cout << "Saved plot: " << outputFilePath << std::endl;

                // Clean up
                delete graphPi0;
            }
        }

        // Generate Mean η Mass vs. Cut Variation Plot
        {
            std::vector<double> xValues;
            std::vector<double> xErrors;
            std::vector<double> yValues;
            std::vector<double> yErrors;
            std::vector<std::string> xLabels;

            for (size_t i = 0; i < etaMeans.size(); ++i) {
                if (std::isnan(etaMeans[i])) continue; // Skip invalid fits
                xValues.push_back(i + 1);
                xErrors.push_back(0);
                yValues.push_back(etaMeans[i]);
                yErrors.push_back(etaErrors[i]);
                xLabels.push_back(cuts[i]);
            }

            if (!xValues.empty()) {
                TGraphErrors* graphEta = new TGraphErrors(xValues.size(), &xValues[0], &yValues[0], &xErrors[0], &yErrors[0]);
                graphEta->SetMarkerStyle(21);
                graphEta->SetMarkerSize(1);
                graphEta->SetLineWidth(2);

                // Determine y-axis range including error bars
                double minY = std::numeric_limits<double>::max();
                double maxY = std::numeric_limits<double>::lowest();
                for (size_t i = 0; i < yValues.size(); ++i) {
                    double yVal = yValues[i];
                    double yErr = yErrors[i];
                    minY = std::min(minY, yVal - yErr);
                    maxY = std::max(maxY, yVal + yErr);
                }
                double yMargin = 0.1 * (maxY - minY);
                minY -= yMargin;
                maxY += yMargin;

                // Create canvas
                TCanvas canvas;
                canvas.SetGrid();

                // Draw graph
                double xMin = 0;
                double xMax = xValues.size() + 1;
                TH1F* hFrame = canvas.DrawFrame(xMin, minY, xMax, maxY, ";Cut Variation;Mean #eta Mass [GeV]");
                graphEta->Draw("P SAME");

                // Set x-axis labels
                TAxis* axis = hFrame->GetXaxis();
                axis->SetNdivisions(xValues.size(), false);
                axis->SetLabelSize(0.025);
                axis->SetTickLength(0);
                axis->SetLabelOffset(999); // Hide default labels

                // Create custom labels
                for (size_t i = 0; i < xValues.size(); ++i) {
                    double x = xValues[i];
                    TLatex label;
                    label.SetTextAlign(33); // Center align
                    label.SetTextSize(0.02);
                    label.SetTextAngle(90); // Vertical text
                    label.DrawLatex(x, minY - (maxY - minY)*0.02, xLabels[i].c_str());
                }

                // Add labels
                TLatex latex;
                latex.SetNDC();
                latex.SetTextSize(0.04);
                latex.DrawLatex(0.6, 0.85, ("Trigger: " + triggerName).c_str());

                // Save plot
                std::string outputDirPath = outputDir + "InvMass/Trigger" + std::to_string(triggerIndex) + "/";
                std::string outputFilePath = outputDirPath + "MeanEta_vs_CutVariation.png";
                gSystem->mkdir(outputDirPath.c_str(), true);
                canvas.SaveAs(outputFilePath.c_str());
                std::cout << "Saved plot: " << outputFilePath << std::endl;

                // Clean up
                delete graphEta;
            }
        }
    }
}


// Main function to run both processes
void plotInvariantMass() {
    gROOT->LoadMacro("sPhenixStyle.C");
    SetsPhenixStyle();
    saveAnnotatedInvariantMassHistograms(inputFilePath);
    std::vector<HistogramData> histogramDataVector = saveAnnotatedInvariantMassHistograms(inputFilePath);
    printHistogramData(histogramDataVector);
    plotMesonMeanVsPt(histogramDataVector);
}
