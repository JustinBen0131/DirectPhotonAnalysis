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
    CutValues cuts;
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
    
    // Resolution parameters
    double pi0FitResolution;
    double pi0FitResolutionError;
    double etaFitResolution;
    double etaFitResolutionError;
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
                    
                    // Calculate fit resolution with safeguards
                    double pi0FitResolution = std::numeric_limits<double>::quiet_NaN();
                    double pi0FitResolutionError = std::numeric_limits<double>::quiet_NaN();

                    // Debugging output for pi0 values
                    std::cout << "Calculating pi0FitResolution:" << std::endl;
                    std::cout << "  meanPi0: " << meanPi0 << ", meanPi0Error: " << meanPi0error << std::endl;
                    std::cout << "  sigmaPi0: " << sigmaPi0 << ", sigmaPi0Error: " << sigmaPi0error << std::endl;

                    if (meanPi0 > 1e-6 && sigmaPi0 > 1e-6) {  // Avoid near-zero values
                        pi0FitResolution = sigmaPi0 / meanPi0;
                        pi0FitResolutionError = pi0FitResolution * sqrt(pow(sigmaPi0error / sigmaPi0, 2) + pow(meanPi0error / meanPi0, 2));
                        std::cout << "  Calculated pi0FitResolution: " << pi0FitResolution << std::endl;
                        std::cout << "  Calculated pi0FitResolutionError: " << pi0FitResolutionError << std::endl;
                    } else {
                        std::cout << "  Skipping pi0FitResolution calculation due to small or zero values." << std::endl;
                    }

                    // Calculate eta fit resolution with safeguards
                    double etaFitResolution = std::numeric_limits<double>::quiet_NaN();
                    double etaFitResolutionError = std::numeric_limits<double>::quiet_NaN();

                    // Debugging output for eta values
                    std::cout << "Calculating etaFitResolution:" << std::endl;
                    std::cout << "  meanEta: " << meanEta << ", meanEtaError: " << meanEtaError << std::endl;
                    std::cout << "  sigmaEta: " << sigmaEta << ", sigmaEtaError: " << sigmaEtaError << std::endl;

                    if (meanEta > 1e-6 && sigmaEta > 1e-6) {  // Avoid near-zero values
                        etaFitResolution = sigmaEta / meanEta;
                        etaFitResolutionError = etaFitResolution * sqrt(pow(sigmaEtaError / sigmaEta, 2) + pow(meanEtaError / meanEta, 2));
                        std::cout << "  Calculated etaFitResolution: " << etaFitResolution << std::endl;
                        std::cout << "  Calculated etaFitResolutionError: " << etaFitResolutionError << std::endl;
                    } else {
                        std::cout << "  Skipping etaFitResolution calculation due to small or zero values." << std::endl;
                    }

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
                    data.pi0FitResolution = pi0FitResolution;
                    data.pi0FitResolutionError = pi0FitResolutionError;
                    data.etaFitResolution = etaFitResolution;
                    data.etaFitResolutionError = etaFitResolutionError;

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
    // Open a CSV file for writing to the specified path
    std::ofstream csvFile("/Users/patsfan753/Desktop/InvariantMassInformation_runnumber46623_runnumber47230.csv");
    
    // Check if the file opened successfully
    if (!csvFile.is_open()) {
        std::cerr << "Error: Could not open the file for writing at /Users/patsfan753/Desktop/InvariantMassInformation_runnumber46623_runnumber47230.csv" << std::endl;
        return;
    }

    // Write the CSV headers with separate columns for Ecore, Chi2, and Asym
    csvFile << "Trigger Index,Ecore,Chi2,Asym,pTMin,pTMax,meanPi0,sigmaPi0,meanEta,sigmaEta\n";
    
    // Organize the data by trigger index, then cut combinations, then pT bins
    std::map<int, std::map<std::string, std::vector<HistogramData>>> dataMap;

    // First, organize the data
    for (const auto& data : histogramDataVector) {
        int triggerIndex = data.cuts.triggerIndex;
        std::ostringstream cutCombinationStream;
        cutCombinationStream << "E" << data.cuts.clusECore
                             << "_Chi" << data.cuts.chi
                             << "_Asym" << data.cuts.asymmetry;
        std::string cutCombination = cutCombinationStream.str();
        
        dataMap[triggerIndex][cutCombination].push_back(data);
    }

    // Print the organized data to console and write it to CSV
    for (const auto& triggerPair : dataMap) {
        int triggerIndex = triggerPair.first;
        std::string triggerName = getTriggerName(triggerIndex);
        std::cout << "Trigger Index: " << triggerIndex << " (" << triggerName << ")\n";
        
        for (const auto& cutPair : triggerPair.second) {
            // Decompose the cut combination
            float Ecore = cutPair.second[0].cuts.clusECore;
            float Chi2 = cutPair.second[0].cuts.chi;
            float Asym = cutPair.second[0].cuts.asymmetry;

            std::cout << "  Cut Combination: E" << Ecore << "_Chi" << Chi2 << "_Asym" << Asym << "\n";
            
            // Print table header to console
            std::cout << "    pTMin - pTMax | meanPi0 ± error | sigmaPi0 ± error | meanEta ± error | sigmaEta ± error\n";
            std::cout << "    ----------------------------------------------------------------------------------------------------------\n";
            
            for (const auto& data : cutPair.second) {
                std::string pTBin = (data.cuts.pTMin == -1) ? "No pT bin" : formatToThreeSigFigs(data.cuts.pTMin) + " - " + formatToThreeSigFigs(data.cuts.pTMax);
                
                // Print to console
                std::cout << "    " << std::setw(14) << pTBin << " | "
                          << formatToThreeSigFigs(data.meanPi0) << " ± " << formatToThreeSigFigs(data.meanPi0Error) << " | "
                          << formatToThreeSigFigs(data.sigmaPi0) << " ± " << formatToThreeSigFigs(data.sigmaPi0Error) << " | "
                          << formatToThreeSigFigs(data.meanEta) << " ± " << formatToThreeSigFigs(data.meanEtaError) << " | "
                          << formatToThreeSigFigs(data.sigmaEta) << " ± " << formatToThreeSigFigs(data.sigmaEtaError) << "\n";

                // Write to CSV file
                csvFile << triggerIndex << ","
                        << Ecore << ","
                        << Chi2 << ","
                        << Asym << ","
                        << data.cuts.pTMin << ","
                        << data.cuts.pTMax << ","
                        << data.meanPi0 << ","
                        << data.sigmaPi0 << ","
                        << data.meanEta << ","
                        << data.sigmaEta << "\n";
            }
            std::cout << "\n"; // Empty line between cut combinations
        }
        std::cout << "\n"; // Empty line between triggers
    }

    // Close the CSV file
    csvFile.close();
}

std::tuple<double, double, double> processDataForTrigger(
    const HistogramData& data,
    std::map<int, std::map<std::string, std::vector<HistogramData>>>& dataMap,
    std::map<int, std::vector<std::string>>& cutVariations,
    std::map<int, std::vector<double>>& meanPi0sPerTrigger,
    std::map<int, std::vector<double>>& meanPi0ErrorsPerTrigger,
    std::map<int, std::vector<double>>& sigmaPi0sPerTrigger,
    std::map<int, std::vector<double>>& sigmaPi0ErrorsPerTrigger,
    std::map<int, std::vector<double>>& signalToBackgroundPi0RatiosPerTrigger,
    std::map<int, std::vector<double>>& signalToBackgroundPi0ErrorsPerTrigger,
    std::map<int, std::vector<double>>& meanEtasPerTrigger,
    std::map<int, std::vector<double>>& meanEtaErrorsPerTrigger,
    std::map<int, std::vector<double>>& sigmaEtasPerTrigger,
    std::map<int, std::vector<double>>& sigmaEtaErrorsPerTrigger,
    std::map<int, std::vector<double>>& signalToBackgroundEtaRatiosPerTrigger,
    std::map<int, std::vector<double>>& signalToBackgroundEtaErrorsPerTrigger,
    std::map<int, std::vector<double>>& signalPi0YieldsPerTrigger,
    std::map<int, std::vector<double>>& signalPi0ErrorsPerTrigger,
    std::map<int, std::vector<double>>& signalEtaYieldsPerTrigger,
    std::map<int, std::vector<double>>& signalEtaErrorsPerTrigger,
    std::map<int, std::vector<double>>& pi0FitResolutionsPerTrigger,
    std::map<int, std::vector<double>>& pi0FitResolutionErrorsPerTrigger,
    std::map<int, std::vector<double>>& etaFitResolutionsPerTrigger,
    std::map<int, std::vector<double>>& etaFitResolutionErrorsPerTrigger) {

    int triggerIndex = data.cuts.triggerIndex;
    double clusECore = data.cuts.clusECore;
    double chi = data.cuts.chi;
    double asymmetry = data.cuts.asymmetry;

    std::ostringstream cutCombinationStream;
    cutCombinationStream << "E" << formatToThreeSigFigs(clusECore)
                         << "_Chi" << formatToThreeSigFigs(chi)
                         << "_Asym" << formatToThreeSigFigs(asymmetry);
    std::string cutCombination = cutCombinationStream.str();

    dataMap[triggerIndex][cutCombination].push_back(data);

    if (data.cuts.pTMin == -1) {
        if (std::find(cutVariations[triggerIndex].begin(), cutVariations[triggerIndex].end(), cutCombination) == cutVariations[triggerIndex].end()) {
            cutVariations[triggerIndex].push_back(cutCombination);

            // Initialize with NaN or default values
            meanPi0sPerTrigger[triggerIndex].push_back(std::numeric_limits<double>::quiet_NaN());
            meanPi0ErrorsPerTrigger[triggerIndex].push_back(0);
            sigmaPi0sPerTrigger[triggerIndex].push_back(std::numeric_limits<double>::quiet_NaN());
            sigmaPi0ErrorsPerTrigger[triggerIndex].push_back(0);
            signalToBackgroundPi0RatiosPerTrigger[triggerIndex].push_back(std::numeric_limits<double>::quiet_NaN());
            signalToBackgroundPi0ErrorsPerTrigger[triggerIndex].push_back(0);

            meanEtasPerTrigger[triggerIndex].push_back(std::numeric_limits<double>::quiet_NaN());
            meanEtaErrorsPerTrigger[triggerIndex].push_back(0);
            sigmaEtasPerTrigger[triggerIndex].push_back(std::numeric_limits<double>::quiet_NaN());
            sigmaEtaErrorsPerTrigger[triggerIndex].push_back(0);
            signalToBackgroundEtaRatiosPerTrigger[triggerIndex].push_back(std::numeric_limits<double>::quiet_NaN());
            signalToBackgroundEtaErrorsPerTrigger[triggerIndex].push_back(0);

            signalPi0YieldsPerTrigger[triggerIndex].push_back(std::numeric_limits<double>::quiet_NaN());
            signalPi0ErrorsPerTrigger[triggerIndex].push_back(0);
            signalEtaYieldsPerTrigger[triggerIndex].push_back(std::numeric_limits<double>::quiet_NaN());
            signalEtaErrorsPerTrigger[triggerIndex].push_back(0);

            // Initialize resolution maps
            pi0FitResolutionsPerTrigger[triggerIndex].push_back(std::numeric_limits<double>::quiet_NaN());
            pi0FitResolutionErrorsPerTrigger[triggerIndex].push_back(0);
            etaFitResolutionsPerTrigger[triggerIndex].push_back(std::numeric_limits<double>::quiet_NaN());
            etaFitResolutionErrorsPerTrigger[triggerIndex].push_back(0);
        }

        size_t index = cutVariations[triggerIndex].size() - 1;

        bool invalidPi0Fit = (data.meanPi0Error > 1.0) || (data.meanPi0Error == 0);
        bool invalidEtaFit = (data.meanEtaError > 1.0) || (data.meanEtaError == 0);

        if (!invalidPi0Fit) {
            meanPi0sPerTrigger[triggerIndex][index] = data.meanPi0;
            meanPi0ErrorsPerTrigger[triggerIndex][index] = data.meanPi0Error;
            sigmaPi0sPerTrigger[triggerIndex][index] = data.sigmaPi0;
            sigmaPi0ErrorsPerTrigger[triggerIndex][index] = data.sigmaPi0Error;
            signalToBackgroundPi0RatiosPerTrigger[triggerIndex][index] = data.signalToBackgroundPi0Ratio;
            signalToBackgroundPi0ErrorsPerTrigger[triggerIndex][index] = data.signalToBackgroundPi0Error;
            signalPi0YieldsPerTrigger[triggerIndex][index] = data.signalPi0Yield;
            signalPi0ErrorsPerTrigger[triggerIndex][index] = data.signalPi0Error;
            pi0FitResolutionsPerTrigger[triggerIndex][index] = data.pi0FitResolution;
            pi0FitResolutionErrorsPerTrigger[triggerIndex][index] = data.pi0FitResolutionError;
        }

        if (!invalidEtaFit) {
            meanEtasPerTrigger[triggerIndex][index] = data.meanEta;
            meanEtaErrorsPerTrigger[triggerIndex][index] = data.meanEtaError;
            sigmaEtasPerTrigger[triggerIndex][index] = data.sigmaEta;
            sigmaEtaErrorsPerTrigger[triggerIndex][index] = data.sigmaEtaError;
            signalToBackgroundEtaRatiosPerTrigger[triggerIndex][index] = data.signalToBackgroundEtaRatio;
            signalToBackgroundEtaErrorsPerTrigger[triggerIndex][index] = data.signalToBackgroundEtaError;
            signalEtaYieldsPerTrigger[triggerIndex][index] = data.signalEtaYield;
            signalEtaErrorsPerTrigger[triggerIndex][index] = data.signalEtaError;
            etaFitResolutionsPerTrigger[triggerIndex][index] = data.etaFitResolution;
            etaFitResolutionErrorsPerTrigger[triggerIndex][index] = data.etaFitResolutionError;
        }
    }

    return std::make_tuple(clusECore, chi, asymmetry);
}


void generateMesonPlotVsPt(
    const std::vector<double>& pTCenters,
    const std::vector<double>& meanValues,
    const std::vector<double>& meanErrors,
    const std::string& yAxisLabel,
    const std::string& outputFilePath,
    const std::string& triggerName,
    const std::string& cutCombination,
    int markerStyle,
    int markerColor,
    double clusECore,
    double chi,
    double asymmetry,
    double yMin = std::numeric_limits<double>::quiet_NaN(),
    double yMax = std::numeric_limits<double>::quiet_NaN(),
    double xMin = std::numeric_limits<double>::quiet_NaN(),
    double xMax = std::numeric_limits<double>::quiet_NaN()) {
    if (pTCenters.empty()) {
        return; // Nothing to plot if pTCenters is empty
    }

    TGraphErrors* graph = new TGraphErrors(pTCenters.size());
    for (size_t i = 0; i < pTCenters.size(); ++i) {
        graph->SetPoint(i, pTCenters[i], meanValues[i]);
        graph->SetPointError(i, 0, meanErrors[i]);
    }
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(1);
    graph->SetLineWidth(2);
    graph->SetMarkerColor(markerColor);
    graph->SetLineColor(markerColor);

    // Automatically calculate y-axis range if yMin or yMax is not provided
    if (std::isnan(yMin) || std::isnan(yMax)) {
        yMin = std::numeric_limits<double>::max();
        yMax = std::numeric_limits<double>::lowest();
        for (size_t i = 0; i < meanValues.size(); ++i) {
            double xVal = pTCenters[i];
            if (xVal >= xMin && xVal <= xMax) {
                double yVal = meanValues[i];
                double yErr = meanErrors[i];
                yMin = std::min(yMin, yVal - yErr);
                yMax = std::max(yMax, yVal + yErr);
            }
        }
        // Add a 5% margin around yMin and yMax
        double yMargin = 0.05 * (yMax - yMin);
        yMin -= yMargin;
        yMax += yMargin;
    }

    // Automatically calculate x-axis range if xMin or xMax is NaN
    if (std::isnan(xMin) || std::isnan(xMax)) {
        xMin = 2;  // Default starting point for x-axis
        xMax = *std::max_element(pTCenters.begin(), pTCenters.end()) * 1.1;  // Adding 10% margin
    }

    // Create canvas
    TCanvas canvas;

    // Draw graph
    TH1F* hFrame = canvas.DrawFrame(xMin, yMin, xMax, yMax, (";p_{T} [GeV];" + yAxisLabel).c_str());
    hFrame->GetXaxis()->SetNdivisions(5);
    hFrame->GetXaxis()->SetLimits(xMin, xMax);
    hFrame->GetXaxis()->CenterLabels(false);

    graph->Draw("P SAME");

    // Add trigger and cut information on the top right of the plot
    TLatex labelText;
    labelText.SetNDC();
    labelText.SetTextSize(0.042);

    TLatex valueText;
    valueText.SetNDC();
    valueText.SetTextSize(0.042);

    labelText.DrawLatex(0.47, 0.87, "#font[62]{Trigger:}");
    valueText.DrawLatex(0.57, 0.87, triggerName.c_str());

    labelText.DrawLatex(0.47, 0.8, "#font[62]{ECore #geq}");
    std::ostringstream eCoreWithUnit;
    eCoreWithUnit << clusECore << "   GeV";
    valueText.DrawLatex(0.62, 0.8, eCoreWithUnit.str().c_str());

    labelText.DrawLatex(0.47, 0.73, "#font[62]{#chi^{2} <}");
    std::ostringstream chiStr;
    chiStr << chi;
    valueText.DrawLatex(0.62, 0.73, chiStr.str().c_str());

    labelText.DrawLatex(0.47, 0.66, "#font[62]{Asymmetry <}");
    std::ostringstream asymmetryStr;
    asymmetryStr << asymmetry;
    valueText.DrawLatex(0.64, 0.66, asymmetryStr.str().c_str());

    // Ensure the directory exists
    std::string outputDirPath = outputFilePath.substr(0, outputFilePath.find_last_of("/"));
    gSystem->mkdir(outputDirPath.c_str(), true);

    // Save plot
    canvas.SaveAs(outputFilePath.c_str());
    std::cout << "Saved plot: " << outputFilePath << std::endl;

    // Clean up
    delete graph;
}

void generateCutVariationPlot(
    const std::vector<double>& xValues,
    const std::vector<double>& xErrors,
    const std::vector<double>& yValues,
    const std::vector<double>& yErrors,
    const std::vector<std::string>& xLabels,
    const std::string& yAxisLabel,
    const std::string& outputFilePath,
    const std::string& triggerName,
    const std::vector<double>& eCoreValues,
    const std::vector<double>& chiValues,
    const std::vector<double>& asymmetryValues,
    int color) {
    if (xValues.empty()) {
        return; // Nothing to plot if xValues is empty
    }

    TGraphErrors* graph = new TGraphErrors(xValues.size(), &xValues[0], &yValues[0], &xErrors[0], &yErrors[0]);
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(color);
    graph->SetMarkerSize(1);
    graph->SetLineWidth(2);
    graph->SetLineColor(color);

    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();
    for (size_t i = 0; i < yValues.size(); ++i) {
        double yVal = yValues[i];
        double yErr = yErrors[i];
        minY = std::min(minY, yVal - yErr);
        maxY = std::max(maxY, yVal + yErr);
    }

    // Add a margin to ensure error bars fit within the range
    double yRange = maxY - minY;
    double yMargin = 0.15 * yRange; // 15% margin
    minY -= yMargin;
    maxY += yMargin;

    // Create canvas
    TCanvas canvas;
    canvas.SetBottomMargin(0.15);  // Increase the bottom margin to accommodate labels

    // Draw graph
    double xMin = 0.5;  // Start at 0.5 to center points between ticks
    double xMax = xValues.size() + 0.5;
    TH1F* hFrame = canvas.DrawFrame(xMin, minY, xMax, maxY, (";Cut Variation;" + yAxisLabel).c_str());
    graph->Draw("P SAME");

    // Set x-axis labels with ticks corresponding to "A", "B", "C", etc.
    TAxis* axis = hFrame->GetXaxis();
    axis->SetNdivisions(xValues.size(), false);
    axis->SetLabelSize(0.03);
    axis->SetTickLength(0.03); // Add visible ticks
    axis->SetLabelOffset(999); // Hide default tick numbers

    // Add simple alphabetical labels between ticks
    for (size_t i = 0; i < xValues.size(); ++i) {
        double x = xValues[i];  // Use xValues directly (centered between ticks)
        TLatex label;
        label.SetTextAlign(22); // Center align
        label.SetTextSize(0.03);
        label.DrawLatex(x, minY - (maxY - minY) * 0.05, xLabels[i].c_str());
    }

    // Add a formatted legend with individual cut values in bold
    TLatex valueText;
    valueText.SetNDC();
    valueText.SetTextSize(0.035);
    double yStart = 0.5; // Starting y-position for the legend

    for (size_t i = 0; i < xLabels.size(); ++i) {
        std::ostringstream entry;
        entry << "#font[62]{" << xLabels[i] << "} : #font[62]{ECore} #geq " << eCoreValues[i]
              << " GeV, #font[62]{#chi^{2}} < " << chiValues[i] << ", #font[62]{Asymmetry} < " << asymmetryValues[i];
        valueText.DrawLatex(0.4, yStart, entry.str().c_str());
        yStart -= 0.05;  // Adjust spacing between legend entries
    }

    // Add trigger name label in bold
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.4, 0.56, ("#font[62]{Trigger:} " + triggerName).c_str());

    // Save plot
    gSystem->mkdir(outputFilePath.substr(0, outputFilePath.find_last_of("/")).c_str(), true);
    canvas.SaveAs(outputFilePath.c_str());
    std::cout << "Saved plot: " << outputFilePath << std::endl;

    // Clean up
    delete graph;
}




void plotMesonInformation(const std::vector<HistogramData>& histogramDataVector) {
    // Organize data by trigger index and cut combination
    std::map<int, std::map<std::string, std::vector<HistogramData>>> dataMap;
    std::map<int, std::vector<std::string>> cutVariations;
    std::map<int, std::vector<double>> meanPi0sPerTrigger, meanPi0ErrorsPerTrigger, sigmaPi0sPerTrigger, sigmaPi0ErrorsPerTrigger;
    std::map<int, std::vector<double>> signalToBackgroundPi0RatiosPerTrigger, signalToBackgroundPi0ErrorsPerTrigger;
    std::map<int, std::vector<double>> meanEtasPerTrigger, meanEtaErrorsPerTrigger, sigmaEtasPerTrigger, sigmaEtaErrorsPerTrigger;
    std::map<int, std::vector<double>> signalToBackgroundEtaRatiosPerTrigger, signalToBackgroundEtaErrorsPerTrigger;
    std::map<int, std::vector<double>> signalPi0YieldsPerTrigger, signalPi0ErrorsPerTrigger;
    std::map<int, std::vector<double>> signalEtaYieldsPerTrigger, signalEtaErrorsPerTrigger;
    std::map<int, std::vector<double>> pi0FitResolutionsPerTrigger, pi0FitResolutionErrorsPerTrigger;
    std::map<int, std::vector<double>> etaFitResolutionsPerTrigger, etaFitResolutionErrorsPerTrigger;


    
    for (const auto& data : histogramDataVector) {
        processDataForTrigger(data, dataMap, cutVariations, meanPi0sPerTrigger, meanPi0ErrorsPerTrigger, sigmaPi0sPerTrigger, sigmaPi0ErrorsPerTrigger, signalToBackgroundPi0RatiosPerTrigger, signalToBackgroundPi0ErrorsPerTrigger, meanEtasPerTrigger, meanEtaErrorsPerTrigger, sigmaEtasPerTrigger, sigmaEtaErrorsPerTrigger, signalToBackgroundEtaRatiosPerTrigger, signalToBackgroundEtaErrorsPerTrigger, signalPi0YieldsPerTrigger, signalPi0ErrorsPerTrigger, signalEtaYieldsPerTrigger, signalEtaErrorsPerTrigger, pi0FitResolutionsPerTrigger, pi0FitResolutionErrorsPerTrigger, etaFitResolutionsPerTrigger, etaFitResolutionErrorsPerTrigger);
    }
    for (const auto& triggerPair : dataMap) {
        int triggerIndex = triggerPair.first;
        std::string triggerName = getTriggerName(triggerIndex);

        for (const auto& cutPair : triggerPair.second) {
            std::string cutCombination = cutPair.first;
            const auto& dataList = cutPair.second;

            // Retrieve cut variables from the first entry in dataList
            double clusECore = dataList[0].cuts.clusECore;
            double chi = dataList[0].cuts.chi;
            double asymmetry = dataList[0].cuts.asymmetry;

            std::vector<double> pTCentersPi0, meanPi0s, meanPi0Errors, sigmaPi0s, sigmaPi0Errors, signalPi0Yields, signalPi0Errors;
            std::vector<double> pTCentersEta, meanEtas, meanEtaErrors, sigmaEtas, sigmaEtaErrors, signalEtaYields, signalEtaErrors;
            std::vector<double> signalToBackgroundPi0Ratios, signalToBackgroundPi0Errors;
            std::vector<double> signalToBackgroundEtaRatios, signalToBackgroundEtaErrors;
            std::vector<double> pi0FitResolutions, pi0FitResolutionErrors;
            std::vector<double> etaFitResolutions, etaFitResolutionErrors;


            for (const auto& data : dataList) {
                if (data.cuts.pTMin == -1) {
                    continue;
                }

                double pTCenter = (data.cuts.pTMin + data.cuts.pTMax) / 2.0;

                bool invalidPi0Fit = (data.meanPi0Error > 1.0) || (data.meanPi0Error == 0);
                bool invalidEtaFit = (data.meanEtaError > 1.0) || (data.meanEtaError == 0);

                if (!invalidPi0Fit) {
                    pTCentersPi0.push_back(pTCenter);
                    meanPi0s.push_back(data.meanPi0);
                    meanPi0Errors.push_back(data.meanPi0Error);
                    sigmaPi0s.push_back(data.sigmaPi0);
                    sigmaPi0Errors.push_back(data.sigmaPi0Error);
                    signalToBackgroundPi0Ratios.push_back(data.signalToBackgroundPi0Ratio);
                    signalToBackgroundPi0Errors.push_back(data.signalToBackgroundPi0Error);
                    signalPi0Yields.push_back(data.signalPi0Yield);
                    signalPi0Errors.push_back(data.signalPi0Error);
                    pi0FitResolutions.push_back(data.pi0FitResolution);
                    pi0FitResolutionErrors.push_back(data.pi0FitResolutionError);
                }

                if (!invalidEtaFit) {
                    pTCentersEta.push_back(pTCenter);
                    meanEtas.push_back(data.meanEta);
                    meanEtaErrors.push_back(data.meanEtaError);
                    sigmaEtas.push_back(data.sigmaEta);
                    sigmaEtaErrors.push_back(data.sigmaEtaError);
                    signalToBackgroundEtaRatios.push_back(data.signalToBackgroundEtaRatio);
                    signalToBackgroundEtaErrors.push_back(data.signalToBackgroundEtaError);
                    signalEtaYields.push_back(data.signalEtaYield);
                    signalEtaErrors.push_back(data.signalEtaError);
                    etaFitResolutions.push_back(data.etaFitResolution);
                    etaFitResolutionErrors.push_back(data.etaFitResolutionError);
                }
            }
            std::string outputDirPath = outputDir + "InvMass/Trigger" + std::to_string(triggerIndex) + "/" + cutCombination + "/";
            generateMesonPlotVsPt(pTCentersPi0, meanPi0s, meanPi0Errors, "#mu_{#pi0} [GeV]",
                                 outputDirPath + "MeanPi0_vs_pT.png", triggerName, cutCombination, 20, kBlue,
                                 clusECore, chi, asymmetry, 0.14, 0.17, 2.0, 7.0);
            generateMesonPlotVsPt(pTCentersEta, meanEtas, meanEtaErrors, "#mu_{#eta} [GeV]", outputDirPath + "MeanEta_vs_pT.png", triggerName, cutCombination, 20, kRed, clusECore, chi, asymmetry, 0.55, 0.62, 2.0, 7.0);
            
            
            generateMesonPlotVsPt(pTCentersPi0, sigmaPi0s, sigmaPi0Errors, "#sigma_{#pi0} [GeV]",
                                  outputDirPath + "SigmaPi0_vs_pT.png", triggerName, cutCombination, 20, kBlue,
                                  clusECore, chi, asymmetry, 0.015, 0.025, 2.0, 7.0);

            generateMesonPlotVsPt(pTCentersEta, sigmaEtas, sigmaEtaErrors, "#sigma_{#eta} [GeV]",
                                  outputDirPath + "SigmaEta_vs_pT.png", triggerName, cutCombination, 20, kRed,
                                  clusECore, chi, asymmetry, 0.045, 0.055, 2.0, 7.0);
            
            generateMesonPlotVsPt(pTCentersPi0, signalToBackgroundPi0Ratios, signalToBackgroundPi0Errors, "Signal-to-Background Ratio (#pi^{0})", outputDirPath + "SignalToBackgroundPi0_vs_pT.png", triggerName, cutCombination, 20, kBlue, clusECore, chi, asymmetry, 0.0, 20.0, 2.0, 7.0);
            
            generateMesonPlotVsPt(pTCentersEta, signalToBackgroundEtaRatios, signalToBackgroundEtaErrors, "Signal-to-Background Ratio (#eta)", outputDirPath + "SignalToBackgroundEta_vs_pT.png", triggerName, cutCombination, 20, kRed, clusECore, chi, asymmetry, 0.0, 2.0, 2.0, 7.0);
            
            generateMesonPlotVsPt(pTCentersPi0, signalPi0Yields, signalPi0Errors, "Signal Yield (#pi^{0})", outputDirPath + "SignalYieldPi0_vs_pT.png", triggerName, cutCombination, 20, kBlue, clusECore, chi, asymmetry, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), 2, 7);
            generateMesonPlotVsPt(pTCentersEta, signalEtaYields, signalEtaErrors, "Signal Yield (#eta)", outputDirPath + "SignalYieldEta_vs_pT.png", triggerName, cutCombination, 20, kRed, clusECore, chi, asymmetry, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), 2, 7);
            generateMesonPlotVsPt(pTCentersPi0, pi0FitResolutions, pi0FitResolutionErrors, "#pi^{0} Fit Resolution", outputDirPath + "Pi0Resolution_vs_pT.png", triggerName, cutCombination, 20, kBlue, clusECore, chi, asymmetry, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), 2.0, 7.0);
            generateMesonPlotVsPt(pTCentersEta, etaFitResolutions, etaFitResolutionErrors, "#eta Fit Resolution", outputDirPath + "EtaResolution_vs_pT.png", triggerName, cutCombination, 20, kRed, clusECore, chi, asymmetry, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), 2.0, 7.0);

        }
    }
    for (const auto& triggerPair : dataMap) {
        int triggerIndex = triggerPair.first;
        std::string triggerName = getTriggerName(triggerIndex);

        const auto& cuts = cutVariations[triggerIndex];
        const auto& meanPi0s = meanPi0sPerTrigger[triggerIndex];
        const auto& meanPi0Errors = meanPi0ErrorsPerTrigger[triggerIndex];
        const auto& sigmaPi0s = sigmaPi0sPerTrigger[triggerIndex];
        const auto& sigmaPi0Errors = sigmaPi0ErrorsPerTrigger[triggerIndex];
        const auto& signalToBackgroundPi0Ratios = signalToBackgroundPi0RatiosPerTrigger[triggerIndex];
        const auto& signalToBackgroundPi0Errors = signalToBackgroundPi0ErrorsPerTrigger[triggerIndex];
        const auto& signalPi0Yields = signalPi0YieldsPerTrigger[triggerIndex];
        const auto& signalPi0Errors = signalPi0ErrorsPerTrigger[triggerIndex];
        const auto& pi0FitResolutions = pi0FitResolutionsPerTrigger[triggerIndex];
        const auto& pi0FitResolutionErrors = pi0FitResolutionErrorsPerTrigger[triggerIndex];

        const auto& meanEtas = meanEtasPerTrigger[triggerIndex];
        const auto& meanEtaErrors = meanEtaErrorsPerTrigger[triggerIndex];
        const auto& sigmaEtas = sigmaEtasPerTrigger[triggerIndex];
        const auto& sigmaEtaErrors = sigmaEtaErrorsPerTrigger[triggerIndex];
        const auto& signalToBackgroundEtaRatios = signalToBackgroundEtaRatiosPerTrigger[triggerIndex];
        const auto& signalToBackgroundEtaErrors = signalToBackgroundEtaErrorsPerTrigger[triggerIndex];
        const auto& signalEtaYields = signalEtaYieldsPerTrigger[triggerIndex];
        const auto& signalEtaErrors = signalEtaErrorsPerTrigger[triggerIndex];
        const auto& etaFitResolutions = etaFitResolutionsPerTrigger[triggerIndex];
        const auto& etaFitResolutionErrors = etaFitResolutionErrorsPerTrigger[triggerIndex];

        // Prepare vectors for cut variation plots
        std::vector<double> xValues, xErrors, eCoreValues, chiValues, asymmetryValues;
        std::vector<std::string> xLabels;
        xErrors.assign(cuts.size(), 0); // Set all xErrors to 0

        for (size_t i = 0; i < cuts.size(); ++i) {
            xValues.push_back(i + 1); // Numerical x-values
            xLabels.push_back(std::string(1, 'A' + i)); // Labels as "A", "B", "C", ...

            // Retrieve cut variables from the first entry in each cut combination
            double clusECore = dataMap[triggerIndex][cuts[i]][0].cuts.clusECore;
            double chi = dataMap[triggerIndex][cuts[i]][0].cuts.chi;
            double asymmetry = dataMap[triggerIndex][cuts[i]][0].cuts.asymmetry;

            eCoreValues.push_back(clusECore);
            chiValues.push_back(chi);
            asymmetryValues.push_back(asymmetry);
        }

        generateCutVariationPlot(xValues, xErrors, meanPi0s, meanPi0Errors, xLabels, "#mu_{#pi^{0}} [GeV]",
                                 outputDir + "InvMass/Trigger" + std::to_string(triggerIndex) + "/MeanPi0_vs_CutVariation.png",
                                 triggerName, eCoreValues, chiValues, asymmetryValues, kBlue);

        generateCutVariationPlot(xValues, xErrors, meanEtas, meanEtaErrors, xLabels, "#mu_{#eta} [GeV]",
                                 outputDir + "InvMass/Trigger" + std::to_string(triggerIndex) + "/MeanEta_vs_CutVariation.png",
                                 triggerName, eCoreValues, chiValues, asymmetryValues, kRed);

        generateCutVariationPlot(xValues, xErrors, sigmaPi0s, sigmaPi0Errors, xLabels, "#sigma_{#pi^{0}} [GeV]",
                                 outputDir + "InvMass/Trigger" + std::to_string(triggerIndex) + "/SigmaPi0_vs_CutVariation.png",
                                 triggerName, eCoreValues, chiValues, asymmetryValues, kBlue);

        generateCutVariationPlot(xValues, xErrors, sigmaEtas, sigmaEtaErrors, xLabels, "#sigma_{#eta} [GeV]",
                                 outputDir + "InvMass/Trigger" + std::to_string(triggerIndex) + "/SigmaEta_vs_CutVariation.png",
                                 triggerName, eCoreValues, chiValues, asymmetryValues, kRed);

        generateCutVariationPlot(xValues, xErrors, signalToBackgroundPi0Ratios, signalToBackgroundPi0Errors, xLabels,
                                 "Signal-to-Background Ratio (#pi^{0})",
                                 outputDir + "InvMass/Trigger" + std::to_string(triggerIndex) + "/SignalToBackgroundPi0_vs_CutVariation.png",
                                 triggerName, eCoreValues, chiValues, asymmetryValues, kBlue);

        generateCutVariationPlot(xValues, xErrors, signalToBackgroundEtaRatios, signalToBackgroundEtaErrors, xLabels,
                                 "Signal-to-Background Ratio (#eta)",
                                 outputDir + "InvMass/Trigger" + std::to_string(triggerIndex) + "/SignalToBackgroundEta_vs_CutVariation.png",
                                 triggerName, eCoreValues, chiValues, asymmetryValues, kRed);

        generateCutVariationPlot(xValues, xErrors, signalPi0Yields, signalPi0Errors, xLabels,
                                 "Signal Yield (#pi^{0})",
                                 outputDir + "InvMass/Trigger" + std::to_string(triggerIndex) + "/SignalYieldPi0_vs_CutVariation.png",
                                 triggerName, eCoreValues, chiValues, asymmetryValues, kBlue);

        generateCutVariationPlot(xValues, xErrors, signalEtaYields, signalEtaErrors, xLabels,
                                 "Signal Yield (#eta)",
                                 outputDir + "InvMass/Trigger" + std::to_string(triggerIndex) + "/SignalYieldEta_vs_CutVariation.png",
                                 triggerName, eCoreValues, chiValues, asymmetryValues, kRed);

        generateCutVariationPlot(xValues, xErrors, pi0FitResolutions, pi0FitResolutionErrors, xLabels,
                                 "#pi^{0} Fit Resolution",
                                 outputDir + "InvMass/Trigger" + std::to_string(triggerIndex) + "/Pi0Resolution_vs_CutVariation.png",
                                 triggerName, eCoreValues, chiValues, asymmetryValues, kBlue);

        generateCutVariationPlot(xValues, xErrors, etaFitResolutions, etaFitResolutionErrors, xLabels,
                                 "#eta Fit Resolution",
                                 outputDir + "InvMass/Trigger" + std::to_string(triggerIndex) + "/EtaResolution_vs_CutVariation.png",
                                 triggerName, eCoreValues, chiValues, asymmetryValues, kRed);
    }
}

// Main function to run both processes
void plotInvariantMass() {
    gROOT->LoadMacro("sPhenixStyle.C");
    SetsPhenixStyle();
    saveAnnotatedInvariantMassHistograms(inputFilePath);
    std::vector<HistogramData> histogramDataVector = saveAnnotatedInvariantMassHistograms(inputFilePath);
    printHistogramData(histogramDataVector);
    plotMesonInformation(histogramDataVector);
}
