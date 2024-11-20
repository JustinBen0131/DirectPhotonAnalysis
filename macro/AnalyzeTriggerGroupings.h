#ifndef ANALYZE_TRIGGER_GROUPINGS_H
#define ANALYZE_TRIGGER_GROUPINGS_H

#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <cctype>

// Namespace for trigger configurations
namespace TriggerConfig {
    // List of all triggers we're interested in
    const std::vector<std::string> allTriggers = {
        "MBD_NandS_geq_1",
        "Photon_2_GeV_plus_MBD_NS_geq_1",
        "Photon_3_GeV_plus_MBD_NS_geq_1",
        "Photon_4_GeV_plus_MBD_NS_geq_1",
        "Photon_5_GeV_plus_MBD_NS_geq_1"
    };

    // List of Photon triggers (excluding MBD_NandS_geq_1)
    const std::vector<std::string> photonTriggers = {
        "Photon_2_GeV_plus_MBD_NS_geq_1",
        "Photon_3_GeV_plus_MBD_NS_geq_1",
        "Photon_4_GeV_plus_MBD_NS_geq_1",
        "Photon_5_GeV_plus_MBD_NS_geq_1"
    };

    // Map of triggers to colors for plotting
    const std::map<std::string, int> triggerColorMap = {
        {"MBD_NandS_geq_1", kBlack},
        {"Photon_2_GeV_plus_MBD_NS_geq_1", kRed},
        {"Photon_3_GeV_plus_MBD_NS_geq_1", kBlue},
        {"Photon_4_GeV_plus_MBD_NS_geq_1", kGreen + 2},
        {"Photon_5_GeV_plus_MBD_NS_geq_1", kMagenta}
    };

    // Map of triggers to human-readable names
    const std::map<std::string, std::string> triggerNameMap = {
        {"MBD_NandS_geq_1", "Minbias"},
        {"Photon_2_GeV_plus_MBD_NS_geq_1", "Photon 2 GeV + Minbias"},
        {"Photon_3_GeV_plus_MBD_NS_geq_1", "Photon 3 GeV + Minbias"},
        {"Photon_4_GeV_plus_MBD_NS_geq_1", "Photon 4 GeV + Minbias"},
        {"Photon_5_GeV_plus_MBD_NS_geq_1", "Photon 5 GeV + Minbias"},
    };

    const std::map<std::pair<float, float>, int> isoEtRangeColorMap = {
        {{-100, 6}, kRed + 1},
        {{-100, 10}, kViolet + 1},
        {{-10, 0}, kGreen + 2},
        {{0, 10}, kMagenta + 2}
    };
}


// Namespace for trigger combination names
namespace TriggerCombinationNames {
    const std::map<std::string, std::string> triggerCombinationNameMap = {
        {"MBD_NandS_geq_1", "Minbias"},
        {"MBD_NandS_geq_1_Photon_2_GeV_plus_MBD_NS_geq_1", "Minbias and Photon 2 GeV"},
        {"MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1", "Minbias and Photon 3 GeV"},
        {"MBD_NandS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1", "Minbias and Photon 4 GeV"},
        {"MBD_NandS_geq_1_Photon_5_GeV_plus_MBD_NS_geq_1", "Minbias and Photon 5 GeV"},
        {"MBD_NandS_geq_1_Photon_2_GeV_plus_MBD_NS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1", "Minbias and Photon 2, 3 GeV"},
        {"MBD_NandS_geq_1_Photon_2_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1", "Minbias and Photon 2, 4 GeV"},
        {"MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1", "Minbias and Photon 3, 4 GeV"},
        {"MBD_NandS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1_Photon_5_GeV_plus_MBD_NS_geq_1", "Minbias and Photon 4, 5 GeV"},
        {"MBD_NandS_geq_1_Photon_2_GeV_plus_MBD_NS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1", "Minbias and Photon 2, 3, 4 GeV"},
        {"MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1_Photon_5_GeV_plus_MBD_NS_geq_1", "Minbias and Photon 3, 4, 5 GeV"},
    };
}

namespace Utils {
    // Helper function to normalize the trigger combination string for case-insensitive and whitespace-insensitive comparison
    std::string normalizeString(const std::string& str) {
        std::string normalized = str;
        std::transform(normalized.begin(), normalized.end(), normalized.begin(), ::tolower);
        normalized.erase(std::remove_if(normalized.begin(), normalized.end(), ::isspace), normalized.end());
        return normalized;
    }

    // Function to retrieve the human-readable combination name
    std::string getTriggerCombinationName(const std::string& dirName, const std::map<std::string, std::string>& nameMap) {
        std::string normalizedDir = normalizeString(dirName);
        for (const auto& entry : nameMap) {
            if (normalizeString(entry.first) == normalizedDir) {
                return entry.second;
            }
        }
        return dirName; // Default to directory name if not found
    }

    // Function to format a double to three significant figures as a string
    std::string formatToThreeSigFigs(double value) {
        std::ostringstream out;
        out << std::fixed << std::setprecision(3) << value;
        return out.str();
    }

    // Function to check if a string ends with another string
    bool EndsWith(const std::string& fullString, const std::string& ending) {
        if (fullString.length() >= ending.length()) {
            return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
        } else {
            return false;
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

        // Set initial parameters
        fitFunc->SetParameter(0, amplitude);  // Amplitude
        fitFunc->SetParameter(1, slope);      // Slope
        fitFunc->SetParameter(2, xOffset);    // XOffset

        // Set parameter limits
        fitFunc->SetParLimits(0, amplitudeMin, amplitudeMax);  // Amplitude limits
        fitFunc->SetParLimits(1, slopeMin, slopeMax);          // Slope limits
        fitFunc->SetParLimits(2, xOffsetMin, xOffsetMax);      // XOffset limits

        return fitFunc;
    }

}

// Define CutValues, FitParameters, and HistogramData within a dedicated namespace
namespace DataStructures {

    struct CutValues {
        float clusECore = 0;
        float asymmetry = 0;
        float chi = 0;
        std::string triggerName;
        float pTMin = -1;  // Default to -1 indicating no pT bin
        float pTMax = -1;  // Default to -1 indicating no pT bin
    };

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

    struct FitParameters {
        double amplitudeEstimate;
        double slopeEstimate;
        double xOffsetEstimate;
        double amplitudeMin;
        double amplitudeMax;
        double slopeMin;
        double slopeMax;
        double xOffsetMin;
        double xOffsetMax;
    };

    struct IsolatedPhotonLog {
        std::string triggerGroupName;
        std::string triggerName;
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
        std::string triggerGroupName;
        std::string triggerName;
        float clusECore;
        float chi;
        float asymmetry;
        float pTMin;
        float pTMax;
        int totalEntries;
        std::string massWindowLabel;
    };

    struct PtWeightingLog {
        std::string triggerGroupName;
        std::string triggerName;
        float clusECore;
        float chi;
        float asymmetry;
        float pTMin;
        float pTMax;
        double weightedAveragePt;
        std::string massWindowLabel;
    };

    struct IsolationData {
        int isolatedCounts;
        int totalCounts;
        double ratio;
        double error;
        double weightedPt;
    };

    struct IsolationDataWithPt {
        float ptMin;
        float ptMax;
        double weightedPt;
        double ratio;
        double error;
        float isoMin;  // Added isoMin
        float isoMax;  // Added isoMax
    };

} // namespace DataStructures

#endif // ANALYZE_TRIGGER_GROUPINGS_H
