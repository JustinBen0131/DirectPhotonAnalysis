#ifndef ANALYZE_TRIGGER_GROUPINGS_H
#define ANALYZE_TRIGGER_GROUPINGS_H

#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <cctype>

namespace ReferenceData {

    // Define first reference dataset
    const std::vector<double> referencePTGamma = {3.36, 4.39, 5.41, 6.42, 7.43, 8.44, 9.80, 11.83, 14.48};
    const std::vector<double> referenceRatio = {0.594, 0.664, 0.626, 0.658, 0.900, 0.715, 0.872, 0.907, 0.802};
    const std::vector<double> referenceStatError = {0.014, 0.028, 0.043, 0.061, 0.113, 0.130, 0.120, 0.190, 0.290};

    // Define second reference data (from the second table - Reference Two)
    const std::vector<double> referenceTwoPTGamma = {3.34, 4.38, 5.40, 6.41, 7.42, 8.43, 9.78, 11.81, 14.41};
    const std::vector<double> referenceTwoRatio = {0.477, 0.455, 0.448, 0.430, 0.338, 0.351, 0.400, 0.286, 0.371};
    const std::vector<double> referenceTwoStatError = {0.0020, 0.0060, 0.012, 0.021, 0.032, 0.053, 0.070, 0.130, 0.180};

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
        double binWidth;
        double binCenter;
        double isolatedYield;
        double isolatedYieldError;
    };


    struct IsolationDataWithPt {
        float ptMin;
        float ptMax;
        double weightedPt;
        double ratio;
        double error;
        double isolatedYield;         // Ensure this member exists
        double isolatedYieldError;    // Ensure this member exists
        float isoMin;
        float isoMax;
        std::string triggerName;
    };


} // namespace DataStructures


// Namespace for trigger configurations
namespace TriggerConfig {
    // List of all triggers we're interested in
    const std::vector<std::string> allTriggers = {
        "MBD_NandS_geq_1",
        "Photon_2_GeV_plus_MBD_NS_geq_1",
        "Photon_3_GeV_plus_MBD_NS_geq_1",
        "Photon_4_GeV_plus_MBD_NS_geq_1",
        "Photon_5_GeV_plus_MBD_NS_geq_1",
//        "Photon_2_GeV",
//        "Photon_3_GeV",
//        "Photon_4_GeV",
//        "Photon_5_GeV"
    };

    // List of Photon triggers (excluding MBD_NandS_geq_1)
    const std::vector<std::string> photonTriggers = {
        "Photon_2_GeV_plus_MBD_NS_geq_1",
        "Photon_3_GeV_plus_MBD_NS_geq_1",
        "Photon_4_GeV_plus_MBD_NS_geq_1",
        "Photon_5_GeV_plus_MBD_NS_geq_1",
//        "Photon_2_GeV",
//        "Photon_3_GeV",
//        "Photon_4_GeV",
//        "Photon_5_GeV"
    };

    // Map of triggers to colors for plotting
    const std::map<std::string, int> triggerColorMap = {
        {"MBD_NandS_geq_1", kBlack},
        {"Photon_2_GeV_plus_MBD_NS_geq_1", kRed},
        {"Photon_3_GeV_plus_MBD_NS_geq_1", kBlue},
        {"Photon_4_GeV_plus_MBD_NS_geq_1", kGreen + 2},
        {"Photon_5_GeV_plus_MBD_NS_geq_1", kMagenta},
//        {"Photon_2_GeV", kRed},
//        {"Photon_3_GeV", kBlue},
//        {"Photon_4_GeV", kGreen + 2},
//        {"Photon_5_GeV", kMagenta}
    };

    // Map of triggers to human-readable names
    const std::map<std::string, std::string> triggerNameMap = {
        {"MBD_NandS_geq_1", "Minbias"},
        {"Photon_2_GeV_plus_MBD_NS_geq_1", "Photon 2 GeV + Minbias"},
        {"Photon_3_GeV_plus_MBD_NS_geq_1", "Photon 3 GeV + Minbias"},
        {"Photon_4_GeV_plus_MBD_NS_geq_1", "Photon 4 GeV + Minbias"},
        {"Photon_5_GeV_plus_MBD_NS_geq_1", "Photon 5 GeV + Minbias"},
//        {"Photon_2_GeV", "Photon 2 GeV"},
//        {"Photon_3_GeV", "Photon 3 GeV"},
//        {"Photon_4_GeV", "Photon 4 GeV"},
//        {"Photon_5_GeV", "Photon 5 GeV"}
        
    };

    const std::map<std::pair<float, float>, int> isoEtRangeColorMap = {
        {{-100, 6}, kRed + 1},
        {{-100, 10}, kRed + 1},
        {{-10, 0}, kGreen + 2},
        {{0, 10}, kMagenta + 2}
    };

    const std::map<std::pair<std::string, std::string>, DataStructures::FitParameters> triggerFitParameters = {
        { {"MBD_NandS_geq_1_Photon_2_GeV_plus_MBD_NS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1", "Photon_2_GeV_plus_MBD_NS_geq_1"}, {
            1.25,   // amplitudeEstimate (closer to 1 for convergence)
            0.72,  // slopeEstimate (increased for sharper rise)
            4.65,   // xOffsetEstimate (shifted left to match flatter region)
            1.0,  // amplitudeMin
            1.3,  // amplitudeMax
            0.71,   // slopeMin
            0.73,   // slopeMax
            4.55,   // xOffsetMin
            4.75    // xOffsetMax
        } },
        { {"MBD_NandS_geq_1_Photon_2_GeV_plus_MBD_NS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1", "Photon_3_GeV_plus_MBD_NS_geq_1"}, {
            1.25,   // amplitudeEstimate (adjusted for y = 1 convergence)
            0.64,   // slopeEstimate (slightly increased for better rise)
            7.2,   // xOffsetEstimate (shifted for alignment)
            1.0,  // amplitudeMin
            1.3,  // amplitudeMax
            0.62,   // slopeMin
            0.66,  // slopeMax
            7.0,   // xOffsetMin
            7.5    // xOffsetMax
        } },
        { {"MBD_NandS_geq_1_Photon_2_GeV_plus_MBD_NS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1", "Photon_4_GeV_plus_MBD_NS_geq_1"}, {
            1.25,   // amplitudeEstimate (adjusted for y = 1 convergence)
            0.64,   // slopeEstimate (slightly increased for better rise)
            7.2,   // xOffsetEstimate (shifted for alignment)
            1.0,  // amplitudeMin
            1.3,  // amplitudeMax
            0.62,   // slopeMin
            0.66,  // slopeMax
            7.0,   // xOffsetMin
            7.5    // xOffsetMax
        } },
        { {"MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1_Photon_5_GeV_plus_MBD_NS_geq_1", "Photon_3_GeV_plus_MBD_NS_geq_1"}, {
            1.25,   // amplitudeEstimate (adjusted for y = 1 convergence)
            0.64,   // slopeEstimate (slightly increased for better rise)
            7.2,   // xOffsetEstimate (shifted for alignment)
            1.0,  // amplitudeMin
            1.3,  // amplitudeMax
            0.62,   // slopeMin
            0.66,  // slopeMax
            7.0,   // xOffsetMin
            7.5    // xOffsetMax
        } },
        { {"MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1_Photon_5_GeV_plus_MBD_NS_geq_1", "Photon_4_GeV_plus_MBD_NS_geq_1"}, {
            1.25,   // amplitudeEstimate (adjusted for y = 1 convergence)
            0.64,   // slopeEstimate (slightly increased for better rise)
            7.2,   // xOffsetEstimate (shifted for alignment)
            1.0,  // amplitudeMin
            1.3,  // amplitudeMax
            0.62,   // slopeMin
            0.66,  // slopeMax
            7.0,   // xOffsetMin
            7.5    // xOffsetMax
        } },
        { {"MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1_Photon_5_GeV_plus_MBD_NS_geq_1", "Photon_5_GeV_plus_MBD_NS_geq_1"}, {
            1.0,   // amplitudeEstimate
            0.48,  // slopeEstimate (reduced for smoother rise)
            10.5,   // xOffsetEstimate (slightly adjusted)
            0.98,  // amplitudeMin
            1.05,  // amplitudeMax
            0.46,   // slopeMin
            0.5,   // slopeMax
            10.4,   // xOffsetMin
            10.6    // xOffsetMax
        } },
        { {"", "Photon_2_GeV_plus_MBD_NS_geq_1"}, {
            1.25,   // amplitudeEstimate (closer to 1 for convergence)
            0.72,  // slopeEstimate (increased for sharper rise)
            4.65,   // xOffsetEstimate (shifted left to match flatter region)
            1.0,  // amplitudeMin
            1.3,  // amplitudeMax
            0.71,   // slopeMin
            0.73,   // slopeMax
            4.55,   // xOffsetMin
            4.75    // xOffsetMax
        } },
        { {"", "Photon_3_GeV_plus_MBD_NS_geq_1"}, {
            1.25,   // amplitudeEstimate (adjusted for y = 1 convergence)
            0.64,   // slopeEstimate (slightly increased for better rise)
            7.2,   // xOffsetEstimate (shifted for alignment)
            1.0,  // amplitudeMin
            1.3,  // amplitudeMax
            0.62,   // slopeMin
            0.66,  // slopeMax
            7.0,   // xOffsetMin
            7.5    // xOffsetMax
        } },
        { {"", "Photon_4_GeV_plus_MBD_NS_geq_1"}, {
            1.25,   // amplitudeEstimate (adjusted for y = 1 convergence)
            0.55,   // slopeEstimate (slightly increased for better rise)
            8.5,   // xOffsetEstimate (shifted for alignment)
            1.0,  // amplitudeMin
            1.3,  // amplitudeMax
            0.53,   // slopeMin
            0.57,  // slopeMax
            8.4,   // xOffsetMin
            8.6    // xOffsetMax
        } },
        { {"", "Photon_5_GeV_plus_MBD_NS_geq_1"}, {
            1.0,   // amplitudeEstimate
            0.48,  // slopeEstimate (reduced for smoother rise)
            10.5,   // xOffsetEstimate (slightly adjusted)
            0.98,  // amplitudeMin
            1.05,  // amplitudeMax
            0.46,   // slopeMin
            0.5,   // slopeMax
            10.4,   // xOffsetMin
            10.6    // xOffsetMax
        } }
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
        {"MBD_NandS_geq_1_Photon_2_GeV_plus_MBD_NS_geq_1_Photon_5_GeV_plus_MBD_NS_geq_1", "Minbias and Photon 2, 5 GeV"},
        {"MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1", "Minbias and Photon 3, 4 GeV"},
        {"MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_5_GeV_plus_MBD_NS_geq_1", "Minbias and Photon 3, 5 GeV"},
        {"MBD_NandS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1_Photon_5_GeV_plus_MBD_NS_geq_1", "Minbias and Photon 4, 5 GeV"},
        {"MBD_NandS_geq_1_Photon_2_GeV_plus_MBD_NS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1", "Minbias and Photon 2, 3, 4 GeV"},
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

#endif // ANALYZE_TRIGGER_GROUPINGS_H
