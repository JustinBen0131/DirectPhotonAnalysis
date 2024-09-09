#include <sPhenixStyle.C>
#include <regex>
#include <map>
#include <sstream>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <vector>
#include <string>

// Structure to store the parsed cut values
struct CutValues {
    float clusEA;
    float clusEB;
    float asymmetry;
    float deltaRMin;
    float deltaRMax;
    float chi;
};

// Function to convert strings with 'point' back to float
float convert(const std::string& input) {
    std::string temp = input;
    size_t pointPos = temp.find("point");
    if (pointPos != std::string::npos) {
        temp.replace(pointPos, 5, ".");
    }
    try {
        return std::stof(temp);
    } catch (const std::exception&) {
        return 0.0f;  // Return 0.0 if conversion fails
    }
}

// Function to parse filename to extract the cut values
CutValues parseFileName(const std::string& filename) {
    CutValues cuts = {0, 0, 0, 0, 0, 0}; // Initialize cut values
    std::regex re("hDiphotonMass_E([0-9]+(?:point[0-9]*)?)_Asym([0-9]+(?:point[0-9]*)?)_DelR([0-9]+(?:point[0-9]*)?)_Chi([0-9]+(?:point[0-9]*)?)");

    std::smatch match; // Object to store regex results
    if (std::regex_search(filename, match, re) && match.size() > 4) {
        cuts.clusEA = convert(match[1].str());  // Energy Cut 1
        cuts.asymmetry = convert(match[2].str());  // Asymmetry Cut
        cuts.deltaRMin = convert(match[3].str());  // Delta R min
        cuts.chi = convert(match[4].str());  // Chi2 cut
    } else {
        std::cerr << "Filename format did not match: " << filename << std::endl;
    }

    return cuts;
}

// Function to format float values for histogram naming (replaces '.' with 'point')
std::string formatFloatForFilename(float value) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(3) << value;
    std::string str = ss.str();
    size_t dotPos = str.find('.');
    if (dotPos != std::string::npos) {
        str = str.substr(0, dotPos) + "point" + str.substr(dotPos + 1);
    }
    return str;
}

void makePhotonPairs(int nEvents = 0, const std::string &fin = "commissioning.root", const std::string &fout = "output.root") {
    SetsPhenixStyle();

    // Open the input ROOT file
    TFile *f = new TFile(fin.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << fin << std::endl;
        return;
    }

    // Get the tree from the input file
    TTree *T = (TTree *)f->Get("T");
    if (!T) {
        std::cerr << "Error: Cannot find tree 'T' in file " << fin << std::endl;
        return;
    }

    std::vector<float> *clusterE{0};
    std::vector<float> *clusterPhi{0};
    std::vector<float> *clusterEta{0};
    std::vector<float> *clusterPt{0};
    std::vector<float> *clusterChi2{0};
    uint64_t b_gl1_scaledvec = 0;

    T->SetBranchAddress("clusterE", &clusterE);
    T->SetBranchAddress("clusterPhi", &clusterPhi);
    T->SetBranchAddress("clusterEta", &clusterEta);
    T->SetBranchAddress("clusterPt", &clusterPt);
    T->SetBranchAddress("clusterChi2", &clusterChi2);
    T->SetBranchAddress("gl1_scaledvec", &b_gl1_scaledvec);

    // Create output file
    TFile *outFile = new TFile(fout.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error: Cannot create output file " << fout << std::endl;
        return;
    }

    // Define cut values
    std::vector<float> maxChi2_values = {3, 4, 5};
    std::vector<float> minClusE_values = {0.5, 0.75, 1, 1.25, 1.5};
    std::vector<float> maxAsym_values = {0.4, 0.5, 0.6};
    std::vector<float> minDeltaR_values = {0, 0.06, 0.07, 0.08};

    // Create histograms for each cut combination
    std::map<std::string, TH1F *> histograms;
    for (float maxChi2 : maxChi2_values) {
        for (float minClusE : minClusE_values) {
            for (float maxAsym : maxAsym_values) {
                for (float minDeltaR : minDeltaR_values) {
                    std::string histName = "hDiphotonMass_E" + formatFloatForFilename(minClusE) +
                                           "_Asym" + formatFloatForFilename(maxAsym) +
                                           "_DelR" + formatFloatForFilename(minDeltaR) +
                                           "_Chi" + formatFloatForFilename(maxChi2);
                    histograms[histName] = new TH1F(histName.c_str(), histName.c_str(), 100, 0, 1);
                    histograms[histName]->SetTitle(";M_{#gamma#gamma};");
                }
            }
        }
    }

    // Get the total number of entries (events)
    int nEntries = (nEvents > 0) ? nEvents : T->GetEntries();
    for (int i = 0; i < nEntries; i++) {
        T->GetEntry(i);
        
        // Vector to store trigger bits
        std::vector<int> trig_bits{};
     
        // Loop over each bit position from 0 to 63
        for (unsigned int bit = 0; bit < 64; bit++) {
            /*
             Check if the bit at the specified position 'bit' in the variable 'b_gl1_scaledvec' is set to 1.
             This is done by shifting 'b_gl1_scaledvec' right by 'bit' positions, which moves the bit of interest to the
             least significant bit position.
           
             Then, it is bitwise AND with 0x1U to isolate this bit. If the result is 1,
             then the bit at position 'bit' was originally set to 1.
           */
            if (((b_gl1_scaledvec >> bit) & 0x1U) == 0x1U) {
              // If the bit is set to 1, add its position to the vector 'trig_bits'
                trig_bits.push_back(bit);
            }
        }

        // Now you have a vector 'trig_bits' containing all set bit positions.
        // You can then check whether specific bits (e.g., 24, 25, 26, 27) are present.
        bool skip = true;
        for (const int& bit : trig_bits) {
            // Check if the bit is either in the range 24-27 or is exactly 10
            if ((bit >= 24 && bit <= 27) || bit == 10) {
                skip = false;  // At least one of the desired bits is set, so do not skip the event
                break;  // No need to check further if a matching bit is found
            }
        }

        // If none of the bits 24, 25, 26, 27, or 10 are set, the event will be skipped
        if (skip) continue;
          

        // Loop over photon candidates
        std::vector<TLorentzVector> photonCandidates;
        for (int clus1 = 0; clus1 < clusterE->size(); clus1++) {
            if (clusterChi2->at(clus1) > *std::max_element(maxChi2_values.begin(), maxChi2_values.end())) continue;
            if (clusterE->at(clus1) < *std::min_element(minClusE_values.begin(), minClusE_values.end())) continue;

            TLorentzVector photon1;
            photon1.SetPtEtaPhiE(clusterPt->at(clus1), clusterEta->at(clus1), clusterPhi->at(clus1), clusterE->at(clus1));
            photonCandidates.push_back(photon1);
        }

        // Process photon pairs
        for (int clus1 = 0; clus1 < photonCandidates.size(); clus1++) {
            for (int clus2 = clus1 + 1; clus2 < photonCandidates.size(); clus2++) {
                TLorentzVector &photon1 = photonCandidates[clus1];
                TLorentzVector &photon2 = photonCandidates[clus2];

                float asym = fabs(photon1.E() - photon2.E()) / (photon1.E() + photon2.E());
                float delta_eta = photon1.Eta() - photon2.Eta();
                float delta_phi = fabs(photon1.Phi() - photon2.Phi());
                if (delta_phi > M_PI) delta_phi = 2 * M_PI - delta_phi;
                float delta_R = TMath::Sqrt(delta_eta * delta_eta + delta_phi * delta_phi);

                // Fill histograms based on cuts
                for (auto &hist : histograms) {
                    CutValues cuts = parseFileName(hist.first);

                    // Apply cuts
                    if (clusterChi2->at(clus1) > cuts.chi || clusterChi2->at(clus2) > cuts.chi) continue;
                    if (clusterE->at(clus1) < cuts.clusEA || clusterE->at(clus2) < cuts.clusEA) continue;
                    if (asym > cuts.asymmetry || delta_R < cuts.deltaRMin) continue;

                    TLorentzVector meson = photon1 + photon2;
                    hist.second->Fill(meson.M());
                }
            }
        }
    }

    // Write histograms
    for (auto &hist : histograms) {
        hist.second->Write();
        delete hist.second;
    }

    outFile->Close();
    delete f;
    delete outFile;

    std::cout << "Processed all events with different cut variations. Output saved to " << fout << std::endl;
}
