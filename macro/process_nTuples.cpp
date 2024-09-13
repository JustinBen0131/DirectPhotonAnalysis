#include <iostream>
#include <fstream>  // To write to log file
#include <vector>
#include <string>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TSystem.h>
#include <TDirectory.h>
#include <TList.h>
#include <TClass.h>

void logCorruptedFile(const std::string& filePath) {
    std::ofstream logFile("/sphenix/tg/tg01/bulk/jbennett/DirectPhotons/outputHists/corrupted_files.txt", std::ios::app);
    if (logFile.is_open()) {
        logFile << "Corrupted file: " << filePath << std::endl;
        logFile.close();
    } else {
        std::cerr << "Error: Unable to open log file to record corrupted file: " << filePath << std::endl;
    }
}
void setTreeBranches(TTree* tree, const std::string& fin,
                     ULong64_t* b_gl1_scaledvec, float* zvertex,
                     std::vector<float>** clusterE, std::vector<float>** clusterECore,
                     std::vector<float>** clusterTowMaxE, std::vector<float>** clusterPhi,
                     std::vector<float>** clusterEta, std::vector<float>** clusterPt,
                     std::vector<float>** clusterChi2, std::vector<float>** clusterNtow,
                     float* totalCaloEEMCal, std::vector<float>** emcTowE, std::vector<int>** emcTowiEta,
                     std::vector<int>** emcTowiPhi, float* totalCaloEOHCal, std::vector<float>** ohcTowE,
                     std::vector<float>** ohcTowiEta, std::vector<float>** ohcTowiPhi, std::vector<float>** ohcChi2,
                     float* totalCaloEIHCal, std::vector<float>** ihcTowE, std::vector<float>** ihcTowiEta,
                     std::vector<float>** ihcTowiPhi, std::vector<float>** ihcChi2) {

    
    if (tree->SetBranchAddress("gl1_scaledvec", b_gl1_scaledvec) < 0 ||
        tree->SetBranchAddress("zvertex", zvertex) < 0 ||
        tree->SetBranchAddress("clusterE", clusterE) < 0 ||
        tree->SetBranchAddress("clusterECore", clusterECore) < 0 ||
        tree->SetBranchAddress("clusterTowMaxE", clusterTowMaxE) < 0 ||
        tree->SetBranchAddress("clusterPhi", clusterPhi) < 0 ||
        tree->SetBranchAddress("clusterEta", clusterEta) < 0 ||
        tree->SetBranchAddress("clusterPt", clusterPt) < 0 ||
        tree->SetBranchAddress("clusterChi2", clusterChi2) < 0 ||
        tree->SetBranchAddress("clusterNtow", clusterNtow) < 0 ||
        tree->SetBranchAddress("totalCaloEEMCal", totalCaloEEMCal) < 0 ||
        tree->SetBranchAddress("emcTowE", emcTowE) < 0 ||
        tree->SetBranchAddress("emcTowiEta", emcTowiEta) < 0 ||
        tree->SetBranchAddress("emcTowiPhi", emcTowiPhi) < 0 ||
        tree->SetBranchAddress("totalCaloEOHCal", totalCaloEOHCal) < 0 ||
        tree->SetBranchAddress("ohcTowE", ohcTowE) < 0 ||
        tree->SetBranchAddress("ohcTowiEta", ohcTowiEta) < 0 ||
        tree->SetBranchAddress("ohcTowiPhi", ohcTowiPhi) < 0 ||
        tree->SetBranchAddress("ohcChi2", ohcChi2) < 0 ||
        tree->SetBranchAddress("totalCaloEIHCal", totalCaloEIHCal) < 0 ||
        tree->SetBranchAddress("ihcTowE", ihcTowE) < 0 ||
        tree->SetBranchAddress("ihcTowiEta", ihcTowiEta) < 0 ||
        tree->SetBranchAddress("ihcTowiPhi", ihcTowiPhi) < 0 ||
        tree->SetBranchAddress("ihcChi2", ihcChi2) < 0) {
        std::cerr << "Error: Corrupted or missing branch in file: " << fin << std::endl;
        logCorruptedFile(fin);  // Log corrupted branch
        throw std::runtime_error("Branch corruption detected, terminating.");
    }
}
void fillQAHistograms(TH2F* h2_EMCal_TowerEtaPhi_2D, TH2F* h2_OHCal_TowerEnergy, TH2F* h2_IHCal_TowerEnergy,
                      TH1F* hTotalCaloEEMCal, TH1F* hTotalCaloEOHCal, TH1F* hTotalCaloEIHCal, TH1F* h8by8TowerEnergySum,
                      TH1F* hClusterECore, TH1F* hClusterPt, TH1F* hClusterChi2, TH1F* h_ohcalChi2, TH1F* h_ihcalChi2, TH1F* hVtxZ,
                      const std::vector<float>* emcTowE, const std::vector<int>* emcTowiEta, const std::vector<int>* emcTowiPhi,
                      const std::vector<float>* ohcTowE, const std::vector<float>* ohcTowiEta, const std::vector<float>* ohcTowiPhi,
                      const std::vector<float>* ihcTowE, const std::vector<float>* ihcTowiEta, const std::vector<float>* ihcTowiPhi,
                      const std::vector<float>* clusterECore, const std::vector<float>* clusterChi2, const std::vector<float>* clusterPt,
                      const std::vector<float>* ohcChi2, const std::vector<float>* ihcChi2, const float* totalCaloEEMCal,
                      const float* totalCaloEOHCal, const float* totalCaloEIHCal, const float* zvertex) {

    // Initialize a 2D array to store energy sums for each 8x8 block
    float energymap[12][35] = {0};
    // Calculate the energy sum for each 8x8 block in the EMCal
    for (size_t j = 0; j < emcTowE->size(); j++) {
        int iEta_emc = emcTowiEta->at(j);
        int iPhi_emc = emcTowiPhi->at(j);
        float energy_emc = emcTowE->at(j);
        // Fill the 2D energy distribution histogram
        h2_EMCal_TowerEtaPhi_2D->Fill(iEta_emc, iPhi_emc, energy_emc);

        // Sum the energy into 8x8 blocks
        int ebin = iEta_emc / 8;
        int pbin = iPhi_emc / 8;
        energymap[ebin][pbin] += energy_emc;
    }
    // Find the maximum energy sum in the 8x8 blocks
    float max_energy_sum = 0.0;
    for (int ebin = 0; ebin < 12; ++ebin) {
        for (int pbin = 0; pbin < 35; ++pbin) {
            if (energymap[ebin][pbin] > max_energy_sum) {
                max_energy_sum = energymap[ebin][pbin];
            }
        }
    }
    // Fill the maximum 8x8 tower energy sum histogram
    h8by8TowerEnergySum->Fill(max_energy_sum);
    // Fill OHCal Tower Energy
    for (size_t j = 0; j < ohcTowE->size(); j++) {
        int iEta_ohc = ohcTowiEta->at(j);
        int iPhi_ohc = ohcTowiPhi->at(j);
        float energy_ohc = ohcTowE->at(j);
        h2_OHCal_TowerEnergy->Fill(iEta_ohc, iPhi_ohc, energy_ohc);
    }
    // Fill IHCal Tower Energy
    for (size_t j = 0; j < ihcTowE->size(); j++) {
        int iEta_ihc = ihcTowiEta->at(j);
        int iPhi_ihc = ihcTowiPhi->at(j);
        float energy_ihc = ihcTowE->at(j);
        h2_IHCal_TowerEnergy->Fill(iEta_ihc, iPhi_ihc, energy_ihc);
    }
    hTotalCaloEEMCal->Fill(*totalCaloEEMCal);
    hTotalCaloEOHCal->Fill(*totalCaloEOHCal);
    hTotalCaloEIHCal->Fill(*totalCaloEIHCal);
    hVtxZ->Fill(*zvertex);

    if (!clusterECore->empty()) {
        float maxClusterECore = *std::max_element(clusterECore->begin(), clusterECore->end());
        hClusterECore->Fill(maxClusterECore);
    }
    if (!clusterPt->empty()) {
        float maxClusterPt = *std::max_element(clusterPt->begin(), clusterPt->end());
        hClusterPt->Fill(maxClusterPt);
    }

    if (!clusterChi2->empty()) {
        float maxClusterChi2 = *std::max_element(clusterChi2->begin(), clusterChi2->end());
        hClusterChi2->Fill(maxClusterChi2);
    }

    if (!ohcChi2->empty()) {
        float maxOhcChi2 = *std::max_element(ohcChi2->begin(), ohcChi2->end());
        h_ohcalChi2->Fill(maxOhcChi2);
    }

    if (!ihcChi2->empty()) {
        float maxIhcChi2 = *std::max_element(ihcChi2->begin(), ihcChi2->end());
        h_ihcalChi2->Fill(maxIhcChi2);
    }
}
std::string formatFloatForFilename(float value) {
    std::ostringstream ss;
    // Increase the precision to handle more decimal places accurately
    ss << std::fixed << std::setprecision(3) << value;
    std::string str = ss.str();
    size_t dotPos = str.find('.');
    if (dotPos != std::string::npos) {
        // Replace '.' with "point"
        str = str.substr(0, dotPos) + "point" + str.substr(dotPos + 1);
    }
    // Remove trailing zeros and 'point' for whole numbers
    if (value == static_cast<int>(value)) {
        size_t pointPos = str.find("point");
        if (pointPos != std::string::npos) {
            str.erase(pointPos);
        }
    } else {
        // Remove trailing zeros for decimal values
        str.erase(str.find_last_not_of('0') + 1, std::string::npos);
    }
    return str;
}

std::vector<int> extractTriggerBits(uint64_t b_gl1_scaledvec, int entry) {
    std::vector<int> trig_bits;
    // Print the 64-bit value as a bitset (binary) for the current entry
    std::bitset<64> bits(b_gl1_scaledvec);
//    std::cout << "Processing entry " << entry << ", gl1_scaledvec (bits): " << bits.to_string() << std::endl;

    
    for (unsigned int bit = 0; bit < 64; bit++) {
        if (((b_gl1_scaledvec >> bit) & 0x1U) == 0x1U) {
            trig_bits.push_back(bit);
//            std::cout << "  Trigger bit set: " << bit << std::endl;
        }
    }
    return trig_bits;
}
bool checkTriggerCondition(const std::vector<int> &trig_bits, int inputBit) {
    for (const int &bit : trig_bits) {
        // Check if the bit is either in the range 24-31 or is exactly 10
        if (bit == inputBit) {
            std::cout << "  Trigger condition met with bit: " << bit << std::endl;
            return true;  // At least one of the desired bits is set
        }
    }
    std::cout << "  No relevant trigger conditions met." << std::endl;
    return false;
}
void createHistogramsForTriggers(TFile* outFile, const std::vector<int>& triggerIndices,
                                 std::map<int, std::map<std::string, TObject*>>& qaHistogramsByTrigger,
                                 std::map<int, std::map<std::string, TH1F*>>& massHistogramsByTrigger,
                                 const std::vector<float>& asymmetry_values,
                                 const std::vector<float>& clus_chi_values,
                                 const std::vector<float>& clus_E_values) {
    // Create directories and histograms for each trigger index
    for (int triggerIndex : triggerIndices) {
        std::cout << "Creating histograms for trigger bit: " << triggerIndex << std::endl;
        
        // Create trigger-specific directories within the output file
        std::string qaDir = "QA/Trigger" + std::to_string(triggerIndex);
        std::string invMassDir = "InvariantMassDistributions/Trigger" + std::to_string(triggerIndex);

        outFile->mkdir(qaDir.c_str());
        outFile->mkdir(invMassDir.c_str());

        // Create QA histograms
        outFile->cd(qaDir.c_str());
        std::map<std::string, TObject*> qaHistograms;
        qaHistograms["h2_EMCal_TowerEtaPhi_2D"] = new TH2F("h2_EMCal_TowerEtaPhi_2D", "EMCal Tower Energy; iEta; iPhi; Energy (GeV)", 96, 0, 96, 256, 0, 256);
        qaHistograms["h2_OHCal_TowerEnergy"] = new TH2F("h2_OHCal_TowerEnergy", "HCal Tower Energy; iEta; iPhi; Energy (GeV)", 24, 0, 24, 64, 0, 64);
        qaHistograms["h2_IHCal_TowerEnergy"] = new TH2F("h2_IHCal_TowerEnergy", "HCal Tower Energy; iEta; iPhi; Energy (GeV)", 24, 0, 24, 64, 0, 64);
        qaHistograms["hTotalCaloEEMCal"] = new TH1F("hTotalCaloEEMCal", "Total EMCal Energy; Energy (GeV)", 100, 0, 500);
        qaHistograms["hTotalCaloEOHCal"] = new TH1F("hTotalCaloEOHCal", "Total OHCal Energy; Energy (GeV)", 100, 0, 500);
        qaHistograms["hTotalCaloEIHCal"] = new TH1F("hTotalCaloEIHCal", "Total IHCal Energy; Energy (GeV)", 100, 0, 500);
        qaHistograms["hClusterChi2"] = new TH1F("hClusterChi2", "Cluster Chi2; Chi2", 100, 0, 100);
        qaHistograms["h_ohcalChi2"] = new TH1F("h_ohcalChi2", "Cluster Chi2; Chi2", 100, 0, 100);
        qaHistograms["h_ihcalChi2"] = new TH1F("h_ihcalChi2", "Cluster Chi2; Chi2", 100, 0, 100);
        qaHistograms["hClusterECore"] = new TH1F("hClusterECore", "Cluster ECore; Cluster ECore [GeV]", 100, 0, 100);
        qaHistograms["hClusterPt"] = new TH1F("hClusterPt", "Cluster pT; Cluster pT [GeV]", 100, 0, 100);
        qaHistograms["hVtxZ"] = new TH1F("hVtxZ", "Z-vertex Distribution; z [cm]", 100, -70, 70);
        qaHistograms["h8by8TowerEnergySum"] = new TH1F("h8by8TowerEnergySum", "Max 8x8 Tower Energy Sum; Energy (GeV); Events", 100, 0, 500);

        qaHistogramsByTrigger[triggerIndex] = qaHistograms;

        // Create invariant mass histograms
        outFile->cd(invMassDir.c_str());
        std::map<std::string, TH1F*> massHistograms;

        for (float maxAsym : asymmetry_values) {
            for (float maxChi2 : clus_chi_values) {
                for (float minClusE : clus_E_values) {
                    std::string histName = "invMass_E" + formatFloatForFilename(minClusE) +
                                           "_Chi" + formatFloatForFilename(maxChi2) +
                                           "_Asym" + formatFloatForFilename(maxAsym);
                    massHistograms[histName] = new TH1F(histName.c_str(), histName.c_str(), 80, 0, 1.0);
                    massHistograms[histName]->SetTitle(";M_{#gamma#gamma};");
                }
            }
        }
        massHistogramsByTrigger[triggerIndex] = massHistograms;
    }
}
void processEvents(TTree* T, const std::vector<int>& triggerIndices,
                   const std::vector<float>& asymmetry_values,
                   const std::vector<float>& clus_chi_values,
                   const std::vector<float>& clus_E_values,
                   std::map<int, std::map<std::string, TObject*>>& qaHistogramsByTrigger,
                   std::map<int, std::map<std::string, TH1F*>>& massHistogramsByTrigger,
                   ULong64_t& b_gl1_scaledvec, float& zvertex,
                   std::vector<float>* clusterE, std::vector<float>* clusterECore,
                   std::vector<float>* clusterTowMaxE, std::vector<float>* clusterPhi,
                   std::vector<float>* clusterEta, std::vector<float>* clusterPt,
                   std::vector<float>* clusterChi2, std::vector<float>* clusterNtow,
                   std::vector<float>* emcTowE, std::vector<int>* emcTowiEta,
                   std::vector<int>* emcTowiPhi, std::vector<float>* ohcTowE,
                   std::vector<float>* ohcTowiEta, std::vector<float>* ohcTowiPhi,
                   std::vector<float>* ihcTowE, std::vector<float>* ihcTowiEta,
                   std::vector<float>* ihcTowiPhi, std::vector<float>* ohcChi2,
                   std::vector<float>* ihcChi2, float& totalCaloEEMCal,
                   float& totalCaloEOHCal, float& totalCaloEIHCal) {

    // Loop over all events once
    for (int i = 0; i < T->GetEntries(); i++) {
        T->GetEntry(i);

        // Extract trigger bits for this event
        std::vector<int> trig_bits = extractTriggerBits(b_gl1_scaledvec, i);

        // Apply event-level cuts (e.g., zvertex)
        if (fabs(zvertex) >= 30) {
            continue;
        }

        // For each trigger bit, fill corresponding histograms if it is active
        for (int triggerIndex : triggerIndices) {
            if (!checkTriggerCondition(trig_bits, triggerIndex)) {
                continue;  // Skip if this trigger bit is not set
            }
            std::cout << "Filling histograms for trigger bit: " << triggerIndex << " in event " << i << std::endl;
            // Fill QA histograms for this trigger
            auto& qaHistograms = qaHistogramsByTrigger[triggerIndex];
            fillQAHistograms(
                (TH2F*)qaHistograms["h2_EMCal_TowerEtaPhi_2D"],
                (TH2F*)qaHistograms["h2_OHCal_TowerEnergy"],
                (TH2F*)qaHistograms["h2_IHCal_TowerEnergy"],
                (TH1F*)qaHistograms["hTotalCaloEEMCal"],
                (TH1F*)qaHistograms["hTotalCaloEOHCal"],
                (TH1F*)qaHistograms["hTotalCaloEIHCal"],
                (TH1F*)qaHistograms["h8by8TowerEnergySum"],
                (TH1F*)qaHistograms["hClusterECore"],
                (TH1F*)qaHistograms["hClusterPt"],
                (TH1F*)qaHistograms["hClusterChi2"],
                (TH1F*)qaHistograms["h_ohcalChi2"],
                (TH1F*)qaHistograms["h_ihcalChi2"],
                (TH1F*)qaHistograms["hVtxZ"],
                emcTowE, emcTowiEta, emcTowiPhi, ohcTowE, ohcTowiEta, ohcTowiPhi,
                ihcTowE, ihcTowiEta, ihcTowiPhi, clusterECore, clusterPt, clusterChi2,
                ohcChi2, ihcChi2, &totalCaloEEMCal, &totalCaloEOHCal, &totalCaloEIHCal, &zvertex
            );

            // Cache cluster properties to avoid repeated access to the vector
            size_t nClusters = clusterE->size();
            std::vector<float> cachedPt(nClusters), cachedE(nClusters), cachedChi2(nClusters);
            for (size_t clus = 0; clus < nClusters; ++clus) {
                cachedPt[clus] = clusterPt->at(clus);
                cachedE[clus] = clusterE->at(clus);
                cachedChi2[clus] = clusterChi2->at(clus);
            }

            // Fill invariant mass histograms for this trigger
            auto& massHistograms = massHistogramsByTrigger[triggerIndex];
            for (size_t clus1 = 0; clus1 < nClusters; ++clus1) {
                for (size_t clus2 = clus1 + 1; clus2 < nClusters; ++clus2) {
                    float pt1 = cachedPt[clus1], pt2 = cachedPt[clus2];
                    float E1 = cachedE[clus1], E2 = cachedE[clus2];
                    float chi1 = cachedChi2[clus1], chi2 = cachedChi2[clus2];

                    if (pt1 < 2 || pt1 >= 20 || pt2 < 2 || pt2 >= 20) {
                        continue;
                    }

                    TLorentzVector photon1, photon2;
                    photon1.SetPtEtaPhiE(pt1, clusterEta->at(clus1), clusterPhi->at(clus1), E1);
                    photon2.SetPtEtaPhiE(pt2, clusterEta->at(clus2), clusterPhi->at(clus2), E2);
                    TLorentzVector meson = photon1 + photon2;
                    float mesonMass = meson.M();
                    float asym = fabs(E1 - E2) / (E1 + E2);

                    // Apply cuts and fill histograms
                    for (float maxAsym : asymmetry_values) {
                        if (asym >= maxAsym) continue;
                        for (float maxChi2 : clus_chi_values) {
                            if (chi1 >= maxChi2 || chi2 >= maxChi2) continue;
                            for (float minClusE : clus_E_values) {
                                if (E1 < minClusE || E2 < minClusE) continue;
                                std::string histName = "invMass_E" + formatFloatForFilename(minClusE) +
                                                       "_Chi" + formatFloatForFilename(maxChi2) +
                                                       "_Asym" + formatFloatForFilename(maxAsym);
                                massHistograms[histName]->Fill(mesonMass);
                            }
                        }
                    }
                }
            }
        }
    }
}

void process_nTuples(int nEvents = 0, const std::string &fin = "/Users/patsfan753/Desktop/DirectPhotonAna/CaloTreeGen_DST_CALO_run2pp_ana430_2024p007-00046649-00084.root", const std::string &fout = "output.root") {
    // Open input file
    TFile* f = new TFile(fin.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << fin << std::endl;
        logCorruptedFile(fin);  // Log the corrupted file
        return;
    }
    
    TTree* T = (TTree*)f->Get("T");
    if (!T) {
        std::cerr << "Error: Cannot find tree 'T' in file " << fin << std::endl;
        logCorruptedFile(fin);  // Log the corrupted file
        f->Close();
        delete f;
        return;
    }

    // Declare variables for branch addresses
    ULong64_t b_gl1_scaledvec = 0;
    float zvertex = 0;
    std::vector<float> *clusterE{0}, *clusterECore{0}, *clusterTowMaxE{0}, *clusterPhi{0}, *clusterEta{0}, *clusterPt{0}, *clusterChi2{0}, *clusterNtow{0};
    float totalCaloEEMCal = 0, totalCaloEOHCal = 0, totalCaloEIHCal = 0;
    std::vector<float> *emcTowE{0}, *ohcTowE{0}, *ihcTowE{0}, *ohcChi2{0}, *ihcChi2{0};
    std::vector<int> *emcTowiEta{0}, *emcTowiPhi{0};
    std::vector<float> *ohcTowiEta{0}, *ohcTowiPhi{0}, *ihcTowiEta{0}, *ihcTowiPhi{0};

    try {
        // Set branch addresses using the utility function carefully in case corruption don't override memory
        setTreeBranches(T, fin, &b_gl1_scaledvec, &zvertex, &clusterE, &clusterECore, &clusterTowMaxE, &clusterPhi,
                        &clusterEta, &clusterPt, &clusterChi2, &clusterNtow, &totalCaloEEMCal, &emcTowE,
                        &emcTowiEta, &emcTowiPhi, &totalCaloEOHCal, &ohcTowE, &ohcTowiEta, &ohcTowiPhi, &ohcChi2,
                        &totalCaloEIHCal, &ihcTowE, &ihcTowiEta, &ihcTowiPhi, &ihcChi2);
    
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        f->Close();
        delete f;
        return;
    }

    // Create output file
    TFile* outFile = new TFile(fout.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error: Cannot create output file " << fout << std::endl;
        logCorruptedFile(fout);  // Log the corrupted output file
        f->Close();
        delete f;
        return;
    }

    // Define the relevant trigger indices , 25, 26, 27, 28, 29, 30, 31
    std::vector<int> triggerIndices = {10, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};  // Add/remove as needed

    // Asymmetry, chi2, and energy cut parameters
    std::vector<float> asymmetry_values = {0.5, 0.7};
    std::vector<float> clus_chi_values = {3, 4};
    std::vector<float> clus_E_values = {0.5, 0.75, 1};

    // Maps to hold histograms for each trigger bit
    std::map<int, std::map<std::string, TObject*>> qaHistogramsByTrigger;
    std::map<int, std::map<std::string, TH1F*>> massHistogramsByTrigger;

    createHistogramsForTriggers(outFile, triggerIndices, qaHistogramsByTrigger, massHistogramsByTrigger, asymmetry_values, clus_chi_values, clus_E_values);

    // Call the function to process all events
    processEvents(T, triggerIndices, asymmetry_values, clus_chi_values, clus_E_values, qaHistogramsByTrigger, massHistogramsByTrigger,
                  b_gl1_scaledvec, zvertex, clusterE, clusterECore, clusterTowMaxE, clusterPhi, clusterEta, clusterPt, clusterChi2, clusterNtow,
                  emcTowE, emcTowiEta, emcTowiPhi, ohcTowE, ohcTowiEta, ohcTowiPhi, ihcTowE, ihcTowiEta, ihcTowiPhi, ohcChi2, ihcChi2, totalCaloEEMCal, totalCaloEOHCal, totalCaloEIHCal);


    // Write all histograms to file and clean up
    for (const auto& [triggerIndex, qaHistograms] : qaHistogramsByTrigger) {
        outFile->cd(("QA/Trigger" + std::to_string(triggerIndex)).c_str());
        for (const auto& [name, hist] : qaHistograms) {
            hist->Write();
            delete hist;
        }
    }

    for (const auto& [triggerIndex, massHistograms] : massHistogramsByTrigger) {
        outFile->cd(("InvariantMassDistributions/Trigger" + std::to_string(triggerIndex)).c_str());
        for (const auto& [name, hist] : massHistograms) {
            hist->Write();
            delete hist;
        }
    }

    outFile->Close();
    delete f;
}
