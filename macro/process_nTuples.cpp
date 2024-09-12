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
                      TH1F* hTotalCaloEEMCal, TH1F* hTotalCaloEOHCal, TH1F* hTotalCaloEIHCal,
                      TH1F* hClusterECore, TH1F* hClusterPt, TH1F* hClusterChi2, TH1F* h_ohcalChi2, TH1F* h_ihcalChi2, TH1F* hVtxZ,
                      const std::vector<float>* emcTowE, const std::vector<int>* emcTowiEta, const std::vector<int>* emcTowiPhi,
                      const std::vector<float>* ohcTowE, const std::vector<float>* ohcTowiEta, const std::vector<float>* ohcTowiPhi,
                      const std::vector<float>* ihcTowE, const std::vector<float>* ihcTowiEta, const std::vector<float>* ihcTowiPhi,
                      const std::vector<float>* clusterECore, const std::vector<float>* clusterChi2, const std::vector<float>* clusterPt,
                      const std::vector<float>* ohcChi2, const std::vector<float>* ihcChi2, const float* totalCaloEEMCal,
                      const float* totalCaloEOHCal, const float* totalCaloEIHCal, const float* zvertex) {


    for (size_t j = 0; j < emcTowE->size(); j++) {
        int iEta_emc = emcTowiEta->at(j);
        int iPhi_emc = emcTowiPhi->at(j);
        float energy_emc = emcTowE->at(j);
        h2_EMCal_TowerEtaPhi_2D->Fill(iEta_emc, iPhi_emc, energy_emc);
    }


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
bool checkTriggerCondition(const std::vector<int> &trig_bits) {
    for (const int &bit : trig_bits) {
        // Check if the bit is either in the range 24-31 or is exactly 10
        if ((bit >= 24 && bit <= 31) || bit == 10 || bit == 12 || bit == 13) {
//            std::cout << "  Trigger condition met with bit: " << bit << std::endl;
            return true;  // At least one of the desired bits is set
        }
    }
//    std::cout << "  No relevant trigger conditions met." << std::endl;
    return false;
}


void compareHistograms(TH1F* h1, TH1F* h2, const std::string& label) {
    bool identical = true;
    // Check if bin contents are identical
    if (h1->GetNbinsX() != h2->GetNbinsX()) {
        identical = false;
    } else {
        for (int i = 1; i <= h1->GetNbinsX(); ++i) {
            if (h1->GetBinContent(i) != h2->GetBinContent(i)) {
                identical = false;
                break;
            }
        }
    }

    if (!identical) {
        std::cout << "\033[1;31m" << "ERROR: The histograms from the two strategies (" << label << ") are NOT identical!" << "\033[0m" << std::endl;
    } else {
        std::cout << "The histograms from the two strategies (" << label << ") are identical." << std::endl;
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
 
  // Create the two directories for the output
  outFile->mkdir("QA");
  outFile->mkdir("InvariantMassDistributions");

    
  outFile->cd("QA");
  // Create 2D histogram for EMC tower energy distribution
  TH2F* h2_EMCal_TowerEtaPhi_2D = new TH2F("h2_EMCal_TowerEtaPhi_2D", "EMCal Tower Energy; iEta; iPhi; Energy (GeV)", 96, 0, 96, 256, 0, 256);  // 96 bins in iEta, 256 bins in iPhi
  TH2F* h2_OHCal_TowerEnergy = new TH2F("h2_OHCal_TowerEtaPhi_2D", "HCal Tower Energy; iEta; iPhi; Energy (GeV)", 24, 0, 24, 64, 0, 64);  // 24 bins in iEta, 64 bins in iPhi
  TH2F* h2_IHCal_TowerEnergy = new TH2F("h2_IHCal_TowerEtaPhi_2D", "HCal Tower Energy; iEta; iPhi; Energy (GeV)", 24, 0, 24, 64, 0, 64);  // 24 bins in iEta, 64 bins in iPhi
  TH1F* hTotalCaloEEMCal = new TH1F("hTotalCaloEEMCal", "Total EMCal Energy; Energy (GeV)", 100, 0, 500);
  TH1F* hTotalCaloEOHCal = new TH1F("hTotalCaloEOHCal", "Total OHCal Energy; Energy (GeV)", 100, 0, 500);
  TH1F* hTotalCaloEIHCal = new TH1F("hTotalCaloEIHCal", "Total IHCal Energy; Energy (GeV)", 100, 0, 500);
  TH1F* hClusterChi2 = new TH1F("hClusterChi2", "Cluster Chi2; Chi2", 100, 0, 100);
  TH1F* h_ohcalChi2 = new TH1F("h_ohcalChi2", "Cluster Chi2; Chi2", 100, 0, 100);
  TH1F* h_ihcalChi2 = new TH1F("h_ihcalChi2", "Cluster Chi2; Chi2", 100, 0, 100);
  TH1F* hClusterECore = new TH1F("hClusterECore", "Cluster ECore; Cluster ECore [GeV]", 100, 0, 100);
  TH1F* hClusterPt = new TH1F("hClusterPt", "Cluster pT; Cluster pT [GeV]", 100, 0, 100);
  TH1F* hVtxZ = new TH1F("hVtxZ", "Z-vertex Distribution; z [cm]", 100, -70, 70);

  std::map<std::string, TObject*> qaHistograms;

  qaHistograms["h2_EMCal_TowerEtaPhi_2D"] = h2_EMCal_TowerEtaPhi_2D;
  qaHistograms["h2_OHCal_TowerEnergy"] = h2_OHCal_TowerEnergy;
  qaHistograms["h2_IHCal_TowerEnergy"] = h2_IHCal_TowerEnergy;
  qaHistograms["hTotalCaloEEMCal"] = hTotalCaloEEMCal;
  qaHistograms["hTotalCaloEOHCal"] = hTotalCaloEOHCal;
  qaHistograms["hTotalCaloEIHCal"] = hTotalCaloEIHCal;
  qaHistograms["hClusterChi2"] = hClusterChi2;
  qaHistograms["h_ohcalChi2"] = h_ohcalChi2;
  qaHistograms["h_ihcalChi2"] = h_ihcalChi2;
  qaHistograms["hClusterECore"] = hClusterECore;
  qaHistograms["hClusterPt"] = hClusterPt;
  qaHistograms["hVtxZ"] = hVtxZ;


  std::map<std::string, TH1F*> massHistograms;
  std::vector<float> asymmetry_values = {0.5, 0.7};
  std::vector<float> clus_chi_values = {3, 4};
  std::vector<float> clus_E_values = {0.5, 0.75, 1};
    
    // Create histograms based on cut parameters
    for (float maxAsym : asymmetry_values) {
        for (float maxChi2 : clus_chi_values) {
            for (float minClusE : clus_E_values) {
                // Add cuts to the histogram name
                std::string histName = "invMass_E" + formatFloatForFilename(minClusE) +
                                        "_Chi" + formatFloatForFilename(maxChi2) +
                                        "_Asym" + formatFloatForFilename(maxAsym);
                TH1F* mass = new TH1F(histName.c_str(), histName.c_str(), 80, 0, 1.0);
                mass->SetTitle(";M_{#gamma#gamma};");
                massHistograms[histName] = mass;  // Store histogram in the map
            }
        }
    }
    /*
     Loop over all events
    */
 
  for (int i = 0; i < T -> GetEntries(); i++) {
      T->GetEntry(i);
      
      std::cout << "Processing event: " << i << std::endl;
      //create QA histograms
      fillQAHistograms(h2_EMCal_TowerEtaPhi_2D, h2_OHCal_TowerEnergy, h2_IHCal_TowerEnergy,
                         hTotalCaloEEMCal, hTotalCaloEOHCal, hTotalCaloEIHCal,
                         hClusterECore, hClusterPt, hClusterChi2, h_ohcalChi2, h_ihcalChi2, hVtxZ,
                         emcTowE, emcTowiEta, emcTowiPhi,
                         ohcTowE, ohcTowiEta, ohcTowiPhi,
                         ihcTowE, ihcTowiEta, ihcTowiPhi,
                         clusterECore, clusterPt, clusterChi2, ohcChi2, ihcChi2,
                         &totalCaloEEMCal, &totalCaloEOHCal, &totalCaloEIHCal, &zvertex);
      
      // Apply event-level cuts (e.g., zvertex, trigger condition)
      if (fabs(zvertex) >= 30) {
          continue;
      }
      std::vector<int> trig_bits = extractTriggerBits(b_gl1_scaledvec, i);

      // Check if the event passes the trigger condition
      if (!checkTriggerCondition(trig_bits)) {
          continue;  // Skip this event if no relevant trigger bits are set
      }
      // Cache cluster properties to avoid repeated access to the vector
      size_t nClusters = clusterE->size();
      std::vector<float> cachedPt(nClusters), cachedE(nClusters), cachedChi2(nClusters);
      for (size_t clus = 0; clus < nClusters; ++clus) {
          cachedPt[clus] = clusterPt->at(clus);
          cachedE[clus] = clusterE->at(clus);
          cachedChi2[clus] = clusterChi2->at(clus);
      }

      // Loop over unique cluster pairs
      for (size_t clus1 = 0; clus1 < nClusters; ++clus1) {
          for (size_t clus2 = clus1 + 1; clus2 < nClusters; ++clus2) {
              float pt1 = cachedPt[clus1], pt2 = cachedPt[clus2];
              float E1 = cachedE[clus1], E2 = cachedE[clus2];
              float chi1 = cachedChi2[clus1], chi2 = cachedChi2[clus2];

              // Apply the pT cut (2-20 GeV)
              if (pt1 < 1 || pt1 >= 50 || pt2 < 1 || pt2 >= 50) {
                  continue;
              }

              TLorentzVector photon1, photon2;
              photon1.SetPtEtaPhiE(pt1, clusterEta->at(clus1), clusterPhi->at(clus1), E1);
              photon2.SetPtEtaPhiE(pt2, clusterEta->at(clus2), clusterPhi->at(clus2), E2);
              TLorentzVector meson = photon1 + photon2;
              float mesonMass = meson.M();
              float asym = fabs(E1 - E2) / (E1 + E2);


              // Now apply all combinations of cuts
              for (float maxAsym : asymmetry_values) {
                  if (asym >= maxAsym) {
                      continue;
                  }

                  for (float maxChi2 : clus_chi_values) {
                      if (chi1 >= maxChi2 || chi2 >= maxChi2) {
                          continue;
                      }

                      for (float minClusE : clus_E_values) {
                          if (E1 < minClusE || E2 < minClusE) {
                              continue;
                          }

                          // Build the histogram name based on cuts
                          std::string histName = "invMass_E" + formatFloatForFilename(minClusE) +
                                                "_Chi" + formatFloatForFilename(maxChi2) +
                                                "_Asym" + formatFloatForFilename(maxAsym);

                          // Fill the corresponding histogram with the meson mass
                          massHistograms[histName]->Fill(mesonMass);
                      }
                  }
              }
          }
      }
  }
  // After the event loop, write histograms to file
  outFile->cd("InvariantMassDistributions");
  for (auto& pair : massHistograms) {
      pair.second->Write();
      delete pair.second;  // Clean up
  }
    
  // Write QA histograms
  outFile->cd("QA");
  for (auto &qaHist : qaHistograms) {
       if (TH1F* hist1D = dynamic_cast<TH1F*>(qaHist.second)) {
           hist1D->Write();
       } else if (TH2F* hist2D = dynamic_cast<TH2F*>(qaHist.second)) {
           hist2D->Write();
       }
       delete qaHist.second;  // Clean up after writing
   }

   outFile->Close();  // Close the output file
   delete f;  // Clean up the input file
}
