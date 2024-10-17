#include "caloTreeGen.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

//Fun4All
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <ffaobjects/EventHeader.h>

//ROOT stuff
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterv1.h>
#include <calobase/RawClusterDefs.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeom.h>

//Tower stuff
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoDefs.h>
#include <chrono>  // Include for timing

#include "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/clusterIsoCopy_src/ClusterIso.h"
//GL1 Information
#include <ffarawobjects/Gl1Packet.h>

//for cluster vertex correction
#include <CLHEP/Geometry/Point3D.h>

//for the vertex
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#define ANSI_COLOR_RED_BOLD "\033[1;31m"
#define ANSI_COLOR_BLUE_BOLD "\033[1;34m"
#define ANSI_COLOR_GREEN_BOLD "\033[1;32m"
#define ANSI_COLOR_RESET "\033[0m"


struct TowerData {
    unsigned int ieta;
    unsigned int iphi;
    double energy;
    int time;
    float chi2;
    float pedestal;
    short good;
    bool isAcceptable;
};

/** \Brief Function to get correct tower eta -- FROM ClusterIso.cc
 *
 * Each calorimeter tower's eta is calculated using the vertex (0,0,0)
 * which is incorrect in many collisions. This function
 * uses geometry to find a given tower's eta using the correct vertex.
 */
double caloTreeGen::getTowerEta(RawTowerGeom *tower_geom, double vx, double vy, double vz) {
  float r;
  if (vx == 0 && vy == 0 && vz == 0) {
    r = tower_geom->get_eta();
  }
    
  else {
    double radius = sqrt((tower_geom->get_center_x() - vx) * (tower_geom->get_center_x() - vx) + (tower_geom->get_center_y() - vy) * (tower_geom->get_center_y() - vy));
    double theta = atan2(radius, tower_geom->get_center_z() - vz);
    r = -log(tan(theta / 2.));
  }
  return r;
}
//____________________________________________________________________________..
caloTreeGen::caloTreeGen(const std::string &name):
SubsysReco("CaloTreeGen")
  ,Outfile(name)
{
  std::cout << "caloTreeGen::caloTreeGen(const std::string &name) Calling ctor" << std::endl;
}



//____________________________________________________________________________..
caloTreeGen::~caloTreeGen() {
  std::cout << "caloTreeGen::~caloTreeGen() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int caloTreeGen::Init(PHCompositeNode *topNode) {
    std::cout << ANSI_COLOR_BLUE_BOLD << "Initializing caloTreeGen..." << ANSI_COLOR_RESET << std::endl;
    out = new TFile(Outfile.c_str(),"RECREATE");
    
    if (verbose) {
        std::cout << ANSI_COLOR_BLUE_BOLD << "Output file created: " << Outfile << ANSI_COLOR_RESET << std::endl;
    }
    
    for (int triggerIndex : triggerIndices) {
        // Create trigger-specific directories within the output file
        std::string qaDir = "QA/Trigger" + std::to_string(triggerIndex);
        std::string invMassDir = "InvariantMassDistributions/Trigger" + std::to_string(triggerIndex);

        std::cout << ANSI_COLOR_BLUE_BOLD << "Creating directories: " << qaDir << ", " << invMassDir << ANSI_COLOR_RESET << std::endl;
        
        out->mkdir(qaDir.c_str());
        out->mkdir(invMassDir.c_str());

        // Create QA histograms
        if (verbose) {
            std::cout << ANSI_COLOR_BLUE_BOLD << "Switching to QA directory: " << qaDir << ANSI_COLOR_RESET << std::endl;
        }
        
        out->cd(qaDir.c_str());
        std::map<std::string, TObject*>& qaHistograms = qaHistogramsByTrigger[triggerIndex];

        
        // Helper functions for creating histograms with logging
        auto createHistogram = [&](const std::string& name, const std::string& title, int bins, double xMin, double xMax) {
            if (verbose) {
                std::cout << ANSI_COLOR_RED_BOLD << "Creating histogram: " << name << " - " << title << ANSI_COLOR_RESET << std::endl;
            }
            
            return new TH1F(name.c_str(), title.c_str(), bins, xMin, xMax);
        };

        auto create2DHistogram = [&](const std::string& name, const std::string& title, int xBins, double xMin, double xMax, int yBins, double yMin, double yMax) {
            if (verbose) {
                std::cout << ANSI_COLOR_RED_BOLD << "Creating 2D histogram: " << name << " - " << title << ANSI_COLOR_RESET << std::endl;
            }
            return new TH2F(name.c_str(), title.c_str(), xBins, xMin, xMax, yBins, yMin, yMax);
        };

        /*
         EMCal QA
         */
        qaHistograms["h2_EMCal_TowerEtaPhi_Weighted_2D_" + std::to_string(triggerIndex)] = create2DHistogram("h2_EMCal_TowerEtaPhi_2D_" + std::to_string(triggerIndex), "EMCal Weighted Tower Energy; #eta; #phi; Energy [GeV]", 96, 0, 96, 256, 0, 256);
        qaHistograms["hTotalCaloEEMCal_" + std::to_string(triggerIndex)] = createHistogram("hTotalCaloEEMCal_" + std::to_string(triggerIndex), "Total EMCal Energy; Energy (GeV)", 100, -50, 100);
        qaHistograms["hTotalCalo_Negative_EEMCal_" + std::to_string(triggerIndex)] = createHistogram("hTotalCalo_Negative_EEMCal_" + std::to_string(triggerIndex), "Total EMCal Energy; Energy (GeV)", 100, 0, 100);
        qaHistograms["h_emcalChi2_" + std::to_string(triggerIndex)] = createHistogram("h_emcalChi2_" + std::to_string(triggerIndex), "Cluster Chi2; Chi2", 100, 0, 100);
        qaHistograms["h_maxTowerEnergy_" + std::to_string(triggerIndex)] = createHistogram("h_maxTowerEnergy_" + std::to_string(triggerIndex), "Max Tower Energy [GeV]; Energy [GeV]", 100, 0, 100);
        

        /*
         Cluster Distributions EMCal
         */
        qaHistograms["hClusterChi2_" + std::to_string(triggerIndex)] = createHistogram("hClusterChi2_" + std::to_string(triggerIndex), "Cluster Chi2; Chi2", 100, 0, 100);
        qaHistograms["hCluster_maxECore_" + std::to_string(triggerIndex)] = createHistogram("hCluster_maxECore_" + std::to_string(triggerIndex), "Max Cluster ECore; Cluster ECore [GeV]", 40, 0, 20);
        qaHistograms["hClusterPt_" + std::to_string(triggerIndex)] = createHistogram("hClusterPt_" + std::to_string(triggerIndex), "Cluster pT; Cluster pT [GeV]", 100, 0, 100);
        qaHistograms["hVtxZ_" + std::to_string(triggerIndex)] = createHistogram("hVtxZ_" + std::to_string(triggerIndex), "Z-vertex Distribution; z [cm]", 100, -70, 70);
        qaHistograms["h_ClusterEnergy_" + std::to_string(triggerIndex)] = createHistogram("h_ClusterEnergy_" + std::to_string(triggerIndex), "Cluster Energy [GeV]; Energy [GeV]", 100, 0, 100);
        qaHistograms["h_ET_" + std::to_string(triggerIndex)] = createHistogram("h_ET_" + std::to_string(triggerIndex), "Cluster Transverse Energy [GeV]; Energy [GeV]", 100, 0, 100);
        qaHistograms["h2_cluster_iso_Ecore_" + std::to_string(triggerIndex)] = create2DHistogram("h2_cluster_iso_Ecore_" + std::to_string(triggerIndex), "Cluster Isolation Energy vs Cluster Energy (getEt function);Cluster ECore [GeV];E_{T}^{iso} [GeV]", 100, 0, 20, 100, -10, 30);
        qaHistograms["h1_isoEt_" + std::to_string(triggerIndex)] = createHistogram("h1_isoEt_" + std::to_string(triggerIndex), "Isolation Energy Distribution;E_{T}^{iso} [GeV];Counts", 100, -10, 30);
        
        /*
         HCal QA
         */
        //inner
        qaHistograms["h2_IHCal_TowerEnergy_" + std::to_string(triggerIndex)] = create2DHistogram("h2_IHCal_TowerEnergy_" + std::to_string(triggerIndex), "HCal Tower Energy; #eta; #phi; Energy [GeV]", 24, 0, 24, 64, 0, 64);
        qaHistograms["hTotalCaloEIHCal_" + std::to_string(triggerIndex)] = createHistogram("hTotalCaloEIHCal_" + std::to_string(triggerIndex), "Total IHCal Energy; Energy (GeV)", 100, 0, 500);
        qaHistograms["h_ihcalChi2_" + std::to_string(triggerIndex)] = createHistogram("h_ihcalChi2_" + std::to_string(triggerIndex), "Cluster Chi2; Chi2", 100, 0, 100);
        
        //outer
        qaHistograms["h2_OHCal_TowerEnergy_" + std::to_string(triggerIndex)] = create2DHistogram("h2_OHCal_TowerEnergy_" + std::to_string(triggerIndex), "HCal Tower Energy; #eta; #phi; Energy [GeV]", 24, 0, 24, 64, 0, 64);
        qaHistograms["hTotalCaloEOHCal_" + std::to_string(triggerIndex)] = createHistogram("hTotalCaloEOHCal_" + std::to_string(triggerIndex), "Total OHCal Energy; Energy (GeV)", 100, 0, 500);
        qaHistograms["h_ohcalChi2_" + std::to_string(triggerIndex)] = createHistogram("h_ohcalChi2_" + std::to_string(triggerIndex), "Cluster Chi2; Chi2", 100, 0, 100);
        
        /*
         Trigger QA Distributions
         */
        //use for turn-on curves
        qaHistograms["h8by8TowerEnergySum_" + std::to_string(triggerIndex)] = createHistogram("h8by8TowerEnergySum_" + std::to_string(triggerIndex), "Max 8x8 Tower Energy Sum; Energy [GeV]; Events", 40, 0, 10); //photon triggers
        qaHistograms["h_jet_energy_" + std::to_string(triggerIndex)] = createHistogram("h_jet_energy_" + std::to_string(triggerIndex), "Maximum 0.8x0.8 Energy Sum (EMCAL + HCAL) [GeV]; Events", 50, 0, 50); //jet triggers
        
        //other possible trigger QA
        qaHistograms["h_hcal_energy_" + std::to_string(triggerIndex)] = createHistogram("h_hcal_energy_" + std::to_string(triggerIndex), "Max HCal Tower Energy Sums; Energy [GeV]; Events", 40, 0, 20);
        qaHistograms["h_jet_emcal_energy_" + std::to_string(triggerIndex)] = createHistogram("h_jet_emcal_energy_" + std::to_string(triggerIndex), "hJetEmcalEnergy [GeV]", 40, 0, 20);
        qaHistograms["h_jet_hcalin_energy_" + std::to_string(triggerIndex)] = createHistogram("h_jet_hcalin_energy_" + std::to_string(triggerIndex), "hJetEmcalEnergy [GeV]", 40, 0, 20);
        qaHistograms["h_jet_hcalout_energy_" + std::to_string(triggerIndex)] = createHistogram("h_jet_hcalout_energy_" + std::to_string(triggerIndex), "hJetEmcalEnergy [GeV]", 40, 0, 20);
        // Initialize the trigger count histogram
        qaHistograms["hTriggerCount_" + std::to_string(triggerIndex)] = new TH1F(
            ("hTriggerCount_" + std::to_string(triggerIndex)).c_str(),
            ("Trigger Count for Index " + std::to_string(triggerIndex) + "; Count; Entries").c_str(),
            1, 0, 1 // Single bin to count occurrences
        );

        qaHistogramsByTrigger[triggerIndex] = qaHistograms;
        
        if (verbose) {
            std::cout << ANSI_COLOR_BLUE_BOLD << "Switching to Invariant Mass directory: " << invMassDir << ANSI_COLOR_RESET << std::endl;
        }
        
        // Create invariant mass histograms
        out->cd(invMassDir.c_str());
        std::map<std::string, TObject*> massHistograms;

        for (float maxAsym : asymmetry_values) {
            for (float maxChi2 : clus_chi_values) {
                for (float minClusE : clus_Ecore_values) {
                    std::string histName = "invMass_E" + formatFloatForFilename(minClusE) +
                                           "_Chi" + formatFloatForFilename(maxChi2) +
                                           "_Asym" + formatFloatForFilename(maxAsym) +
                                           "_" + std::to_string(triggerIndex);
                    if (verbose) {
                        std::cout << ANSI_COLOR_RED_BOLD << "Creating invariant mass histogram: " << histName << ANSI_COLOR_RESET << std::endl;
                    }
                    
                    TH1F* hist = new TH1F(histName.c_str(), histName.c_str(), 80, 0, 1.0);
                    hist->SetTitle(";M_{#gamma#gamma};");
                    massHistograms[histName] = hist; // Store it as TObject*
                }
            }
        }
        massHistogramsByTrigger[triggerIndex] = massHistograms;
    }

    //so that the histos actually get written out
    Fun4AllServer *se = Fun4AllServer::instance();
    se -> Print("NODETREE");
    std::cout << "caloTreeGen::Init(PHCompositeNode *topNode) Initializing" << std::endl;

    return Fun4AllReturnCodes::EVENT_OK;
}

void caloTreeGen::collectTowerData(TowerInfoContainer* towerContainer,
                                   std::vector<TowerData>& towerDataList) {
    if (!towerContainer) {
        std::cout << ANSI_COLOR_RED_BOLD << "Error: Tower container is null." << ANSI_COLOR_RESET << std::endl;
        return;
    }

    unsigned int tower_range = towerContainer->size();
    towerDataList.reserve(tower_range);
    if (verbose) {
        std::cout << ANSI_COLOR_BLUE_BOLD << "Collecting tower data from container, expecting " << tower_range << " towers..." << ANSI_COLOR_RESET << std::endl;
    }
    
    unsigned int validTowerCount = 0;
    for (unsigned int iter = 0; iter < tower_range; ++iter) {
        TowerInfo* tower = towerContainer->get_tower_at_channel(iter);
        if (!tower) {
            std::cout << "Warning: Null tower at index " << iter << std::endl;
            continue;
        }

        unsigned int towerkey = towerContainer->encode_key(iter);
        unsigned int ieta = TowerInfoDefs::getCaloTowerEtaBin(towerkey);
        unsigned int iphi = TowerInfoDefs::getCaloTowerPhiBin(towerkey);

        TowerData data = {
            .ieta = ieta,
            .iphi = iphi,
            .energy = tower->get_energy(),
            .time = tower->get_time(),
            .chi2 = tower->get_chi2(),
            .pedestal = tower->get_pedestal(),
            .good = static_cast<short>(tower->get_isGood() ? 1 : 0),
            .isAcceptable = IsAcceptableTower(tower) // For calculateIsoEt
        };
        towerDataList.push_back(data);
        validTowerCount++;
    }
    if (verbose) {
        std::cout << ANSI_COLOR_BLUE_BOLD << "Collected data for " << validTowerCount << " valid towers out of " << tower_range << " total." << ANSI_COLOR_RESET << std::endl;
    }
    
}

void caloTreeGen::processTowers(TowerInfoContainer* towerContainer,
                                float& totalCaloE,
                                std::vector<float>& towEta,
                                std::vector<float>& towPhi,
                                std::vector<float>& towE,
                                std::vector<int>& towTime,
                                std::vector<float>& towChi2,
                                std::vector<float>& towPed,
                                std::vector<short>& towGood) {
    std::vector<TowerData> towerDataList;
    collectTowerData(towerContainer, towerDataList);

    totalCaloE = 0;
    towEta.reserve(towerDataList.size());
    towPhi.reserve(towerDataList.size());
    towE.reserve(towerDataList.size());
    towTime.reserve(towerDataList.size());
    towChi2.reserve(towerDataList.size());
    towPed.reserve(towerDataList.size());
    towGood.reserve(towerDataList.size());
    if (verbose) {
        std::cout << ANSI_COLOR_BLUE_BOLD << "Processing " << towerDataList.size() << " towers..." << ANSI_COLOR_RESET << std::endl;
    }
    

    for (const auto& data : towerDataList) {
        totalCaloE += data.energy;
        towEta.push_back(data.ieta);
        towPhi.push_back(data.iphi);
        towE.push_back(data.energy);
        towTime.push_back(data.time);
        towChi2.push_back(data.chi2);
        towPed.push_back(data.pedestal);
        towGood.push_back(data.good);
    }
    if (verbose) {
        std::cout << ANSI_COLOR_BLUE_BOLD
                  << "Completed tower processing: "
                  << "Total energy = " << totalCaloE << ", "
                  << "towEta size = " << towEta.size() << ", "
                  << "towPhi size = " << towPhi.size() << ", "
                  << "towE size = " << towE.size() << ", "
                  << "towTime size = " << towTime.size() << ", "
                  << "towChi2 size = " << towChi2.size() << ", "
                  << "towPed size = " << towPed.size() << ", "
                  << "towGood size = " << towGood.size() << "."
                  << ANSI_COLOR_RESET << std::endl;
    }

}
void caloTreeGen::processEnergyMaps(const std::vector<float>* m_emcTowE, const std::vector<float>* m_emciEta, const std::vector<float>* m_emciPhi, const std::vector<float>* m_ohcTowE, const std::vector<float>* m_ohciTowEta, const std::vector<float>* m_ohciTowPhi, const std::vector<float>* m_ihcTowE, const std::vector<float>* m_ihciTowEta, const std::vector<float>* m_ihciTowPhi, std::vector<short>* m_emcal_good, std::vector<short>* m_ohc_good, std::vector<short>* m_ihc_good, std::map<std::string, TObject*>& qaHistograms, int triggerIndex) {
    
    if (verbose) {
        std::cout << ANSI_COLOR_BLUE_BOLD << "Processing tower energy maps for Trigger Index: " << triggerIndex << ANSI_COLOR_RESET << std::endl;
        std::cout << "EMCal Towers: " << m_emcTowE->size()
                  << ", IHCal Towers: " << m_ihcTowE->size()
                  << ", OHCal Towers: " << m_ohcTowE->size() << std::endl;
    }
    
    float energymap[12][32] = {0};
    float energymap_emcal[12][35] = {0};
    float energymap_hcalin[12][35] = {0};
    float energymap_hcalout[12][35] = {0};
    float energymap_extend[12][35] = {0};
    
    float energymap_jet[9][32] = {0};
    float energymap_jet_emcal[9][32] = {0};
    float energymap_jet_hcalin[9][32] = {0};
    float energymap_jet_hcalout[9][32] = {0};
    
    if (verbose) {
        std::cout << "Reset energy maps for processing..." << std::endl;
    }
    for (int j = 0; j < 35; j++) {
        for (int k =0 ; k < 12; k++) {
            if (j < 32) {
                energymap[k][j] = 0.0;
            }
            energymap_hcalin[k][j] = 0.0;
            energymap_hcalout[k][j] = 0.0;
            energymap_emcal[k][j] = 0.0;
            energymap_extend[k][j] = 0.0;
        }
    }

    if (verbose) {
        std::cout << "Processing EMCal Tower Maps..." << std::endl;
    }
    for (size_t ie = 0; ie < m_emcTowE->size(); ie++) {
        int ebin = m_emciEta->at(ie) / 8;
        int pbin = m_emciPhi->at(ie) / 8;
        
        
        if (!m_emcal_good->at(ie)) continue;
        
        energymap[ebin][pbin] += m_emcTowE->at(ie);
        energymap_emcal[ebin][pbin] += m_emcTowE->at(ie);
        energymap_extend[ebin][pbin] += m_emcTowE->at(ie);
        
        if (pbin < 3) {
            energymap_emcal[ebin][pbin+32] += m_emcTowE->at(ie);
            energymap_extend[ebin][pbin+32] += m_emcTowE->at(ie);
        }
    }
    
    if (verbose) {
        std::cout << "Processing IHCal Tower Maps..." << std::endl;
    }
    for (size_t ie = 0; ie < m_ihcTowE->size(); ie++) {
        int ebin_ihc = m_ihciTowEta->at(ie) / 8;
        int pbin_ihc = m_ihciTowPhi->at(ie) / 8;
        
        if (!m_ihc_good->at(ie)) continue;
        energymap_extend[ebin_ihc][pbin_ihc] += m_ihcTowE->at(ie);
        energymap_hcalin[ebin_ihc][pbin_ihc] += m_ihcTowE->at(ie);
        if (pbin_ihc < 3) {
            energymap_hcalin[ebin_ihc][pbin_ihc+32] += m_ihcTowE->at(ie);
            energymap_hcalin[ebin_ihc][pbin_ihc+32] += m_ihcTowE->at(ie);
        }
    }
    
    if (verbose) {
        std::cout << "Processing OHCal Tower Maps..." << std::endl;
    }
    for (size_t ie = 0; ie < m_ohcTowE->size(); ie++) {
        int ebin_ohc = m_ohciTowEta->at(ie) / 8;
        int pbin_ohc = m_ohciTowPhi->at(ie) / 8;

        if (!m_ohc_good->at(ie)) continue;
        
        energymap_extend[ebin_ohc][pbin_ohc] += m_ohcTowE->at(ie);
        energymap_hcalout[ebin_ohc][pbin_ohc] += m_ohcTowE->at(ie);
        if (pbin_ohc < 3) {
            energymap_hcalout[ebin_ohc][pbin_ohc+32] += m_ohcTowE->at(ie);
            energymap_extend[ebin_ohc][pbin_ohc+32] += m_ohcTowE->at(ie);
        }
    }
    for (int ie = 0; ie< 9; ie++) {
        for (int ip = 0 ; ip < 32; ip++) {
            energymap_jet[ie][ip] = 0.0;
            energymap_jet_emcal[ie][ip] = 0.0;
            energymap_jet_hcalin[ie][ip] = 0.0;
            energymap_jet_hcalout[ie][ip] = 0.0;

            for (int is = 0; is < 16; is++) {
                energymap_jet[ie][ip] += energymap_extend[ie + is%4][ip + is/4];
                energymap_jet_emcal[ie][ip] += energymap_emcal[ie + is%4][ip + is/4];
                energymap_jet_hcalin[ie][ip] += energymap_hcalin[ie + is%4][ip + is/4];
                energymap_jet_hcalout[ie][ip] += energymap_hcalout[ie + is%4][ip + is/4];
            }
        }
    }
    //in case want good trigger map for EMCal
//    int ebin = 0;
//    int pbin = 0;

    int jet_ebin = 0;
    int jet_pbin = 0;
    
    float max_8by8energy_emcal = 0.0;
    float max_energy_hcal = 0.0;
    float max_energy_jet = 0.0;

    
    if (verbose) {
        std::cout << "Reset max energy values for new trigger index:  " << triggerIndex << std::endl;
    }


    // Loop over eta bins and phi bins to find the maximum energy in EMCal
    for (int j = 0; j < 32; j++) {
        for (int k = 0; k < 12; k++) {
            if (k < 9) {
                if (energymap_jet[k][j] > max_energy_jet) {
                    max_energy_jet = energymap_jet[k][j];
                    jet_ebin = k;
                    jet_pbin = j;
                }
            }
            if (energymap[k][j] > max_8by8energy_emcal) {
                max_8by8energy_emcal = energymap[k][j]; // Update maximum EMCal energy
                //in case want good trigger map
//                ebin = k;
//                pbin = j;
            }
            if (energymap_hcalout[k][j] + energymap_hcalin[k][j] > max_energy_hcal) {
                max_energy_hcal = energymap_hcalin[k][j] + energymap_hcalout[k][j];
            }
        }
    }
    
    if (verbose) {
        std::cout << "Filling histograms with maximum energy values for Trigger Index: " << triggerIndex << std::endl;
    }
    // Fill histograms for overall maximum values
    ((TH1F*)qaHistograms["h8by8TowerEnergySum_" + std::to_string(triggerIndex)])->Fill(max_8by8energy_emcal);
    ((TH1F*)qaHistograms["h_hcal_energy_" + std::to_string(triggerIndex)])->Fill(max_energy_hcal);
    ((TH1F*)qaHistograms["h_jet_energy_" + std::to_string(triggerIndex)])->Fill(max_energy_jet);

    // Fill histograms for contributions from each calorimeter component
    ((TH1F*)qaHistograms["h_jet_emcal_energy_" + std::to_string(triggerIndex)])->Fill(energymap_jet_emcal[jet_ebin][jet_pbin]);
    ((TH1F*)qaHistograms["h_jet_hcalin_energy_" + std::to_string(triggerIndex)])->Fill(energymap_jet_hcalin[jet_ebin][jet_pbin]);
    ((TH1F*)qaHistograms["h_jet_hcalout_energy_" + std::to_string(triggerIndex)])->Fill(energymap_jet_hcalout[jet_ebin][jet_pbin]);
}

void caloTreeGen::processClusterInvariantMass(
    const std::vector<float>& clusterE,
    const std::vector<float>& clusterPt,
    const std::vector<float>& clusterChi2,
    const std::vector<float>& clusterEta,
    const std::vector<float>& clusterPhi,
    int triggerIndex,
    std::map<std::string, TObject*>& massHistograms,
    std::map<std::string, TObject*>& qaHistograms,
    const std::vector<float>& asymmetry_values,
    const std::vector<float>& clus_chi_values,
    const std::vector<float>& clus_Ecore_values)
{
    if (verbose) {
        std::cout << "Processing invariant mass for Trigger Index: " << triggerIndex << std::endl;
        std::cout << "Cluster sizes - Energies: " << clusterE.size()
                  << ", Pt: " << clusterPt.size()
                  << ", Chi2: " << clusterChi2.size()
                  << ", Eta: " << clusterEta.size()
                  << ", Phi: " << clusterPhi.size() << std::endl;
    }

    // Check for mismatched vector sizes, which could indicate an issue
    if (clusterE.size() != clusterPt.size() || clusterE.size() != clusterChi2.size() ||
        clusterE.size() != clusterEta.size() || clusterE.size() != clusterPhi.size()) {
        std::cerr << "Error: Mismatched cluster vector sizes for Trigger Index: " << triggerIndex << std::endl;
        return;
    }

    // Counters for summary
    size_t totalPairs = 0;
    size_t skippedPairsDueToPt = 0;
    size_t skippedPairsDueToAsymmetry = 0;
    size_t skippedPairsDueToChi2 = 0;
    size_t skippedPairsDueToEcore = 0;
    size_t filledHistogramCount = 0;
    size_t noHistogramsFilledCount = 0;

    // Loop over all pairs of clusters
    for (size_t clus1 = 0; clus1 < clusterE.size(); clus1++) {
        for (size_t clus2 = clus1 + 1; clus2 < clusterE.size(); clus2++) {
            totalPairs++;

            // Extract cluster properties
            float E1 = clusterE[clus1], E2 = clusterE[clus2];
            float pt1 = clusterPt[clus1], pt2 = clusterPt[clus2];
            float chi1 = clusterChi2[clus1], chi2 = clusterChi2[clus2];
            float eta1 = clusterEta[clus1], eta2 = clusterEta[clus2];
            float phi1 = clusterPhi[clus1], phi2 = clusterPhi[clus2];

//            // Skip combinations based on basic pT cuts
//            if (pt1 < 2 || pt1 >= 30 || pt2 < 2 || pt2 >= 30) {
//                skippedPairsDueToPt++;
//                continue;
//            }

            // Create Lorentz vectors for the clusters
            TLorentzVector photon1, photon2;
            photon1.SetPtEtaPhiE(pt1, eta1, phi1, E1);
            photon2.SetPtEtaPhiE(pt2, eta2, phi2, E2);
            TLorentzVector meson = photon1 + photon2;
            float mesonMass = meson.M();
            float asym = (E1 + E2 != 0) ? fabs(E1 - E2) / (E1 + E2) : 0;

            // Apply all combinations of cuts and fill histograms
            bool filledHistogram = false; // Track if any histogram was filled
            bool pairSkipped = false;     // Track if the pair has been skipped by any cut

            for (float maxAsym : asymmetry_values) {
                for (float maxChi2 : clus_chi_values) {
                    for (float minClusEcore : clus_Ecore_values) {
                        // Apply selection criteria and count specific cut skips, ensuring each pair is counted once
                        if (!pairSkipped && asym >= maxAsym) {
                            skippedPairsDueToAsymmetry++;
                            pairSkipped = true;
                        } else if (!pairSkipped && (chi1 >= maxChi2 || chi2 >= maxChi2)) {
                            skippedPairsDueToChi2++;
                            pairSkipped = true;
                        } else if (!pairSkipped && (E1 < minClusEcore || E2 < minClusEcore)) {
                            skippedPairsDueToEcore++;
                            pairSkipped = true;
                        }

                        // If already skipped, skip the rest of checks
                        if (pairSkipped) continue;

                        // Format the histogram name based on the cut values
                        std::string histName = "invMass_E" + formatFloatForFilename(minClusEcore) +
                                               "_Chi" + formatFloatForFilename(maxChi2) +
                                               "_Asym" + formatFloatForFilename(maxAsym) +
                                               "_" + std::to_string(triggerIndex);

                        // Attempt to find and fill the histogram
                        TH1F* hist = (TH1F*)massHistograms[histName];
                        if (!hist) {
                            std::cerr << "Error: Histogram " << histName
                                      << " not found when trying to fill for Trigger Index: "
                                      << triggerIndex << std::endl;
                            continue;
                        }

                        hist->Fill(mesonMass);
                        filledHistogram = true;
                        filledHistogramCount++;
                    }
                }
            }

            // Increment count if no histograms were filled for this pair
            if (!filledHistogram) {
                noHistogramsFilledCount++;
            }
        }
    }

    // Print summary of processing
    if (verbose) {
        std::cout << "Completed processing invariant mass for Trigger Index: " << triggerIndex << std::endl;
        std::cout << "Total pairs processed: " << totalPairs << std::endl;
        std::cout << "Pairs skipped due to pT cuts: " << skippedPairsDueToPt
                  << " (" << (100.0 * skippedPairsDueToPt / totalPairs) << "%)" << std::endl;
        std::cout << "Pairs skipped due to asymmetry cuts: " << skippedPairsDueToAsymmetry
                  << " (" << (100.0 * skippedPairsDueToAsymmetry / totalPairs) << "%)" << std::endl;
        std::cout << "Pairs skipped due to chi2 cuts: " << skippedPairsDueToChi2
                  << " (" << (100.0 * skippedPairsDueToChi2 / totalPairs) << "%)" << std::endl;
        std::cout << "Pairs skipped due to Ecore cuts: " << skippedPairsDueToEcore
                  << " (" << (100.0 * skippedPairsDueToEcore / totalPairs) << "%)" << std::endl;
        std::cout << "Pairs with no histograms filled: " << noHistogramsFilledCount
                  << " (" << (100.0 * noHistogramsFilledCount / totalPairs) << "%)" << std::endl;
        std::cout << "Histograms filled: " << filledHistogramCount << std::endl;
    }
}

//____________________________________________________________________________..
int caloTreeGen::process_event(PHCompositeNode *topNode) {
    event_count++;
    // If event limit is enabled and the limit is exceeded, stop processing
    if (m_limitEvents && event_count > m_eventLimit) {
        if (verbose) {
            std::cout << "Event limit of " << m_eventLimit << " reached, skipping further events." << std::endl;
        }
        return Fun4AllReturnCodes::EVENT_OK;  // Skip further events
    }
    
    std::cout << "\n========== Processing CALOTREEGEN -- Event " << event_count << " ==========\n";

    _gl1_packet = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    if (_gl1_packet) {
        b_gl1_scaledvec = _gl1_packet->lValue(0, "ScaledVector");
    }
    std::vector<int> activeTriggerBits = extractTriggerBits(b_gl1_scaledvec, event_count);
    if (verbose) {
        for (int bit : activeTriggerBits) {
            std::cout << bit << " ";
        }
    }

    std::cout << std::endl;

    for (int triggerIndex : triggerIndices) {
        // Check if the current trigger bit is active
        if (!checkTriggerCondition(activeTriggerBits, triggerIndex)) {
            continue;  // Skip if this trigger bit is not active
        }
        // Fill the trigger count histogram to count each time this trigger fires
        TH1F* hTriggerCount = (TH1F*)qaHistogramsByTrigger[triggerIndex]["hTriggerCount_" + std::to_string(triggerIndex)];
        if (hTriggerCount) {
            hTriggerCount->Fill(0.5); // Fill the single bin to count the trigger occurrence
        }
        
        if (verbose) {
            std::cout << "Processing Trigger Bit: " << triggerIndex << std::endl;
        }

        
        if (verbose) {
            std::cout << "Vectors cleared for Trigger Bit: " << triggerIndex << std::endl;
        }
        
        auto& qaHistograms = qaHistogramsByTrigger[triggerIndex];
        
        GlobalVertexMap* vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
        m_vx = m_vy = m_vz = 0; // Initialize vertex coordinates to zero

        // Check if the GlobalVertexMap node is missing
        if (!vertexmap) {
            std::cout << "Error: GlobalVertexMap node is missing." << std::endl;
        } else {
            if (verbose) {
                std::cout << "GlobalVertexMap node found." << std::endl;
            }
            // Check if the vertex map is empty
            if (vertexmap->empty()) {
                if (verbose) {
                    std::cout << "Warning: GlobalVertexMap is empty." << std::endl;
                }
                return Fun4AllReturnCodes::ABORTEVENT;
            }

            // Access the first vertex in the map
            GlobalVertex* vtx = vertexmap->begin()->second;

            if (vtx) {
                // Retrieve vertex coordinates
                m_vx = vtx->get_x();
                m_vy = vtx->get_y();
                m_vz = vtx->get_z();
                if (verbose) {
                    std::cout << "Vertex coordinates retrieved: "
                              << "x = " << m_vx << ", "
                              << "y = " << m_vy << ", "
                              << "z = " << m_vz << std::endl;
                }

                // Check if the histogram exists before filling it
                std::string histName = "hVtxZ_" + std::to_string(triggerIndex);
                if (qaHistograms.find(histName) != qaHistograms.end()) {
                    // Fill the histogram with the z-coordinate of the vertex
                    ((TH1F*)qaHistograms[histName])->Fill(m_vz);
                    if (verbose) {
                        std::cout << "Filled histogram " << histName << " with m_vz = " << m_vz << std::endl;
                    }
                } else if (verbose) {
                    std::cerr << "Error: Histogram " << histName << " not found in qaHistograms!" << std::endl;
                }

                // Apply a cut on the absolute value of the z vertex
                if (std::abs(m_vz) >= 60) {
                    if (verbose) {
                        std::cout << "Skipping event: |m_vz| = " << std::abs(m_vz) << " is outside the allowed range of |30 cm|." << std::endl;
                    }
                    return Fun4AllReturnCodes::ABORTEVENT; // Skip the rest of the event processing
                }

            } else if (verbose) {
                std::cout << "Warning: Vertex object is null." << std::endl;
            }
        }
        //switched to UE subtracted
        TowerInfoContainer* emcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER_SUB1");
        TowerInfoContainer* ohcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN_SUB1");
        TowerInfoContainer* ihcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT_SUB1");
        if (!emcTowerContainer && !ihcTowerContainer && !ohcTowerContainer) {
            std::cout << ANSI_COLOR_RED_BOLD << "No tower containers found, skipping tower processing." << ANSI_COLOR_RESET << std::endl;
            continue;
        }

        // Declare geometry container pointers
        RawTowerGeomContainer* geomEM = nullptr;
        RawTowerGeomContainer* geomIH = nullptr;
        RawTowerGeomContainer* geomOH = nullptr;

        // Fetch geometry containers
        std::vector<std::tuple<RawTowerGeomContainer**, const char*, const char*>> geomContainers = {
            {&geomEM, "TOWERGEOM_CEMC", "EMC"}, //switched to HCAL geometry with UE subtracted
            {&geomIH, "TOWERGEOM_HCALIN", "Inner HCal"},
            {&geomOH, "TOWERGEOM_HCALOUT", "Outer HCal"}
        };

        for (auto& [geomContainer, geomName, geomLabel] : geomContainers) {
            *geomContainer = findNode::getClass<RawTowerGeomContainer>(topNode, geomName);
            if (!*geomContainer) {
                std::cout << ANSI_COLOR_RED_BOLD << "Error: Missing " << geomLabel << " Geometry Container: " << geomName << ANSI_COLOR_RESET << std::endl;
                return Fun4AllReturnCodes::ABORTEVENT;
            }
            if (verbose) {
                std::cout << ANSI_COLOR_BLUE_BOLD << "Loaded " << geomLabel << " Geometry Container." << ANSI_COLOR_RESET << std::endl;
            }
        }
        
        // Process towers and fill histograms
        if (emcTowerContainer) {
            if (verbose) {
                std::cout << ANSI_COLOR_BLUE_BOLD << "Processing EMCal Towers..." << ANSI_COLOR_RESET << std::endl;
            }

            processTowers(emcTowerContainer, totalCaloEEMCal, m_emciEta, m_emciPhi, m_emcTowE, m_emcTime, m_emcChi2, m_emcPed, m_emcal_good);
            for (size_t i = 0; i < m_emcTowE.size(); ++i) {
                ((TH2F*)qaHistograms["h2_EMCal_TowerEtaPhi_Weighted_2D_" + std::to_string(triggerIndex)])->Fill(m_emciEta[i], m_emciPhi[i], m_emcTowE[i]);
            }
            // Check if the total EMCal energy is positive or negative and fill the appropriate histogram
            if (totalCaloEEMCal >= 0) {
                ((TH1F*)qaHistograms["hTotalCaloEEMCal_" + std::to_string(triggerIndex)])->Fill(totalCaloEEMCal);
            } else {
                ((TH1F*)qaHistograms["hTotalCalo_Negative_EEMCal_" + std::to_string(triggerIndex)])->Fill(totalCaloEEMCal);
            }
            //maximum chi2 of all the towers fill the histograms for EMCAL chi2 for this event -- similar for other calos
            if (!m_emcChi2.empty()) {
                float max_emcalChi2 = *std::max_element(m_emcChi2.begin(), m_emcChi2.end());
                ((TH1F*)qaHistograms["h_emcalChi2_" + std::to_string(triggerIndex)])->Fill(max_emcalChi2);
            }

        }
        if (ihcTowerContainer) {
            if (verbose) {
                std::cout << ANSI_COLOR_BLUE_BOLD << "Processing IHCal Towers..." << ANSI_COLOR_RESET << std::endl;
            }
            
            processTowers(ihcTowerContainer, totalCaloEIHCal, m_ihciTowEta, m_ihciTowPhi, m_ihcTowE, m_ihcTime, m_ihcChi2, m_ihcPed, m_ihc_good);
            for (size_t i = 0; i < m_ihcTowE.size(); ++i) {
                ((TH2F*)qaHistograms["h2_IHCal_TowerEnergy_" + std::to_string(triggerIndex)])->Fill(m_ihciTowEta[i], m_ihciTowPhi[i], m_ihcTowE[i]);
            }
            ((TH1F*)qaHistograms["hTotalCaloEIHCal_" + std::to_string(triggerIndex)])->Fill(totalCaloEIHCal);
            if (!m_ihcChi2.empty()) {
                float max_ihcalChi2 = *std::max_element(m_ihcChi2.begin(), m_ihcChi2.end());
                ((TH1F*)qaHistograms["h_ihcalChi2_" + std::to_string(triggerIndex)])->Fill(max_ihcalChi2);
            }
        }
        if (ohcTowerContainer) {
            if (verbose) {
                std::cout << ANSI_COLOR_BLUE_BOLD << "Processing OHCal Towers..." << ANSI_COLOR_RESET << std::endl;
            }
            
            processTowers(ohcTowerContainer, totalCaloEOHCal, m_ohciTowEta, m_ohciTowPhi, m_ohcTowE, m_ohcTime, m_ohcChi2, m_ohcPed, m_ohc_good);
            for (size_t i = 0; i < m_ohcTowE.size(); ++i) {
                ((TH2F*)qaHistograms["h2_OHCal_TowerEnergy_" + std::to_string(triggerIndex)])->Fill(m_ohciTowEta[i], m_ohciTowPhi[i], m_ohcTowE[i]);
            }
            ((TH1F*)qaHistograms["hTotalCaloEOHCal_" + std::to_string(triggerIndex)])->Fill(totalCaloEOHCal);
            if (!m_ohcChi2.empty()) {
                float max_ohcalChi2 = *std::max_element(m_ohcChi2.begin(), m_ohcChi2.end());
                ((TH1F*)qaHistograms["h_ohcalChi2_" + std::to_string(triggerIndex)])->Fill(max_ohcalChi2);
            }
        }

        // Display vector sizes and confirm entry into the function call if verbose is enabled
        if (verbose) {
            std::cout << ANSI_COLOR_RED_BOLD << "Preparing to call processEnergyMaps for Trigger Index: "
                      << triggerIndex << ANSI_COLOR_RESET << std::endl;
            std::cout << ANSI_COLOR_RED_BOLD << "Vector Sizes: " << ANSI_COLOR_RESET << std::endl;
            std::cout << ANSI_COLOR_RED_BOLD << "m_emcTowE size: " << m_emcTowE.size()
                      << ", m_emciEta size: " << m_emciEta.size()
                      << ", m_emciPhi size: " << m_emciPhi.size() << ANSI_COLOR_RESET << std::endl;
            std::cout << ANSI_COLOR_RED_BOLD << "m_ohcTowE size: " << m_ohcTowE.size()
                      << ", m_ohciTowEta size: " << m_ohciTowEta.size()
                      << ", m_ohciTowPhi size: " << m_ohciTowPhi.size() << ANSI_COLOR_RESET << std::endl;
            std::cout << ANSI_COLOR_RED_BOLD << "m_ihcTowE size: " << m_ihcTowE.size()
                      << ", m_ihciTowEta size: " << m_ihciTowEta.size()
                      << ", m_ihciTowPhi size: " << m_ihciTowPhi.size() << ANSI_COLOR_RESET << std::endl;
            std::cout << ANSI_COLOR_RED_BOLD << "m_emcal_good size: " << m_emcal_good.size()
                      << ", m_ohc_good size: " << m_ohc_good.size()
                      << ", m_ihc_good size: " << m_ihc_good.size() << ANSI_COLOR_RESET << std::endl;
            std::cout << ANSI_COLOR_RED_BOLD << "Entering processEnergyMaps function for Trigger Index: "
                      << triggerIndex << ANSI_COLOR_RESET << std::endl;
        }

        // Pass the addresses of the pointers to the function
        processEnergyMaps(&m_emcTowE, &m_emciEta, &m_emciPhi, &m_ohcTowE, &m_ohciTowEta, &m_ohciTowPhi,
                          &m_ihcTowE, &m_ihciTowEta, &m_ihciTowPhi,
                          &m_emcal_good, &m_ohc_good, &m_ihc_good, qaHistograms, triggerIndex);

        if (verbose) {
            std::cout << ANSI_COLOR_RED_BOLD << "Exited processEnergyMaps function for Trigger Index: "
                      << triggerIndex << ANSI_COLOR_RESET << std::endl;
        }
        

        if (verbose) {
            std::cout << "Finished processing energy maps for Trigger Index: " << triggerIndex << std::endl;
        }
        
        
        if (verbose) {
            std::cout << "Fetching RawClusterContainer..." << std::endl;
        }

        RawClusterContainer *clusterContainer = findNode::getClass<RawClusterContainer>(topNode, "CLUSTERINFO_CEMC");
        if(!clusterContainer) {
            std::cout << PHWHERE << "Cluster node is missing...ruh roh scooby doo" << std::endl;
            return 0;
        }
        
        if (verbose) {
            std::cout << "Processing clusters..." << std::endl;
        }
        
        RawClusterContainer::ConstRange clusterRange = clusterContainer->getClusters();
        int nan_count = 0;
        float max_energy_clus = 0.0;
        float max_isoEt = std::numeric_limits<float>::lowest();
        float min_isoEt = std::numeric_limits<float>::max();
        
        for (auto clusterIter = clusterRange.first; clusterIter != clusterRange.second; ++clusterIter) {
            RawCluster* cluster = clusterIter->second;
            if (!cluster) {
                std::cout << "Warning: Null cluster found." << std::endl;
                continue;
            }
            /*
             the first 1k events in every segment are used for MBD calibration and are guaranteed to not have a vertex -- can see after the first 1000 events its found
             */
            CLHEP::Hep3Vector vertex(m_vx, m_vy, m_vz);
            
            CLHEP::Hep3Vector Ecore_vec_cluster = RawClusterUtility::GetECoreVec(*cluster, vertex);
            CLHEP::Hep3Vector E_vec_cluster_Full = RawClusterUtility::GetEVec(*cluster, vertex);

            float clusE = E_vec_cluster_Full.mag(); //only vartiable that uses not GetECoreVec--NOT ECORE
            // Check if histograms exist and fill them
            TH1F* h_clusE = (TH1F*)qaHistograms["h_ClusterEnergy_" + std::to_string(triggerIndex)];
            if (!h_clusE) {
                std::cerr << "Error: h_ClusterEnergy_" << triggerIndex << " histogram is null!" << std::endl;
            } else {
                h_clusE->Fill(clusE);
            }
            float clusEcore = Ecore_vec_cluster.mag();

            if (clusEcore < 1) { // cut on ecore
                continue;
            }
            
            
            m_clusterECore.push_back(clusEcore);
            float clus_eta = Ecore_vec_cluster.pseudoRapidity();
            float clus_phi = Ecore_vec_cluster.phi();
            m_clusterPhi.push_back(clus_phi);
            m_clusterEta.push_back(clus_eta);
            
            float clus_eT = clusE / std::cosh(clus_eta);

            // Fill the eT histogram
            TH1F* h_ET = (TH1F*)qaHistograms["h_ET_" + std::to_string(triggerIndex)];
            if (!h_ET) {
                std::cerr << "Error: h_ET_" << triggerIndex << " histogram is null!" << std::endl;
            } else {
                h_ET->Fill(clus_eT);
            }


            float clus_pt = Ecore_vec_cluster.perp();
            m_clusterPt.push_back(clus_pt);
        
            float clus_chi = cluster -> get_chi2();
            m_clusterChi.push_back(clus_chi);
            
            
            float maxTowerEnergy = getMaxTowerE(cluster,emcTowerContainer);
            // Check if histograms exist and fill them
            TH1F* h_maxTowE = (TH1F*)qaHistograms["h_maxTowerEnergy_" + std::to_string(triggerIndex)];
            if (!h_maxTowE) {
                std::cerr << "Error: h_maxTowerEnergy_" << triggerIndex << " histogram is null!" << std::endl;
            } else {
                h_maxTowE->Fill(maxTowerEnergy);
            }
            
            m_clusTowPhi.push_back(returnClusterTowPhi(cluster,emcTowerContainer));
            m_clusTowEta.push_back(returnClusterTowEta(cluster,emcTowerContainer));
            m_clusTowE.push_back(returnClusterTowE(cluster,emcTowerContainer));
            
            float et_iso = cluster->get_et_iso(3, 1, 1);
            // Check for NaN and update nan_count
            if (std::isnan(et_iso)) {
                std::cout << "Warning: Isolation energy is NaN for cluster ID: " << cluster->get_id() << std::endl;
                nan_count++;  // Increment the NaN count
            } else {
                // Update max and min isoEt
                if (et_iso > max_isoEt) {
                    max_isoEt = et_iso;
                }
                if (et_iso < min_isoEt) {
                    min_isoEt = et_iso;
                }
            }
            
            m_clusterEtIso.push_back(et_iso);
            
            // Access histograms from the map
            TH1F* hPt = (TH1F*)qaHistograms["hClusterPt_" + std::to_string(triggerIndex)];
            TH1F* hChi2 = (TH1F*)qaHistograms["hClusterChi2_" + std::to_string(triggerIndex)];

            // Check if histograms exist and fill them
            if (!hPt) {
                std::cerr << "Error: hClusterPt_" << triggerIndex << " histogram is null!" << std::endl;
            } else {
                hPt->Fill(clus_pt);
            }

            if (!hChi2) {
                std::cerr << "Error: hClusterChi2_" << triggerIndex << " histogram is null!" << std::endl;
            } else {
                hChi2->Fill(clus_chi);
            }
            
            if (clusEcore > max_energy_clus) {
                max_energy_clus = clusEcore; // Update the maximum cluster energy
            }

        }
        
        if (verbose) {
            std::cout << "Finished processing clusters for Trigger Index: " << triggerIndex << std::endl;
            // After processing all clusters, print the total NaN count and max/min isoEt
            std::cout << "Trigger Index " << triggerIndex << " had " << nan_count << " clusters with NaN isolation energy out of "
                      << m_clusterEtIso.size() << " total clusters." << std::endl;
            std::cout << "Max isolation energy (isoEt) for Trigger Index " << triggerIndex << ": " << max_isoEt << std::endl;
            std::cout << "Min isolation energy (isoEt) for Trigger Index " << triggerIndex << ": " << min_isoEt << std::endl;
        }


        try {
            for (size_t j = 0; j < m_clusterEtIso.size(); ++j) {
                // Error checking for accessing the vectors
                try {
                    // Get the Ecore and isoEt from the ClusterIso calculation
                    float ecore = m_clusterECore.at(j);  // This can throw out_of_range
                    float isoEt_CoreSoft = m_clusterEtIso.at(j);  // This can throw out_of_range

                    // Error checking for histogram existence before filling
                    auto hist2D = qaHistograms["h2_cluster_iso_Ecore_" + std::to_string(triggerIndex)];
                    auto hist1D = qaHistograms["h1_isoEt_" + std::to_string(triggerIndex)];

                    if (!hist2D) {
                        throw std::runtime_error("Error: h2_cluster_iso_Ecore_" + std::to_string(triggerIndex) + " is null.");
                    }
                    if (!hist1D) {
                        throw std::runtime_error("Error: h1_isoEt_" + std::to_string(triggerIndex) + " is null.");
                    }

                    ((TH1F*)hist2D)->Fill(ecore, isoEt_CoreSoft);
                    ((TH1F*)hist1D)->Fill(isoEt_CoreSoft);
                }
                catch (const std::out_of_range& e) {
                    if (verbose) {
                        std::cerr << "Out of range error: " << e.what() << " while accessing cluster data at index " << j << std::endl;
                    }
                }
                catch (const std::runtime_error& e) {
                    if (verbose) {
                        std::cerr << e.what() << " for trigger index " << triggerIndex << " at cluster index " << j << std::endl;
                    }
                }
            }
        }
        catch (const std::exception& e) {
            if (verbose) {
                std::cerr << "Unexpected error during the cluster et ecore histogram filling: " << e.what() << std::endl;
            }
        }


        
        // Fill the histogram with the maximum cluster energy core value
        TH1F* h_maxECore = (TH1F*)qaHistograms["hCluster_maxECore_" + std::to_string(triggerIndex)];

        if (!h_maxECore) {
            std::cerr << "Error: Histogram hCluster_maxECore_" << triggerIndex << " is null and cannot be filled." << std::endl;
        } else {
            h_maxECore->Fill(max_energy_clus);
            
            if (verbose) {
                std::cout << "Filled histogram hCluster_maxECore_" << triggerIndex
                          << " with value: " << max_energy_clus << std::endl;
            }
        }

        
        // Before calling the invariant mass processing, check the size of the vectors
        if (verbose) {
            std::cout << ANSI_COLOR_BLUE_BOLD << "Trigger Index: " << triggerIndex << ANSI_COLOR_RESET << std::endl;
            std::cout << "Cluster vectors sizes before processing invariant mass:" << std::endl;
            std::cout << "  m_clusterECore size: " << m_clusterECore.size() << std::endl;
            std::cout << "  m_clusterPt size: " << m_clusterPt.size() << std::endl;
            std::cout << "  m_clusterChi size: " << m_clusterChi.size() << std::endl;
            std::cout << "  m_clusterEta size: " << m_clusterEta.size() << std::endl;
            std::cout << "  m_clusterPhi size: " << m_clusterPhi.size() << std::endl;
        }




        // Call the invariant mass processing function
        processClusterInvariantMass(m_clusterECore, m_clusterPt, m_clusterChi, m_clusterEta, m_clusterPhi,
            triggerIndex, massHistogramsByTrigger[triggerIndex], qaHistogramsByTrigger[triggerIndex],
            asymmetry_values, clus_chi_values, clus_Ecore_values);


        const auto& massHistograms = massHistogramsByTrigger[triggerIndex];
        if (verbose) {
            if (massHistograms.empty()) {
                std::cerr << "Warning: No histograms found in massHistogramsByTrigger[" << triggerIndex << "] after processing invariant mass." << std::endl;
            } else {
                for (const auto& [name, hist] : massHistograms) {
                    if (!hist) {
                        std::cerr << "Error: Histogram " << name << " is null for Trigger Index: " << triggerIndex << "." << std::endl;
                    } else {
                        std::cout << "Histogram " << name << " found in massHistogramsByTrigger[" << triggerIndex << "]." << std::endl;
                    }
                }
            }
        }

        
        ResetEvent(topNode);// Reset after each trigger bit
    }
    return Fun4AllReturnCodes::EVENT_OK;
}


//____________________________________________________________________________..
int caloTreeGen::ResetEvent(PHCompositeNode *topNode) {
    if (verbose) {
        std::cout << ANSI_COLOR_BLUE_BOLD << "Resetting event..." << ANSI_COLOR_RESET << std::endl;
    }

    m_clusterE.clear();
    m_clusterPhi.clear();
    m_clusterEta.clear();
    m_clusterPt.clear();
    m_clusterChi.clear();
    m_clusterTowMaxE.clear();
    m_clusterECore.clear();
    m_clusterEtIso.clear();
      
    m_emcTowE.clear();
    m_emciEta.clear();
    m_emciPhi.clear();
    m_emcTime.clear();
    m_emcChi2.clear();
    m_emcPed.clear();
    m_emcal_good.clear();

    m_ihcTowE.clear();
    m_ihciTowEta.clear();
    m_ihciTowPhi.clear();
    m_ihcTime.clear();
    m_ihcChi2.clear();
    m_ihcPed.clear();
    m_ihc_good.clear();

    m_ohcTowE.clear();
    m_ohciTowEta.clear();
    m_ohciTowPhi.clear();
    m_ohcTime.clear();
    m_ohcChi2.clear();
    m_ohcPed.clear();
    m_ohc_good.clear();

    m_clusTowPhi.clear();
    m_clusTowEta.clear();
    m_clusTowE.clear();
    
    return Fun4AllReturnCodes::EVENT_OK;

    if (verbose) {
        std::cout << ANSI_COLOR_GREEN_BOLD << "Event reset complete." << ANSI_COLOR_RESET << std::endl;
    }
}

//____________________________________________________________________________..
int caloTreeGen::EndRun(const int runnumber) {
  std::cout << "caloTreeGen::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int caloTreeGen::End(PHCompositeNode *topNode) {
    std::cout << ANSI_COLOR_BLUE_BOLD << "caloTreeGen::End(PHCompositeNode *topNode) All events have been processed. Beginning final analysis steps..." << ANSI_COLOR_RESET << std::endl;

    // Write all QA histograms to file and clean up
    for (const auto& [triggerIndex, qaHistograms] : qaHistogramsByTrigger) {
        std::string qaDir = "QA/Trigger" + std::to_string(triggerIndex);
        if (!out->cd(qaDir.c_str())) {
            std::cerr << ANSI_COLOR_RED_BOLD << "Error: Failed to change directory to " << qaDir << " when writing QA histograms." << ANSI_COLOR_RESET << std::endl;
            continue;  // Skip this trigger index since the directory change failed
        } else if (verbose) {
            std::cout << ANSI_COLOR_BLUE_BOLD << "Changed to QA directory: " << qaDir << ANSI_COLOR_RESET << std::endl;
        }

        for (const auto& [name, hist] : qaHistograms) {
            if (!hist) {
                std::cerr << ANSI_COLOR_RED_BOLD << "Warning: QA Histogram " << name << " is null, skipping write." << ANSI_COLOR_RESET << std::endl;
                continue;
            }

            int writeResult = hist->Write();
            if (writeResult == 0) {
                std::cerr << ANSI_COLOR_RED_BOLD << "Error: Failed to write QA histogram " << name << " to directory " << qaDir << "." << ANSI_COLOR_RESET << std::endl;
            } else if (verbose) {
                std::cout << ANSI_COLOR_GREEN_BOLD << "Successfully wrote QA histogram " << name << " to directory " << qaDir << "." << ANSI_COLOR_RESET << std::endl;
            }

            try {
                delete hist;
                if (verbose) {
                    std::cout << ANSI_COLOR_GREEN_BOLD << "Successfully deleted QA histogram " << name << " from memory." << ANSI_COLOR_RESET << std::endl;
                }
            } catch (const std::exception& e) {
                std::cerr << ANSI_COLOR_RED_BOLD << "Exception occurred while deleting QA histogram " << name << ": " << e.what() << ANSI_COLOR_RESET << std::endl;
            } catch (...) {
                std::cerr << ANSI_COLOR_RED_BOLD << "Unknown error occurred while deleting QA histogram " << name << "." << ANSI_COLOR_RESET << std::endl;
            }
        }
    }

    // Write all invariant mass histograms to file and clean up
    for (const auto& [triggerIndex, massHistograms] : massHistogramsByTrigger) {
        std::string invMassDir = "InvariantMassDistributions/Trigger" + std::to_string(triggerIndex);
        
        if (!out->cd(invMassDir.c_str())) {
            std::cerr << ANSI_COLOR_RED_BOLD << "Error: Failed to change directory to " << invMassDir << " when writing invariant mass histograms." << ANSI_COLOR_RESET << std::endl;
            continue;  // Skip this trigger index since the directory change failed
        } else if (verbose) {
            std::cout << ANSI_COLOR_BLUE_BOLD << "Changed to directory: " << invMassDir << ANSI_COLOR_RESET << std::endl;
        }

        for (const auto& [name, hist] : massHistograms) {
            if (!hist) {
                std::cerr << ANSI_COLOR_RED_BOLD << "Warning: Invariant Mass Histogram " << name << " is null, skipping write." << ANSI_COLOR_RESET << std::endl;
                continue;
            }

            int writeResult = hist->Write();
            if (writeResult == 0) {
                std::cerr << ANSI_COLOR_RED_BOLD << "Error: Failed to write invariant mass histogram " << name << " to directory " << invMassDir << "." << ANSI_COLOR_RESET << std::endl;
            } else if (verbose) {
                std::cout << ANSI_COLOR_GREEN_BOLD << "Successfully wrote invariant mass histogram " << name << " to directory " << invMassDir << "." << ANSI_COLOR_RESET << std::endl;
            }

            try {
                delete hist;
                if (verbose) {
                    std::cout << ANSI_COLOR_GREEN_BOLD << "Successfully deleted invariant mass histogram " << name << " from memory." << ANSI_COLOR_RESET << std::endl;
                }
            } catch (const std::exception& e) {
                std::cerr << ANSI_COLOR_RED_BOLD << "Exception occurred while deleting invariant mass histogram " << name << ": " << e.what() << ANSI_COLOR_RESET << std::endl;
            } catch (...) {
                std::cerr << ANSI_COLOR_RED_BOLD << "Unknown error occurred while deleting invariant mass histogram " << name << "." << ANSI_COLOR_RESET << std::endl;
            }
        }
    }

    // Close the output file and clean up
    std::cout << ANSI_COLOR_BLUE_BOLD << "Closing output file and cleaning up..." << ANSI_COLOR_RESET << std::endl;
    if (out) {
        out->Close();
        delete out;
        out = nullptr;
        std::cout << ANSI_COLOR_GREEN_BOLD << "Output file successfully closed and deleted." << ANSI_COLOR_RESET << std::endl;
    } else {
        std::cerr << ANSI_COLOR_RED_BOLD << "Warning: Output file was already null, possibly already closed or deleted." << ANSI_COLOR_RESET << std::endl;
    }

    std::cout << ANSI_COLOR_GREEN_BOLD << "End of caloTreeGen::End(PHCompositeNode *topNode). Exiting smoothly." << ANSI_COLOR_RESET << std::endl;

    return Fun4AllReturnCodes::EVENT_OK;
}


//____________________________________________________________________________..
int caloTreeGen::Reset(PHCompositeNode *topNode) {
 std::cout << "caloTreeGen::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void caloTreeGen::Print(const std::string &what) const {
  std::cout << "caloTreeGen::Print(const std::string &what) const Printing info for " << what << std::endl;
}
//____________________________________________________________________________..
float caloTreeGen::getMaxTowerE(RawCluster *cluster, TowerInfoContainer *towerContainer) {
  RawCluster::TowerConstRange towers = cluster -> get_towers();
  RawCluster::TowerConstIterator toweriter;
  
  float maxEnergy = 0;
  for(toweriter = towers.first; toweriter != towers.second; toweriter++)
    {
      float towE = toweriter -> second;
   
      if( towE > maxEnergy)  maxEnergy = towE;
    }
  return maxEnergy;
}
//____________________________________________________________________________..
std::vector<int> caloTreeGen::returnClusterTowEta(RawCluster *cluster, TowerInfoContainer *towerContainer) {
  RawCluster::TowerConstRange towers = cluster -> get_towers();
  RawCluster::TowerConstIterator toweriter;
  
  std::vector<int> towerIDsEta;
  for(toweriter = towers.first; toweriter != towers.second; toweriter++) towerIDsEta.push_back(RawTowerDefs::decode_index1(toweriter -> first));

  return towerIDsEta;
}
//____________________________________________________________________________..
std::vector<int> caloTreeGen::returnClusterTowPhi(RawCluster *cluster, TowerInfoContainer *towerContainer) {
  RawCluster::TowerConstRange towers = cluster -> get_towers();
  RawCluster::TowerConstIterator toweriter;
  
  std::vector<int> towerIDsPhi;
  for(toweriter = towers.first; toweriter != towers.second; toweriter++) towerIDsPhi.push_back(RawTowerDefs::decode_index2(toweriter -> first));
  return towerIDsPhi;
}
//____________________________________________________________________________..
std::vector<float> caloTreeGen::returnClusterTowE(RawCluster *cluster, TowerInfoContainer *towerContainer) {
  RawCluster::TowerConstRange towers = cluster -> get_towers();
  RawCluster::TowerConstIterator toweriter;
  
  std::vector<float> towerE;
  for(toweriter = towers.first; toweriter != towers.second; toweriter++) towerE.push_back(toweriter -> second);
  
  return towerE;
}

bool caloTreeGen::IsAcceptableTower(TowerInfo *tower) {
  if (tower->get_isBadTime()) {
    return false;
  }
  if (tower->get_isHot()) {
    return false;
  }

  if (tower->get_isBadChi2()) {
    return false;
  }

  if (tower->get_isNotInstr()) {
    return false;
  }
  return true;
}
