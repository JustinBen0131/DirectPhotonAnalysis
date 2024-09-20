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


//GL1 Information
#include <ffarawobjects/Gl1Packet.h>

//for cluster vertex correction
#include <CLHEP/Geometry/Point3D.h>

//for the vertex
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#define ANSI_COLOR_RED_BOLD "\033[1;31m"
#define ANSI_COLOR_BLUE_BOLD "\033[1;34m"
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
caloTreeGen::~caloTreeGen()
{
  std::cout << "caloTreeGen::~caloTreeGen() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int caloTreeGen::Init(PHCompositeNode *topNode) {
    out = new TFile(Outfile.c_str(),"RECREATE");
    
    for (int triggerIndex : triggerIndices) {
        // Create trigger-specific directories within the output file
        std::string qaDir = "QA/Trigger" + std::to_string(triggerIndex);
        std::string invMassDir = "InvariantMassDistributions/Trigger" + std::to_string(triggerIndex);

        out->mkdir(qaDir.c_str());
        out->mkdir(invMassDir.c_str());

        // Create QA histograms
        out->cd(qaDir.c_str());
        std::map<std::string, TObject*> qaHistograms;
        qaHistograms["h2_EMCal_TowerEtaPhi_2D"] = new TH2F("h2_EMCal_TowerEtaPhi_2D", "EMCal Tower Energy; iEta; iPhi; Energy (GeV)", 96, 0, 96, 256, 0, 256);
        qaHistograms["h2_OHCal_TowerEnergy"] = new TH2F("h2_OHCal_TowerEnergy", "HCal Tower Energy; iEta; iPhi; Energy (GeV)", 24, 0, 24, 64, 0, 64);
        qaHistograms["h2_IHCal_TowerEnergy"] = new TH2F("h2_IHCal_TowerEnergy", "HCal Tower Energy; iEta; iPhi; Energy (GeV)", 24, 0, 24, 64, 0, 64);
        qaHistograms["hTotalCaloEEMCal"] = new TH1F("hTotalCaloEEMCal", "Total EMCal Energy; Energy (GeV)", 100, 0, 500);
        qaHistograms["hTotalCaloEOHCal"] = new TH1F("hTotalCaloEOHCal", "Total OHCal Energy; Energy (GeV)", 100, 0, 500);
        qaHistograms["hTotalCaloEIHCal"] = new TH1F("hTotalCaloEIHCal", "Total IHCal Energy; Energy (GeV)", 100, 0, 500);
        qaHistograms["hClusterChi2"] = new TH1F("hClusterChi2", "Cluster Chi2; Chi2", 100, 0, 100);
        qaHistograms["h_emcalChi2"] = new TH1F("h_ihcalChi2", "Cluster Chi2; Chi2", 100, 0, 100);
        qaHistograms["h_ohcalChi2"] = new TH1F("h_ohcalChi2", "Cluster Chi2; Chi2", 100, 0, 100);
        qaHistograms["h_ihcalChi2"] = new TH1F("h_ihcalChi2", "Cluster Chi2; Chi2", 100, 0, 100);
        
        qaHistograms["hCluster_maxECore"] = new TH1F("hCluster_maxECore", "Max Cluster ECore; Cluster ECore [GeV]", 40, 0, 20);
        qaHistograms["hClusterECore"] = new TH1F("hClusterECore", "Cluster ECore; Cluster ECore [GeV]", 40, 0, 20);
        qaHistograms["hClusterPt"] = new TH1F("hClusterPt", "Cluster pT; Cluster pT [GeV]", 100, 0, 100);
        qaHistograms["hVtxZ"] = new TH1F("hVtxZ", "Z-vertex Distribution; z [cm]", 100, -70, 70);
        qaHistograms["h8by8TowerEnergySum"] = new TH1F("h8by8TowerEnergySum", "Max 8x8 Tower Energy Sum; Energy [GeV]; Events", 40, 0, 20);
        qaHistograms["h_hcal_energy"] = new TH1F("h_hcal_energy", "Max HCal Tower Energy Sums; Energy [GeV]; Events", 40, 0, 20);
        qaHistograms["h_jet_energy"] = new TH1F("h_jet_energy", "Maximum 0.8x0.8 Energy Sum (EMCAL + HCAL) [GeV]; Events", 50, 0, 50);
        qaHistograms["h2_isoEtE"] = new TH2F("h2_isoEtE", "Cluster Isolation Energy vs Cluster Energy;Cluster E [GeV];E_{T}^{iso} [GeV]", 500, 20, 100, 500, -10, 10);

        qaHistogramsByTrigger[triggerIndex] = qaHistograms;

        // Create invariant mass histograms
        out->cd(invMassDir.c_str());
        std::map<std::string, TObject*> massHistograms;

        for (float maxAsym : asymmetry_values) {
            for (float maxChi2 : clus_chi_values) {
                for (float minClusE : clus_Ecore_values) {
                    std::string histName = "invMass_E" + formatFloatForFilename(minClusE) +
                                           "_Chi" + formatFloatForFilename(maxChi2) +
                                           "_Asym" + formatFloatForFilename(maxAsym);
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
    // Initialize counters
    positive_isoEt_count = 0;
    negative_isoEt_count = 0;
    skipped_tower_count = 0;
    towers_in_cone_count = 0;
    min_isoEt = std::numeric_limits<double>::max();
    max_isoEt = std::numeric_limits<double>::lowest();

    return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int caloTreeGen::InitRun(PHCompositeNode *topNode) {
  std::cout << "caloTreeGen::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

void caloTreeGen::collectTowerData(TowerInfoContainer* towerContainer,
                                   std::vector<TowerData>& towerDataList) {
    if (!towerContainer) {
        std::cout << "Error: Tower container is null." << std::endl;
        return;
    }

    unsigned int tower_range = towerContainer->size();
    towerDataList.reserve(tower_range);

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
}
void caloTreeGen::processEnergyMaps(const std::vector<float>* m_emcTowE, const std::vector<float>* m_emciEta, const std::vector<float>* m_emciPhi, const std::vector<float>* m_ohcTowE, const std::vector<float>* m_ohciTowEta, const std::vector<float>* m_ohciTowPhi, const std::vector<float>* m_ihcTowE, const std::vector<float>* m_ihciTowEta, const std::vector<float>* m_ihciTowPhi, std::vector<short>** m_emcal_good, std::vector<short>** m_ohc_good, std::vector<short>** m_ihc_good, float& max_8by8energy_emcal, float& max_energy_hcal, float& max_energy_jet) {
    
    float energymap[12][32] = {0};
    float energymap_emcal[12][35] = {0};
    float energymap_hcalin[12][35] = {0};
    float energymap_hcalout[12][35] = {0};
    float energymap_extend[12][35] = {0};
    
    float energymap_jet[9][32] = {0};
    float energymap_jet_emcal[9][32] = {0};
    float energymap_jet_hcalin[9][32] = {0};
    float energymap_jet_hcalout[9][32] = {0};
    
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

    for (size_t ie = 0; ie < m_emcTowE->size(); ie++) {
        float eta = static_cast<float>(m_emciEta->at(ie));
        float phi = static_cast<float>(m_emciPhi->at(ie));

        int ebin = static_cast<int>(eta / 8.0);
        int pbin = static_cast<int>(phi / 8.0);
        
        // Skip entries that are not marked as good
        if (!(*m_emcal_good)->at(ie)) {
            continue;
        }
        energymap[ebin][pbin] += m_emcTowE->at(ie);
        energymap_emcal[ebin][pbin] += m_emcTowE->at(ie);
        energymap_extend[ebin][pbin] += m_emcTowE->at(ie);
        
        if (pbin < 3) {
            energymap_emcal[ebin][pbin+32] += m_emcTowE->at(ie);
            energymap_extend[ebin][pbin+32] += m_emcTowE->at(ie);
        }
    }
    
    for (size_t ie = 0; ie < m_ihcTowE->size(); ie++) {
        float eta_ihc = static_cast<float>(m_ihciTowEta->at(ie));
        float phi_ihc = static_cast<float>(m_ihciTowPhi->at(ie));
        
        int ebin_ihc = static_cast<int>(eta_ihc / 8.0);
        int pbin_ihc = static_cast<int>(phi_ihc / 8.0);
        
        if (!(*m_ihc_good)->at(ie)) continue;
        energymap_extend[ebin_ihc][pbin_ihc] += m_ihcTowE->at(ie);
        energymap_hcalin[ebin_ihc][pbin_ihc] += m_ihcTowE->at(ie);
        if (pbin_ihc < 3) {
            energymap_hcalin[ebin_ihc][pbin_ihc+32] += m_ihcTowE->at(ie);
            energymap_hcalin[ebin_ihc][pbin_ihc+32] += m_ihcTowE->at(ie);
        }
    }
    for (size_t ie = 0; ie < m_ohcTowE->size(); ie++) {
        float eta_ohc = static_cast<float>(m_ohciTowEta->at(ie));
        float phi_ohc = static_cast<float>(m_ohciTowPhi->at(ie));
        
        int ebin_ohc = static_cast<int>(eta_ohc / 8.0);
        int pbin_ohc = static_cast<int>(phi_ohc / 8.0);

        if (!(*m_ohc_good)->at(ie)) continue;
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
    
//    int ebin = 0;
//    int pbin = 0;
//    int hcal_ebin = 0;
//    int hcal_pbin = 0;
//    int jet_ebin = 0;
//    int jet_pbin = 0;


    // Loop over eta bins and phi bins to find the maximum energy in EMCal
    for (int j = 0; j < 32; j++) {
        for (int k = 0; k < 12; k++) {
            if (k < 9) {
                if (energymap_jet[k][j] > max_energy_jet) {
                    max_energy_jet = energymap_jet[k][j];
//                    jet_ebin = k;
//                    jet_pbin = j;
                }
            }
            if (energymap[k][j] > max_8by8energy_emcal) {
                max_8by8energy_emcal = energymap[k][j]; // Update maximum EMCal energy
//                ebin = k;
//                pbin = j;
            }
            if (energymap_hcalout[k][j] + energymap_hcalin[k][j] > max_energy_hcal) {
                max_energy_hcal = energymap_hcalin[k][j] + energymap_hcalout[k][j];
//                hcal_ebin = k;
//                hcal_pbin = j;
            }
        }
    }
}



void caloTreeGen::calculateIsoEt(TowerInfoContainer* towerContainer,
                                 RawTowerGeomContainer* geomContainer,
                                 double& isoEt,
                                 const double clus_eta,
                                 const double clus_phi,
                                 const double m_vx,
                                 const double m_vy,
                                 const double m_vz,
                                 int& skipped_tower_count,
                                 int& towers_in_cone_count,
                                 const std::string& geomContainerName,
                                 double dR_cut) {
    std::vector<TowerData> towerDataList;
    collectTowerData(towerContainer, towerDataList);

    // Determine the correct CalorimeterId based on the geometry container name
    RawTowerDefs::CalorimeterId caloId;
    if (geomContainerName == "TOWERGEOM_CEMC") {
        caloId = RawTowerDefs::CalorimeterId::CEMC;
    } else if (geomContainerName == "TOWERGEOM_HCALIN") {
        caloId = RawTowerDefs::CalorimeterId::HCALIN;
    } else if (geomContainerName == "TOWERGEOM_HCALOUT") {
        caloId = RawTowerDefs::CalorimeterId::HCALOUT;
    } else {
        std::cout << "Error: Unknown geometry container name: " << geomContainerName << std::endl;
        return;
    }

    for (const auto& data : towerDataList) {
        if (!data.isAcceptable) {
            skipped_tower_count++;
            continue;
        }

        RawTowerGeom* tower_geom = geomContainer->get_tower_geometry(
            RawTowerDefs::encode_towerid(caloId, data.ieta, data.iphi)
        );

        if (!tower_geom) {
            std::cout << "Warning: Tower geometry is null." << std::endl;
            continue;
        }

        double tower_phi = tower_geom->get_phi();
        double tower_eta = getTowerEta(tower_geom, m_vx, m_vy, m_vz);
        double dR = deltaR(clus_eta, tower_eta, clus_phi, tower_phi);

        if (dR < dR_cut) {
            isoEt += data.energy / cosh(tower_eta);
            towers_in_cone_count++;
        }
    }
}


//____________________________________________________________________________..
int caloTreeGen::process_event(PHCompositeNode *topNode) {
    event_count++;
    std::cout << "\n========== Processing Event " << event_count << " ==========\n";
    
    
    // Initialize vertex coordinates to zero early
    m_vertex = -9999;
    m_vx = m_vy = m_vz = 0;

    // Fetch the GlobalVertexMap node
    GlobalVertexMap* vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
    if (vertexmap && !vertexmap->empty()) {
        GlobalVertex* vtx = vertexmap->begin()->second;  // Fetch the first vertex
        if (vtx) {
            m_vx = vtx->get_x();
            m_vy = vtx->get_y();
            m_vz = vtx->get_z();
            m_vertex = vtx->get_z();
        } else {
            m_vx = 0;
            m_vy = 0;
            m_vz = 0;
        }
    } else {
        m_vx = 0;
        m_vy = 0;
        m_vz = 0;
    }
    
    _gl1_packet = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    if (_gl1_packet) {
        b_gl1_scaledvec = _gl1_packet->lValue(0, "ScaledVector");
    }
    std::vector<int> activeTriggerBits = extractTriggerBits(b_gl1_scaledvec, event_count);
    std::cout << "We are now looping over trigger bits: ";
    for (int bit : activeTriggerBits) {
        std::cout << bit << " ";
    }
    std::cout << std::endl;

    for (int triggerIndex : triggerIndices) {
        // Check if the current trigger bit is active
        if (!checkTriggerCondition(activeTriggerBits, triggerIndex)) {
            continue;  // Skip if this trigger bit is not active
        }
        std::cout << "Processing Trigger Bit: " << triggerIndex << std::endl;
        auto& qaHistograms = qaHistogramsByTrigger[triggerIndex];
        if (fabs(m_vertex) >= 30) {
            continue;
        }
        ((TH1F*)qaHistograms["hVtxZ"])->Fill(m_vertex);

        
        TowerInfoContainer* emcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
        TowerInfoContainer* ohcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
        TowerInfoContainer* ihcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");

        // Declare geometry container pointers
        RawTowerGeomContainer* geomEM = nullptr;
        RawTowerGeomContainer* geomIH = nullptr;
        RawTowerGeomContainer* geomOH = nullptr;

        // Fetch geometry containers
        std::vector<std::tuple<RawTowerGeomContainer**, const char*, const char*>> geomContainers = {
            {&geomEM, "TOWERGEOM_CEMC", "EMC"},
            {&geomIH, "TOWERGEOM_HCALIN", "Inner HCal"},
            {&geomOH, "TOWERGEOM_HCALOUT", "Outer HCal"}
        };

        for (auto& [geomContainer, geomName, geomLabel] : geomContainers) {
            *geomContainer = findNode::getClass<RawTowerGeomContainer>(topNode, geomName);
            if (!*geomContainer) {
                std::cout << "Error: Missing " << geomLabel << " Geometry Container: " << geomName << std::endl;
                return Fun4AllReturnCodes::ABORTEVENT;
            }
        }
        
        // Process towers and fill histograms
        if (emcTowerContainer) {
            std::cout << "Processing EMCal Towers..." << std::endl;
            processTowers(emcTowerContainer, totalCaloEEMCal, m_emciEta, m_emciPhi, m_emcTowE, m_emcTime, m_emcChi2, m_emcPed, m_emcal_good);
            for (size_t i = 0; i < m_emcTowE.size(); ++i) {
                ((TH2F*)qaHistograms["h2_EMCal_TowerEtaPhi_2D"])->Fill(m_emciEta[i], m_emciPhi[i], m_emcTowE[i]);
            }
            ((TH1F*)qaHistograms["hTotalCaloEEMCal"])->Fill(totalCaloEEMCal);
            for (const auto& chi2 : m_emcChi2) {
                ((TH1F*)qaHistograms["h_emcalChi2"])->Fill(chi2);
            }
        }
        if (ihcTowerContainer) {
            std::cout << "Processing IHCal Towers..." << std::endl;
            processTowers(ihcTowerContainer, totalCaloEIHCal, m_ihciTowEta, m_ihciTowPhi, m_ihcTowE, m_ihcTime, m_ihcChi2, m_ihcPed, m_ihc_good);
            for (size_t i = 0; i < m_ihcTowE.size(); ++i) {
                ((TH2F*)qaHistograms["h2_IHCal_TowerEnergy"])->Fill(m_ihciTowEta[i], m_ihciTowPhi[i], m_ihcTowE[i]);
            }
            ((TH1F*)qaHistograms["hTotalCaloEIHCal"])->Fill(totalCaloEIHCal);
            for (const auto& chi2 : m_ihcChi2) {
                ((TH1F*)qaHistograms["h_ihcalChi2"])->Fill(chi2);
            }
        }
        if (ohcTowerContainer) {
            std::cout << "Processing OHCal Towers..." << std::endl;
            processTowers(ohcTowerContainer, totalCaloEOHCal, m_ohciTowEta, m_ohciTowPhi, m_ohcTowE, m_ohcTime, m_ohcChi2, m_ohcPed, m_ohc_good);
            for (size_t i = 0; i < m_ohcTowE.size(); ++i) {
                ((TH2F*)qaHistograms["h2_OHCal_TowerEnergy"])->Fill(m_ohciTowEta[i], m_ohciTowPhi[i], m_ohcTowE[i]);
            }
            ((TH1F*)qaHistograms["hTotalCaloEOHCal"])->Fill(totalCaloEOHCal);
            for (const auto& chi2 : m_ohcChi2) {
                ((TH1F*)qaHistograms["h_ohcalChi2"])->Fill(chi2);
            }
        }
        

        RawClusterContainer *clusterContainer = findNode::getClass<RawClusterContainer>(topNode, "CLUSTERINFO_CEMC");
        if(!clusterContainer) {
            std::cout << PHWHERE << "Cluster node is missing. Output related to this node will be empty" << std::endl;
            return 0;
        }
        
        RawClusterContainer::ConstRange clusterRange = clusterContainer->getClusters();
        float max_energy_clus = 0.0;
        
        for (auto clusterIter = clusterRange.first; clusterIter != clusterRange.second; ++clusterIter) {
            RawCluster* cluster = clusterIter->second;
            if (!cluster) {
                std::cout << "Warning: Null cluster found." << std::endl;
                continue;
            }
            // Use vertex coordinates for cluster calculations
            // CLHEP::Hep3Vector vertex(0,0,0); --before in caloTreeGen -- IT WILL DO THIS if vertex map missing
            // Use vertex coordinates for cluster calculations
            CLHEP::Hep3Vector vertex(m_vx, m_vy, m_vz);
            if (m_vertex != -9999) {
                vertex.setZ(m_vertex);
            }
            CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*cluster, vertex);
            CLHEP::Hep3Vector E_vec_cluster_Full = RawClusterUtility::GetEVec(*cluster, vertex);

            float clusE = E_vec_cluster_Full.mag(); //only vartiable that uses not GetECoreVec--NOT ECORE
            float clusEcore = E_vec_cluster.mag();
            
            if (clusEcore > max_energy_clus) {
                max_energy_clus = clusEcore; // Update the maximum cluster energy
            }
            
            float clus_eta = E_vec_cluster.pseudoRapidity();
            float clus_phi = E_vec_cluster.phi();
            float clus_pt = E_vec_cluster.perp();
            float clus_chi = cluster -> get_chi2();
            float nTowers = cluster ->getNTowers();
            float maxTowerEnergy = getMaxTowerE(cluster,emcTowerContainer);

            m_clusterE.push_back(clusE);
            m_clusterECore.push_back(clusEcore);
            m_clusterPhi.push_back(clus_phi);
            m_clusterEta.push_back(clus_eta);
            m_clusterPt.push_back(clus_pt);
            m_clusterChi.push_back(clus_chi);
            m_clusterNtow.push_back(nTowers);
            m_clusterTowMaxE.push_back(maxTowerEnergy);
            m_clusTowPhi.push_back(returnClusterTowPhi(cluster,emcTowerContainer));
            m_clusTowEta.push_back(returnClusterTowEta(cluster,emcTowerContainer));
            m_clusTowE.push_back(returnClusterTowE(cluster,emcTowerContainer));
            
            ((TH1F*)qaHistograms["hClusterPt"])->Fill(clus_pt);
            ((TH1F*)qaHistograms["hClusterChi2"])->Fill(clus_chi);
            ((TH1F*)qaHistograms["hClusterEcore"])->Fill(clusEcore);
            
            if (clusEcore < 1) { //cut on ecore
                continue;
            }
            double et = clusEcore / cosh(clus_eta);
            double isoEt = 0;
            if (et < 2) { //cut on transverse energy
                continue;
            }
    

            for (double dR_cut : dR_values) {
                std::string dR_cut_str = formatFloatForFilename(dR_cut);
                
                calculateIsoEt(emcTowerContainer, geomEM, isoEt, clus_eta, clus_phi, m_vx, m_vy, m_vz,
                               skipped_tower_count, towers_in_cone_count, "TOWERGEOM_CEMC", dR_cut);
                calculateIsoEt(ihcTowerContainer, geomIH, isoEt, clus_eta, clus_phi, m_vx, m_vy, m_vz,
                               skipped_tower_count, towers_in_cone_count, "TOWERGEOM_HCALIN", dR_cut);
                calculateIsoEt(ohcTowerContainer, geomOH, isoEt, clus_eta, clus_phi, m_vx, m_vy, m_vz,
                               skipped_tower_count, towers_in_cone_count, "TOWERGEOM_HCALOUT", dR_cut);
                
                // Subtract cluster energy from isolation energy
                isoEt -= et;

                if (isoEt < 0) {
                    negative_isoEt_count++;
                } else {
                    positive_isoEt_count++;
                }
                // Update the min and max isoEt values
                if (isoEt < min_isoEt) {
                    min_isoEt = isoEt;
                }
                if (isoEt > max_isoEt) {
                    max_isoEt = isoEt;
                }

                m_clusterEtIso.push_back(isoEt);
                ((TH1F*)qaHistograms["h2_isoEtE_" + dR_cut_str])->Fill(clusEcore, isoEt);
                
            }
        }
        ((TH1F*)qaHistograms["hCluster_maxECore"])->Fill(max_energy_clus);
    }
    return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int caloTreeGen::ResetEvent(PHCompositeNode *topNode) {
  m_clusterE.clear();
  m_clusterPhi.clear();
  m_clusterEta.clear();
  m_clusterPt.clear();
  m_clusterChi.clear();
  m_clusterTowMaxE.clear();
  m_clusterNtow.clear();
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
}

//____________________________________________________________________________..
int caloTreeGen::EndRun(const int runnumber) {
  std::cout << "caloTreeGen::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int caloTreeGen::End(PHCompositeNode *topNode) {
    std::cout << "caloTreeGen::End(PHCompositeNode *topNode) All events have been processed. Beginning final analysis steps......" << std::endl;
    for (int triggerIndex : triggerIndices) {
        auto& qaHistograms = qaHistogramsByTrigger[triggerIndex];
        
        float max_8by8energy_emcal = 0.0;
        float max_energy_hcal = 0.0;
        float max_energy_jet = 0.0;

        // Create pointers to the vectors
        std::vector<short>* emcal_good_ptr = &m_emcal_good;
        std::vector<short>* ohc_good_ptr = &m_ohc_good;
        std::vector<short>* ihc_good_ptr = &m_ihc_good;

        // Pass the addresses of the pointers to the function
        processEnergyMaps(&m_emcTowE, &m_emciEta, &m_emciPhi, &m_ohcTowE, &m_ohciTowEta, &m_ohciTowPhi,
                          &m_ihcTowE, &m_ihciTowEta, &m_ihciTowPhi,
                          &emcal_good_ptr, &ohc_good_ptr, &ihc_good_ptr,
                          max_8by8energy_emcal, max_energy_hcal, max_energy_jet);



        std::cout << "Finished processing energy maps for Trigger Index: " << triggerIndex << std::endl;
        
        ((TH1F*)qaHistograms["h8by8TowerEnergySum"])->Fill(max_8by8energy_emcal);
        ((TH1F*)qaHistograms["h_hcal_energy"])->Fill(max_energy_hcal);
        ((TH1F*)qaHistograms["h_jet_energy"])->Fill(max_energy_jet);
        
        
        std::cout << "Processing Invariant Mass Distributions for Trigger Index: " << triggerIndex << std::endl;
        
        // Cache cluster properties to avoid repeated access to the vector
        size_t nClusters = m_clusterECore.size();
        std::vector<float> cachedPt(nClusters), cachedEcore(nClusters), cachedChi2(nClusters);
        for (size_t clus = 0; clus < nClusters; ++clus) {
            cachedPt[clus] = m_clusterPt.at(clus);
            cachedEcore[clus] = m_clusterECore.at(clus);
            cachedChi2[clus] = m_clusterChi.at(clus);
        }

        // Fill invariant mass histograms for this trigger
        auto& massHistograms = massHistogramsByTrigger[triggerIndex];
        for (size_t clus1 = 0; clus1 < nClusters; ++clus1) {
            for (size_t clus2 = clus1 + 1; clus2 < nClusters; ++clus2) {
                float pt1 = cachedPt[clus1], pt2 = cachedPt[clus2];
                float E1 = cachedEcore[clus1], E2 = cachedEcore[clus2];
                float chi1 = cachedChi2[clus1], chi2 = cachedChi2[clus2];

                if (pt1 < 2 || pt1 >= 10 || pt2 < 2 || pt2 >= 10) {
                    continue;
                }

                TLorentzVector photon1, photon2;
                photon1.SetPtEtaPhiE(pt1, m_clusterEta.at(clus1), m_clusterPhi.at(clus1), E1);
                photon2.SetPtEtaPhiE(pt2, m_clusterEta.at(clus2), m_clusterPhi.at(clus2), E2);
                TLorentzVector meson = photon1 + photon2;
                float mesonMass = meson.M();
                float asym = fabs(E1 - E2) / (E1 + E2);

                // Apply cuts and fill histograms
                for (float maxAsym : asymmetry_values) {
                    for (float maxChi2 : clus_chi_values) {
                        for (float minClusEcore : clus_Ecore_values) {
                            std::cout << "Trigger Index: " << triggerIndex
                                      << ", Cut Values -> Asymmetry: " << maxAsym
                                      << ", Chi2: " << maxChi2
                                      << ", Ecore: " << minClusEcore << std::endl;

                            if (asym >= maxAsym || chi1 >= maxChi2 || chi2 >= maxChi2 || E1 < minClusEcore || E2 < minClusEcore) {
                                continue;
                            }
                            std::string histName = "invMass_E" + formatFloatForFilename(minClusEcore) +
                                                   "_Chi" + formatFloatForFilename(maxChi2) +
                                                   "_Asym" + formatFloatForFilename(maxAsym);
                            // When filling the histogram, cast to TH1F*
                            ((TH1F*)massHistograms[histName])->Fill(mesonMass);
                        }
                    }
                }
            }
        }
    }

    
    // Write all histograms to file and clean up
    for (const auto& [triggerIndex, qaHistograms] : qaHistogramsByTrigger) {
        out->cd(("QA/Trigger" + std::to_string(triggerIndex)).c_str());
        for (const auto& [name, hist] : qaHistograms) {
            hist->Write();
            delete hist;
        }
    }

    for (const auto& [triggerIndex, massHistograms] : massHistogramsByTrigger) {
        out->cd(("InvariantMassDistributions/Trigger" + std::to_string(triggerIndex)).c_str());
        for (const auto& [name, hist] : massHistograms) {
            hist->Write();
            delete hist;
        }
    }

    out -> Close();
    delete out;

    // At the end of the event loop, print out the summary
    std::cout << "Summary of Isolation Energy Calculation:" << std::endl;
    std::cout << "Positive isolation energy events: " << positive_isoEt_count << std::endl;
    std::cout << "Negative isolation energy events: " << negative_isoEt_count << std::endl;
    std::cout << "Skipped towers due to IsAcceptableTower flag: " << skipped_tower_count << std::endl;
    std::cout << "Towers contributing to isolation energy (inside cone): " << towers_in_cone_count << std::endl;
    std::cout << "Minimum isolation energy: " << min_isoEt << std::endl;
    std::cout << "Maximum isolation energy: " << max_isoEt << std::endl;

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
