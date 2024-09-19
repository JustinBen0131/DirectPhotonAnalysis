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
  ,T(nullptr)
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
int caloTreeGen::Init(PHCompositeNode *topNode)
{
  
  out = new TFile(Outfile.c_str(),"RECREATE");

  
  T = new TTree("T","T");
  
  //Electromagnetic Calorimeter
  if(storeEMCal)
    {
      T -> Branch("emcTowE",&m_emcTowE);
      T -> Branch("emcTowiEta",&m_emciEta);
      T -> Branch("emcTowiPhi",&m_emciPhi);
      T -> Branch("emcTime",&m_emcTime);
      T -> Branch("emcChi2",&m_emcChi2);
      T -> Branch("emcPed",&m_emcPed);
      T -> Branch("emcal_good",&m_emcal_good);
  
      //EMCal Cluster information
      if(storeEMCal && storeClusters) {
      T -> Branch("clusterE",&m_clusterE);
      T -> Branch("clusterPhi",&m_clusterPhi);
      T -> Branch("clusterEta", &m_clusterEta);
      T -> Branch("clusterPt", &m_clusterPt);
      T -> Branch("clusterChi2", &m_clusterChi);
      T -> Branch("clusterNtow",&m_clusterNtow);
      T -> Branch("clusterTowMaxE",&m_clusterTowMaxE);
      T -> Branch("clusterECore",&m_clusterECore);

      //Information for towers within clusters
      //Enabled by setting "DoFineClusters" in the macro
      if(storeEMCal && storeClusters && storeClusterDetails)
        {
          T -> Branch("clusTowPhi","vector<vector<int> >",&m_clusTowPhi);
          T -> Branch("clusTowEta","vector<vector<int> >",&m_clusTowEta);
          T -> Branch("clusTowE","vector<vector<float> >",&m_clusTowE);
        }
    }
  }
  //Outer Hadronic Calorimeter
  if(storeHCals)
    {
      T -> Branch("ohcTowE",&m_ohcTowE);
      T -> Branch("ohcTowiEta",&m_ohciTowEta);
      T -> Branch("ohcTowiPhi",&m_ohciTowPhi);
      T -> Branch("ohcTime",&m_ohcTime);
      T -> Branch("ohcChi2",&m_ohcChi2);
      T -> Branch("ohcPed",&m_ohcPed);
      T -> Branch("ohc_good",&m_ohc_good);
  
      //Inner Hadronic Calorimeter
      T -> Branch("ihcTowE",&m_ihcTowE);
      T -> Branch("ihcTowiEta",&m_ihciTowEta);
      T -> Branch("ihcTowiPhi",&m_ihciTowPhi);
      T -> Branch("ihcTime",&m_ihcTime);
      T -> Branch("ihcChi2",&m_ihcChi2);
      T -> Branch("ihcPed",&m_ihcPed);
      T -> Branch("ihc_good",&m_ihc_good);
        
    }
  //ZDC information
  if(storeZDC)
    {
      T -> Branch("zdcTowE",&m_zdcTowE);
      T -> Branch("zdcTowside",&m_zdcSide);
  
      //SMD information
      T -> Branch("smdE",&m_smdE);
      T -> Branch("smdSide",&m_smdSide);
    }
    

  T->Branch("clusterEtIso", &m_clusterEtIso);


  //Total
  T -> Branch("totalCaloEEMCal",&totalCaloEEMCal);
  T -> Branch("totalCaloEOHCal",&totalCaloEOHCal);
  T -> Branch("totalCaloEIHCal",&totalCaloEIHCal);
  T -> Branch("totalCaloEZDC",&totalCaloEZDC);
  T -> Branch("zvertex",&m_vertex);
  

  T->Branch("gl1_clock",&b_gl1_clock, "gl1_clock/l");
  T->Branch("gl1_scaled",b_gl1_scaled, "gl1_scaled[64]/l");
  T->Branch("gl1_live",b_gl1_live, "gl1_live[64]/l");
  T->Branch("gl1_raw",b_gl1_raw, "gl1_raw[64]/l");
  T->Branch("gl1_rawvec",&b_gl1_rawvec, "gl1_rawvec/l");
  T->Branch("gl1_livevec",&b_gl1_livevec, "gl1_livevec/l");
  T->Branch("gl1_scaledvec",&b_gl1_scaledvec, "gl1_scaledvec/l");
        
  

  
  
  zVertex = new TH1F("zVertex","zVertex",200,-100,100);

 //so that the histos actually get written out
  Fun4AllServer *se = Fun4AllServer::instance();
  se -> Print("NODETREE");
  
  std::cout << "caloTreeGen::Init(PHCompositeNode *topNode) Initializing" << std::endl;
    
  // Initialize counters
  positive_isoEt_count = 0;
  negative_isoEt_count = 0;
  skipped_tower_count = 0;
  towers_in_cone_count = 0;

    
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int caloTreeGen::InitRun(PHCompositeNode *topNode) {
  std::cout << "caloTreeGen::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
// Helper structure to hold tower data
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

// Helper function to collect tower data
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
                                 const std::string& geomContainerName) {
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

        if (dR < 0.3) {
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
            zVertex->Fill(m_vertex);
        } else {
//            std::cout << "GlobalVertex pointer is null, setting vertex to (0,0,0)" << std::endl;
            m_vx = 0;
            m_vy = 0;
            m_vz = 0;
        }
    } else {
//        std::cout << "GlobalVertexMap node is missing or empty, setting vertex to (0,0,0)" << std::endl;
        m_vx = 0;
        m_vy = 0;
        m_vz = 0;
    }


    // Fetch tower containers and geometry once
    std::cout << "Fetching tower containers and geometry..." << std::endl;

    // Declare tower container pointers
    TowerInfoContainer* emcTowerContainer = nullptr;
    TowerInfoContainer* ohcTowerContainer = nullptr;
    TowerInfoContainer* ihcTowerContainer = nullptr;
    TowerInfoContainer* zdcTowerContainer = nullptr;

    // Fetch tower containers
    std::vector<std::pair<TowerInfoContainer**, std::string>> towerContainers = {
        {&emcTowerContainer, m_emcTowerNode},
        {&ohcTowerContainer, m_ohcTowerNode},
        {&ihcTowerContainer, m_ihcTowerNode},
        {&zdcTowerContainer, m_zdcTowerNode}
    };

    for (auto& [container, nodeName] : towerContainers) {
        *container = findNode::getClass<TowerInfoContainer>(topNode, nodeName.c_str());
        if (!*container) {
            std::cout << "Error: Missing Tower Container: " << nodeName << std::endl;
            return Fun4AllReturnCodes::ABORTEVENT;
        }
    }

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
        } else {
            std::cout << "Loaded " << geomLabel << " Geometry Container: " << geomName << std::endl;
        }
    }
    
    if (storeEMCal && emcTowerContainer) {
        processTowers(emcTowerContainer, totalCaloEEMCal, m_emciEta, m_emciPhi, m_emcTowE, m_emcTime, m_emcChi2, m_emcPed, m_emcal_good);
    }

    // Process Inner and Outer HCal towers
    if (storeHCals) {
        if (ohcTowerContainer) {
            processTowers(ohcTowerContainer, totalCaloEOHCal, m_ohciTowEta, m_ohciTowPhi, m_ohcTowE, m_ohcTime, m_ohcChi2, m_ohcPed, m_ohc_good);
        }
        if (ihcTowerContainer) {
            processTowers(ihcTowerContainer, totalCaloEIHCal, m_ihciTowEta, m_ihciTowPhi, m_ihcTowE, m_ihcTime, m_ihcChi2, m_ihcPed, m_ihc_good);
        }
    }
   
    if(storeClusters && storeEMCal && storeHCals) {
        RawClusterContainer *clusterContainer = findNode::getClass<RawClusterContainer>(topNode, m_clusterNode.c_str());
        if(!clusterContainer && storeClusters) {
            std::cout << PHWHERE << "caloTreeGen::process_event: "<<  m_clusterNode << " node is missing. Output related to this node will be empty" << std::endl;
            return 0;
        }
        
        RawClusterContainer::ConstRange clusterRange = clusterContainer->getClusters();
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

            float clusE = E_vec_cluster_Full.mag(); //only vartiable that uses not GetECoreVec
            
            float clusEcore = E_vec_cluster.mag();
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
            
            if(storeClusterDetails) {
                m_clusTowPhi.push_back(returnClusterTowPhi(cluster,emcTowerContainer));
                m_clusTowEta.push_back(returnClusterTowEta(cluster,emcTowerContainer));
                m_clusTowE.push_back(returnClusterTowE(cluster,emcTowerContainer));
            }
            // Calculate isolation energy
            double et = clusEcore / cosh(clus_eta);
            double isoEt = 0;
            if (et < 0) {
                std::cout << "Skipping cluster with et < 0 in event " << event_count << std::endl;
                continue;  // Skip clusters below eT cut of 0.5 GeV
            }

            // Example call with enhanced logging
            calculateIsoEt(emcTowerContainer, geomEM, isoEt, clus_eta, clus_phi, m_vx, m_vy, m_vz,
                           skipped_tower_count, towers_in_cone_count, "TOWERGEOM_CEMC");
            calculateIsoEt(ihcTowerContainer, geomIH, isoEt, clus_eta, clus_phi, m_vx, m_vy, m_vz,
                           skipped_tower_count, towers_in_cone_count, "TOWERGEOM_HCALIN");
            calculateIsoEt(ohcTowerContainer, geomOH, isoEt, clus_eta, clus_phi, m_vx, m_vy, m_vz,
                           skipped_tower_count, towers_in_cone_count, "TOWERGEOM_HCALOUT");

            // Subtract cluster energy from isolation energy
            isoEt -= et;

            if (isoEt < 0) {
                negative_isoEt_count++;
            } else {
                positive_isoEt_count++;
            }

            m_clusterEtIso.push_back(isoEt);
        }
    }
  

    

    if(storeZDC && zdcTowerContainer) {
        unsigned int tower_range = zdcTowerContainer->size();
        totalCaloEZDC = 0;
        // Pre-allocate vector memory
        m_zdcTowE.reserve(16);
        m_zdcSide.reserve(16);
        m_smdE.reserve(32);
        m_smdSide.reserve(32);
        
        for(unsigned int iter = 0; iter < tower_range; iter++) {
            TowerInfo* tower = zdcTowerContainer->get_tower_at_channel(iter);
            if (!tower) continue;

            if(iter < 16) {
                float energy = tower->get_energy();
                unsigned int towerkey = zdcTowerContainer->encode_key(iter);
                unsigned int side = TowerInfoDefs::get_zdc_side(towerkey);
                
                totalCaloEZDC += energy;
                m_zdcTowE.push_back(energy);
                m_zdcSide.push_back(side);
            }
            if(iter > 15 && iter < 48) {
                //smd north stuff
                float energy = tower->get_energy();
                unsigned int towerkey = zdcTowerContainer->encode_key(iter);
                unsigned int side = TowerInfoDefs::get_zdc_side(towerkey);
                
                m_smdE.push_back(energy);
                m_smdSide.push_back(side);
            }
         }
      }
  
    _gl1_packet = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    if (_gl1_packet) {
        b_gl1_clock = _gl1_packet->lValue(0, "BCO");
        b_gl1_rawvec = _gl1_packet->lValue(0, "TriggerInput");
        b_gl1_livevec = _gl1_packet->lValue(0, "TriggerVector");
        b_gl1_scaledvec = _gl1_packet->lValue(0, "ScaledVector");
        
        for (int i = 0; i < 64; i++) {
          b_gl1_scaled[i] = _gl1_packet->lValue(i, 2);
          b_gl1_raw[i] = _gl1_packet->lValue(i, 0);
          b_gl1_live[i] = _gl1_packet->lValue(i, 1);
        }
    }
    T -> Fill();
  
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
  m_clusterEtIso.clear(); //to use isolation energy
    

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
  std::cout << "caloTreeGen::End(PHCompositeNode *topNode) This is the End..." << std::endl;

  out -> cd();
  T -> Write();
  zVertex -> Write();
  out -> Close();
  delete out;
    
    
  // At the end of the event loop, print out the summary
  std::cout << "Summary of Isolation Energy Calculation:" << std::endl;
  std::cout << "Positive isolation energy events: " << positive_isoEt_count << std::endl;
  std::cout << "Negative isolation energy events: " << negative_isoEt_count << std::endl;
  std::cout << "Skipped towers due to IsAcceptableTower flag: " << skipped_tower_count << std::endl;
  std::cout << "Towers contributing to isolation energy (inside cone): " << towers_in_cone_count << std::endl;

  //hm -> dumpHistos(Outfile.c_str(), "UPDATE");
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
std::vector<int> caloTreeGen::returnClusterTowEta(RawCluster *cluster, TowerInfoContainer *towerContainer)
{
  RawCluster::TowerConstRange towers = cluster -> get_towers();
  RawCluster::TowerConstIterator toweriter;
  
  std::vector<int> towerIDsEta;
  for(toweriter = towers.first; toweriter != towers.second; toweriter++) towerIDsEta.push_back(RawTowerDefs::decode_index1(toweriter -> first));

  return towerIDsEta;
}
//____________________________________________________________________________..
std::vector<int> caloTreeGen::returnClusterTowPhi(RawCluster *cluster, TowerInfoContainer *towerContainer)
{
  RawCluster::TowerConstRange towers = cluster -> get_towers();
  RawCluster::TowerConstIterator toweriter;
  
  std::vector<int> towerIDsPhi;
  for(toweriter = towers.first; toweriter != towers.second; toweriter++) towerIDsPhi.push_back(RawTowerDefs::decode_index2(toweriter -> first));
  return towerIDsPhi;
}
//____________________________________________________________________________..
std::vector<float> caloTreeGen::returnClusterTowE(RawCluster *cluster, TowerInfoContainer *towerContainer)
{
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
