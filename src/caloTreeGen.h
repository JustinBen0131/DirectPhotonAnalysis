// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOTREEGEN_H
#define CALOTREEGEN_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>
#include <TTree.h>
#include <ffarawobjects/Gl1Packet.h>
#include <ffarawobjects/Gl1Packetv1.h>
#include <ffarawobjects/Gl1Packetv2.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoDefs.h>
#include <calobase/RawTowerGeomContainer.h>


//for the vertex
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>



class TTree;
class PHCompositeNode;
class Fun4AllHistoManager;
class TFile;
class RawCluster;
class TowerInfoContainer;
class TH1F;

class caloTreeGen : public SubsysReco
{
 public:

  caloTreeGen(const std::string &name = "caloTreeGen");

  ~caloTreeGen() override;

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
  */
  int Init(PHCompositeNode *topNode) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
  */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
  */
  int process_event(PHCompositeNode *topNode) override;

  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode *topNode) override;

  /// Called at the end of each run.
  int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

  void doClusters(int clusters, std::string clusterNode)  {storeClusters = clusters; m_clusterNode = clusterNode;}

  void doClusterDetails(int fineCluster) {storeClusterDetails = fineCluster;}
  
  void doEMCal(int emcalOn, std::string emcNode) {storeEMCal = emcalOn; m_emcTowerNode = emcNode;}

  void doHCals(int hcalsOn, std::string ohcNode, std::string ihcNode) {storeHCals = hcalsOn; m_ohcTowerNode = ohcNode; m_ihcTowerNode = ihcNode;}

  void doZDC(int zdcOn, std::string zdcNode) {storeZDC = zdcOn; m_zdcTowerNode = zdcNode;}
  


 private:

  TTree *T;
    // Declare counters to track across events
  int positive_isoEt_count;
  int negative_isoEt_count;
  int skipped_tower_count;
  int towers_in_cone_count;
  int event_count = 0;
    
  //EMCal
  std::vector<float> m_emcTowE;
  std::vector<float> m_emciEta;
  std::vector<float> m_emciPhi;
  std::vector<int> m_emcTime;
  std::vector<float> m_emcChi2;
  std::vector<float> m_emcPed;
  std::vector<short> m_emcal_good;
  
  //OHCal
  std::vector<float> m_ohciTowPhi;
  std::vector<float> m_ohciTowEta;
  std::vector<float> m_ohcTowE;
  std::vector<int> m_ohcTime;
  std::vector<float> m_ohcChi2;
  std::vector<float> m_ohcPed;
  std::vector<short> m_ohc_good;
    
  //IHCal
  std::vector<float> m_ihciTowPhi;
  std::vector<float> m_ihciTowEta;
  std::vector<float> m_ihcTowE;
  std::vector<int> m_ihcTime;
  std::vector<float> m_ihcChi2;
  std::vector<float> m_ihcPed;
  std::vector<short> m_ihc_good;

  //ZDC
  std::vector<float> m_zdcTowE;
  std::vector<int> m_zdcSide;
  
  //SMD
  std::vector<float> m_smdE;
  std::vector<int> m_smdSide;
  
  //Clusters
  std::vector<float> m_clusterE;
  std::vector<float> m_clusterPhi;
  std::vector<float> m_clusterEta;
  std::vector<float> m_clusterPt;
  std::vector<float> m_clusterChi;
  std::vector<float> m_clusterNtow;
  std::vector<float> m_clusterTowMaxE;
  std::vector<float> m_clusterECore;
  std::vector<float> m_clusterEtIso;
  
  std::vector<std::vector<int> > m_clusTowEta;
  std::vector<std::vector<int> > m_clusTowPhi;
  std::vector<std::vector<float> > m_clusTowE;

  //GL1 information
  Gl1Packet *_gl1_packet;
  uint64_t b_gl1_rawvec;
  uint64_t b_gl1_livevec;
  uint64_t b_gl1_scaledvec;
  uint64_t b_gl1_clock;
  uint64_t b_gl1_raw[64];
  uint64_t b_gl1_live[64];
  uint64_t b_gl1_scaled[64];

  float m_vertex;
  double m_vx, m_vy, m_vz;
    
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

  TFile *out;
  //Fun4AllHistoManager *hm = nullptr;
  std::string Outfile = "commissioning.root";

  TH1F *zVertex;
  void collectTowerData(TowerInfoContainer* towerContainer, std::vector<TowerData>& towerDataList);
  void processTowers(TowerInfoContainer* towerContainer, float& totalCaloE, std::vector<float>& towEta, std::vector<float>& towPhi, std::vector<float>& towE, std::vector<int>& towTime, std::vector<float>& towChi2, std::vector<float>& towPed, std::vector<short>& towGood);
    
  void calculateIsoEt(TowerInfoContainer* towerContainer, RawTowerGeomContainer* geomContainer, double& isoEt, const double clus_eta, const double clus_phi, const double m_vx, const double m_vy, const double m_vz, int& skipped_tower_count, int& towers_in_cone_count, const std::string& geomContainerName);

  float getMaxTowerE(RawCluster *cluster, TowerInfoContainer *towerContainer);
  std::vector<float> returnClusterTowE(RawCluster *cluster, TowerInfoContainer *towerContainer);
  std::vector<int> returnClusterTowPhi(RawCluster *cluster, TowerInfoContainer *towerContainer);
  std::vector<int> returnClusterTowEta(RawCluster *cluster, TowerInfoContainer *towerContainer);
    

  
  float totalCaloEEMCal;
  float totalCaloEOHCal;
  float totalCaloEIHCal;
  float totalCaloEZDC;
  float totalChargeMBD;

  int storeClusters = 1;
  int storeClusterDetails = 1;
  int storeEMCal = 1;
  int storeHCals = 1;
  int storeZDC = 1;

  
  std::string m_emcTowerNode;
  std::string m_ohcTowerNode;
  std::string m_ihcTowerNode;
  std::string m_zdcTowerNode;
  std::string m_clusterNode;
  
  // Inline deltaR function for calculating distance between points in η-φ space
  inline float deltaR(float eta1, float eta2, float phi1, float phi2) {
    float deta = eta1 - eta2;
    float dphi = phi1 - phi2;
    if (dphi > M_PI) dphi -= 2 * M_PI;
    if (dphi < -M_PI) dphi += 2 * M_PI;
    return sqrt(deta * deta + dphi * dphi);
  }
  double getTowerEta(RawTowerGeom* tower_geom, double vx, double vy, double vz);
  bool IsAcceptableTower(TowerInfo* tower);

};

#endif  // CALOTREEGEN_H
