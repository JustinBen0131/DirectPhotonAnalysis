// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOTREEGEN_H
#define CALOTREEGEN_H

#include <fun4all/SubsysReco.h>
#include <calotrigger/TriggerAnalyzer.h>

#include <string>
#include <chrono>
#include <vector>
#include <sstream>
#include <iomanip>
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



class PHCompositeNode;
class Fun4AllHistoManager;
class TFile;
class RawCluster;
class TowerInfoContainer;
class TH1F;
class TH2F;

class caloTreeGen : public SubsysReco{
public:
    
    // A constructor that takes two std::string arguments
    caloTreeGen(const std::string &dataOutFile = "caloTreeData.root",
                const std::string &simOutFile  = "caloTreeSim.root");
    
    ~caloTreeGen() override;

    int Init(PHCompositeNode *topNode) override;
    
    int process_event(PHCompositeNode *topNode) override;
    
    int ResetEvent(PHCompositeNode *topNode) override;
    
    int End(PHCompositeNode *topNode) override;
    
    int Reset(PHCompositeNode * /*topNode*/) override;
    
    void Print(const std::string &what = "ALL") const override;
    
    void setGenEvent(int eventGet)     {getEvent = eventGet;}
    
    void setVerbose(bool v) { verbose = v; }

    
    struct EnergyMaps {
        float max_8by8energy_emcal;
        float max_energy_hcal;
        float max_energy_jet;
        int jet_ebin;
        int jet_pbin;
        float energymap_jet_emcal[9][32];
        float energymap_jet_hcalin[9][32];
        float energymap_jet_hcalout[9][32];
    };
    
    void setWantSim(bool sim) { wantSim = sim; }
    void setWantData(bool data) { wantData = data; }
    
    bool getWantSim()  const    { return wantSim;  }
    bool getWantData() const    { return wantData; }
    
private:
    
    bool wantSim = false;   // default false
    bool wantData = true;   // default true
    
    // 2) File pointers
    //    - one for data output
    //    - one for simulation output
    TFile *out     = nullptr;  // data
    TFile *outSim  = nullptr;  // sim

    // 3) Filenames:
    std::string Outfile;     // data output file
    std::string SimOutfile;  // sim output file
    
    int getEvent;
    TriggerAnalyzer* trigAna{nullptr};
    
    std::map<std::string, std::map<std::string, TObject*>> qaHistogramsByTrigger;
    std::map<std::string,
             std::map<std::tuple<float, float, float>,
                      std::map<std::pair<float, float>,
                               std::map<std::string, TObject*>>>> massAndIsolationHistograms;

    std::map<std::string, std::map<std::string, TH1F*>> massAndIsolationHistogramsNoPtBins;
    std::map<std::string, std::map<std::pair<float, float>, std::map<std::string, TObject*>>> qaIsolationHistogramsByTriggerAndPt;
    
    static const std::string IN_MASS_WINDOW_LABEL;
    static const std::string OUTSIDE_MASS_WINDOW_LABEL;
    
    struct MesonMassWindow {
        std::string triggerName;
        float clusEnergy;
        float Chi2;
        float Asym;
        float pTMin;
        float pTMax;
        float meanPi0;
        float sigmaPi0;
        float meanEta;
        float sigmaEta;
    };
    std::map<std::tuple<std::string, float, float, float, float, float>, MesonMassWindow> mesonMassWindowsMap;

    
    bool verbose = true;
    bool m_limitEvents = true;   // Enable event limiting by default
    int m_eventLimit = 50000;    // Maximum number of events to process (10,000 by default)
    
    std::vector<float> asymmetry_values = {0.5, 0.7};
    std::vector<float> clus_chi_values = {3, 4, 5};
    std::vector<float> clus_Energy_values = {1.0, 1.5, 3.0, 5.0, 8.0};
    std::vector<std::pair<float, float>> pT_bins = {
        {2.0, 3.0}, {3.0, 4.0}, {4.0, 5.0}, {5.0, 6.0}, {6.0, 7.0}, {7.0, 8.0}, {8.0, 9.0}, {9.0, 10.0}, {10.0, 12.0}, {12.0, 15.0}, {15, 20}, {20, 30}
    };
    
    std::vector<std::pair<float, float>> isoEtRanges = {
        {-100, 6},
        {-100, 10},
        {-10, 0},
        {0, 10}
    };
    
    std::map<std::string, std::string> triggerNameMap = {
        {"MBD N&S >= 1",          "MBD_NandS_geq_1"},
        {"Jet 6 GeV + MBD NS >=1","Jet_6_GeV_plus_MBD_NS_geq_1"},
        {"Jet 8 GeV + MBD NS >= 1","Jet_8_GeV_plus_MBD_NS_geq_1"},
        {"Jet 10 GeV + MBD NS >= 1","Jet_10_GeV_plus_MBD_NS_geq_1"},
        {"Jet 12 GeV + MBD NS >= 1","Jet_12_GeV_plus_MBD_NS_geq_1"},
        {"Photon 2 GeV+ MBD NS >= 1","Photon_2_GeV_plus_MBD_NS_geq_1"},
        {"Photon 3 GeV + MBD NS >= 1","Photon_3_GeV_plus_MBD_NS_geq_1"},
        {"Photon 4 GeV + MBD NS >= 1","Photon_4_GeV_plus_MBD_NS_geq_1"},
        {"Photon 5 GeV + MBD NS >= 1","Photon_5_GeV_plus_MBD_NS_geq_1"}
    };
    
    // Pointer to the active trigger name map for the current run
    std::map<int, std::string>* activeTriggerNameMap = nullptr;

    int event_count = 0;

    //EMCal
    std::vector<float> m_emcTowE;
    std::vector<float> m_emciEta;
    std::vector<float> m_emciPhi;
    std::vector<int> m_emcTime;
    std::vector<float> m_emcChi2;
    std::vector<float> m_emcPed;
    std::vector<short> m_emcal_good;
    std::vector<float> m_maxTowEnergy;
    
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

    //Clusters
    std::vector<float> m_clusterE;
    std::vector<float> m_clusterEt;
    std::vector<float> m_clusterPhi;
    std::vector<float> m_clusterEta;
    std::vector<float> m_clusterPt;
    std::vector<float> m_clusterChi;
    std::vector<float> m_clusterTowMaxE;
    std::vector<float> m_clusterEnergy;
    std::vector<float> m_clusterEtIso;
    std::vector<int> m_clusterIds;
    
    std::vector<std::vector<int> > m_clusTowEta;
    std::vector<std::vector<int> > m_clusTowPhi;
    std::vector<std::vector<float> > m_clusTowE;
    
    std::map<int, std::pair<float, float>> clusterEtIsoMap_unsubtracted;
    std::map<int, std::pair<float, float>> clusterEtIsoMap_subtracted;
    
    // 4) Photon classification histograms for SIM
    TH1F* hIsoFromPi0Eta  = nullptr;  // e.g. iso <= 6, from pi0/eta
    TH1F* hIsoNotPi0Eta   = nullptr;  // e.g. iso <= 6, not from pi0/eta
    
    //GL1 information
    Gl1Packet *_gl1_packet;
    uint64_t b_gl1_scaledvec;

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
    bool loadMesonMassWindows(const std::string& csvFilePath);
    
    void createHistos_Data();
    void createHistos_ForSimulation();
    
    int process_event_Sim(PHCompositeNode *topNode);
    int process_event_Data(PHCompositeNode *topNode);
    
    int resetEvent_Data(PHCompositeNode* topNode);
    int resetEvent_Sim(PHCompositeNode* topNode);
    
    int endData(PHCompositeNode *topNode);
    int endSim(PHCompositeNode *topNode);
    
    void collectTowerData(TowerInfoContainer* towerContainer, std::vector<TowerData>& towerDataList);

    void processTowers(TowerInfoContainer* towerContainer, float& totalCaloE, std::vector<float>& towEta, std::vector<float>& towPhi, std::vector<float>& towE, std::vector<int>& towTime, std::vector<float>& towChi2, std::vector<float>& towPed, std::vector<short>& towGood);
    
    EnergyMaps processEnergyMaps(
        const std::vector<float>* m_emcTowE, const std::vector<float>* m_emciEta, const std::vector<float>* m_emciPhi,
        const std::vector<float>* m_ohcTowE, const std::vector<float>* m_ohciTowEta, const std::vector<float>* m_ohciTowPhi,
        const std::vector<float>* m_ihcTowE, const std::vector<float>* m_ihciTowEta, const std::vector<float>* m_ihciTowPhi,
        std::vector<short>* m_emcal_good, std::vector<short>* m_ohc_good, std::vector<short>* m_ihc_good, std::vector<std::string>& activeTriggerNames);
    
    void processClusterIsolationHistograms(
        int clusterID,
        float mesonMass,
        float minClusEnergy,
        float maxChi2,
        float maxAsym,
        const std::string& massWindowLabel,
        float pT_min,
        float pT_max,
        const std::string& triggerName,
        const std::map<int, std::pair<float, float>>& clusterEtIsoMap,
        std::map<std::pair<float, float>, std::map<std::string, TObject*>>& cutHistMap,
        size_t& filledHistogramCount,
        bool& filledHistogram,
        bool verbose,
        float pionMass,
        float pionMassWindow,
        float etaMass,
        float etaMassWindow,
        const std::pair<float, float>& pT_bin
    );
    
    void processIsolationRanges(
        const std::vector<std::pair<float, float>>& isoEtRanges,
        const std::vector<int>& clusterIDs,
        size_t clus1,
        size_t clus2,
        float minClusEnergy,
        float maxChi2,
        float maxAsym,
        const std::string& massWindowLabel,
        float pT_min,
        float pT_max,
        const std::string& triggerName,
        const std::map<int, std::pair<float, float>>& clusterEtIsoMap,
        std::map<std::pair<float, float>, std::map<std::string, TObject*>>& cutHistMap,
        bool& filledHistogram,
        bool verbose,
        const std::pair<float, float>& pT_bin
    );
    
    void fillHistogramsForTriggers(
        float mesonMass,
        size_t clus1,
        size_t clus2,
        float pt1,
        float pt2,
        float E1,
        float E2,
        float minClusEnergy,
        float maxChi2,
        float maxAsym,
        size_t& filledHistogramCount,
        const std::vector<int>& clusterIDs,
        const std::map<int, std::pair<float, float>>& clusterEtIsoMap,
        const std::vector<std::string>& activeTriggerNames,
        bool& filledHistogram);

    void processClusterInvariantMass(
        const std::vector<float>& clusterE,
        const std::vector<float>& clusterPt,
        const std::vector<float>& clusterChi2,
        const std::vector<float>& clusterEta,
        const std::vector<float>& clusterPhi,
        const std::vector<int>& clusterIDs,
        const std::map<int, std::pair<float, float>>& clusterEtIsoMap,
        const std::vector<std::string>& activeTriggerNames);



    float getMaxTowerE(RawCluster *cluster, TowerInfoContainer *towerContainer);
    std::vector<float> returnClusterTowE(RawCluster *cluster, TowerInfoContainer *towerContainer);
    std::vector<int> returnClusterTowPhi(RawCluster *cluster, TowerInfoContainer *towerContainer);
    std::vector<int> returnClusterTowEta(RawCluster *cluster, TowerInfoContainer *towerContainer);


    float totalCaloEEMCal;
    float totalCaloEOHCal;
    float totalCaloEIHCal;
    float totalCaloEZDC;
    float totalChargeMBD;

  
    // Inline deltaR function for calculating distance between points in η-φ space
    inline float deltaR(float eta1, float eta2, float phi1, float phi2) {
        float deta = eta1 - eta2;
        float dphi = phi1 - phi2;
        if (dphi > M_PI) dphi -= 2 * M_PI;
        if (dphi < -M_PI) dphi += 2 * M_PI;
        return sqrt(deta * deta + dphi * dphi);
    }
    
    inline std::string formatFloatForFilename(float value) {
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
    
    bool IsAcceptableTower(TowerInfo* tower);

};

#endif  // CALOTREEGEN_H
