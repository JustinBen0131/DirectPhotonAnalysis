// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOTREEGEN_H
#define CALOTREEGEN_H

#include <fun4all/SubsysReco.h>

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

    caloTreeGen(const std::string &name = "caloTreeGen");

    ~caloTreeGen() override;

    /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
    */
    int Init(PHCompositeNode *topNode) override;

    /** Called for each event.
      This is where you do the real work.
    */
    int process_event(PHCompositeNode *topNode) override;

    /// Clean up internals after each event.
    int ResetEvent(PHCompositeNode *topNode) override;

    /// Called at the end of all processing.
    int End(PHCompositeNode *topNode) override;

    /// Reset
    int Reset(PHCompositeNode * /*topNode*/) override;

    void Print(const std::string &what = "ALL") const override;
    
    void setGenEvent(int eventGet)     {getEvent = eventGet;}
    
    void setVerbose(bool v) { verbose = v; }

 private:

    TFile *out;
    std::string Outfile = "commissioning.root";
    int getEvent;
    std::map<int, std::map<std::string, TObject*>> qaHistogramsByTrigger;
    // Declare the map to hold histograms for each trigger, cut combination, and pT bin
    std::map<int, std::map<std::tuple<float, float, float>, std::map<std::pair<float, float>, std::map<std::string, TObject*>>>> massAndIsolationHistograms;
    std::map<int, std::map<std::string, TH1F*>> massAndIsolationHistogramsNoPtBins;
    std::map<int, std::map<std::pair<float, float>, std::map<std::string, TObject*>>> qaIsolationHistogramsByTriggerAndPt;

    static const std::string IN_MASS_WINDOW_LABEL;
    static const std::string OUTSIDE_MASS_WINDOW_LABEL;
    
    struct MesonMassWindow {
        int triggerIndex;
        float Ecore;
        float Chi2;
        float Asym;
        float pTMin;
        float pTMax;
        float meanPi0;
        float sigmaPi0;
        float meanEta;
        float sigmaEta;
    };

    // Define the map to store mass windows by trigger, Ecore, Chi2, Asym, pTMin, and pTMax
    std::map<std::tuple<int, float, float, float, float, float>, MesonMassWindow> mesonMassWindowsMap;
    
    bool verbose = true;
    bool m_limitEvents = true;   // Enable event limiting by default
    int m_eventLimit = 5000;    // Maximum number of events to process (10,000 by default)

    std::vector<int> triggerIndices = {10, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};
    std::vector<float> asymmetry_values = {0.5, 0.7};
    std::vector<float> clus_chi_values = {4};
    std::vector<float> clus_Ecore_values = {1.0, 1.5};
    std::vector<std::pair<float, float>> pT_bins = {
        {2.0, 3.0}, {3.0, 4.0}, {4.0, 5.0}, {5.0, 6.0}, {6.0, 7.0}, {7.0, 8.0}, {8.0, 9.0}, {9.0, 10.0}, {10.0, 12.0}, {12.0, 15.0}, {15, 20}, {20, 30}
    };
    
    std::vector<std::pair<float, float>> isoEtRanges = {
        {-100, 6},
        {-100, 6},
        {-10, 0},
        {0, 10}
    };
    
    triggerNameMap1 = {
        {0, "Clock"},
        {1, "ZDC_South"},
        {2, "ZDC_North"},
        {3, "ZDC_Coincidence"},
        {4, "HCAL_Singles"},
        {5, "HCAL_Coincidence"},
        {8, "MBD_S_>=_1"},
        {9, "MBD_N_>=_1"},
        {10, "MBD_N&S_>=_1"},
        {11, "MBD_N&S_>=_2"},
        {12, "MBD_N&S_>=_1_vtx_<_10_cm"},
        {13, "MBD_N&S_>=_1_vtx_<_30_cm"},
        {14, "MBD_N&S_>=_1_vtx_<_60_cm"},
        {15, "HCAL_Singles_+_MBD_NS_>=_1"},
        {16, "Jet_4_GeV_+_MBD_NS_>=_1"},
        {17, "Jet_6_GeV_+_MBD_NS_>=_1"},
        {18, "Jet_8_GeV_+_MBD_NS_>=_1"},
        {19, "Jet_10_GeV_+_MBD_NS_>=_1"},
        {20, "Jet_4_GeV"},
        {21, "Jet_6_GeV"},
        {22, "Jet_8_GeV"},
        {23, "Jet_10_GeV"},
        {24, "Photon_1_GeV_+_MBD_NS_>=_1"},
        {25, "Photon_2_GeV_+_MBD_NS_>=_1"},
        {26, "Photon_3_GeV_+_MBD_NS_>=_1"},
        {27, "Photon_4_GeV_+_MBD_NS_>=_1"},
        {28, "Photon_1_GeV"},
        {29, "Photon_2_GeV"},
        {30, "Photon_3_GeV"},
        {31, "Photon_4_GeV"}
    };

    triggerNameMap2 = {
        {0, "Clock"},
        {1, "ZDC_South"},
        {2, "ZDC_North"},
        {3, "ZDC_Coincidence"},
        {4, "HCAL_Singles"},
        {5, "HCAL_Coincidence"},
        {8, "MBD_S_>=_1"},
        {9, "MBD_N_>=_1"},
        {10, "MBD_N&S_>=_1"},
        {11, "MBD_N&S_>=_2"},
        {12, "MBD_N&S_>=_1_vtx_<_10_cm"},
        {13, "MBD_N&S_>=_1_vtx_<_30_cm"},
        {14, "MBD_N&S_>=_1_vtx_<_60_cm"},
        {15, "HCAL_Singles_+_MBD_NS_>=_1"},
        {16, "Jet_6_GeV_+_MBD_NS_>=_1"},
        {17, "Jet_8_GeV_+_MBD_NS_>=_1"},
        {18, "Jet_10_GeV_+_MBD_NS_>=_1"},
        {19, "Jet_12_GeV_+_MBD_NS_>=_1"},
        {20, "Jet_6_GeV"},
        {21, "Jet_8_GeV"},
        {22, "Jet_10_GeV"},
        {23, "Jet_12_GeV"},
        {24, "Photon_2_GeV_+_MBD_NS_>=_1"},
        {25, "Photon_3_GeV_+_MBD_NS_>=_1"},
        {26, "Photon_4_GeV_+_MBD_NS_>=_1"},
        {27, "Photon_5_GeV_+_MBD_NS_>=_1"},
        {28, "Photon_2_GeV"},
        {29, "Photon_3_GeV"},
        {30, "Photon_4_GeV"},
        {31, "Photon_5_GeV"}
    };
    runNumbersForMap1 = {44477, 44478, 44482, 44483, 44495, 44498, 44499, 44503, 44505, 44506, 44507, 44509, 44510, 44511, 44512, 44513, 44533
        44534, 44604, 44608, 44611, 44616, 44618, 44619, 44621, 44631, 44638, 44642, 45034, 45035, 45036, 45038, 45041, 45048
        45051, 45052, 45090, 45100, 45103, 45105, 45106, 45107, 45150, 45151, 45153, 45154, 45155, 45157, 45159, 45160, 45161
        45162, 45164, 45166, 45167, 45170, 45172, 45176, 45177, 45178, 45181, 45183, 45186, 45189, 45190, 45191, 45196, 45199
        45201, 45203, 45246, 45248, 45249, 45252, 45255, 45256, 45258, 45274, 45288, 45290, 45291, 45292, 45315, 45316, 45318
        45325, 45390, 45391, 45393, 45394, 45401, 45402, 45414, 45443, 45485, 45486, 45487, 45489, 45490, 45491, 45493, 45494
        45495, 45531, 45540, 45541, 45547, 45548, 45551, 45552, 45620, 45624, 45627, 45628, 45633, 45637, 45645, 45724, 45807
        45816, 45837, 45841, 45842, 45851, 45852, 45856, 45858, 45872, 45883, 46011, 46019, 46022, 46023, 46025, 46029, 46036};
    
    runNumbersForMap2 = {44477, 44478, 44482, 44483, 44495, 44498, 44499, 44503, 44505, 44506, 44507, 44509, 44510, 44511, 44512, 44513, 44533
        44534, 44604, 44608, 44611, 44616, 44618, 44619, 44621, 44631, 44638, 44642, 45034, 45035, 45036, 45038, 45041, 45048
        45051, 45052, 45090, 45100, 45103, 45105, 45106, 45107, 45150, 45151, 45153, 45154, 45155, 45157, 45159, 45160, 45161
        45162, 45164, 45166, 45167, 45170, 45172, 45176, 45177, 45178, 45181, 45183, 45186, 45189, 45190, 45191, 45196, 45199
        45201, 45203, 45246, 45248, 45249, 45252, 45255, 45256, 45258, 45274, 45288, 45290, 45291, 45292, 45315, 45316, 45318
        45325, 45390, 45391, 45393, 45394, 45401, 45402, 45414, 45443, 45485, 45486, 45487, 45489, 45490, 45491, 45493, 45494
        45495, 45531, 45540, 45541, 45547, 45548, 45551, 45552, 45620, 45624, 45627, 45628, 45633, 45637, 45645, 45724, 45807
        45816, 45837, 45841, 45842, 45851, 45852, 45856, 45858, 45872, 45883, 46011, 46019, 46022, 46023, 46025, 46029, 46036};
    
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
    std::vector<float> m_clusterECore;
    std::vector<float> m_clusterEtIso;
    std::vector<int> m_clusterIds;
    
    std::vector<std::vector<int> > m_clusTowEta;
    std::vector<std::vector<int> > m_clusTowPhi;
    std::vector<std::vector<float> > m_clusTowE;
    
    std::map<int, std::pair<float, float>> clusterEtIsoMap_unsubtracted;
    std::map<int, std::pair<float, float>> clusterEtIsoMap_subtracted;
    
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
    
    void collectTowerData(TowerInfoContainer* towerContainer, std::vector<TowerData>& towerDataList);

    void processTowers(TowerInfoContainer* towerContainer, float& totalCaloE, std::vector<float>& towEta, std::vector<float>& towPhi, std::vector<float>& towE, std::vector<int>& towTime, std::vector<float>& towChi2, std::vector<float>& towPed, std::vector<short>& towGood);
    
    EnergyMaps processEnergyMaps(
        const std::vector<float>* m_emcTowE, const std::vector<float>* m_emciEta, const std::vector<float>* m_emciPhi,
        const std::vector<float>* m_ohcTowE, const std::vector<float>* m_ohciTowEta, const std::vector<float>* m_ohciTowPhi,
        const std::vector<float>* m_ihcTowE, const std::vector<float>* m_ihciTowEta, const std::vector<float>* m_ihciTowPhi,
        std::vector<short>* m_emcal_good, std::vector<short>* m_ohc_good, std::vector<short>* m_ihc_good, std::vector<int> activeTriggerBits);
    
    void processClusterIsolationHistograms(
        int clusterID,
        float mesonMass,
        float minClusEcore,
        float maxChi2,
        float maxAsym,
        const std::string& massWindowLabel,
        float pT_min,
        float pT_max,
        int triggerIndex,
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
        float minClusEcore,
        float maxChi2,
        float maxAsym,
        const std::string& massWindowLabel,
        float pT_min,
        float pT_max,
        int triggerIndex,
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
        float minClusEcore,
        float maxChi2,
        float maxAsym,
        size_t& filledHistogramCount,
        const std::vector<int>& clusterIDs,
        const std::map<int, std::pair<float, float>>& clusterEtIsoMap,
        const std::vector<int>& activeTriggerBits,
        bool& filledHistogram);

    void processClusterInvariantMass(
        const std::vector<float>& clusterE,
        const std::vector<float>& clusterPt,
        const std::vector<float>& clusterChi2,
        const std::vector<float>& clusterEta,
        const std::vector<float>& clusterPhi,
        const std::vector<int>& clusterIDs,
        const std::map<int, std::pair<float, float>>& clusterEtIsoMap, std::vector<int> activeTriggerBits);



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

    // Inline function to extract trigger bits from GL1 scaled vector
    inline std::vector<int> extractTriggerBits(uint64_t b_gl1_scaledvec, int entry) {
        std::vector<int> trig_bits;
        std::bitset<64> bits(b_gl1_scaledvec);
        if (verbose) {
            std::cout << "Processing entry " << entry << ", gl1_scaledvec (bits): " << bits.to_string() << std::endl;
        }
        
        for (unsigned int bit = 0; bit < 64; bit++) {
            if (((b_gl1_scaledvec >> bit) & 0x1U) == 0x1U) {
                trig_bits.push_back(bit);
            }
        }
        return trig_bits;
    }

    // Inline function to check trigger condition
    inline bool checkTriggerCondition(const std::vector<int> &trig_bits, int inputBit) {
        for (const int &bit : trig_bits) {
            if (bit == inputBit) {
                if (verbose) {
                    std::cout << "  Trigger condition met with bit: " << bit << std::endl;
                }
                
                return true;
            }
        }
        if (verbose) {
            std::cout << "  No relevant trigger conditions met." << std::endl;
        }
        
        return false;
    }
    // Function to set the correct trigger map based on the run number
    void setTriggerNameMapForRun(int runNumber);
    
    bool IsAcceptableTower(TowerInfo* tower);

};

#endif  // CALOTREEGEN_H
