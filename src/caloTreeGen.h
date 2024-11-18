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
    
    // Function to set the correct trigger map based on the run number
    void setTriggerNameMapForRun(int runNumber);
    
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
    
    bool verbose = false;
    bool m_limitEvents = false;   // Enable event limiting by default
    int m_eventLimit = 2000;    // Maximum number of events to process (10,000 by default)
    
    std::vector<int> triggerIndices = {10, 24, 25, 26, 27, 28, 29, 30, 31};
    std::vector<float> asymmetry_values = {0.5, 0.7};
    std::vector<float> clus_chi_values = {4};
    std::vector<float> clus_Ecore_values = {1.0, 1.5};
    std::vector<std::pair<float, float>> pT_bins = {
        {2.0, 3.0}, {3.0, 4.0}, {4.0, 5.0}, {5.0, 6.0}, {6.0, 7.0}, {7.0, 8.0}, {8.0, 9.0}, {9.0, 10.0}, {10.0, 12.0}, {12.0, 15.0}, {15, 20}, {20, 30}
    };
    
    std::vector<std::pair<float, float>> isoEtRanges = {
        {-100, 6},
        {-100, 10},
        {-10, 0},
        {0, 10}
    };
    
    std::map<int, std::string> triggerNameMap1 = {
        {0, "Clock"},
        {1, "ZDC_South"},
        {2, "ZDC_North"},
        {3, "ZDC_Coincidence"},
        {4, "HCAL_Singles"},
        {5, "HCAL_Coincidence"},
        {8, "MBD_S_geq_1"},
        {9, "MBD_N_geq_1"},
        {10, "MBD_NandS_geq_1"},
        {11, "MBD_NandS_geq_2"},
        {12, "MBD_NandS_geq_1_vtx_lessThen_10_cm"},
        {13, "MBD_NandS_geq_1_vtx_lessThen_30_cm"},
        {14, "MBD_NandS_geq_1_vtx_lessThen_60_cm"},
        {15, "HCAL_Singles_plus_MBD_NS_geq_1"},
        {16, "Jet_4_GeV_plus_MBD_NS_geq_1"},
        {17, "Jet_6_GeV_plus_MBD_NS_geq_1"},
        {18, "Jet_8_GeV_plus_MBD_NS_geq_1"},
        {19, "Jet_10_GeV_plus_MBD_NS_geq_1"},
        {20, "Jet_4_GeV"},
        {21, "Jet_6_GeV"},
        {22, "Jet_8_GeV"},
        {23, "Jet_10_GeV"},
        {24, "Photon_1_GeV_plus_MBD_NS_geq_1"},
        {25, "Photon_2_GeV_plus_MBD_NS_geq_1"},
        {26, "Photon_3_GeV_plus_MBD_NS_geq_1"},
        {27, "Photon_4_GeV_plus_MBD_NS_geq_1"},
        {28, "Photon_1_GeV"},
        {29, "Photon_2_GeV"},
        {30, "Photon_3_GeV"},
        {31, "Photon_4_GeV"}
    };
    
    std::map<int, std::string> triggerNameMap2 = {
        {0, "Clock"},
        {1, "ZDC_South"},
        {2, "ZDC_North"},
        {3, "ZDC_Coincidence"},
        {4, "HCAL_Singles"},
        {5, "HCAL_Coincidence"},
        {8, "MBD_S_geq_1"},
        {9, "MBD_N_geq_1"},
        {10, "MBD_NandS_geq_1"},
        {11, "MBD_NandS_geq_2"},
        {12, "MBD_NandS_geq_1_vtx_lessThen_10_cm"},
        {13, "MBD_NandS_geq_1_vtx_lessThen_30_cm"},
        {14, "MBD_NandS_geq_1_vtx_lessThen_60_cm"},
        {15, "HCAL_Singles_plus_MBD_NS_geq_1"},
        {16, "Jet_6_GeV_plus_MBD_NS_geq_1"},
        {17, "Jet_8_GeV_plus_MBD_NS_geq_1"},
        {18, "Jet_10_GeV_plus_MBD_NS_geq_1"},
        {19, "Jet_12_GeV_plus_MBD_NS_geq_1"},
        {20, "Jet_6_GeV"},
        {21, "Jet_8_GeV"},
        {22, "Jet_10_GeV"},
        {23, "Jet_12_GeV"},
        {24, "Photon_2_GeV_plus_MBD_NS_geq_1"},
        {25, "Photon_3_GeV_plus_MBD_NS_geq_1"},
        {26, "Photon_4_GeV_plus_MBD_NS_geq_1"},
        {27, "Photon_5_GeV_plus_MBD_NS_geq_1"},
        {28, "Photon_2_GeV"},
        {29, "Photon_3_GeV"},
        {30, "Photon_4_GeV"},
        {31, "Photon_5_GeV"}
    };
    
    
    std::vector<int>  runNumbersForMap1 = { 44477,44478,44482,44483,44495,44498,44499,44503,44505,44506,44507,44509,44510,44511,44512,44513,44533,44534,44604,44608,44611,
        44616,44618,44619,44621,44631,44638,44642,45034,45035,45036,45038,45041,45048,45051,45052,45090,45100,45103,45105,45106,45107,
        45150,45151,45153,45154,45155,45157,45159,45160,45161,45162,45164,45166,45167,45170,45172,45176,45177,45178,45181,45183,45186,
        45189,45190,45191,45196,45199,45201,45203,45246,45248,45249,45252,45255,45256,45258,45274,45288,45290,45291,45292,45315,45316,
        45318,45325,45390,45391,45393,45394,45401,45402,45414,45443,45485,45486,45487,45489,45490,45491,45493,45494,45495,45531,45540,
        45541,45547,45548,45551,45552,45620,45624,45627,45628,45633,45637,45645,45724,45807,45816,45837,45841,45842,45851,45852,45856,
        45858,45872,45883,46011,46019,46022,46023,46025,46029,46036};
    
    
     std::vector<int> runNumbersForMap2 = {
        46065,46068,46105,46107,46133,46135,46137,46427,46429,46430,46431,46432,46434,46436,46437,46438,46451,46452,46455,46473,46477,
        46480,46523,46524,46529,46530,46531,46537,46539,46541,46545,46546,46553,46559,46563,46567,46569,46577,46587,46588,46593,46595,
        46598,46603,46605,46619,46623,46640,46649,46655,46656,46663,46665,46668,46669,46676,46697,46699,46704,46721,46735,46750,46753,
        46755,46756,46759,46772,46775,46912,46917,46918,46941,46943,46944,46945,46946,46947,46949,46950,46951,46952,46965,46967,46968,
        47002,47006,47007,47009,47014,47017,47019,47022,47032,47033,47034,47036,47038,47040,47043,47051,47052,47053,47055,47056,47058,
        47060,47061,47064,47066,47068,47089,47098,47101,47102,47114,47115,47116,47124,47125,47129,47131,47135,47137,47138,47139,47140,
        47141,47143,47146,47155,47156,47158,47160,47161,47162,47201,47202,47203,47204,47211,47216,47219,47229,47230,47289,47293,47303,
        47306,47310,47315,47316,47323,47330,47332,47334,47360,47375,47376,47377,47378,47381,47382,47391,47393,47395,47396,47399,47443,
        47451,47455,47457,47458,47459,47464,47474,47476,47480,47484,47485,47491,47492,47494,47495,47497,47502,47503,47505,47506,47507,
        47513,47514,47516,47522,47524,47525,47538,47540,47548,47552,47557,47568,47634,47636,47638,47657,47658,47659,47661,47662,47666,
        47667,47698,47715,47716,47720,47722,47723,47724,47725,47727,47729,47730,47732,47733,47769,47777,47778,47783,47807,47831,47846,
        47848,47867,47893,47939,47946,47962,47966,47982,48073,48080,48085,48086,48100,48166,48180,48181,48231,48233,48234,48237,48239,
        48240,48244,48245,48253,48255,48256,48257,48258,48260,48261,48262,48263,48265,48287,48291,48293,48294,48295,48307,48313,48318,
        48320,48323,48325,48326,48327,48335,48337,48338,48341,48342,48343,48346,48347,48348,48349,48352,48356,48357,48358,48359,48366,
        48367,48369,48409,48410,48412,48416,48417,48418,48421,48422,48423,48454,48455,48456,48459,48461,48462,48469,48536,48635,48636,
        48638,48645,48656,48657,48658,48660,48701,48720,48721,48722,48725,48726,48727,48730,48731,48732,48734,48736,48742,48743,48745,
        48746,48801,48803,48805,48806,48807,48810,48811,48813,48824,48826,48828,48829,48832,48836,48838,48839,48859,48861,48863,48864,
        48865,48867,48868,48869,48870,48872,48873,48874,48877,48883,48884,48885,48891,48892,48893,48894,48895,48896,48897,48899,48900,
        48901,48902,48903,48918,48935,48936,48938,48940,48943,48946,48949,48951,48971,48976,48981,48982,48983,48984,48985,48986,48987,
        48988,48990,48991,49017,49023,49026,49027,49028,49029,49030,49031,49035,49042,49044,49047,49048,49050,49052,49053,49054,49060,
        49061,49062,49063,49066,49067,49069,49070,49071,49072,49073,49098,49125,49128,49133,49138,49219,49224,49226,49227,49228,49229,
        49230,49233,49240,49241,49244,49247,49248,49249,49250,49251,49254,49263,49264,49266,49267,49268,49269,49270,49307,49308,49309,
        49310,49311,49312,49313,49314,49316,49317,49324,49329,49330,49331,49332,49333,49336,49337,49338,49339,49340,49343,49345,49346,
        49347,49348,49349,49350,49351,49352,49356,49357,49358,49359,49361,49362,49363,49365,49366,49367,49368,49372,49374,49376,49377,
        49378,49379,49380,49381,49382,49383,49384,49385,49386,49389,49390,49433,49434,49435,49437,49438,49439,49440,49445,49446,49447,
        49448,49449,49451,49452,49453,49454,49455,49456,49457,49458,49464,49466,49467,49650,49651,49652,49653,49655,49656,49658,49660,
        49661,49662,49663,49664,49736,49737,49742,49743,49748,49749,49750,49751,49752,49758,49760,49761,49762,50343,50507,50546,50554,
        50595,50599,50600,50601,50602,50603,50605,50606,50607,50613,50615,50650,50655,50658,50662,50663,50664,50666,50668,50670,50671,
        50673,50674,50676,50677,50678};
    
    
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
        float minClusEcore,
        float maxChi2,
        float maxAsym,
        const std::string& massWindowLabel,
        float pT_min,
        float pT_max,
        int triggerIndex,
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
    
    bool IsAcceptableTower(TowerInfo* tower);

};

#endif  // CALOTREEGEN_H
