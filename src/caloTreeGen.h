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

 private:

    TFile *out;
    std::string Outfile = "commissioning.root";
    std::map<int, std::map<std::string, TObject*>> qaHistogramsByTrigger;
    std::map<int, std::map<std::string, TObject*>> massHistogramsByTrigger;
    std::map<int, TTree*> treesByTrigger; // Map of TTrees for each trigger bit


    std::vector<int> triggerIndices = {10, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};
    std::vector<double> dR_values = {0.3, 0.4, 0.5};
    std::vector<float> asymmetry_values = {0.7};
    std::vector<float> clus_chi_values = {4};
    std::vector<float> clus_Ecore_values = {0.5, 1};

    // Declare counters to track across events
    double min_isoEt;
    double max_isoEt;
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
    void clearVectors();
    void collectTowerData(TowerInfoContainer* towerContainer, std::vector<TowerData>& towerDataList);

    void processTowers(TowerInfoContainer* towerContainer, float& totalCaloE, std::vector<float>& towEta, std::vector<float>& towPhi, std::vector<float>& towE, std::vector<int>& towTime, std::vector<float>& towChi2, std::vector<float>& towPed, std::vector<short>& towGood);
    
    void processEnergyMaps(
        const std::vector<float>* m_emcTowE, const std::vector<float>* m_emciEta, const std::vector<float>* m_emciPhi,
        const std::vector<float>* m_ohcTowE, const std::vector<float>* m_ohciTowEta, const std::vector<float>* m_ohciTowPhi,
        const std::vector<float>* m_ihcTowE, const std::vector<float>* m_ihciTowEta, const std::vector<float>* m_ihciTowPhi,
        std::vector<short>* m_emcal_good, std::vector<short>* m_ohc_good, std::vector<short>* m_ihc_good,
        float& max_8by8energy_emcal, float& max_energy_hcal, float& max_energy_jet,
        std::map<std::string, TObject*>& qaHistograms, int triggerIndex);
    
    void calculateIsoEt(TowerInfoContainer* towerContainer, RawTowerGeomContainer* geomContainer, double& isoEt, const double clus_eta, const double clus_phi, const double m_vx, const double m_vy, const double m_vz, int& skipped_tower_count, int& towers_in_cone_count, const std::string& geomContainerName, double dR_cut);

    void processIsolationEnergy(
        TowerInfoContainer* emcTowerContainer,
        TowerInfoContainer* ihcTowerContainer,
        TowerInfoContainer* ohcTowerContainer,
        RawTowerGeomContainer* geomEM,
        RawTowerGeomContainer* geomIH,
        RawTowerGeomContainer* geomOH,
        float clusEcore,
        float clus_eta,
        float clus_phi,
        float m_vx,
        float m_vy,
        float m_vz,
        std::map<std::string, TObject*>& qaHistograms,
        int triggerIndex,
        std::vector<double> dR_values,
        int& skipped_tower_count,
        int& towers_in_cone_count,
        int& positive_isoEt_count,
        int& negative_isoEt_count,
        double& min_isoEt,
        double& max_isoEt,
        std::vector<float>& m_clusterEtIso
    );
    void processClusterInvariantMass(
        const std::vector<float>& clusterE,
        const std::vector<float>& clusterPt,
        const std::vector<float>& clusterChi2,
        const std::vector<float>& clusterEta,
        const std::vector<float>& clusterPhi,
        int triggerIndex,
        std::map<std::string, TObject*>& massHistograms,
        const std::vector<float>& asymmetry_values,
        const std::vector<float>& clus_chi_values,
        const std::vector<float>& clus_Ecore_values);

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
//        std::cout << "Processing entry " << entry << ", gl1_scaledvec (bits): " << bits.to_string() << std::endl;
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
//                std::cout << "  Trigger condition met with bit: " << bit << std::endl;
                return true;
            }
        }
//        std::cout << "  No relevant trigger conditions met." << std::endl;
        return false;
    }
    
    double getTowerEta(RawTowerGeom* tower_geom, double vx, double vy, double vz);
    bool IsAcceptableTower(TowerInfo* tower);

};

#endif  // CALOTREEGEN_H
