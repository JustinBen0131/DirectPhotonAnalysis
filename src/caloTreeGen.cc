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

#include <filesystem>
#include <fstream>
#include <locale>

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

const std::string caloTreeGen::IN_MASS_WINDOW_LABEL = "_inMassWindow";
const std::string caloTreeGen::OUTSIDE_MASS_WINDOW_LABEL = "_outsideMassWindow";


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

bool caloTreeGen::loadMesonMassWindows(const std::string& csvFilePath) {
    // Set the global locale to classic (C locale)
    std::locale::global(std::locale::classic());

    std::ifstream csvFile(csvFilePath);
    if (!csvFile.is_open()) {
        std::cerr << "Error: Could not open CSV file " << csvFilePath << std::endl;
        return false;
    }

    if (verbose) {
        std::cout << "Opening CSV file: " << csvFilePath << std::endl;
    }

    std::string line;
    // Skip the header line
    std::getline(csvFile, line);

    int lineNumber = 2; // Start from line 2 since we skipped the header
    int entriesLoaded = 0;

    while (std::getline(csvFile, line)) {
        std::istringstream lineStream(line);
        std::string cell;

        try {
            // Read and parse each field from the CSV line
            std::string triggerName;
            float Ecore, Chi2, Asym, pTMin, pTMax;
            float meanPi0, sigmaPi0, meanEta, sigmaEta;

            // Trim function to remove whitespace
            auto trim = [](std::string& s) {
                s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
                    return !std::isspace(ch);
                }));
                s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
                    return !std::isspace(ch);
                }).base(), s.end());
            };

            // Read the columns using ',' as delimiter
            std::getline(lineStream, cell, ','); triggerName = cell;
            std::getline(lineStream, cell, ','); trim(cell); Ecore = std::stof(cell);
            std::getline(lineStream, cell, ','); trim(cell); Chi2 = std::stof(cell);
            std::getline(lineStream, cell, ','); trim(cell); Asym = std::stof(cell);
            std::getline(lineStream, cell, ','); trim(cell); pTMin = std::stof(cell);
            std::getline(lineStream, cell, ','); trim(cell); pTMax = std::stof(cell);
            std::getline(lineStream, cell, ','); trim(cell); meanPi0 = std::stof(cell);
            std::getline(lineStream, cell, ','); /* Skip meanPi0Error */
            std::getline(lineStream, cell, ','); trim(cell); sigmaPi0 = std::stof(cell);
            std::getline(lineStream, cell, ','); /* Skip sigmaPi0Error */
            std::getline(lineStream, cell, ','); trim(cell); meanEta = std::stof(cell);
            std::getline(lineStream, cell, ','); /* Skip meanEtaError */
            std::getline(lineStream, cell, ','); trim(cell); sigmaEta = std::stof(cell);
            std::getline(lineStream, cell, ','); /* Skip sigmaEtaError */

            // Create the MesonMassWindow struct and tuple key
            MesonMassWindow massWindow = {triggerName, Ecore, Chi2, Asym, pTMin, pTMax, meanPi0, sigmaPi0, meanEta, sigmaEta};
            auto key = std::make_tuple(triggerName, Ecore, Chi2, Asym, pTMin, pTMax);

            // Insert into the map
            mesonMassWindowsMap[key] = massWindow;
            entriesLoaded++;

            if (verbose) {
                std::cout << "Loaded entry " << entriesLoaded << " (Line " << lineNumber << "):" << std::endl;
                std::cout << " - Trigger Name: " << triggerName << std::endl;
                std::cout << " - Ecore: " << Ecore << ", Chi2: " << Chi2 << ", Asym: " << Asym << std::endl;
                std::cout << " - pT Range: [" << pTMin << ", " << pTMax << "]" << std::endl;
                std::cout << " - Mean Pi0: " << meanPi0 << " ± " << sigmaPi0 << std::endl;
                std::cout << " - Mean Eta: " << meanEta << " ± " << sigmaEta << std::endl;
            }
        } catch (const std::exception& e) {
            std::cerr << "Error parsing line " << lineNumber << " in CSV file: " << e.what() << std::endl;
            std::cerr << "Line content: " << line << std::endl;
            std::cerr << "Exception occurred when parsing cell: '" << cell << "'" << std::endl;
            continue;  // Skip this line and continue
        }

        lineNumber++;
    }

    csvFile.close();

    if (verbose) {
        std::cout << "Finished loading meson mass windows from CSV." << std::endl;
        std::cout << "Total entries loaded: " << entriesLoaded << std::endl;
        std::cout << "MesonMassWindow map size: " << mesonMassWindowsMap.size() << std::endl;
    }

    return true;
}


void caloTreeGen::setTriggerNameMapForRun(int runNumber) {
    if (std::find(runNumbersForMap1.begin(), runNumbersForMap1.end(), runNumber) != runNumbersForMap1.end()) {
        activeTriggerNameMap = &triggerNameMap1;
        if (verbose) {
            std::cout << "[INFO] Using triggerNameMap1 for run number " << runNumber << std::endl;
        }
    } else if (std::find(runNumbersForMap2.begin(), runNumbersForMap2.end(), runNumber) != runNumbersForMap2.end()) {
        activeTriggerNameMap = &triggerNameMap2;
        if (verbose) {
            std::cout << "[INFO] Using triggerNameMap2 for run number " << runNumber << std::endl;
        }
    } else {
        std::cerr << "[ERROR] Run number " << runNumber << " does not match any trigger map." << std::endl;
        activeTriggerNameMap = nullptr;
    }
}

//____________________________________________________________________________..
int caloTreeGen::Init(PHCompositeNode *topNode) {
    std::cout << ANSI_COLOR_BLUE_BOLD << "Initializing caloTreeGen..." << ANSI_COLOR_RESET << std::endl;
    out = new TFile(Outfile.c_str(),"RECREATE");
    
    if (verbose) {
        std::cout << ANSI_COLOR_BLUE_BOLD << "Output file created: " << Outfile << ANSI_COLOR_RESET << std::endl;
    }
    
    // Load meson mass windows from the CSV file if it exists
    const std::string csvFilePath = "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/InvMassCsvFiles/InvariantMassInformation_MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1_Photon_5_GeV_plus_MBD_NS_geq_1_afterTriggerFirmwareUpdate.csv";

    if (std::filesystem::exists(csvFilePath)) {
        if (!loadMesonMassWindows(csvFilePath)) {
            std::cerr << "Warning: Failed to load meson mass windows from CSV. Continuing without it." << std::endl;
        }
    } else {
        std::cout << "No CSV file found at " << csvFilePath << ". Skipping meson mass windows loading." << std::endl;
    }
    
    for (int triggerIndex : triggerIndices) {
        // Get the trigger name from the active trigger name map
        std::string triggerName;
        if (activeTriggerNameMap && activeTriggerNameMap->find(triggerIndex) != activeTriggerNameMap->end()) {
            triggerName = activeTriggerNameMap->at(triggerIndex);
        } else {
            std::cerr << "[ERROR] Trigger index " << triggerIndex << " not found in active trigger map." << std::endl;
            continue; // Skip to next triggerIndex if not found
        }

        if (verbose) {
            std::cout << ANSI_COLOR_BLUE_BOLD << "Creating histograms for trigger: " << triggerName << " (index " << triggerIndex << ")" << ANSI_COLOR_RESET << std::endl;
        }

        // Create a directory for the current trigger
        TDirectory* triggerDir = out->mkdir(triggerName.c_str());
        if (!triggerDir) {
            std::cerr << "[ERROR] Failed to create directory for trigger: " << triggerName << std::endl;
            exit(EXIT_FAILURE);
        }
        triggerDir->cd(); // Set the current directory to the trigger directory


        std::map<std::string, TObject*>& qaHistograms = qaHistogramsByTrigger[triggerIndex];
        // Helper functions for creating histograms with logging
        auto createHistogram = [&](const std::string& name, const std::string& title, int bins, double xMin, double xMax) {
            // Check if a histogram with the same name already exists
            if (out->Get(name.c_str()) != nullptr) {
                std::cerr << "\n[ERROR] Duplicate histogram detected: " << name << std::endl;
                std::cerr << "A histogram with this name already exists in the output file." << std::endl;
                std::cerr << "Aborting Fun4All macro to avoid further conflicts." << std::endl;
                exit(EXIT_FAILURE);
            }
            
            if (verbose) {
                std::cout << ANSI_COLOR_RED_BOLD << "Creating histogram: " << name << " - " << title << ANSI_COLOR_RESET << std::endl;
            }
            
            TH1F* hist = new TH1F(name.c_str(), title.c_str(), bins, xMin, xMax);
            hist->SetDirectory(out); // Ensure it is linked to the output file
            return hist;
        };

        auto create2DHistogram = [&](const std::string& name, const std::string& title, int xBins, double xMin, double xMax, int yBins, double yMin, double yMax) {
            // Check if a histogram with the same name already exists
            if (out->Get(name.c_str()) != nullptr) {
                std::cerr << "\n[ERROR] Duplicate 2D histogram detected: " << name << std::endl;
                std::cerr << "A histogram with this name already exists in the output file." << std::endl;
                std::cerr << "Aborting Fun4All macro to avoid further conflicts." << std::endl;
                exit(EXIT_FAILURE);
            }
            
            if (verbose) {
                std::cout << ANSI_COLOR_RED_BOLD << "Creating 2D histogram: " << name << " - " << title << ANSI_COLOR_RESET << std::endl;
            }

            TH2F* hist = new TH2F(name.c_str(), title.c_str(), xBins, xMin, xMax, yBins, yMin, yMax);
            hist->SetDirectory(out); // Ensure it is linked to the output file
            return hist;
        };

        /*
         EMCal QA
         */
        qaHistograms["h2_EMCal_TowerEtaPhi_Weighted_2D_" + triggerName] = create2DHistogram("h2_EMCal_TowerEtaPhi_2D_" + triggerName, "EMCal Weighted Tower Energy; #eta; #phi; Energy [GeV]", 96, 0, 96, 256, 0, 256);
        qaHistograms["hTotalCaloEEMCal_" + triggerName] = createHistogram("hTotalCaloEEMCal_" + triggerName, "Total EMCal Energy; Energy (GeV)", 100, -50, 100);
        qaHistograms["hTotalCalo_Negative_EEMCal_" + triggerName] = createHistogram("hTotalCalo_Negative_EEMCal_" + triggerName, "Total EMCal Energy; Energy (GeV)", 100, 0, 100);
        qaHistograms["h_emcalChi2_" + triggerName] = createHistogram("h_emcalChi2_" + triggerName, "Cluster Chi2; Chi2", 100, 0, 100);
        qaHistograms["h_maxTowerEnergy_" + triggerName] = createHistogram("h_maxTowerEnergy_" + triggerName, "Max Tower Energy [GeV]; Energy [GeV]", 100, 0, 100);
        

        /*
         Cluster Distributions EMCal
         */
        qaHistograms["hClusterChi2_" + triggerName] = createHistogram("hClusterChi2_" + triggerName, "Cluster Chi2; Chi2", 100, 0, 100);
        qaHistograms["hCluster_maxECore_" + triggerName] = createHistogram("hCluster_maxECore_" + triggerName, "Max Cluster ECore; Cluster ECore [GeV]", 40, 0, 20);
        qaHistograms["hClusterPt_" + triggerName] = createHistogram("hClusterPt_" + triggerName, "Cluster pT; Cluster pT [GeV]", 100, 0, 100);
        qaHistograms["hVtxZ_" + triggerName] = createHistogram("hVtxZ_" + triggerName, "Z-vertex Distribution; z [cm]", 100, -70, 70);
        qaHistograms["h_ET_" + triggerName] = createHistogram("h_ET_" + triggerName, "Cluster Transverse Energy [GeV]; Energy [GeV]", 100, 0, 100);
        qaHistograms["h2_cluster_iso_Et_unsubtracted_" + triggerName] =
            create2DHistogram("h2_cluster_iso_Et_unsubtracted_" + triggerName,
                              "Cluster Isolation Energy vs Cluster Et (Unsubtracted);Cluster Et [GeV];E_{T}^{iso} [GeV]",
                              100, 0, 20, 100, -20, 20);

        qaHistograms["h1_isoEt_unsubtracted_" + triggerName] =
            createHistogram("h1_isoEt_unsubtracted_" + triggerName,
                            "Isolation Energy Distribution (Unsubtracted);E_{T}^{iso} [GeV];Counts",
                            100, -20, 20);

        // Initialize histograms for the subtracted version
        qaHistograms["h2_cluster_iso_Et_subtracted_" + triggerName] =
            create2DHistogram("h2_cluster_iso_Et_subtracted_" + triggerName,
                              "Cluster Isolation Energy vs Cluster Et (Subtracted);Cluster Et [GeV];E_{T}^{iso} [GeV]",
                              100, 0, 20, 100, -20, 20);

        qaHistograms["h1_isoEt_subtracted_" + triggerName] =
            createHistogram("h1_isoEt_subtracted_" + triggerName,
                            "Isolation Energy Distribution (Subtracted);E_{T}^{iso} [GeV];Counts",
                            100, -20, 20);

        
        /*
         HCal QA
         */
        //inner
        qaHistograms["h2_IHCal_TowerEnergy_" + triggerName] = create2DHistogram("h2_IHCal_TowerEnergy_" + triggerName, "HCal Tower Energy; #eta; #phi; Energy [GeV]", 24, 0, 24, 64, 0, 64);
        qaHistograms["hTotalCaloEIHCal_" + triggerName] = createHistogram("hTotalCaloEIHCal_" + triggerName, "Total IHCal Energy; Energy (GeV)", 100, 0, 500);
        qaHistograms["h_ihcalChi2_" + triggerName] = createHistogram("h_ihcalChi2_" + triggerName, "Cluster Chi2; Chi2", 100, 0, 100);
        
        //outer
        qaHistograms["h2_OHCal_TowerEnergy_" + triggerName] = create2DHistogram("h2_OHCal_TowerEnergy_" + triggerName, "HCal Tower Energy; #eta; #phi; Energy [GeV]", 24, 0, 24, 64, 0, 64);
        qaHistograms["hTotalCaloEOHCal_" + triggerName] = createHistogram("hTotalCaloEOHCal_" + triggerName, "Total OHCal Energy; Energy (GeV)", 100, 0, 500);
        qaHistograms["h_ohcalChi2_" + triggerName] = createHistogram("h_ohcalChi2_" + triggerName, "Cluster Chi2; Chi2", 100, 0, 100);
        
        /*
         Trigger QA Distributions
         */
        //use for turn-on curves
        qaHistograms["h8by8TowerEnergySum_" + triggerName] = createHistogram("h8by8TowerEnergySum_" + triggerName, "Max 8x8 Tower Energy Sum; Energy [GeV]; Events", 40, 0, 20); //photon triggers
        qaHistograms["h_jet_energy_" + triggerName] = createHistogram("h_jet_energy_" + triggerName, "Maximum 0.8x0.8 Energy Sum (EMCAL + HCAL) [GeV]; Events", 50, 0, 50); //jet triggers
        
        //other possible trigger QA
        qaHistograms["h_hcal_energy_" + triggerName] = createHistogram("h_hcal_energy_" + triggerName, "Max HCal Tower Energy Sums; Energy [GeV]; Events", 40, 0, 20);
        qaHistograms["h_jet_emcal_energy_" + triggerName] = createHistogram("h_jet_emcal_energy_" + triggerName, "hJetEmcalEnergy [GeV]", 40, 0, 20);
        qaHistograms["h_jet_hcalin_energy_" + triggerName] = createHistogram("h_jet_hcalin_energy_" + triggerName, "hJetEmcalEnergy [GeV]", 40, 0, 20);
        qaHistograms["h_jet_hcalout_energy_" + triggerName] = createHistogram("h_jet_hcalout_energy_" + triggerName, "hJetEmcalEnergy [GeV]", 40, 0, 20);
        // Initialize the trigger count histogram
        qaHistograms["hTriggerCount_" + triggerName] = new TH1F(
            ("hTriggerCount_" + triggerName).c_str(),
            ("Trigger Count for Index " + triggerName + "; Count; Entries").c_str(),
            1, 0, 1 // Single bin to count occurrences
        );
        
        std::map<std::pair<float, float>, std::map<std::string, TObject*>>& qaIsolationHistograms = qaIsolationHistogramsByTriggerAndPt[triggerIndex];

        for (const auto& pT_bin : pT_bins) {
            float pT_min = pT_bin.first;
            float pT_max = pT_bin.second;
            std::pair<float, float> pT_range = {pT_min, pT_max};

            // Create unique names for the unsubtracted histograms based on the pT range
            std::string hist2DName_unsubtracted = "h2_cluster_iso_Et_unsubtracted_pT_" +
                formatFloatForFilename(pT_min) + "to" +
                formatFloatForFilename(pT_max) + "_" + triggerName;
            std::string hist1DName_unsubtracted = "h1_isoEt_unsubtracted_pT_" +
                formatFloatForFilename(pT_min) + "to" +
                formatFloatForFilename(pT_max) + "_" + triggerName;

            // Create and store the unsubtracted histograms in the existing map
            qaIsolationHistograms[pT_range][hist2DName_unsubtracted] = create2DHistogram(
                hist2DName_unsubtracted,
                "Cluster Isolation Energy vs Cluster Et (Unsubtracted);Cluster Et [GeV];E_{T}^{iso} [GeV]",
                100, 0, 20, 100, -20, 20);
            qaIsolationHistograms[pT_range][hist1DName_unsubtracted] = createHistogram(
                hist1DName_unsubtracted,
                "Isolation Energy Distribution (Unsubtracted);E_{T}^{iso} [GeV];Counts",
                100, -20, 20);

            // Create unique names for the subtracted histograms
            std::string hist2DName_subtracted = "h2_cluster_iso_Et_subtracted_pT_" +
                formatFloatForFilename(pT_min) + "to" +
                formatFloatForFilename(pT_max) + "_" + triggerName;
            std::string hist1DName_subtracted = "h1_isoEt_subtracted_pT_" +
                formatFloatForFilename(pT_min) + "to" +
                formatFloatForFilename(pT_max) + "_" + triggerName;

            // Create and store the subtracted histograms in the same map
            qaIsolationHistograms[pT_range][hist2DName_subtracted] = create2DHistogram(
                hist2DName_subtracted,
                "Cluster Isolation Energy vs Cluster Et (Subtracted);Cluster Et [GeV];E_{T}^{iso} [GeV]",
                100, 0, 20, 100, -20, 20);
            qaIsolationHistograms[pT_range][hist1DName_subtracted] = createHistogram(
                hist1DName_subtracted,
                "Isolation Energy Distribution (Subtracted);E_{T}^{iso} [GeV];Counts",
                100, -20, 20);
        }


        qaHistogramsByTrigger[triggerIndex] = qaHistograms;
        qaIsolationHistogramsByTriggerAndPt[triggerIndex] = qaIsolationHistograms;

        for (float maxAsym : asymmetry_values) {
            for (float maxChi2 : clus_chi_values) {
                for (float minClusE : clus_Ecore_values) {
                    
                    std::string invMassHistName_noBinsOfPt_name = "invMass_noPtBins_E" + formatFloatForFilename(minClusE) +
                                                  "_Chi" + formatFloatForFilename(maxChi2) +
                                                  "_Asym" + formatFloatForFilename(maxAsym) +
                                                  "_" + triggerName;
                    TH1F* invMassHistName_noBinsOfPt = new TH1F(invMassHistName_noBinsOfPt_name.c_str(), invMassHistName_noBinsOfPt_name.c_str(), 80, 0, 1.0);
                    invMassHistName_noBinsOfPt->SetTitle(";M_{#gamma#gamma};");
                    massAndIsolationHistogramsNoPtBins[triggerIndex][invMassHistName_noBinsOfPt_name] = invMassHistName_noBinsOfPt;
                    

                    for (const auto& pT_bin : pT_bins) {
                        std::string invMassHistName = "invMass_E" + formatFloatForFilename(minClusE) +
                                                      "_Chi" + formatFloatForFilename(maxChi2) +
                                                      "_Asym" + formatFloatForFilename(maxAsym) +
                                                      "_pT_" + formatFloatForFilename(pT_bin.first) + "to" + formatFloatForFilename(pT_bin.second) +
                                                      "_" + triggerName;
                        
                        if (verbose) {
                            std::cout << ANSI_COLOR_RED_BOLD << "Creating invariant mass histogram: " << invMassHistName << ANSI_COLOR_RESET << std::endl;
                        }
                        
                        TH1F* invMassHist = new TH1F(invMassHistName.c_str(), invMassHistName.c_str(), 80, 0, 1.0);
                        invMassHist->SetTitle(";M_{#gamma#gamma};");
                        
                        massAndIsolationHistograms[triggerIndex][std::make_tuple(maxAsym, maxChi2, minClusE)][pT_bin][invMassHistName] = invMassHist;
                        
                        // Define pion mass vs isolation energy histogram
                        std::string pionHistName = "pionMass_vs_isoEt_E" + formatFloatForFilename(minClusE) +
                                                   "_Chi" + formatFloatForFilename(maxChi2) +
                                                   "_Asym" + formatFloatForFilename(maxAsym) +
                                                   "_pT_" + formatFloatForFilename(pT_bin.first) + "to" + formatFloatForFilename(pT_bin.second) +
                                                   "_" + triggerName;

                        TH2F* pionMassVsIsoHist = new TH2F(pionHistName.c_str(),
                                                           "Pion Mass vs Isolation Energy;M_{#pi^{0}} [GeV];E_{T}^{iso} [GeV]",
                                                           80, 0, 1, 100, -20, 20);
                        massAndIsolationHistograms[triggerIndex][std::make_tuple(maxAsym, maxChi2, minClusE)][pT_bin][pionHistName] = pionMassVsIsoHist;

                        // Define eta mass vs isolation energy histogram
                        std::string etaHistName = "etaMass_vs_isoEt_E" + formatFloatForFilename(minClusE) +
                                                  "_Chi" + formatFloatForFilename(maxChi2) +
                                                  "_Asym" + formatFloatForFilename(maxAsym) +
                                                  "_pT_" + formatFloatForFilename(pT_bin.first) + "to" + formatFloatForFilename(pT_bin.second) +
                                                  "_" + triggerName;

                        TH2F* etaMassVsIsoHist = new TH2F(etaHistName.c_str(),
                                                          "Eta Mass vs Isolation Energy;M_{#eta} [GeV];E_{T}^{iso} [GeV]",
                                                          80, 0, 1, 100, -20, 20);
                        massAndIsolationHistograms[triggerIndex][std::make_tuple(maxAsym, maxChi2, minClusE)][pT_bin][etaHistName] = etaMassVsIsoHist;
                        

                        for (const std::string& massWindowLabel : {IN_MASS_WINDOW_LABEL, OUTSIDE_MASS_WINDOW_LABEL}) {
                            // Adjust histogram names to include massWindowLabel
                            std::string hist2DName = "h2_cluster_iso_Et_E" + formatFloatForFilename(minClusE) +
                                                      "_Chi" + formatFloatForFilename(maxChi2) +
                                                      "_Asym" + formatFloatForFilename(maxAsym) +
                                                      massWindowLabel +
                                                      "_pT_" + formatFloatForFilename(pT_bin.first) + "to" + formatFloatForFilename(pT_bin.second) +
                                                      "_" + triggerName;
                             
                             std::string hist1DName = "h1_isoEt_E" + formatFloatForFilename(minClusE) +
                                                      "_Chi" + formatFloatForFilename(maxChi2) +
                                                      "_Asym" + formatFloatForFilename(maxAsym) +
                                                      massWindowLabel +
                                                      "_pT_" + formatFloatForFilename(pT_bin.first) + "to" + formatFloatForFilename(pT_bin.second) +
                                                      "_" + triggerName;

                            // Create and store the 2D histogram (Isolation energy vs Ecore)
                            TH2F* hist2D = new TH2F(hist2DName.c_str(),
                                                  "Cluster Isolation Energy vs Cluster Et;Cluster Et [GeV];E_{T}^{iso} [GeV]",
                                                  100, 0, 20, 100, -20, 20);
                            massAndIsolationHistograms[triggerIndex][std::make_tuple(maxAsym, maxChi2, minClusE)][pT_bin][hist2DName] = hist2D;

                            // Create and store the 1D histogram (Isolation energy)
                            TH1F* hist1D = new TH1F(hist1DName.c_str(),
                                                  "Isolation Energy Distribution;E_{T}^{iso} [GeV];Counts",
                                                  100, -20, 20);
                            massAndIsolationHistograms[triggerIndex][std::make_tuple(maxAsym, maxChi2, minClusE)][pT_bin][hist1DName] = hist1D;

                            
                            for (const auto& isoRange : isoEtRanges) {
                                float isoMin = isoRange.first;
                                float isoMax = isoRange.second;
                                
                                std::string isolatedPhotonHistName = "isolatedPhotonCount_E" + formatFloatForFilename(minClusE) +
                                                                     "_Chi" + formatFloatForFilename(maxChi2) +
                                                                     "_Asym" + formatFloatForFilename(maxAsym) +
                                                                     massWindowLabel +
                                                                     "_isoEt_" + formatFloatForFilename(isoMin) + "to" + formatFloatForFilename(isoMax) +
                                                                     "_pT_" + formatFloatForFilename(pT_bin.first) + "to" + formatFloatForFilename(pT_bin.second) +
                                                                     "_" + triggerName;
                                
                                TH1F* isolatedPhotonHist = new TH1F(isolatedPhotonHistName.c_str(), "Isolated Photon Count; pT [GeV]; Count", 100, pT_bin.first, pT_bin.second);
                                
                                massAndIsolationHistograms[triggerIndex][std::make_tuple(maxAsym, maxChi2, minClusE)][pT_bin][isolatedPhotonHistName] = isolatedPhotonHist;
                            }
                            
                            // Create the histogram for all photons from pi0 decays
                            std::string allPhotonHistName = "allPhotonCount_E" + formatFloatForFilename(minClusE) +
                                                            "_Chi" + formatFloatForFilename(maxChi2) +
                                                            "_Asym" + formatFloatForFilename(maxAsym) +
                                                            massWindowLabel +
                                                            "_pT_" + formatFloatForFilename(pT_bin.first) + "to" + formatFloatForFilename(pT_bin.second) +
                                                            "_" + triggerName;
                            
                            TH1F* allPhotonHist = new TH1F(allPhotonHistName.c_str(), "All Photon Count; pT [GeV]; Count", 100, pT_bin.first, pT_bin.second);

                            // Create histogram for the pT distribution of all photons
                            std::string ptPhotonHistName = "ptPhoton_E" + formatFloatForFilename(minClusE) +
                                                           "_Chi" + formatFloatForFilename(maxChi2) +
                                                           "_Asym" + formatFloatForFilename(maxAsym) +
                                                           massWindowLabel +
                                                           "_pT_" + formatFloatForFilename(pT_bin.first) + "to" + formatFloatForFilename(pT_bin.second) +
                                                           "_" + triggerName;

                            TH1F* ptPhotonHist = new TH1F(ptPhotonHistName.c_str(), "pT of Photons; pT [GeV]; Count", 100, pT_bin.first, pT_bin.second);
                            
                            // Store these histograms in the massAndIsolationHistograms structure
                            massAndIsolationHistograms[triggerIndex][std::make_tuple(maxAsym, maxChi2, minClusE)][pT_bin][allPhotonHistName] = allPhotonHist;
                            massAndIsolationHistograms[triggerIndex][std::make_tuple(maxAsym, maxChi2, minClusE)][pT_bin][ptPhotonHistName] = ptPhotonHist;
                            
                        }
                    }
                }
            }
        }
    }
    //so that the histos actually get written out
    Fun4AllServer *se = Fun4AllServer::instance();
    if (verbose) {
        se -> Print("NODETREE");
    }
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


caloTreeGen::EnergyMaps caloTreeGen::processEnergyMaps(const std::vector<float>* m_emcTowE, const std::vector<float>* m_emciEta, const std::vector<float>* m_emciPhi, const std::vector<float>* m_ohcTowE, const std::vector<float>* m_ohciTowEta, const std::vector<float>* m_ohciTowPhi, const std::vector<float>* m_ihcTowE, const std::vector<float>* m_ihciTowEta, const std::vector<float>* m_ihciTowPhi, std::vector<short>* m_emcal_good, std::vector<short>* m_ohc_good, std::vector<short>* m_ihc_good, std::vector<int> activeTriggerBits) {
    
    if (verbose) {
        std::cout << ANSI_COLOR_BLUE_BOLD << "Processing tower energy maps: " << ANSI_COLOR_RESET << std::endl;
        std::cout << "EMCal Towers: " << m_emcTowE->size()
                  << ", IHCal Towers: " << m_ihcTowE->size()
                  << ", OHCal Towers: " << m_ohcTowE->size() << std::endl;
    }
    EnergyMaps result;
    // Initialize energy maps
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
    result.jet_ebin = 0;
    result.jet_pbin = 0;
    
    result.max_8by8energy_emcal = 0.0;
    result.max_energy_hcal = 0.0;
    result.max_energy_jet = 0.0;

    // Loop over eta bins and phi bins to find the maximum energy in EMCal
    for (int j = 0; j < 32; j++) {
        for (int k = 0; k < 12; k++) {
            if (k < 9) {
                if (energymap_jet[k][j] > result.max_energy_jet) {
                    result.max_energy_jet = energymap_jet[k][j];
                    result.jet_ebin = k;  // Store jet eta bin
                    result.jet_pbin = j;  // Store jet phi bin
                }
            }
            if (energymap[k][j] > result.max_8by8energy_emcal) {
                result.max_8by8energy_emcal = energymap[k][j];
                //in case want good trigger map
//                ebin = k;
//                pbin = j;
            }
            if (energymap_hcalout[k][j] + energymap_hcalin[k][j] > result.max_energy_hcal) {
                result.max_energy_hcal = energymap_hcalin[k][j] + energymap_hcalout[k][j];
            }
        }
    }
    // Copy the energy maps for jets
    std::copy(&energymap_jet_emcal[0][0], &energymap_jet_emcal[0][0] + 9 * 32, &result.energymap_jet_emcal[0][0]);
    std::copy(&energymap_jet_hcalin[0][0], &energymap_jet_hcalin[0][0] + 9 * 32, &result.energymap_jet_hcalin[0][0]);
    std::copy(&energymap_jet_hcalout[0][0], &energymap_jet_hcalout[0][0] + 9 * 32, &result.energymap_jet_hcalout[0][0]);

    return result;
}

void caloTreeGen::processClusterIsolationHistograms(
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
    const std::map<int, std::pair<float, float>>& clusterEtIsoMap_unsubtracted,
    std::map<std::pair<float, float>, std::map<std::string, TObject*>>& cutHistMap,
    size_t& filledHistogramCount,
    bool& filledHistogram,
    bool verbose,
    float pionMass,
    float pionMassWindow,
    float etaMass,
    float etaMassWindow,
    const std::pair<float, float>& pT_bin) {
    
    if (clusterEtIsoMap_unsubtracted.count(clusterID)) {
        float cluster_et_fromMap = clusterEtIsoMap_unsubtracted.at(clusterID).first;
        float isoEt_FromMap = clusterEtIsoMap_unsubtracted.at(clusterID).second;
        
        // Check pion mass window
        if (fabs(mesonMass - pionMass) <= pionMassWindow) {
            std::string pionHistName = "pionMass_vs_isoEt_E" + formatFloatForFilename(minClusEcore) +
            "_Chi" + formatFloatForFilename(maxChi2) +
            "_Asym" + formatFloatForFilename(maxAsym) +
            "_pT_" + formatFloatForFilename(pT_min) + "to" + formatFloatForFilename(pT_max) +
            "_" + triggerName;
            TH2F* pionMassVsIsoHist = dynamic_cast<TH2F*>(cutHistMap[pT_bin][pionHistName]);
            if (pionMassVsIsoHist) {
                pionMassVsIsoHist->Fill(mesonMass, isoEt_FromMap);
                if (verbose) {
                    std::cout << "Filled pion mass vs isolation energy histogram for pion with mass " << mesonMass << " and isolation energy " << isoEt_FromMap << std::endl;
                }
            }
        }
        
        if (fabs(mesonMass - etaMass) <= etaMassWindow) {
            std::string etaHistName = "etaMass_vs_isoEt_E" + formatFloatForFilename(minClusEcore) +
            "_Chi" + formatFloatForFilename(maxChi2) +
            "_Asym" + formatFloatForFilename(maxAsym) +
            "_pT_" + formatFloatForFilename(pT_min) + "to" + formatFloatForFilename(pT_max) +
            "_" + triggerName;
            TH2F* etaMassVsIsoHist = dynamic_cast<TH2F*>(cutHistMap[pT_bin][etaHistName]);
            if (etaMassVsIsoHist) {
                etaMassVsIsoHist->Fill(mesonMass, isoEt_FromMap);
                if (verbose) {
                    std::cout << "Filled eta mass vs isolation energy histogram for eta with mass " << mesonMass << " and isolation energy " << isoEt_FromMap << std::endl;
                }
            }
        }
        
        if (verbose) {
            std::cout << "Cluster Et: " << cluster_et_fromMap << ", IsoEt: " << isoEt_FromMap << std::endl;
        }
        
        std::string hist2DName = "h2_cluster_iso_Et_E" + formatFloatForFilename(minClusEcore) +
                                 "_Chi" + formatFloatForFilename(maxChi2) +
                                 "_Asym" + formatFloatForFilename(maxAsym) +
                                 massWindowLabel +
                                 "_pT_" + formatFloatForFilename(pT_min) + "to" + formatFloatForFilename(pT_max) +
                                 "_" + triggerName;

        std::string hist1DName = "h1_isoEt_E" + formatFloatForFilename(minClusEcore) +
                                 "_Chi" + formatFloatForFilename(maxChi2) +
                                 "_Asym" + formatFloatForFilename(maxAsym) +
                                 massWindowLabel +
                                 "_pT_" + formatFloatForFilename(pT_min) + "to" + formatFloatForFilename(pT_max) +
                                 "_" + triggerName;

        
        if (verbose) {
            std::cout << "Attempting to fill histograms: " << hist2DName << " and " << hist1DName << std::endl;
        }
        
        
        // Retrieve or create the histograms
        auto hist2D = dynamic_cast<TH2F*>(cutHistMap[pT_bin][hist2DName]);
        auto hist1D = dynamic_cast<TH1F*>(cutHistMap[pT_bin][hist1DName]);
        
        if (hist2D && hist1D) {
            // Fill the histograms
            hist2D->Fill(cluster_et_fromMap, isoEt_FromMap);
            hist1D->Fill(isoEt_FromMap);
            filledHistogramCount++;
            filledHistogram = true;

            if (verbose) {
                std::cout << "Filled histograms " << hist2DName << " and " << hist1DName << std::endl;
            }
        } else {
            std::cerr << "Error: Histograms for isolation energy are null." << std::endl;
        }
    }
}


void caloTreeGen::processIsolationRanges(
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
    const std::map<int, std::pair<float, float>>& clusterEtIsoMap_unsubtracted,
    std::map<std::pair<float, float>, std::map<std::string, TObject*>>& cutHistMap,
    bool& filledHistogram,
    bool verbose,
    const std::pair<float, float>& pT_bin) {
    
    for (const auto& isoRange : isoEtRanges) {
        float isoMin = isoRange.first;
        float isoMax = isoRange.second;
        
        if (verbose) {
            std::cout << "Processing isolation Et range: " << isoMin << " to " << isoMax << std::endl;
        }
        // Get the cluster IDs and check for isolation energy in the defined ranges
        bool clus1_isolated = clusterEtIsoMap_unsubtracted.count(clusterIDs[clus1]) &&
                              (clusterEtIsoMap_unsubtracted.at(clusterIDs[clus1]).second >= isoMin &&
                               clusterEtIsoMap_unsubtracted.at(clusterIDs[clus1]).second < isoMax);

        bool clus2_isolated = clusterEtIsoMap_unsubtracted.count(clusterIDs[clus2]) &&
                              (clusterEtIsoMap_unsubtracted.at(clusterIDs[clus2]).second >= isoMin &&
                               clusterEtIsoMap_unsubtracted.at(clusterIDs[clus2]).second < isoMax);
        if (verbose) {
            std::cout << "clus1_isolated: " << clus1_isolated << ", clus2_isolated: " << clus2_isolated << std::endl;
        }
        
        if (clus1_isolated || clus2_isolated) {
            if (verbose) {
                std::cout << "At least one cluster is isolated in this isoEt range." << std::endl;
            }
            std::string isolatedPhotonHistName = "isolatedPhotonCount_E" + formatFloatForFilename(minClusEcore) +
                                                 "_Chi" + formatFloatForFilename(maxChi2) +
                                                 "_Asym" + formatFloatForFilename(maxAsym) +
                                                 massWindowLabel +
                                                 "_isoEt_" + formatFloatForFilename(isoMin) + "to" + formatFloatForFilename(isoMax) +
                                                 "_pT_" + formatFloatForFilename(pT_min) + "to" + formatFloatForFilename(pT_max) +
                                                 "_" + triggerName;
            if (verbose) {
                std::cout << "Attempting to fill isolated photon histogram: " << isolatedPhotonHistName << std::endl;
            }
            TH1F* isolatedPhotonHist = dynamic_cast<TH1F*>(cutHistMap[pT_bin][isolatedPhotonHistName]);
            if (isolatedPhotonHist) {
                if (clus1_isolated) {
                    isolatedPhotonHist->Fill(1);
                    filledHistogram = true;
                    if (verbose) {
                        std::cout << "Filled isolatedPhotonHist for clus1." << std::endl;
                    }
                }
                if (clus2_isolated) {
                    isolatedPhotonHist->Fill(1);
                    filledHistogram = true;
                    if (verbose) {
                        std::cout << "Filled isolatedPhotonHist for clus2." << std::endl;
                    }
                }
            }
        }
    }
}

void caloTreeGen::fillHistogramsForTriggers(
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
    const std::map<int, std::pair<float, float>>& clusterEtIsoMap_unsubtracted,
    const std::vector<int>& activeTriggerBits,
    bool& filledHistogram) {
    
    /*
     temporary -- next should do the invariant mass analysis in first passs go back through using proper mass window for each pT bin
     */
    const float defaultPionMass = 0.15;
    const float defaultPionMassWindow = 0.02;
    const float defaultEtaMass = 0.59;
    const float defaultEtaMassWindow = 0.05;
    
    for (int triggerIndex : triggerIndices) {
        if (!checkTriggerCondition(activeTriggerBits, triggerIndex)) {
            if (verbose) {
                std::cout << "Trigger index " << triggerIndex << " is not active. Skipping." << std::endl;
            }
            continue;  // Skip if this trigger bit is not active
        }
        // Get the trigger name from the active trigger name map
        std::string triggerName;
        if (activeTriggerNameMap && activeTriggerNameMap->find(triggerIndex) != activeTriggerNameMap->end()) {
            triggerName = activeTriggerNameMap->at(triggerIndex);
        } else {
            std::cerr << "[ERROR] Trigger index " << triggerIndex << " not found in active trigger map." << std::endl;
            continue; // Skip to next triggerIndex
        }

        
        if (verbose) {
            std::cout << "Processing trigger index: " << triggerIndex << std::endl;
        }
        /*
         Filling Histograms outside of pT binning
         */
        auto& noPtBinHistMap = massAndIsolationHistogramsNoPtBins[triggerIndex];
        
        std::string invMassHistName_noBinsOfPt_name = "invMass_noPtBins_E" + formatFloatForFilename(minClusEcore) +
                               "_Chi" + formatFloatForFilename(maxChi2) +
                               "_Asym" + formatFloatForFilename(maxAsym) +
                               "_" + triggerName;
        
        if (verbose) {
            std::cout << "Attempting to fill histogram: " << invMassHistName_noBinsOfPt_name << std::endl;
        }

        // Attempt to find and fill the histogram
        TH1F* invMassHistName_noBinsOfPt = dynamic_cast<TH1F*>(noPtBinHistMap[invMassHistName_noBinsOfPt_name]);

        if (!invMassHistName_noBinsOfPt) {
            std::cerr << "Error: Histogram " << invMassHistName_noBinsOfPt_name
                      << " not found when trying to fill for Trigger Index: "
                      << triggerIndex << std::endl;
            continue;
        }
        
        invMassHistName_noBinsOfPt->Fill(mesonMass);
        filledHistogramCount++;
        
        if (verbose) {
            std::cout << "Filled histogram " << invMassHistName_noBinsOfPt_name << " with meson mass " << mesonMass << std::endl;
        }
        
        auto& cutHistMap = massAndIsolationHistograms[triggerIndex][std::make_tuple(maxAsym, maxChi2, minClusEcore)];

        for (const auto& pT_bin : pT_bins) {
            float pT_min = pT_bin.first;
            float pT_max = pT_bin.second;

            if (verbose) {
                std::cout << "Processing pT bin: " << pT_min << " to " << pT_max << std::endl;
            }
            
            if ((pt1 >= pT_min && pt1 < pT_max) || (pt2 >= pT_min && pt2 < pT_max)) {
                if (verbose) {
                    std::cout << "Cluster pt1: " << pt1 << ", pt2: " << pt2 << std::endl;
                    std::cout << "At least one cluster is within the pT bin." << std::endl;
                }
                // Fill invariant mass histogram for the current pT bin
                std::string invMassHistName = "invMass_E" + formatFloatForFilename(minClusEcore) +
                                              "_Chi" + formatFloatForFilename(maxChi2) +
                                              "_Asym" + formatFloatForFilename(maxAsym) +
                                              "_pT_" + formatFloatForFilename(pT_min) + "to" + formatFloatForFilename(pT_max) +
                                              "_" + triggerName;

                TH1F* invMassHist = dynamic_cast<TH1F*>(cutHistMap[pT_bin][invMassHistName]);
                if (invMassHist) {
                    invMassHist->Fill(mesonMass);
                    if (verbose) {
                        std::cout << "Filled invariant mass histogram " << invMassHistName << " with meson mass " << mesonMass << std::endl;
                    }
                }
                
                // Initialize mass values with defaults, adjusting windows to ±3σ
                float pionMass = defaultPionMass;
                float pionMassWindow = defaultPionMassWindow * 3;
                float etaMass = defaultEtaMass;
                float etaMassWindow = defaultEtaMassWindow * 3;

                if (verbose) {
                    std::cout << "Initialized with default mass values (using ±3σ range):" << std::endl;
                    std::cout << " - Pion Mass: " << pionMass << " ± " << pionMassWindow << std::endl;
                    std::cout << " - Eta Mass: " << etaMass << " ± " << etaMassWindow << std::endl;
                }

                auto massWindowKey = std::make_tuple(triggerName, minClusEcore, maxChi2, maxAsym, pT_min, pT_max);
                if (!mesonMassWindowsMap.empty() && mesonMassWindowsMap.count(massWindowKey)) {
                    const MesonMassWindow& massWindow = mesonMassWindowsMap[massWindowKey];

                    if (verbose) {
                        std::cout << "Mass window found in map for trigger index " << triggerIndex << ":" << std::endl;
                        std::cout << " - Mean Pi0: " << massWindow.meanPi0 << ", Sigma Pi0: " << massWindow.sigmaPi0 << std::endl;
                        std::cout << " - Mean Eta: " << massWindow.meanEta << ", Sigma Eta: " << massWindow.sigmaEta << std::endl;
                    }

                    // Apply custom values if within valid range
                    if ((massWindow.meanPi0 >= 0.12 && massWindow.meanPi0 <= 0.3) &&
                        (massWindow.meanEta >= 0.4 && massWindow.meanEta <= 0.7)) {
                        pionMass = massWindow.meanPi0;
                        pionMassWindow = massWindow.sigmaPi0 * 3;  // Use ±3σ for custom values
                        etaMass = massWindow.meanEta;
                        etaMassWindow = massWindow.sigmaEta * 3;

                        if (verbose) {
                            std::cout << "Custom mass window within valid range; updated mass values to ±3σ:" << std::endl;
                            std::cout << " - Pion Mass: " << pionMass << " ± " << pionMassWindow << std::endl;
                            std::cout << " - Eta Mass: " << etaMass << " ± " << etaMassWindow << std::endl;
                        }
                    } else if (verbose) {
                        std::cout << "Mass window outside valid range; retaining default ±3σ values." << std::endl;
                    }
                } else if (verbose) {
                    std::cout << "No valid mass window found in map; using default ±3σ values." << std::endl;
                }

                // Determine if meson mass is within the ±3σ mass window
                bool isInMassWindow = (fabs(mesonMass - pionMass) <= pionMassWindow) ||
                                      (fabs(mesonMass - etaMass) <= etaMassWindow);

                if (verbose) {
                    std::cout << "Meson mass " << mesonMass << (isInMassWindow ? " is " : " is not ")
                              << "within the ±3σ mass window." << std::endl;
                }
                
                std::string massWindowLabel = isInMassWindow ? "_inMassWindow" : "_outsideMassWindow";

                if (verbose) {
                    if (isInMassWindow) {
                        std::cout << "Meson mass " << mesonMass << " is within pion or eta mass window." << std::endl;
                    } else {
                        std::cout << "Meson mass " << mesonMass << " is outside pion and eta mass window." << std::endl;
                    }
                }
                for (size_t clusterIndex : {clus1, clus2}) {
                    int clusterID = clusterIDs[clusterIndex];
                    if (verbose) {
                        std::cout << "Processing cluster ID: " << clusterID << std::endl;
                    }
                    processClusterIsolationHistograms(
                         clusterID,
                         mesonMass,
                         minClusEcore,
                         maxChi2,
                         maxAsym,
                         massWindowLabel,
                         pT_min,
                         pT_max,
                         triggerIndex,
                         triggerName,
                         clusterEtIsoMap_unsubtracted,
                         cutHistMap,
                         filledHistogramCount,
                         filledHistogram,
                         verbose,
                         pionMass,
                         pionMassWindow,
                         etaMass,
                         etaMassWindow,
                         pT_bin
                     );
                }
                processIsolationRanges(
                    isoEtRanges,
                    clusterIDs,
                    clus1,
                    clus2,
                    minClusEcore,
                    maxChi2,
                    maxAsym,
                    massWindowLabel,
                    pT_min,
                    pT_max,
                    triggerIndex,
                    triggerName,
                    clusterEtIsoMap_unsubtracted,
                    cutHistMap,
                    filledHistogram,
                    verbose,
                    pT_bin
                );
                // Fill all photons histogram
                std::string allPhotonHistName = "allPhotonCount_E" + formatFloatForFilename(minClusEcore) +
                                                "_Chi" + formatFloatForFilename(maxChi2) +
                                                "_Asym" + formatFloatForFilename(maxAsym) +
                                                massWindowLabel +
                                                "_pT_" + formatFloatForFilename(pT_min) + "to" + formatFloatForFilename(pT_max) +
                                                "_" + triggerName;

                if (verbose) {
                    std::cout << "Attempting to fill all photons histogram: " << allPhotonHistName << std::endl;
                }

                
                TH1F* allPhotonHist = dynamic_cast<TH1F*>(cutHistMap[pT_bin][allPhotonHistName]);
                if (allPhotonHist) {
                    allPhotonHist->Fill(1);
                    allPhotonHist->Fill(1); // 1 for each photon
                    filledHistogram = true;
                    if (verbose) {
                        std::cout << "Filled allPhotonHist with 2 entries." << std::endl;
                    }
                }
                
                // Fill pT distribution histograms
                std::string ptPhotonHistName = "ptPhoton_E" + formatFloatForFilename(minClusEcore) +
                                               "_Chi" + formatFloatForFilename(maxChi2) +
                                               "_Asym" + formatFloatForFilename(maxAsym) +
                                               massWindowLabel +
                                               "_pT_" + formatFloatForFilename(pT_min) + "to" + formatFloatForFilename(pT_max) +
                                               "_" + triggerName;

                if (verbose) {
                    std::cout << "Attempting to fill pT photon histogram: " << ptPhotonHistName << std::endl;
                }

                TH1F* ptPhotonHist = dynamic_cast<TH1F*>(cutHistMap[pT_bin][ptPhotonHistName]);
                if (ptPhotonHist) {
                    if (pt1 >= pT_min && pt1 < pT_max) {
                        ptPhotonHist->Fill(pt1);
                        filledHistogram = true;
                        if (verbose) {
                            std::cout << "Filled ptPhotonHist with pt1: " << pt1 << std::endl;
                        }
                    }
                    if (pt2 >= pT_min && pt2 < pT_max) {
                        ptPhotonHist->Fill(pt2);
                        filledHistogram = true;
                        if (verbose) {
                            std::cout << "Filled ptPhotonHist with pt2: " << pt2 << std::endl;
                        }
                    }
                }
            }
        }
    }
}


void caloTreeGen::processClusterInvariantMass(
    const std::vector<float>& clusterE,
    const std::vector<float>& clusterPt,
    const std::vector<float>& clusterChi2,
    const std::vector<float>& clusterEta,
    const std::vector<float>& clusterPhi,
    const std::vector<int>& clusterIDs,
    const std::map<int, std::pair<float, float>>& clusterEtIsoMap_unsubtracted,
    std::vector<int> activeTriggerBits)
{
    
    std::map<std::pair<float, float>, int> totalPhotonCountInBin;
    std::map<std::pair<float, float>, int> isolatedPhotonCountInBin;
    
    // Initialize counters for each pT bin
    for (const auto& bin : pT_bins) {
        totalPhotonCountInBin[bin] = 0;
        isolatedPhotonCountInBin[bin] = 0;
    }
    
    // Check for mismatched vector sizes, which could indicate an issue
    if (clusterE.size() != clusterPt.size() || clusterE.size() != clusterChi2.size() ||
        clusterE.size() != clusterEta.size() || clusterE.size() != clusterPhi.size()) {
        std::cerr << "Error: Mismatched cluster vector sizes" << std::endl;
        return;
    }

    // Counters for summary
    size_t totalPairs = 0;
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

            // Create Lorentz vectors for the clusters
            TLorentzVector photon1, photon2;
            photon1.SetPtEtaPhiE(pt1, eta1, phi1, E1);
            photon2.SetPtEtaPhiE(pt2, eta2, phi2, E2);
            TLorentzVector meson = photon1 + photon2;
            float mesonMass = meson.M();
            float asym = (E1 + E2 != 0) ? fabs(E1 - E2) / (E1 + E2) : 0;

            // Apply all combinations of cuts and fill histograms
            bool filledHistogram = false;
            bool pairSkipped = false;

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

                        fillHistogramsForTriggers(
                            mesonMass,
                            clus1,
                            clus2,
                            pt1,
                            pt2,
                            E1,
                            E2,
                            minClusEcore,
                            maxChi2,
                            maxAsym,
                            filledHistogramCount,
                            clusterIDs,
                            clusterEtIsoMap_unsubtracted,
                            activeTriggerBits,
                            filledHistogram
                        );
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
        std::cout << "Total pairs processed: " << totalPairs << std::endl;
        std::cout << "Pairs skipped due to asymmetry cuts: " << skippedPairsDueToAsymmetry
                  << " (" << (100.0 * skippedPairsDueToAsymmetry / totalPairs) << "%)" << std::endl;
        std::cout << "Pairs skipped due to chi2 cuts: " << skippedPairsDueToChi2
                  << " (" << (100.0 * skippedPairsDueToChi2 / totalPairs) << "%)" << std::endl;
        std::cout << "Pairs skipped due to Ecore cuts: " << skippedPairsDueToEcore
                  << " (" << (100.0 * skippedPairsDueToEcore / totalPairs) << "%)" << std::endl;
        std::cout << "Pairs with no histograms filled: " << noHistogramsFilledCount
                  << " (" << (100.0 * noHistogramsFilledCount / totalPairs) << "%)" << std::endl;
        std::cout << "Histograms filled: " << filledHistogramCount << std::endl;
        std::cout << "Cluster pairs that passed all cuts with corresponding meson masses and cuts:" << std::endl;

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
        return Fun4AllReturnCodes::ABORTRUN;
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

    /*
     Process vertex, tower information, cluster information, outside of loop
     */
    
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
    
    /*switched to UE subtracted
     use:
     -TOWERINFO_CALIB_CEMC_RETOWER_SUB1
     -TOWERINFO_CALIB_HCALIN_SUB1
     -TOWERINFO_CALIB_HCALOUT_SUB1
     */
    /*
     switch to non retowered emcal -- TOWERINFO_CALIB_CEMC
     */
    TowerInfoContainer* emcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
    TowerInfoContainer* ohcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
    TowerInfoContainer* ihcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
    if (!emcTowerContainer && !ihcTowerContainer && !ohcTowerContainer) {
        std::cout << ANSI_COLOR_RED_BOLD << "No tower containers found, skipping tower processing." << ANSI_COLOR_RESET << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;  // Return if no tower containers
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
    }
    if (ihcTowerContainer) {
        if (verbose) {
            std::cout << ANSI_COLOR_BLUE_BOLD << "Processing IHCal Towers..." << ANSI_COLOR_RESET << std::endl;
        }
        processTowers(ihcTowerContainer, totalCaloEIHCal, m_ihciTowEta, m_ihciTowPhi, m_ihcTowE, m_ihcTime, m_ihcChi2, m_ihcPed, m_ihc_good);
    }
    if (ohcTowerContainer) {
        if (verbose) {
            std::cout << ANSI_COLOR_BLUE_BOLD << "Processing OHCal Towers..." << ANSI_COLOR_RESET << std::endl;
        }
        
        processTowers(ohcTowerContainer, totalCaloEOHCal, m_ohciTowEta, m_ohciTowPhi, m_ohcTowE, m_ohcTime, m_ohcChi2, m_ohcPed, m_ohc_good);
    }
    
    // Display vector sizes and confirm entry into the function call if verbose is enabled
    if (verbose) {
        std::cout << ANSI_COLOR_RED_BOLD << "Calling processEnergyMaps with Vector Sizes: " << ANSI_COLOR_RESET << std::endl;
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
    }

    // Pass the addresses of the pointers to the function
    EnergyMaps energyMaps = processEnergyMaps(&m_emcTowE, &m_emciEta, &m_emciPhi, &m_ohcTowE, &m_ohciTowEta, &m_ohciTowPhi,
                      &m_ihcTowE, &m_ihciTowEta, &m_ihciTowPhi, &m_emcal_good, &m_ohc_good, &m_ihc_good, activeTriggerBits);
    
    if (verbose) {
        std::cout << "Fetching RawClusterContainer..." << std::endl;
    }

    RawClusterContainer *clusterContainer = findNode::getClass<RawClusterContainer>(topNode, "CLUSTERINFO_CEMC");
    if(!clusterContainer) {
        std::cout << PHWHERE << "Cluster node is missing...ruh roh scooby doo" << std::endl;
        return 0;
    }
    
    RawClusterContainer::ConstRange clusterRange = clusterContainer->getClusters();
    if (verbose) {
        std::cout << "Number of clusters to be processed: " << " = " << std::distance(clusterRange.first, clusterRange.second) << std::endl;

    }
    
    int nan_count = 0;
    int skippedEcoreCount = 0; // Counter for clusters skipped due to ECore cut
    float max_energy_clus = 0.0;
    float max_isoEt = std::numeric_limits<float>::lowest();
    float min_isoEt = std::numeric_limits<float>::max();
    std::map<int, std::pair<float, float>> clusterEtIsoMap_unsubtracted; //to store cluster ID and corresponding et iso value
    std::map<int, std::pair<float, float>> clusterEtIsoMap_subtracted;
    std::vector<int> m_clusterIds; // Store cluster IDs
    
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
        
        
        float clusEcore = Ecore_vec_cluster.mag();
        
        if (clusEcore < 1.0) { // cut on ecore
            skippedEcoreCount++;
            continue;
        }
        float clus_eta = Ecore_vec_cluster.pseudoRapidity();
        float clus_phi = Ecore_vec_cluster.phi();
        float clus_eT = clusEcore / std::cosh(clus_eta);
        float clus_pt = Ecore_vec_cluster.perp();
        float clus_chi = cluster -> get_chi2();
        float maxTowerEnergy = getMaxTowerE(cluster,emcTowerContainer);
        
        m_clusterIds.push_back(cluster->get_id());
        m_clusterECore.push_back(clusEcore);
        m_clusterEt.push_back(clus_eT);
        m_clusterPt.push_back(clus_pt);
        m_clusterChi.push_back(clus_chi);
        m_clusterPhi.push_back(clus_phi);
        m_clusterEta.push_back(clus_eta);
        m_clusTowPhi.push_back(returnClusterTowPhi(cluster,emcTowerContainer));
        m_clusTowEta.push_back(returnClusterTowEta(cluster,emcTowerContainer));
        m_clusTowE.push_back(returnClusterTowE(cluster,emcTowerContainer));
        m_maxTowEnergy.push_back(maxTowerEnergy);
        
        float et_iso_unsubtracted = cluster->get_et_iso(3, /*unsubtracted=*/0, /*clusterTower=*/1);

        // Retrieve subtracted isolation energy
        float et_iso_subtracted = cluster->get_et_iso(3, /*subtracted=*/1, /*clusterTower=*/1);


        // Check if the isolation energy is NaN
        if (!std::isnan(et_iso_unsubtracted)) {
            clusterEtIsoMap_unsubtracted[cluster->get_id()] = std::make_pair(clus_eT, et_iso_unsubtracted);
            if (et_iso_unsubtracted > max_isoEt) {
                max_isoEt = et_iso_unsubtracted;
            }
            if (et_iso_unsubtracted < min_isoEt) {
                min_isoEt = et_iso_unsubtracted;
            }
            if (verbose) {
                std::cout << "Cluster passed isolation cut: ID " << cluster->get_id()
                          << ", Et = " << clus_eT << ", Unsubtracted isoEt = " << et_iso_unsubtracted << std::endl;
            }
        } else {
            if (verbose) {
                std::cout << "Warning: Unsubtracted Isolation energy is NaN for cluster ID: " << cluster->get_id() << std::endl;
            }
            nan_count++;
        }

        if (!std::isnan(et_iso_subtracted)) {
            clusterEtIsoMap_subtracted[cluster->get_id()] = std::make_pair(clus_eT, et_iso_subtracted);
            if (verbose) {
                std::cout << "Cluster passed isolation cut: ID " << cluster->get_id()
                          << ", Et = " << clus_eT << ", Subtracted isoEt = " << et_iso_subtracted << std::endl;
            }
        } else {
            if (verbose) {
                std::cout << "Warning: UE Subtracted Isolation energy is NaN for cluster ID: " << cluster->get_id() << std::endl;
            }
            nan_count++;
        }
        
        
        if (clusEcore > max_energy_clus) {
            max_energy_clus = clusEcore; // Update the maximum cluster energy
        }
    }
    
    if (verbose) {
        std::cout << "\n--- Cluster Processing Summary ---" << std::endl;
        std::cout << "Clusters processed: " << m_clusterIds.size() << std::endl;
        std::cout << "Clusters skipped due to ECore < 1: " << skippedEcoreCount << std::endl;

        std::cout << "\nVector Sizes:" << std::endl;
        std::cout << "  m_clusterIds size: " << m_clusterIds.size() << std::endl;
        std::cout << "  m_clusterECore size: " << m_clusterECore.size() << std::endl;
        std::cout << "  m_clusterEta size: " << m_clusterEta.size() << std::endl;
        std::cout << "  m_clusterPhi size: " << m_clusterPhi.size() << std::endl;
        std::cout << "  m_clusterPt size: " << m_clusterPt.size() << std::endl;
        std::cout << "  m_clusterChi size: " << m_clusterChi.size() << std::endl;

        std::cout << "\nCluster Isolation Summary:\n";
        std::cout << "Clusters with NaN isolation energy: " << nan_count << std::endl;
        std::cout << "Size of clusterEtIsoMap_unsubtracted: " << clusterEtIsoMap_unsubtracted.size() << std::endl;
        std::cout << "Max isolation energy (isoEt): " << max_isoEt << std::endl;
        std::cout << "Min isolation energy (isoEt): " << min_isoEt << std::endl;

        std::cout << "\nCluster Isolation Energy Table:" << std::endl;
        std::cout << "------------------------------------------------" << std::endl;
        std::cout << std::setw(12) << "Cluster ID" << std::setw(20) << "Ecore" << std::setw(20) << "Isolation Energy (et_iso)" << std::endl;
        std::cout << "------------------------------------------------" << std::endl;
        for (const auto& entry : clusterEtIsoMap_unsubtracted) {
            std::cout << std::setw(12) << entry.first << std::setw(20) << entry.second.first << std::setw(20) << entry.second.second << std::endl;
        }
        std::cout << "------------------------------------------------\n";
    }
    processClusterInvariantMass(m_clusterECore, m_clusterPt, m_clusterChi, m_clusterEta, m_clusterPhi, m_clusterIds, clusterEtIsoMap_unsubtracted, activeTriggerBits);

    
    for (int triggerIndex : triggerIndices) {
        // Check if the current trigger bit is active
        if (!checkTriggerCondition(activeTriggerBits, triggerIndex)) {
            continue;  // Skip if this trigger bit is not active
        }
        std::string triggerName;
        if (activeTriggerNameMap && activeTriggerNameMap->find(triggerIndex) != activeTriggerNameMap->end()) {
            triggerName = activeTriggerNameMap->at(triggerIndex);
        } else {
            std::cerr << "[ERROR] Trigger index " << triggerIndex << " not found in active trigger map." << std::endl;
            continue; // Skip to next triggerIndex
        }
        
        auto& qaHistograms = qaHistogramsByTrigger[triggerIndex];
        std::map<std::pair<float, float>, std::map<std::string, TObject*>>& qaIsolationHistograms = qaIsolationHistogramsByTriggerAndPt[triggerIndex];

        // Check if the histogram exists before filling it
        std::string histName = "hVtxZ_" + triggerName;
        if (qaHistograms.find(histName) != qaHistograms.end()) {
            // Fill the histogram with the z-coordinate of the vertex
            ((TH1F*)qaHistograms[histName])->Fill(m_vz);
            if (verbose) {
                std::cout << "Filled histogram " << histName << " with m_vz = " << m_vz << std::endl;
            }
        } else if (verbose) {
            std::cerr << "Error: Histogram " << histName << " not found in qaHistograms!" << std::endl;
        }

        // Fill the trigger count histogram to count each time this trigger fires
        TH1F* hTriggerCount = (TH1F*)qaHistogramsByTrigger[triggerIndex]["hTriggerCount_" + triggerName];
        if (hTriggerCount) {
            hTriggerCount->Fill(0.5); // Fill the single bin to count the trigger occurrence
        }
        
        if (verbose) {
            std::cout << "Processing Trigger Bit: " << triggerIndex << std::endl;
        }
        
        // Process towers and fill histograms
        if (emcTowerContainer) {
            for (size_t i = 0; i < m_emcTowE.size(); ++i) {
                ((TH2F*)qaHistograms["h2_EMCal_TowerEtaPhi_Weighted_2D_" + triggerName])->Fill(m_emciEta[i], m_emciPhi[i], m_emcTowE[i]);
            }
            // Check if the total EMCal energy is positive or negative and fill the appropriate histogram
            if (totalCaloEEMCal >= 0) {
                ((TH1F*)qaHistograms["hTotalCaloEEMCal_" + triggerName])->Fill(totalCaloEEMCal);
            } else {
                ((TH1F*)qaHistograms["hTotalCalo_Negative_EEMCal_" + triggerName])->Fill(totalCaloEEMCal);
            }
            //maximum chi2 of all the towers fill the histograms for EMCAL chi2 for this event -- similar for other calos
            if (!m_emcChi2.empty()) {
                float max_emcalChi2 = *std::max_element(m_emcChi2.begin(), m_emcChi2.end());
                ((TH1F*)qaHistograms["h_emcalChi2_" + triggerName])->Fill(max_emcalChi2);
            }

        }
        if (ihcTowerContainer) {
            for (size_t i = 0; i < m_ihcTowE.size(); ++i) {
                ((TH2F*)qaHistograms["h2_IHCal_TowerEnergy_" + triggerName])->Fill(m_ihciTowEta[i], m_ihciTowPhi[i], m_ihcTowE[i]);
            }
            ((TH1F*)qaHistograms["hTotalCaloEIHCal_" + triggerName])->Fill(totalCaloEIHCal);
            if (!m_ihcChi2.empty()) {
                float max_ihcalChi2 = *std::max_element(m_ihcChi2.begin(), m_ihcChi2.end());
                ((TH1F*)qaHistograms["h_ihcalChi2_" + triggerName])->Fill(max_ihcalChi2);
            }
        }
        
        if (ohcTowerContainer) {
            for (size_t i = 0; i < m_ohcTowE.size(); ++i) {
                ((TH2F*)qaHistograms["h2_OHCal_TowerEnergy_" + triggerName])->Fill(m_ohciTowEta[i], m_ohciTowPhi[i], m_ohcTowE[i]);
            }
            ((TH1F*)qaHistograms["hTotalCaloEOHCal_" + triggerName])->Fill(totalCaloEOHCal);
            if (!m_ohcChi2.empty()) {
                float max_ohcalChi2 = *std::max_element(m_ohcChi2.begin(), m_ohcChi2.end());
                ((TH1F*)qaHistograms["h_ohcalChi2_" + triggerName])->Fill(max_ohcalChi2);
            }
        }
        
        ((TH1F*)qaHistograms["h8by8TowerEnergySum_" + triggerName])->Fill(energyMaps.max_8by8energy_emcal);
        ((TH1F*)qaHistograms["h_hcal_energy_" + triggerName])->Fill(energyMaps.max_energy_hcal);
        ((TH1F*)qaHistograms["h_jet_energy_" + triggerName])->Fill(energyMaps.max_energy_jet);

        ((TH1F*)qaHistograms["h_jet_emcal_energy_" + triggerName])->Fill(energyMaps.energymap_jet_emcal[energyMaps.jet_ebin][energyMaps.jet_pbin]);
        ((TH1F*)qaHistograms["h_jet_hcalin_energy_" + triggerName])->Fill(energyMaps.energymap_jet_hcalin[energyMaps.jet_ebin][energyMaps.jet_pbin]);
        ((TH1F*)qaHistograms["h_jet_hcalout_energy_" + triggerName])->Fill(energyMaps.energymap_jet_hcalout[energyMaps.jet_ebin][energyMaps.jet_pbin]);
        
        
        for (size_t i = 0; i < m_clusterIds.size(); ++i) {
            // Fill the eT histogram
            TH1F* h_ET = (TH1F*)qaHistograms["h_ET_" + triggerName];
            if (!h_ET) {
                std::cerr << "Error: h_ET_" << triggerIndex << " histogram is null!" << std::endl;
            } else {
                h_ET->Fill(m_clusterEt[i]);
            }
            
            // Check if histograms exist and fill them
            TH1F* h_maxTowE = (TH1F*)qaHistograms["h_maxTowerEnergy_" + triggerName];
            if (!h_maxTowE) {
                std::cerr << "Error: h_maxTowerEnergy_" << triggerIndex << " histogram is null!" << std::endl;
            } else {
                h_maxTowE->Fill(m_maxTowEnergy[i]);
            }
            
            // Access histograms from the map
            TH1F* hPt = (TH1F*)qaHistograms["hClusterPt_" + triggerName];
            TH1F* hChi2 = (TH1F*)qaHistograms["hClusterChi2_" + triggerName];

            // Check if histograms exist and fill them
            if (!hPt) {
                std::cerr << "Error: hClusterPt_" << triggerIndex << " histogram is null!" << std::endl;
            } else {
                hPt->Fill(m_clusterPt[i]);
            }

            if (!hChi2) {
                std::cerr << "Error: hClusterChi2_" << triggerIndex << " histogram is null!" << std::endl;
            } else {
                hChi2->Fill(m_clusterChi[i]);
            }
        }
        
        try {
            // Filling unsubtracted histograms
            for (const auto& entry : clusterEtIsoMap_unsubtracted) {
                // Get Et and isoEt from the map
                float cluster_et_fromMap = entry.second.first;
                float isoEt_FromMap = entry.second.second;

                // Unsubtracted histogram names
                std::string hist2DName = "h2_cluster_iso_Et_unsubtracted_" + triggerName;
                std::string hist1DName = "h1_isoEt_unsubtracted_" + triggerName;

                // Error checking for histogram existence before filling
                auto hist2D = qaHistograms[hist2DName];
                auto hist1D = qaHistograms[hist1DName];

                if (!hist2D) {
                    throw std::runtime_error("Error: " + hist2DName + " is null.");
                }
                if (!hist1D) {
                    throw std::runtime_error("Error: " + hist1DName + " is null.");
                }

                // Fill the 1D and 2D histograms
                ((TH2F*)hist2D)->Fill(cluster_et_fromMap, isoEt_FromMap);
                ((TH1F*)hist1D)->Fill(isoEt_FromMap);
            }

            // Loop over clusters and pT bins to fill additional unsubtracted histograms
            for (size_t i = 0; i < m_clusterIds.size(); ++i) {
                int clusterId = m_clusterIds[i];
                float clusterPt = m_clusterPt[i];

                // Check if this clusterId exists in the clusterEtIsoMap_unsubtracted
                if (clusterEtIsoMap_unsubtracted.count(clusterId)) {
                    float cluster_et_fromMap = clusterEtIsoMap_unsubtracted[clusterId].first;
                    float isoEt_FromMap = clusterEtIsoMap_unsubtracted[clusterId].second;

                    // Loop over the predefined pT bins
                    for (const auto& pT_bin : pT_bins) {
                        float pT_min = pT_bin.first;
                        float pT_max = pT_bin.second;
                        std::pair<float, float> pT_range = {pT_min, pT_max};

                        // Check if the cluster pT falls within the current bin
                        if (clusterPt >= pT_min && clusterPt < pT_max) {
                            // Unsubtracted histogram names
                            std::string hist2DName = "h2_cluster_iso_Et_unsubtracted_pT_" +
                                formatFloatForFilename(pT_min) + "to" +
                                formatFloatForFilename(pT_max) + "_" + triggerName;
                            std::string hist1DName = "h1_isoEt_unsubtracted_pT_" +
                                formatFloatForFilename(pT_min) + "to" +
                                formatFloatForFilename(pT_max) + "_" + triggerName;

                            // Retrieve the histograms from the existing map
                            auto& hist2D_pT = qaIsolationHistograms[pT_range][hist2DName];
                            auto& hist1D_pT = qaIsolationHistograms[pT_range][hist1DName];

                            if (!hist2D_pT) {
                                throw std::runtime_error("Error: " + hist2DName + " is null.");
                            }
                            if (!hist1D_pT) {
                                throw std::runtime_error("Error: " + hist1DName + " is null.");
                            }

                            // Fill the 2D and 1D histograms for the current pT bin
                            ((TH2F*)hist2D_pT)->Fill(cluster_et_fromMap, isoEt_FromMap);
                            ((TH1F*)hist1D_pT)->Fill(isoEt_FromMap);

                            if (verbose) {
                                std::cout << "Filled unsubtracted histograms for cluster ID " << clusterId
                                          << " in pT range [" << pT_min << ", " << pT_max << "] for trigger "
                                          << triggerIndex << std::endl;
                            }
                        }
                    }
                }
            }
        }
        catch (const std::exception& ex) {
            std::cerr << "Exception occurred in unsubtracted histogram filling: " << ex.what() << std::endl;
        }

        try {
            // Filling subtracted histograms
            for (const auto& entry : clusterEtIsoMap_subtracted) {
                // Get Et and isoEt from the map
                float cluster_et_fromMap = entry.second.first;
                float isoEt_FromMap = entry.second.second;

                // Subtracted histogram names
                std::string hist2DName = "h2_cluster_iso_Et_subtracted_" + triggerName;
                std::string hist1DName = "h1_isoEt_subtracted_" + triggerName;

                // Error checking for histogram existence before filling
                auto hist2D = qaHistograms[hist2DName];
                auto hist1D = qaHistograms[hist1DName];

                if (!hist2D) {
                    throw std::runtime_error("Error: " + hist2DName + " is null.");
                }
                if (!hist1D) {
                    throw std::runtime_error("Error: " + hist1DName + " is null.");
                }

                // Fill the 1D and 2D histograms
                ((TH2F*)hist2D)->Fill(cluster_et_fromMap, isoEt_FromMap);
                ((TH1F*)hist1D)->Fill(isoEt_FromMap);
            }

            // Loop over clusters and pT bins to fill additional subtracted histograms
            for (size_t i = 0; i < m_clusterIds.size(); ++i) {
                int clusterId = m_clusterIds[i];
                float clusterPt = m_clusterPt[i];

                // Check if this clusterId exists in the clusterEtIsoMap_subtracted
                if (clusterEtIsoMap_subtracted.count(clusterId)) {
                    float cluster_et_fromMap = clusterEtIsoMap_subtracted[clusterId].first;
                    float isoEt_FromMap = clusterEtIsoMap_subtracted[clusterId].second;

                    // Loop over the predefined pT bins
                    for (const auto& pT_bin : pT_bins) {
                        float pT_min = pT_bin.first;
                        float pT_max = pT_bin.second;
                        std::pair<float, float> pT_range = {pT_min, pT_max};

                        // Check if the cluster pT falls within the current bin
                        if (clusterPt >= pT_min && clusterPt < pT_max) {
                            // Subtracted histogram names
                            std::string hist2DName = "h2_cluster_iso_Et_subtracted_pT_" +
                                formatFloatForFilename(pT_min) + "to" +
                                formatFloatForFilename(pT_max) + "_" + triggerName;
                            std::string hist1DName = "h1_isoEt_subtracted_pT_" +
                                formatFloatForFilename(pT_min) + "to" +
                                formatFloatForFilename(pT_max) + "_" + triggerName;

                            // Retrieve the histograms from the existing map
                            auto& hist2D_pT = qaIsolationHistograms[pT_range][hist2DName];
                            auto& hist1D_pT = qaIsolationHistograms[pT_range][hist1DName];

                            if (!hist2D_pT) {
                                throw std::runtime_error("Error: " + hist2DName + " is null.");
                            }
                            if (!hist1D_pT) {
                                throw std::runtime_error("Error: " + hist1DName + " is null.");
                            }

                            // Fill the 2D and 1D histograms for the current pT bin
                            ((TH2F*)hist2D_pT)->Fill(cluster_et_fromMap, isoEt_FromMap);
                            ((TH1F*)hist1D_pT)->Fill(isoEt_FromMap);

                            if (verbose) {
                                std::cout << "Filled subtracted histograms for cluster ID " << clusterId
                                          << " in pT range [" << pT_min << ", " << pT_max << "] for trigger "
                                          << triggerIndex << std::endl;
                            }
                        }
                    }
                }
            }
        }
        catch (const std::exception& ex) {
            std::cerr << "Exception occurred in subtracted histogram filling: " << ex.what() << std::endl;
        }
        
        // Fill the histogram with the maximum cluster energy core value
        TH1F* h_maxECore = (TH1F*)qaHistograms["hCluster_maxECore_" + triggerName];

        if (!h_maxECore) {
            std::cerr << "Error: Histogram hCluster_maxECore_" << triggerIndex << " is null and cannot be filled." << std::endl;
        } else {
            h_maxECore->Fill(max_energy_clus);
            
            if (verbose) {
                std::cout << "Filled histogram hCluster_maxECore_" << triggerIndex
                          << " with value: " << max_energy_clus << std::endl;
            }
        }
    }
    return Fun4AllReturnCodes::EVENT_OK;
}


//____________________________________________________________________________..
int caloTreeGen::ResetEvent(PHCompositeNode *topNode) {
    if (verbose) {
        std::cout << ANSI_COLOR_BLUE_BOLD << "Resetting event..." << ANSI_COLOR_RESET << std::endl;
    }
    m_clusterE.clear();
    m_clusterEt.clear();
    m_clusterPhi.clear();
    m_clusterEta.clear();
    m_clusterPt.clear();
    m_clusterChi.clear();
    m_clusterTowMaxE.clear();
    m_clusterECore.clear();
    m_clusterEtIso.clear();
    clusterEtIsoMap_unsubtracted.clear();
    clusterEtIsoMap_subtracted.clear();
    m_clusterIds.clear();
    
    m_emcTowE.clear();
    m_emciEta.clear();
    m_emciPhi.clear();
    m_emcTime.clear();
    m_emcChi2.clear();
    m_emcPed.clear();
    m_emcal_good.clear();
    m_maxTowEnergy.clear();

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


int caloTreeGen::End(PHCompositeNode *topNode) {
    std::cout << ANSI_COLOR_BLUE_BOLD << "caloTreeGen::End(PHCompositeNode *topNode) All events have been processed. Beginning final analysis steps..." << ANSI_COLOR_RESET << std::endl;
    // Ensure the output file is open and set as the current directory
    if (out && out->IsOpen()) {
        out->cd();
        gDirectory = out; // Explicitly set the global ROOT directory
    } else {
        std::cerr << ANSI_COLOR_RED_BOLD << "Error: Output file is not open. Exiting End() without writing histograms." << ANSI_COLOR_RESET << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
    }

    auto writeHistogram = [&](TObject* hist, const std::string& name) {
        if (!hist) {
            std::cerr << ANSI_COLOR_RED_BOLD << "Warning: Histogram " << name << " is null, skipping write." << ANSI_COLOR_RESET << std::endl;
            return;
        }

        if (hist->Write() == 0) {
            std::cerr << ANSI_COLOR_RED_BOLD << "Error: Failed to write histogram " << name << "." << ANSI_COLOR_RESET << std::endl;
        } else if (verbose) {
            std::cout << ANSI_COLOR_GREEN_BOLD << "Successfully wrote histogram " << name << "." << ANSI_COLOR_RESET << std::endl;
        }

        delete hist; // Clean up after writing
    };
    
    
    // Iterate over each trigger and write the QA histograms to the correct directory
    for (const auto& [triggerIndex, qaHistograms] : qaHistogramsByTrigger) {
        // Get trigger name for the directory
        std::string triggerName;
        if (activeTriggerNameMap && activeTriggerNameMap->find(triggerIndex) != activeTriggerNameMap->end()) {
            triggerName = activeTriggerNameMap->at(triggerIndex);
        } else {
            std::cerr << "[ERROR] Trigger index " << triggerIndex << " not found in active trigger map." << std::endl;
            continue;
        }

        // Navigate to the specific trigger directory
        TDirectory* triggerDir = out->GetDirectory(triggerName.c_str());
        if (triggerDir) {
            triggerDir->cd();
        } else {
            std::cerr << ANSI_COLOR_RED_BOLD << "Warning: Directory for trigger " << triggerName << " does not exist." << ANSI_COLOR_RESET << std::endl;
            continue;
        }

        // Write QA histograms for this trigger
        for (const auto& [name, hist] : qaHistograms) {
            writeHistogram(hist, name);
        }

        // Write isolation histograms for this trigger
        if (qaIsolationHistogramsByTriggerAndPt.count(triggerIndex)) {
            for (const auto& [pT_bin, histMap] : qaIsolationHistogramsByTriggerAndPt.at(triggerIndex)) {
                for (const auto& [name, hist] : histMap) {
                    writeHistogram(hist, name);
                }
            }
        }

        // Write mass and isolation histograms with pT bins
        if (massAndIsolationHistograms.count(triggerIndex)) {
            for (const auto& [cutCombination, pTHistMap] : massAndIsolationHistograms.at(triggerIndex)) {
                for (const auto& [pT_bin, histMap] : pTHistMap) {
                    for (const auto& [histName, hist] : histMap) {
                        writeHistogram(hist, histName);
                    }
                }
            }
        }

        // Write mass histograms without pT bins
        if (massAndIsolationHistogramsNoPtBins.count(triggerIndex)) {
            for (const auto& [histName, hist] : massAndIsolationHistogramsNoPtBins.at(triggerIndex)) {
                writeHistogram(hist, histName);
            }
        }

        out->cd(); // Return to the root directory
    }
    // Close the output file and clean up
    std::cout << ANSI_COLOR_BLUE_BOLD << "Closing output file and cleaning up..." << ANSI_COLOR_RESET << std::endl;
    if (out) {
        out->Close();
        delete out;
        out = nullptr;
        if (verbose) {
            std::cout << ANSI_COLOR_GREEN_BOLD << "Output file successfully closed and deleted." << ANSI_COLOR_RESET << std::endl;
        }

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
//____________________________________________________________________________..
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
