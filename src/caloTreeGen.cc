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
#include <unordered_map>

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
caloTreeGen::caloTreeGen(const std::string &dataOutFile,
                         const std::string &simOutFile)
  : SubsysReco("CaloTreeGen")
  , Outfile(dataOutFile)
  , SimOutfile(simOutFile)
{
  std::cout << "[DEBUG] caloTreeGen::caloTreeGen() constructor called." << std::endl;
  std::cout << "    Data output will go to: " << Outfile << std::endl;
  std::cout << "    Sim  output will go to: " << SimOutfile << std::endl;
}

//____________________________________________________________________________..
caloTreeGen::~caloTreeGen() {
    std::cout << "[DEBUG] caloTreeGen::~caloTreeGen() destructor called." << std::endl;
}



//____________________________________________________________________________..
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
            float Energy, Chi2, Asym, pTMin, pTMax;
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
            std::getline(lineStream, cell, ','); trim(cell); Energy = std::stof(cell);
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
            MesonMassWindow massWindow = {triggerName, Energy, Chi2, Asym, pTMin, pTMax, meanPi0, sigmaPi0, meanEta, sigmaEta};
            auto key = std::make_tuple(triggerName, Energy, Chi2, Asym, pTMin, pTMax);

            // Insert into the map
            mesonMassWindowsMap[key] = massWindow;
            entriesLoaded++;

            if (verbose) {
                std::cout << "Loaded entry " << entriesLoaded << " (Line " << lineNumber << "):" << std::endl;
                std::cout << " - Trigger Name: " << triggerName << std::endl;
                std::cout << " - Energy: " << Energy << ", Chi2: " << Chi2 << ", Asym: " << Asym << std::endl;
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

//____________________________________________________________________________..
int caloTreeGen::Init(PHCompositeNode *topNode) {
    if (verbose) {
        std::cout << ANSI_COLOR_BLUE_BOLD << "Initializing caloTreeGen -- RUNNING Init" << ANSI_COLOR_RESET << std::endl;
    }
    
    //====================================//
    //          DATA MODE
    //====================================//
    if (wantData) {
        if (verbose) {
            std::cout << "[INFO] Running in DATA mode." << std::endl;
        }
        out = new TFile(Outfile.c_str(),"RECREATE");
        
        if (verbose) {
            std::cout << ANSI_COLOR_BLUE_BOLD << "Output file created: " << Outfile << ANSI_COLOR_RESET << std::endl;
        }
        
        // Load meson mass windows from the CSV file if it exists
        
        /*
         NEED TO UPDATE to ana540_p009 Data still
         */
        const std::string csvFilePath = "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/InvMassCsvFiles/InvariantMassInformation_MBD_NandS_geq_1_Photon_3_GeV_plus_MBD_NS_geq_1_Photon_4_GeV_plus_MBD_NS_geq_1_Photon_5_GeV_plus_MBD_NS_geq_1_afterTriggerFirmwareUpdate.csv";

        if (std::filesystem::exists(csvFilePath)) {
            if (!loadMesonMassWindows(csvFilePath)) {
                std::cerr << "Warning: Failed to load meson mass windows from CSV. Continuing without it." << std::endl;
            }
        } else {
            std::cout << "No CSV file found at " << csvFilePath << ". Skipping meson mass windows loading." << std::endl;
        }
        trigAna = new TriggerAnalyzer();
        
        createHistos_Data();
        
        //so that the histos actually get written out
        Fun4AllServer *se = Fun4AllServer::instance();
        if (verbose) {
            se -> Print("NODETREE");
        }
        std::cout << "caloTreeGen::Init(PHCompositeNode *topNode) Initializing" << std::endl;
        
        // Notice we do NOT 'return' here if wantSim is also true;
        // instead, we keep going so we can also do the sim part below.
    }
    
    // -----------------------------------------------------
    // 2) SIMULATION MODE
    // -----------------------------------------------------
    if (wantSim)
    {
        if (verbose)
        {
            std::cout << "[INFO] Running in SIMULATION mode." << std::endl;
        }

        // Create a separate output file for simulation histograms
        outSim = new TFile(SimOutfile.c_str(), "RECREATE");
        if (!outSim || outSim->IsZombie())
        {
            std::cerr << "[ERROR] Could not open simulation output file: "
                      << SimOutfile << std::endl;
            return Fun4AllReturnCodes::ABORTEVENT;
        }
        // Build simulation-specific histograms
        createHistos_ForSimulation();
        
        // **Open** the simulation input file "simInputFileName"
        simInFile = TFile::Open(simInputFileName.c_str(), "READ");
        if (!simInFile || simInFile->IsZombie())
        {
          std::cerr << "[ERROR] Could not open simulation input file: "
                    << simInputFileName << std::endl;
          return Fun4AllReturnCodes::ABORTEVENT;
        }
        if (verbose)
          std::cout << "[SIM] Reading from sim input file: " << simInputFileName << std::endl;

        // Retrieve TTree "slimtree"
        slimTree = dynamic_cast<TTree*>(simInFile->Get("slimtree"));
        if (!slimTree)
        {
          std::cerr << "[ERROR] TTree 'slimtree' not found in " << simInputFileName << std::endl;
          return Fun4AllReturnCodes::ABORTEVENT;
        }
        // Set branch addresses (cluster-level)
        slimTree->SetBranchAddress("ncluster_CLUSTERINFO_CEMC", &ncluster_CEMC_SIM);
        slimTree->SetBranchAddress("cluster_pid_CEMC",          cluster_pid_CEMC_SIM);
        slimTree->SetBranchAddress("cluster_iso_04_CEMC",       cluster_iso_04_CEMC_SIM);
        slimTree->SetBranchAddress("cluster_Et_CEMC",           cluster_Et_CEMC_SIM);
        slimTree->SetBranchAddress("cluster_Eta_CEMC",          cluster_Eta_CEMC_SIM);
        slimTree->SetBranchAddress("cluster_Phi_CEMC",          cluster_Phi_CEMC_SIM);

        // Set branch addresses (truth-level)
        slimTree->SetBranchAddress("nparticles",         &nparticles_SIM);
        slimTree->SetBranchAddress("particle_pid",       particle_pid_SIM);
        slimTree->SetBranchAddress("particle_photonclass", particle_photonclass_SIM);
        slimTree->SetBranchAddress("particle_Pt",        particle_Pt_SIM);
        slimTree->SetBranchAddress("particle_Eta",       particle_Eta_SIM);
        slimTree->SetBranchAddress("particle_Phi",       particle_Phi_SIM);

        std::cout << "[SIM] Done setting TTree branch addresses.\n";
      }

      // If neither data nor sim was requested => error
      if (!wantData && !wantSim)
      {
        std::cerr << "[ERROR] Neither wantData nor wantSim is set.\n";
        return Fun4AllReturnCodes::ABORTEVENT;
      }

    // If we get here, we have at least one of them set to true
    return Fun4AllReturnCodes::EVENT_OK;
}




//____________________________________________________________________________..
void caloTreeGen::createHistos_Data() {
    std::cout << "[DEBUG] Entering caloTreeGen::createHistos()..." << std::endl;

    for (const auto& kv : triggerNameMap)
    {
        // kv.first  -> DB name  (unused here except for logging)
        // kv.second -> shortName
        const std::string& triggerName = kv.second;

        if (verbose)
        {
            std::cout << ANSI_COLOR_BLUE_BOLD
                      << "Creating histograms for trigger: " << triggerName
                      << ANSI_COLOR_RESET << std::endl;
        }

        // Create a directory for the current trigger
        TDirectory* triggerDir = out->mkdir(triggerName.c_str());
        if (!triggerDir) {
            std::cerr << "[ERROR] Failed to create directory for trigger: " << triggerName << std::endl;
            exit(EXIT_FAILURE);
        }
        triggerDir->cd(); // Set the current directory to the trigger directory


        std::map<std::string, TObject*>& qaHistograms = qaHistogramsByTrigger[triggerName];
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
        qaHistograms["h_maxEnergyClus_" + triggerName] = createHistogram("h_maxEnergyClus_" + triggerName, "Max Cluster Energy; Cluster Energy [GeV]", 40, 0, 20);
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

        // E3x3 / E_cluster
        qaHistograms["E3x3_over_ClusterE_NoShowerShapeCuts_" + triggerName] =
            createHistogram("E3x3_over_ClusterE_NoShowerShapeCuts_" + triggerName,
                            "E_{3x3}/E_{cluster} (No Cuts);E_{3x3}/E_{cluster};Counts",
                            100, 0.0, 1.2);

        qaHistograms["E3x3_over_ClusterE_withShowerShapeCuts_" + triggerName] =
            createHistogram("E3x3_over_ClusterE_withShowerShapeCuts_" + triggerName,
                            "E_{3x3}/E_{cluster} (With Cuts);E_{3x3}/E_{cluster};Counts",
                            100, 0.0, 1.2);

        // E3x3 / E3x7
        qaHistograms["E3x3_over_E3x7_NoShowerShapeCuts_" + triggerName] =
            createHistogram("E3x3_over_E3x7_NoShowerShapeCuts_" + triggerName,
                            "E_{3x3}/E_{3x7} (No Cuts);E_{3x3}/E_{3x7};Counts",
                            100, 0.0, 1.2);

        qaHistograms["E3x3_over_E3x7_withShowerShapeCuts_" + triggerName] =
            createHistogram("E3x3_over_E3x7_withShowerShapeCuts_" + triggerName,
                            "E_{3x3}/E_{3x7} (With Cuts);E_{3x3}/E_{3x7};Counts",
                            100, 0.0, 1.2);
        
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
        
        /*
         Shower Shape QA
         */
        
        //E3x7/E7x7
        qaHistograms["E3by7_over_E7by7_NoShowerShapeCuts_" + triggerName] = createHistogram("E3by7_over_E7by7_NoShowerShapeCuts_" + triggerName,
                            "E3x7/E7x7 (No Cuts);E_{3x7}/E_{7x7};Counts",
                             100, 0.0, 1.2);

        qaHistograms["E3by7_over_E7by7_withShowerShapeCuts_" + triggerName] = createHistogram("E3by7_over_E7by7_withShowerShapeCuts_" + triggerName,
                            "E3x7/E7x7 (With ShowerShape Cuts);E_{3x7}/E_{7x7};Counts",
                             100, 0.0, 1.2);
        
        //w7x2
        qaHistograms["w72_NoShowerShapeCuts_" + triggerName] =
            createHistogram("w72_NoShowerShapeCuts_" + triggerName,
                            "w_{72} (No Cuts);w_{72};Counts",
                             100, 0.0, 2.0);  // choose x-range as needed

        qaHistograms["w72_withShowerShapeCuts_" + triggerName] =
            createHistogram("w72_withShowerShapeCuts_" + triggerName,
                            "w_{72} (With ShowerShape Cuts);w_{72};Counts",
                             100, 0.0, 2.0);
        
        // E1x1 / Ecluster
        qaHistograms["E1x1_over_ClusterE_NoShowerShapeCuts_" + triggerName] =
            createHistogram("E1x1_over_ClusterE_NoShowerShapeCuts_" + triggerName,
                            "E_{1x1}/E_{cluster} (No Cuts);E_{1x1}/E;Counts",
                             100, 0.0, 1.2);

        qaHistograms["E1x1_over_ClusterE_withShowerShapeCuts_" + triggerName] =
            createHistogram("E1x1_over_ClusterE_withShowerShapeCuts_" + triggerName,
                            "E_{1x1}/E_{cluster} (With Cuts);E_{1x1}/E;Counts",
                             100, 0.0, 1.2);

        // E1x1 / E3x3
        qaHistograms["E1x1_over_E3x3_NoShowerShapeCuts_" + triggerName] =
            createHistogram("E1x1_over_E3x3_NoShowerShapeCuts_" + triggerName,
                            "E_{1x1}/E_{3x3} (No Cuts);E_{1x1}/E_{3x3};Counts",
                             100, 0.0, 1.2);

        qaHistograms["E1x1_over_E3x3_withShowerShapeCuts_" + triggerName] =
            createHistogram("E1x1_over_E3x3_withShowerShapeCuts_" + triggerName,
                            "E_{1x1}/E_{3x3} (With Cuts);E_{1x1}/E_{3x3};Counts",
                             100, 0.0, 1.2);

        
        // E1x7 / E7x7
        qaHistograms["E1by7_over_E7by7_NoShowerShapeCuts_" + triggerName] =
            createHistogram("E1by7_over_E7by7_NoShowerShapeCuts_" + triggerName,
                            "E_{1x7}/E_{7x7} (No Cuts);E_{1x7}/E_{7x7};Counts",
                            100, 0.0, 1.2);

        qaHistograms["E1by7_over_E7by7_withShowerShapeCuts_" + triggerName] =
            createHistogram("E1by7_over_E7by7_withShowerShapeCuts_" + triggerName,
                            "E_{1x7}/E_{7x7} (With Cuts);E_{1x7}/E_{7x7};Counts",
                            100, 0.0, 1.2);
        
        // E3x2 / E3x5
        qaHistograms["E3x2_over_E3x5_NoShowerShapeCuts_" + triggerName] =
            createHistogram("E3x2_over_E3x5_NoShowerShapeCuts_" + triggerName,
                            "E_{3x2}/E_{3x5} (No Cuts);E_{3x2}/E_{3x5};Counts",
                            100, 0.0, 1.2);

        qaHistograms["E3x2_over_E3x5_withShowerShapeCuts_" + triggerName] =
            createHistogram("E3x2_over_E3x5_withShowerShapeCuts_" + triggerName,
                            "E_{3x2}/E_{3x5} (With Cuts);E_{3x2}/E_{3x5};Counts",
                            100, 0.0, 1.2);
        
        // w_eta histograms
        qaHistograms["weta_NoShowerShapeCuts_" + triggerName] =
            createHistogram("weta_NoShowerShapeCuts_" + triggerName,
                            "w_{#eta} (No Cuts);w_{#eta};Counts",
                            100, 0.0, 0.1);  // Example bin range: 0.0 to 0.1

        qaHistograms["weta_withShowerShapeCuts_" + triggerName] =
            createHistogram("weta_withShowerShapeCuts_" + triggerName,
                            "w_{#eta} (With Cuts);w_{#eta};Counts",
                            100, 0.0, 0.1);

        // w_phi histograms
        qaHistograms["wphi_NoShowerShapeCuts_" + triggerName] =
            createHistogram("wphi_NoShowerShapeCuts_" + triggerName,
                            "w_{#phi} (No Cuts);w_{#phi};Counts",
                            100, 0.0, 0.1);

        qaHistograms["wphi_withShowerShapeCuts_" + triggerName] =
            createHistogram("wphi_withShowerShapeCuts_" + triggerName,
                            "w_{#phi} (With Cuts);w_{#phi};Counts",
                            100, 0.0, 0.1);

        
        std::map<std::pair<float, float>, std::map<std::string, TObject*>>& qaIsolationHistograms = qaIsolationHistogramsByTriggerAndPt[triggerName];

        for (const auto& pT_bin : pT_bins)
        {
            float pT_min = pT_bin.first;
            float pT_max = pT_bin.second;
            std::pair<float, float> pT_range = {pT_min, pT_max};

            // ---------- 1) "Unsubtracted" iso hist for this pT bin
            {
                std::string hist2DName_unsub = "h2_cluster_iso_Et_unsubtracted_pT_" +
                    formatFloatForFilename(pT_min) + "to" +
                    formatFloatForFilename(pT_max) + "_" + triggerName;

                std::string hist1DName_unsub = "h1_isoEt_unsubtracted_pT_" +
                    formatFloatForFilename(pT_min) + "to" +
                    formatFloatForFilename(pT_max) + "_" + triggerName;

                qaIsolationHistograms[pT_range][hist2DName_unsub] = create2DHistogram(
                    hist2DName_unsub,
                    "Cluster Isolation Energy vs Cluster Et (Unsub);Cluster Et [GeV];E_{T}^{iso} [GeV]",
                    100, 0, 20, 100, -20, 20);

                qaIsolationHistograms[pT_range][hist1DName_unsub] = createHistogram(
                    hist1DName_unsub,
                    "Isolation Energy Distribution (Unsub);E_{T}^{iso} [GeV];Counts",
                    100, -20, 20);
            }

            // ---------- 2) "Subtracted" iso hist for this pT bin
            {
                std::string hist2DName_sub = "h2_cluster_iso_Et_subtracted_pT_" +
                    formatFloatForFilename(pT_min) + "to" +
                    formatFloatForFilename(pT_max) + "_" + triggerName;

                std::string hist1DName_sub = "h1_isoEt_subtracted_pT_" +
                    formatFloatForFilename(pT_min) + "to" +
                    formatFloatForFilename(pT_max) + "_" + triggerName;

                qaIsolationHistograms[pT_range][hist2DName_sub] = create2DHistogram(
                    hist2DName_sub,
                    "Cluster Isolation Energy vs Cluster Et (Sub);Cluster Et [GeV];E_{T}^{iso} [GeV]",
                    100, 0, 20, 100, -20, 20);

                qaIsolationHistograms[pT_range][hist1DName_sub] = createHistogram(
                    hist1DName_sub,
                    "Isolation Energy Distribution (Sub);E_{T}^{iso} [GeV];Counts",
                    100, -20, 20);
            }

            // ---------- 3) Create the nClustersInEvent_pT_... histogram here
            {
                std::string nClusHistName = "nClustersInEvent_pT_" +
                    formatFloatForFilename(pT_min) + "to" +
                    formatFloatForFilename(pT_max) + "_" + triggerName;

                TH1F* nClusHist = new TH1F(
                    nClusHistName.c_str(),
                    (nClusHistName + ";N_{clusters};Events").c_str(),
                    200, 0, 200
                );
                nClusHist->SetDirectory(out);

                // Store it in the same map
                qaIsolationHistograms[pT_range][nClusHistName] = nClusHist;
            }
        }


        qaHistogramsByTrigger[triggerName] = qaHistograms;
        qaIsolationHistogramsByTriggerAndPt[triggerName] = qaIsolationHistograms;


        for (float maxAsym : asymmetry_values) {
            for (float maxChi2 : clus_chi_values) {
                for (float minClusE : clus_Energy_values) {
                    
                    std::string invMassHistName_noBinsOfPt_name = "invMass_noPtBins_E" + formatFloatForFilename(minClusE) +
                                                  "_Chi" + formatFloatForFilename(maxChi2) +
                                                  "_Asym" + formatFloatForFilename(maxAsym) +
                                                  "_" + triggerName;
                    TH1F* invMassHistName_noBinsOfPt = new TH1F(invMassHistName_noBinsOfPt_name.c_str(), invMassHistName_noBinsOfPt_name.c_str(), 80, 0, 1.0);
                    invMassHistName_noBinsOfPt->SetTitle(";M_{#gamma#gamma};");
                    massAndIsolationHistogramsNoPtBins[triggerName][invMassHistName_noBinsOfPt_name] = invMassHistName_noBinsOfPt;


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
                        
                        massAndIsolationHistograms[triggerName][std::make_tuple(maxAsym, maxChi2, minClusE)][pT_bin][invMassHistName] = invMassHist;
                        
                        // Define pion mass vs isolation energy histogram
                        std::string pionHistName = "pionMass_vs_isoEt_E" + formatFloatForFilename(minClusE) +
                                                   "_Chi" + formatFloatForFilename(maxChi2) +
                                                   "_Asym" + formatFloatForFilename(maxAsym) +
                                                   "_pT_" + formatFloatForFilename(pT_bin.first) + "to" + formatFloatForFilename(pT_bin.second) +
                                                   "_" + triggerName;

                        TH2F* pionMassVsIsoHist = new TH2F(pionHistName.c_str(),
                                                           "Pion Mass vs Isolation Energy;M_{#pi^{0}} [GeV];E_{T}^{iso} [GeV]",
                                                           80, 0, 1, 100, -20, 20);
                        massAndIsolationHistograms[triggerName][std::make_tuple(maxAsym, maxChi2, minClusE)][pT_bin][pionHistName] = pionMassVsIsoHist;

                        // Define eta mass vs isolation energy histogram
                        std::string etaHistName = "etaMass_vs_isoEt_E" + formatFloatForFilename(minClusE) +
                                                  "_Chi" + formatFloatForFilename(maxChi2) +
                                                  "_Asym" + formatFloatForFilename(maxAsym) +
                                                  "_pT_" + formatFloatForFilename(pT_bin.first) + "to" + formatFloatForFilename(pT_bin.second) +
                                                  "_" + triggerName;

                        TH2F* etaMassVsIsoHist = new TH2F(etaHistName.c_str(),
                                                          "Eta Mass vs Isolation Energy;M_{#eta} [GeV];E_{T}^{iso} [GeV]",
                                                          80, 0, 1, 100, -20, 20);
                        massAndIsolationHistograms[triggerName][std::make_tuple(maxAsym, maxChi2, minClusE)][pT_bin][etaHistName] = etaMassVsIsoHist;
                        

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

                            TH2F* hist2D = new TH2F(hist2DName.c_str(),
                                                  "Cluster Isolation Energy vs Cluster Et;Cluster Et [GeV];E_{T}^{iso} [GeV]",
                                                  100, 0, 20, 100, -20, 20);
                            massAndIsolationHistograms[triggerName][std::make_tuple(maxAsym, maxChi2, minClusE)][pT_bin][hist2DName] = hist2D;

                            // Create and store the 1D histogram (Isolation energy)
                            TH1F* hist1D = new TH1F(hist1DName.c_str(),
                                                  "Isolation Energy Distribution;E_{T}^{iso} [GeV];Counts",
                                                  100, -20, 20);
                            massAndIsolationHistograms[triggerName][std::make_tuple(maxAsym, maxChi2, minClusE)][pT_bin][hist1DName] = hist1D;

                            
                            for (const auto& isoRange : isoEtRanges)
                            {
                                float isoMin = isoRange.first;
                                float isoMax = isoRange.second;

                                // Build a reusable base for the histogram name
                                std::string baseName = "isolatedPhotonCount_E" + formatFloatForFilename(minClusE) +
                                                       "_Chi" + formatFloatForFilename(maxChi2) +
                                                       "_Asym" + formatFloatForFilename(maxAsym) +
                                                       massWindowLabel +
                                                       "_isoEt_" + formatFloatForFilename(isoMin) + "to" + formatFloatForFilename(isoMax) +
                                                       "_pT_" + formatFloatForFilename(pT_bin.first) + "to" + formatFloatForFilename(pT_bin.second) +
                                                       "_"; // trailing underscore

                                // 1) Create the "withoutShowerShapeCuts" version
                                std::string noCutHistName = baseName + "withoutShowerShapeCuts_" + triggerName;
                                TH1F* hist_noCuts = new TH1F(
                                    noCutHistName.c_str(),
                                    "Isolated Photon Count (No ShowerShapeCuts); pT [GeV]; Count",
                                    /* nbins */ 100, pT_bin.first, pT_bin.second
                                );

                                // 2) Create the "withShowerShapeCuts" version
                                std::string withCutHistName = baseName + "withShowerShapeCuts_" + triggerName;
                                TH1F* hist_withCuts = new TH1F(
                                    withCutHistName.c_str(),
                                    "Isolated Photon Count (With ShowerShapeCuts); pT [GeV]; Count",
                                    /* nbins */ 100, pT_bin.first, pT_bin.second
                                );

                                // Store them in your existing map data structure
                                massAndIsolationHistograms[triggerName][std::make_tuple(maxAsym, maxChi2, minClusE)]
                                                           [pT_bin][noCutHistName] = hist_noCuts;
                                massAndIsolationHistograms[triggerName][std::make_tuple(maxAsym, maxChi2, minClusE)]
                                                           [pT_bin][withCutHistName] = hist_withCuts;
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
                            massAndIsolationHistograms[triggerName][std::make_tuple(maxAsym, maxChi2, minClusE)][pT_bin][allPhotonHistName] = allPhotonHist;
                            massAndIsolationHistograms[triggerName][std::make_tuple(maxAsym, maxChi2, minClusE)][pT_bin][ptPhotonHistName] = ptPhotonHist;
                            
                        }
                    }
                }
            }
        }
    }
}


void caloTreeGen::createHistos_ForSimulation()
{
  // [1] If verbose, print a message that we are entering this method
  if (verbose)
    std::cout << "[DEBUG] Entering caloTreeGen::createHistos_ForSimulation()...\n";

  // [2] Basic error check: ensure we have a valid TFile pointer for SIM
  if (!outSim)
  {
    std::cerr << "[ERROR] 'outSim' is null! Did you forget to open the sim output file?\n";
    return;
  }

  // [3] Create or retrieve the sub-directory "SimulationQA" in the outSim TFile
  //     so we can keep simulation histograms separate from data histos
  TDirectory* simDir = outSim->mkdir("SimulationQA");
  if (!simDir)
  {
    std::cerr << "[ERROR] Failed to create 'SimulationQA' directory in outSim.\n";
    exit(EXIT_FAILURE);
  }
  simDir->cd(); // Move into that directory, so new histograms go there

  // [4] Define a small lambda that helps create 1D histograms and checks for duplicates
  auto createHist1D = [&](const std::string &hname,
                          const std::string &htitle,
                          int nbins, double xmin, double xmax)
  {
    // Check for duplicates
    if (outSim->Get(hname.c_str()))
    {
      std::cerr << "[ERROR] Duplicate histogram name in SIM: " << hname << std::endl;
      exit(EXIT_FAILURE);
    }
    // Actually create the histogram
    TH1F* h = new TH1F(hname.c_str(), htitle.c_str(), nbins, xmin, xmax);
    // Assign it to the simDir so it is owned by the TFile
    h->SetDirectory(simDir);

    // If we want extra debug info
    if (verbose)
    {
      std::cout << "[DEBUG] Creating 1D hist: " << hname
                << ", title='" << htitle << "', nbins=" << nbins
                << ", range=[" << xmin << "," << xmax << "]\n";
    }

    return h;
  };

  // [5] We build classification & purity histograms for each pT bin
  //     1) classification: "hClusterTruthClass_pT_<min>to<max>"
  //     2) purity:         "hPromptPurity_pT_<min>to<max>"
  for (auto& bin : pT_bins)
  {
    float pTmin = bin.first;
    float pTmax = bin.second;

    // Construct a safe name for classification histogram
    std::string hname_class = "hClusterTruthClass_pT_" +
       formatFloatForFilename(pTmin) + "to" + formatFloatForFilename(pTmax);

    // Title might mention that we store 1=prompt,2=frag,3=decay,4=other
    std::string htitle_class = hname_class + ";Truth Class Code;Counts";
    
    // Create the classification histogram
    TH1F* hClass = createHist1D(hname_class, htitle_class, /*nbins=*/4, /*xmin=*/0.5, /*xmax=*/4.5);
    
    // Save the pointer in our map so we can fill it later
    hClusterTruthClass_pTbin[bin] = hClass;

    // Now define the purity histogram name
    std::string hname_purity = "hPromptPurity_pT_" +
       formatFloatForFilename(pTmin) + "to" + formatFloatForFilename(pTmax);

    // Title for that histogram
    std::string htitle_purity = hname_purity + ";Purity;Counts";
    // We'll fill it with a ratio from 0..1 => 50 bins from 0..1 is typical
    TH1F* hPurity = createHist1D(hname_purity, htitle_purity, /*nbins=*/50, /*xmin=*/0., /*xmax=*/1.);
    hPromptPurity_pTbin[bin] = hPurity;
  }

  // [6] Create “global” or “unbinned” histograms for isolation from pi0/eta, etc.
  hIsoFromPi0Eta = createHist1D(
    "hIsoFromPi0Eta",
    "Isolated Clusters from #pi^{0}/#eta (iso <= 6 GeV);Cluster E_{T} [GeV];Counts",
    50, 0.0, 50.0
  );

  hIsoNotPi0Eta  = createHist1D(
    "hIsoNotPi0Eta",
    "Isolated Clusters not from #pi^{0}/#eta (iso <= 6 GeV);Cluster E_{T} [GeV];Counts",
    50, 0.0, 50.0
  );

  // Perhaps a global histogram for prompt photon pT
  hPhotonPtPrompt = createHist1D(
    "hPhotonPtPrompt",
    "Prompt Photon pT; pT [GeV];Counts",
    60, 0, 30
  );

  // [7] Return to top-level directory so we don’t nest future objects incorrectly
  outSim->cd();

  // [8] If verbose, print a final message that we created the histograms
  if (verbose)
  {
    std::cout << "[INFO] createHistos_ForSimulation: pT-binned and global histograms created.\n"
              << "       Histograms stored in 'SimulationQA' directory of outSim.\n";
  }
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


caloTreeGen::EnergyMaps caloTreeGen::processEnergyMaps(const std::vector<float>* m_emcTowE, const std::vector<float>* m_emciEta, const std::vector<float>* m_emciPhi, const std::vector<float>* m_ohcTowE, const std::vector<float>* m_ohciTowEta, const std::vector<float>* m_ohciTowPhi, const std::vector<float>* m_ihcTowE, const std::vector<float>* m_ihciTowEta, const std::vector<float>* m_ihciTowPhi, std::vector<short>* m_emcal_good, std::vector<short>* m_ohc_good, std::vector<short>* m_ihc_good, std::vector<std::string>& activeTriggerNames) {
    
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

std::vector<int> caloTreeGen::find_closest_hcal_tower(
    float eta, float phi,
    RawTowerGeomContainer* geom,
    TowerInfoContainer* towerContainer,
    float vertex_z,
    bool isihcal
)
{
  if (verbose)
  {
    std::cout << "[DEBUG] find_closest_hcal_tower(...) called with: "
              << "eta=" << eta << ", phi=" << phi
              << ", vertex_z=" << vertex_z
              << ", isihcal=" << isihcal << std::endl
              << "[DEBUG]   towerContainer->size()=" << towerContainer->size()
              << ", geom->size()=" << geom->size() << std::endl;
  }
    
  // Initialize "best match" results
  int    matchedieta = -1;
  int    matchediphi = -1;
  double matchedeta  = -999.;
  double matchedphi  = -999.;
    
  // We'll search for the tower with the smallest deltaR in (eta,phi)
  float  minR = 999.f;
  unsigned int ntowers = towerContainer->size();

  // Loop over all towers in TowerInfoContainer
  for (unsigned int channel = 0; channel < ntowers; channel++)
  {
    TowerInfo* tower = towerContainer->get_tower_at_channel(channel);
    if (!tower)
    {
      if (verbose)
      {
        std::cout << "[DEBUG] tower pointer for channel=" << channel
                  << " is NULL, skipping.\n";
      }
      continue;
    }

    // Decode the tower key => ieta, iphi
    unsigned int towerkey = towerContainer->encode_key(channel);
    int ieta = towerContainer->getTowerEtaBin(towerkey);
    int iphi = towerContainer->getTowerPhiBin(towerkey);

    //  geometry has 24 ieta bins [0..23] & 64 iphi bins [0..63],
    // skip anything out of that range:
    if (ieta < 0 || ieta >= 24)
    {
      if (verbose)
      {
        std::cout << "[DEBUG] ieta=" << ieta
                  << " is out of range => skipping.\n";
      }
      continue;
    }
    if (iphi < 0 || iphi >= 64)
    {
      if (verbose)
      {
        std::cout << "[DEBUG] iphi=" << iphi
                  << " is out of range => skipping.\n";
      }
      continue;
    }

    // Build the raw tower geometry key
    RawTowerDefs::CalorimeterId calID = isihcal ? RawTowerDefs::HCALIN
                                                : RawTowerDefs::HCALOUT;
    RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(calID, ieta, iphi);

    // Grab geometry for that tower
    RawTowerGeom* tower_geom = geom->get_tower_geometry(key);
    if (!tower_geom)
    {
      // If geometry is missing, skip
      if (verbose)
      {
        std::cout << "[DEBUG] tower_geom is NULL for channel=" << channel
                  << " ieta=" << ieta << ", iphi=" << iphi
                  << " => skipping.\n";
      }
      continue;
    }

    // Get this tower’s (eta, phi), compute deltaR to the cluster
    double this_phi = tower_geom->get_phi();
    double this_eta = getTowerEta(tower_geom, 0, 0, vertex_z);
    double dR       = deltaR(eta, this_eta, phi, this_phi);

    // Update "closest tower" if this tower is closer
    if (dR < minR)
    {
      minR       = dR;
      matchedieta= ieta;
      matchediphi= iphi;
      matchedeta = this_eta;
      matchedphi = this_phi;
    }
  } // end for(channel)

  // Now we have the best-match tower => compute sign of (dEta, dPhi)
  float deta = eta - matchedeta;
  float dphi = phi - matchedphi;
  float dphiwrap = 2 * M_PI - std::abs(dphi);
  if (std::abs(dphiwrap) < std::abs(dphi))
  {
    dphi = (dphi > 0) ? -dphiwrap : dphiwrap;
  }
  int dphisign = (dphi > 0) ? 1 : -1;
  int detasign = (deta > 0) ? 1 : -1;

  // Return results in a small integer vector
  std::vector<int> result = {matchedieta, matchediphi, detasign, dphisign};
  return result;
}



caloTreeGen::ShowerShapeVars caloTreeGen::computeShowerShapesForCluster(
    RawCluster*            cluster,
    TowerInfoContainer*    emcTowerContainer,
    RawTowerGeomContainer* geomEM,
    TowerInfoContainer*    ihcalTowerContainer,
    RawTowerGeomContainer* geomIH,
    TowerInfoContainer*    ohcalTowerContainer,
    RawTowerGeomContainer* geomOH,
    float                  vtxz,         // z-vertex
    float                  clusterEta,   // cluster's eta
    float                  clusterPhi    // cluster's phi
)
{
  /*
    Calculates various energy-based “shower shape” variables.
    If something fails, returns a zero-initialized ShowerShapeVars struct.
  */

  ShowerShapeVars svars;  // default-constructed (all zero)

  if (verbose)
  {
    std::cout << "[DEBUG] computeShowerShapesForCluster() entered."
              << " cluster=" << cluster
              << ", emcTowerContainer=" << emcTowerContainer
              << ", geomEM=" << geomEM
              << std::endl;
  }

  // -----------------------------------------------------------------
  // Step 0) Basic checks
  // -----------------------------------------------------------------
  if (!cluster || !emcTowerContainer || !geomEM)
  {
    if (verbose)
    {
      std::cout << "[DEBUG] computeShowerShapes: missing pointer(s)."
                << " cluster=" << cluster
                << ", emcTowerContainer=" << emcTowerContainer
                << ", geomEM=" << geomEM
                << " => returning default svars.\n";
    }
    return svars;
  }

  if (verbose)
  {
    std::cout << "[DEBUG] Step 0 passed. Attempting to retrieve shower shape from cluster...\n";
  }

  // -----------------------------------------------------------------
  // Step 1) Retrieve cluster shower shapes
  // -----------------------------------------------------------------
  std::vector<float> showershape = cluster->get_shower_shapes(0.070);
  if (showershape.size() < 6)
  {
    if (verbose)
    {
      std::cout << "[DEBUG] cluster->get_shower_shapes(0.070) returned only "
                << showershape.size() << " elements < 6 => cannot proceed.\n";
    }
    return svars;
  }

  float avg_eta = showershape[4] + 0.5f;
  float avg_phi = showershape[5] + 0.5f;
  int maxieta   = static_cast<int>(std::floor(avg_eta));
  int maxiphi   = static_cast<int>(std::floor(avg_phi));

  if (verbose)
  {
    std::cout << "[DEBUG] Step 1 done. maxieta=" << maxieta
              << ", maxiphi=" << maxiphi
              << " (clusterEta=" << clusterEta
              << ", clusterPhi=" << clusterPhi << ")\n";
  }

  // Check boundary
  if (maxieta < 3 || maxieta > 92)
  {
    if (verbose)
    {
      std::cout << "[DEBUG] maxieta=" << maxieta
                << " is too close to boundary for a 7x7 => returning.\n";
    }
    return svars;
  }

  // -----------------------------------------------------------------
  // Step 2) Build local 7×7 array => E77
  // -----------------------------------------------------------------
  float E77[7][7];
  for (int i = 0; i < 7; i++)
  {
    for (int j = 0; j < 7; j++)
    {
      E77[i][j] = 0.f;
    }
  }

  for (int di = -3; di <= 3; di++)
  {
    for (int dj = -3; dj <= 3; dj++)
    {
      int ieta = maxieta + di;
      int iphi = maxiphi + dj;
      shift_tower_index(ieta, iphi, 96, 256);

      if (ieta < 0 || ieta > 95) continue;

      unsigned int twKey = TowerInfoDefs::encode_emcal(ieta, iphi);
      TowerInfo* tower   = emcTowerContainer->get_tower_at_key(twKey);
      if (!tower) continue;

      float tE = tower->get_energy();
      if (tE < 0.f) tE = 0.f;

      int arrI = di + 3;
      int arrJ = dj + 3;
      E77[arrI][arrJ] = tE;
    }
  }

  if (verbose)
  {
    std::cout << "[DEBUG] Step 2 done. E77 local 7x7 constructed around ("
              << maxieta << "," << maxiphi << ").\n";
  }

  // -----------------------------------------------------------------
  // Step 3) NxN sums
  // -----------------------------------------------------------------
  auto blockSum = [&](int i1, int i2, int j1, int j2)
  {
    float sum = 0.f;
    for (int ii = i1; ii <= i2; ii++)
    {
      for (int jj = j1; jj <= j2; jj++)
      {
        sum += E77[ii][jj];
      }
    }
    return sum;
  };
  auto rowSum = [&](int row, int j1, int j2)
  {
    float sum = 0.f;
    for (int jj = j1; jj <= j2; jj++)
    {
      sum += E77[row][jj];
    }
    return sum;
  };
  auto colSum = [&](int col, int i1, int i2)
  {
    float sum = 0.f;
    for (int ii = i1; ii <= i2; ii++)
    {
      sum += E77[ii][col];
    }
    return sum;
  };

  float e77  = blockSum(0, 6, 0, 6);
  svars.e7x7 = e77;
  svars.e77  = e77;

  svars.e3x7 = blockSum(0, 6, 2, 4);
  svars.e3x3 = blockSum(2, 4, 2, 4);
  svars.e1x1 = E77[3][3];
  svars.e3x2 = blockSum(2, 4, 2, 3);
  svars.e5x5 = blockSum(1, 5, 1, 5);

  // partial sums
  svars.e11 = svars.e1x1;
  svars.e22 = blockSum(2, 3, 2, 3);
  svars.e13 = rowSum(3, 2, 4);
  svars.e15 = rowSum(3, 1, 5);
  svars.e17 = rowSum(3, 0, 6);

  svars.e31 = colSum(3, 2, 4);
  svars.e51 = colSum(3, 1, 5);
  svars.e71 = colSum(3, 0, 6);

  svars.e33 = blockSum(2, 4, 2, 4);
  svars.e35 = blockSum(2, 4, 1, 5);
  svars.e37 = blockSum(2, 4, 0, 6);
  svars.e53 = blockSum(1, 5, 2, 4);
  svars.e73 = blockSum(0, 6, 2, 4);

  svars.e55 = blockSum(1, 5, 1, 5);
  svars.e57 = blockSum(1, 5, 0, 6);
  svars.e75 = blockSum(0, 6, 1, 5);

  if (verbose)
  {
    std::cout << "[DEBUG] Step 3 done. e7x7=" << svars.e7x7
              << ", e3x3=" << svars.e3x3
              << ", e5x5=" << svars.e5x5
              << "\n";
  }

  // -----------------------------------------------------------------
  // Step 4) Weighted widths (w72, w32, w52)
  // -----------------------------------------------------------------
  int signphi = ((avg_phi - std::floor(avg_phi)) > 0.5f) ? 1 : -1;

  // w72
  {
    float sumE72 = 0.f;
    float sumW72 = 0.f;
    for (int i = 0; i < 7; i++)
    {
      for (int j = 2; j <= 4; j++)
      {
        float cE = E77[i][j];
        sumE72   += cE;
        float di  = (float)(i - 3);
        sumW72   += cE * di * di;
      }
    }
    svars.e72 = sumE72;
    svars.w72 = (sumE72 > 1e-6f) ? std::sqrt(sumW72 / sumE72) : 0.f;
  }

  // w32
  {
    float sumE32 = 0.f;
    float sumW32 = 0.f;
    for (int i = 0; i < 7; i++)
    {
      for (int j = 2; j <= 4; j++)
      {
        int di = std::abs(i - 3);
        int dj = std::abs(j - 3);
        // picks 2 columns in the cross‐section
        if (di <= 1 && (dj == 0 || j == (3 + signphi)))
        {
          float cE = E77[i][j];
          sumE32   += cE;
          float dI  = (float)(i - 3);
          sumW32   += cE * (dI * dI);
        }
      }
    }
    svars.e32 = sumE32;
    svars.w32 = (sumE32 > 1e-6f) ? std::sqrt(sumW32 / sumE32) : 0.f;
  }

  // w52
  {
    float sumE52 = 0.f;
    float sumW52 = 0.f;
    for (int i = 0; i < 7; i++)
    {
      for (int j = 2; j <= 4; j++)
      {
        int di = std::abs(i - 3);
        int dj = std::abs(j - 3);
        if (di <= 2 && (dj == 0 || j == (3 + signphi)))
        {
          float cE = E77[i][j];
          sumE52   += cE;
          float dI  = (float)(i - 3);
          sumW52   += cE * (dI * dI);
        }
      }
    }
    svars.e52 = sumE52;
    svars.w52 = (sumE52 > 1e-6f) ? std::sqrt(sumW52 / sumE52) : 0.f;
  }

  if (verbose)
  {
    std::cout << "[DEBUG] Step 4 done. w72=" << svars.w72
              << ", w32=" << svars.w32
              << ", w52=" << svars.w52
              << "\n";
  }

  // -----------------------------------------------------------------
  // Step 5) detamax, dphimax from cluster's towermap
  // -----------------------------------------------------------------
  {
    int dEtaMax = 0;
    int dPhiMax = 0;
    const auto &twmap = cluster->get_towermap();
    if (verbose)
    {
      std::cout << "[DEBUG] cluster->get_towermap() size="
                << twmap.size() << "\n";
    }

    for (const auto &it : twmap)
    {
      auto rawkey = it.first;
      int ieta2   = RawTowerDefs::decode_index1(rawkey);
      int iphi2   = RawTowerDefs::decode_index2(rawkey);

      int ddEta   = std::abs(ieta2 - maxieta);
      int raw_dPhi = iphi2 - maxiphi;
      if (raw_dPhi > 128)  raw_dPhi -= 256;
      if (raw_dPhi < -128) raw_dPhi += 256;
      int ddPhi   = std::abs(raw_dPhi);

      if (ddEta > dEtaMax) dEtaMax = ddEta;
      if (ddPhi > dPhiMax) dPhiMax = ddPhi;
    }
    svars.detamax = dEtaMax;
    svars.dphimax = dPhiMax;
  }

  if (verbose)
  {
    std::cout << "[DEBUG] Step 5 done. detamax=" << svars.detamax
              << ", dphimax=" << svars.dphimax << "\n";
  }

  // -----------------------------------------------------------------
  // Step 6) iHCal / oHCal sums
  // -----------------------------------------------------------------
  if (ihcalTowerContainer && geomIH && ohcalTowerContainer && geomOH)
  {
    if (verbose)
    {
      std::cout << "[DEBUG] Step 6: Summing iHCal/oHCal behind cluster.\n";
    }

    auto ihcalVec = find_closest_hcal_tower(
      clusterEta, clusterPhi,
      geomIH, ihcalTowerContainer, vtxz, true
    );
    svars.ihcal_ieta = ihcalVec[0];
    svars.ihcal_iphi = ihcalVec[1];
      
    // Add debug prints to confirm the returned tower indices
    if (verbose)
    {
      std::cout << "[DEBUG] iHCal tower indices from find_closest_hcal_tower: "
              << "ieta=" << svars.ihcal_ieta << ", iphi=" << svars.ihcal_iphi
              << std::endl;
    }


    auto ohcalVec = find_closest_hcal_tower(
      clusterEta, clusterPhi,
      geomOH, ohcalTowerContainer, vtxz, false
    );
    svars.ohcal_ieta = ohcalVec[0];
    svars.ohcal_iphi = ohcalVec[1];
      
    if (verbose)
      {
        std::cout << "[DEBUG] oHCal tower indices from find_closest_hcal_tower: "
                  << "ieta=" << svars.ohcal_ieta << ", iphi=" << svars.ohcal_iphi
                  << std::endl;
    }
      

      auto sumHcalBlock = [&](TowerInfoContainer* cTC,
                              RawTowerGeomContainer* cGeom,
                              bool inOrOut,
                              int centerEta, int centerPhi,
                              int dEtaMin, int dEtaMax,
                              int dPhiMin, int dPhiMax)
      {
        if (verbose)
        {
          std::cout << "[DEBUG] sumHcalBlock() called with:\n"
                    << "   centerEta=" << centerEta
                    << ", centerPhi=" << centerPhi
                    << ", dEtaMin=" << dEtaMin
                    << ", dEtaMax=" << dEtaMax
                    << ", dPhiMin=" << dPhiMin
                    << ", dPhiMax=" << dPhiMax
                    << ", inOrOut=" << (inOrOut ? "HCALIN" : "HCALOUT")
                    << std::endl;
        }

        float accum = 0.f;
        for (int di = dEtaMin; di <= dEtaMax; di++)
        {
          for (int dj = dPhiMin; dj <= dPhiMax; dj++)
          {
            int ieta_ = centerEta + di;
            int iphi_ = centerPhi + dj;

            // Wrap phi only; ieta_ remains as-is if negative or beyond
            shift_tower_index(ieta_, iphi_, 24, 64);

            // --- NEW RANGE CHECK: ieta_ must be in [0..23], iphi_ in [0..63] ---
            if (ieta_ < 0 || ieta_ >= 24)
            {
              if (verbose)
              {
                std::cout << "   -> ieta_=" << ieta_
                          << " out of range [0..23], skipping.\n";
              }
              continue;
            }
            if (iphi_ < 0 || iphi_ >= 64)
            {
              if (verbose)
              {
                std::cout << "   -> iphi_=" << iphi_
                          << " out of range [0..63], skipping.\n";
              }
              continue;
            }

            // Now safe to encode
            unsigned int tKey = TowerInfoDefs::encode_hcal(ieta_, iphi_);
            TowerInfo* tw = cTC->get_tower_at_key(tKey);

            if (verbose)
            {
              std::cout << "   -> Summation tower: (ieta_=" << ieta_
                        << ", iphi_=" << iphi_
                        << "), TowerInfoKey=" << tKey
                        << ". Tower pointer=" << tw << std::endl;
            }

            if (!tw)
            {
              if (verbose)
              {
                std::cout << "      [DEBUG] Tower pointer is NULL, skipping.\n";
              }
              continue;
            }
            if (!tw->get_isGood())
            {
              if (verbose)
              {
                std::cout << "      [DEBUG] Tower is not 'good', skipping.\n";
              }
              continue;
            }

            RawTowerDefs::CalorimeterId calID =
                inOrOut ? RawTowerDefs::HCALIN : RawTowerDefs::HCALOUT;
            auto rkey = RawTowerDefs::encode_towerid(calID, ieta_, iphi_);
            auto tg   = cGeom->get_tower_geometry(rkey);
            if (!tg)
            {
              if (verbose)
              {
                std::cout << "      [DEBUG] Could not find tower geometry (tg is NULL). "
                          << "Skipping.\n";
              }
              continue;
            }

            float e_   = tw->get_energy();
            float eta_ = getTowerEta(tg, 0, 0, vtxz);
            float et_  = e_ / std::cosh(eta_);
            accum     += et_;

            if (verbose)
            {
              std::cout << "      [DEBUG] Tower is good. E=" << e_
                        << ", eta=" << eta_ << ", ET=" << et_
                        << ", new accum=" << accum << std::endl;
            }
          }
        }
        return accum;
      };
      
      
    // single tower
    svars.ihcal_et = sumHcalBlock(
      ihcalTowerContainer, geomIH, true,
      svars.ihcal_ieta, svars.ihcal_iphi,
      0, 0, 0, 0
    );
    svars.ohcal_et = sumHcalBlock(
      ohcalTowerContainer, geomOH, false,
      svars.ohcal_ieta, svars.ohcal_iphi,
      0, 0, 0, 0
    );

    // 2×2
    svars.ihcal_et22 = sumHcalBlock(
      ihcalTowerContainer, geomIH, true,
      svars.ihcal_ieta, svars.ihcal_iphi,
      -1, 0, -1, 0
    );
    svars.ohcal_et22 = sumHcalBlock(
      ohcalTowerContainer, geomOH, false,
      svars.ohcal_ieta, svars.ohcal_iphi,
      -1, 0, -1, 0
    );

    // 3×3
    svars.ihcal_et33 = sumHcalBlock(
      ihcalTowerContainer, geomIH, true,
      svars.ihcal_ieta, svars.ihcal_iphi,
      -1, 1, -1, 1
    );
    svars.ohcal_et33 = sumHcalBlock(
      ohcalTowerContainer, geomOH, false,
      svars.ohcal_ieta, svars.ohcal_iphi,
      -1, 1, -1, 1
    );
  }
  if (verbose)
  {
    std::cout << "[DEBUG] Step 6 done. ihcal_et=" << svars.ihcal_et
              << ", ohcal_et=" << svars.ohcal_et
              << ", ihcal_et33=" << svars.ihcal_et33
              << ", ohcal_et33=" << svars.ohcal_et33
              << "\n";
  }

  // -----------------------------------------------------------------
  // Step 7) RMS in real (eta,phi)
  // -----------------------------------------------------------------
  {
    float sumEeta      = 0.f;
    float sumWetaLocal = 0.f;
    float sumEphi      = 0.f;
    float sumWphiLocal = 0.f;

    auto wrapDeltaPhi = [](float dphi)
    {
      while (dphi >  M_PI)  dphi -= 2.f*M_PI;
      while (dphi < -M_PI)  dphi += 2.f*M_PI;
      return dphi;
    };

    for (int di = -3; di <= 3; di++)
    {
      for (int dj = -3; dj <= 3; dj++)
      {
        int ieta = maxieta + di;
        int iphi = maxiphi + dj;
        shift_tower_index(ieta, iphi, 96, 256);

        if (ieta < 0 || ieta >= 96)
          {
            if (verbose)
            {
              std::cout << "[DEBUG] ieta=" << ieta
                        << " out of [0..95] => skipping.\n";
            }
            continue;
          }
          if (iphi < 0 || iphi >= 256)
          {
            if (verbose)
            {
              std::cout << "[DEBUG] iphi=" << iphi
                        << " out of [0..255] => skipping.\n";
            }
            continue;
        }


        float tE = E77[di + 3][dj + 3];
        if (tE < 1e-9f) continue;

        RawTowerDefs::CalorimeterId caloID = RawTowerDefs::CEMC;
        RawTowerDefs::keytype twKey = RawTowerDefs::encode_towerid(caloID, ieta, iphi);
        RawTowerGeom* towerGeom = geomEM->get_tower_geometry(twKey);
        if (!towerGeom) continue;


        float towerEta = getTowerEta(towerGeom, 0, 0, vtxz);
        float towerPhi = towerGeom->get_phi();

        float dEta = towerEta - clusterEta;
        float dPhi = wrapDeltaPhi(towerPhi - clusterPhi);

        sumEeta      += tE;
        sumWetaLocal += tE * (dEta * dEta);

        sumEphi      += tE;
        sumWphiLocal += tE * (dPhi * dPhi);
      }
    }

    if (sumEeta > 1e-9f)
      svars.weta = std::sqrt(sumWetaLocal / sumEeta);
    if (sumEphi > 1e-9f)
      svars.wphi = std::sqrt(sumWphiLocal / sumEphi);
  }
  if (verbose)
  {
    std::cout << "[DEBUG] Step 7 done. weta=" << svars.weta
              << ", wphi=" << svars.wphi << "\n";
  }

  if (verbose)
  {
    std::cout << "[DEBUG] computeShowerShapesForCluster finished successfully.\n";
  }
  return svars;
}




void caloTreeGen::processClusterIsolationHistograms(
    int clusterID,
    float mesonMass,
    float minClusEnergy,
    float maxChi2,
    float maxAsym,
    const std::string& massWindowLabel,
    float pT_min,
    float pT_max,
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
    const std::pair<float, float>& pT_bin)
{
    // 1) Must have a valid iso entry for this cluster
    if (!clusterEtIsoMap_unsubtracted.count(clusterID)) return;

    float cluster_et_fromMap = clusterEtIsoMap_unsubtracted.at(clusterID).first;
    float isoEt_FromMap      = clusterEtIsoMap_unsubtracted.at(clusterID).second;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // NOTE: mesonMass is from the *pair*, but we fill histograms per *cluster*.
    // We'll do "once per cluster ID" approach to avoid double-counting
    // the *same cluster* if it reappears in multiple pairs in the same event.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // ------------------------------------------------------------------
    // (A) 2D histogram: "pionMass_vs_isoEt" for clusters from a pi0-candidate pair
    // ------------------------------------------------------------------
    if (std::fabs(mesonMass - pionMass) <= pionMassWindow)
    {
        std::string pionHistName =
            "pionMass_vs_isoEt_E" + formatFloatForFilename(minClusEnergy) +
            "_Chi"  + formatFloatForFilename(maxChi2) +
            "_Asym" + formatFloatForFilename(maxAsym) +
            "_pT_"  + formatFloatForFilename(pT_min) + "to" + formatFloatForFilename(pT_max) +
            "_"     + triggerName;

        // If this cluster hasn't been used yet for the 'pionMass_vs_isoEt' histogram:
        if (!m_usedClustersThisEvent[pionHistName].count(clusterID))
        {
            TH2F* pionMassVsIsoHist = dynamic_cast<TH2F*>(cutHistMap[pT_bin][pionHistName]);
            if (pionMassVsIsoHist)
            {
                // Fill with (M_{gamma-gamma}, cluster iso)
                pionMassVsIsoHist->Fill(mesonMass, isoEt_FromMap);

                if (verbose)
                {
                    std::cout << "[INFO] Filled " << pionHistName
                              << " for clusterID=" << clusterID
                              << " with (mesonMass=" << mesonMass
                              << ", isoEt=" << isoEt_FromMap << ")\n";
                }
                // Mark cluster as used for this hist, so we don't double-count
                m_usedClustersThisEvent[pionHistName].insert(clusterID);

                // Update counters
                filledHistogramCount++;
                filledHistogram = true;
            }
        }
    }

    // ------------------------------------------------------------------
    // (B) 2D histogram: "etaMass_vs_isoEt" for clusters from an eta-candidate pair
    // ------------------------------------------------------------------
    if (std::fabs(mesonMass - etaMass) <= etaMassWindow)
    {
        std::string etaHistName =
            "etaMass_vs_isoEt_E" + formatFloatForFilename(minClusEnergy) +
            "_Chi"  + formatFloatForFilename(maxChi2) +
            "_Asym" + formatFloatForFilename(maxAsym) +
            "_pT_"  + formatFloatForFilename(pT_min) + "to" + formatFloatForFilename(pT_max) +
            "_"     + triggerName;

        if (!m_usedClustersThisEvent[etaHistName].count(clusterID))
        {
            TH2F* etaMassVsIsoHist = dynamic_cast<TH2F*>(cutHistMap[pT_bin][etaHistName]);
            if (etaMassVsIsoHist)
            {
                // Fill with (M_{gamma-gamma}, cluster iso)
                etaMassVsIsoHist->Fill(mesonMass, isoEt_FromMap);

                if (verbose)
                {
                    std::cout << "[INFO] Filled " << etaHistName
                              << " for clusterID=" << clusterID
                              << " with (mesonMass=" << mesonMass
                              << ", isoEt=" << isoEt_FromMap << ")\n";
                }
                m_usedClustersThisEvent[etaHistName].insert(clusterID);

                filledHistogramCount++;
                filledHistogram = true;
            }
        }
    }

    // ------------------------------------------------------------------
    // (C) "Main" iso‐Et histograms => cluster-level hist, once per cluster
    // ------------------------------------------------------------------
    if (verbose)
    {
        std::cout << "[DEBUG] clusterID=" << clusterID
                  << ", clusterEt=" << cluster_et_fromMap
                  << ", isoEt="     << isoEt_FromMap << std::endl;
    }

    // 2D
    std::string hist2DName =
        "h2_cluster_iso_Et_E" + formatFloatForFilename(minClusEnergy) +
        "_Chi"  + formatFloatForFilename(maxChi2) +
        "_Asym" + formatFloatForFilename(maxAsym) +
        massWindowLabel +
        "_pT_"  + formatFloatForFilename(pT_min) + "to" + formatFloatForFilename(pT_max) +
        "_"     + triggerName;

    // 1D
    std::string hist1DName =
        "h1_isoEt_E" + formatFloatForFilename(minClusEnergy) +
        "_Chi"  + formatFloatForFilename(maxChi2) +
        "_Asym" + formatFloatForFilename(maxAsym) +
        massWindowLabel +
        "_pT_"  + formatFloatForFilename(pT_min) + "to" + formatFloatForFilename(pT_max) +
        "_"     + triggerName;

    if (verbose)
    {
        std::cout << "[DEBUG] Attempting to fill:\n  " << hist2DName
                  << "\n  " << hist1DName
                  << " for clusterID=" << clusterID << std::endl;
    }

    // 2D hist fill
    if (!m_usedClustersThisEvent[hist2DName].count(clusterID))
    {
        TH2F* hist2D = dynamic_cast<TH2F*>(cutHistMap[pT_bin][hist2DName]);
        if (hist2D)
        {
            hist2D->Fill(cluster_et_fromMap, isoEt_FromMap);
            m_usedClustersThisEvent[hist2DName].insert(clusterID);

            filledHistogramCount++;
            filledHistogram = true;

            if (verbose)
            {
                std::cout << "[DEBUG] Filled " << hist2DName
                          << " for clusterID=" << clusterID << std::endl;
            }
        }
        else
        {
            std::cerr << "[WARN] " << hist2DName << " is null.\n";
        }
    }

    // 1D hist fill
    if (!m_usedClustersThisEvent[hist1DName].count(clusterID))
    {
        TH1F* hist1D = dynamic_cast<TH1F*>(cutHistMap[pT_bin][hist1DName]);
        if (hist1D)
        {
            hist1D->Fill(isoEt_FromMap);
            m_usedClustersThisEvent[hist1DName].insert(clusterID);

            filledHistogramCount++;
            filledHistogram = true;

            if (verbose)
            {
                std::cout << "[DEBUG] Filled " << hist1DName
                          << " for clusterID=" << clusterID << std::endl;
            }
        }
        else
        {
            std::cerr << "[WARN] " << hist1DName << " is null.\n";
        }
    }
}



void caloTreeGen::processIsolationRanges(
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
    const std::map<int, std::pair<float, float>>& clusterEtIsoMap_unsubtracted,
    std::map<std::pair<float, float>, std::map<std::string, TObject*>>& cutHistMap,
    bool& filledHistogram,
    bool verbose,
    const std::pair<float, float>& pT_bin,
    const std::unordered_map<int,bool> &clusterPassedShowerCuts){

    for (const auto& isoRange : isoEtRanges)
    {
        float isoMin = isoRange.first;
        float isoMax = isoRange.second;

        if (verbose)
        {
            std::cout << "Processing isolation Et range: " << isoMin
                      << " to " << isoMax << std::endl;
        }

        // Check if cluster #1 is isolated
        bool clus1_isolated = false;
        if (clusterEtIsoMap_unsubtracted.count(clusterIDs[clus1]))
        {
            float isoVal_1 = clusterEtIsoMap_unsubtracted.at(clusterIDs[clus1]).second;
            clus1_isolated = (isoVal_1 >= isoMin && isoVal_1 < isoMax);
        }

        // Check if cluster #2 is isolated
        bool clus2_isolated = false;
        if (clusterEtIsoMap_unsubtracted.count(clusterIDs[clus2]))
        {
            float isoVal_2 = clusterEtIsoMap_unsubtracted.at(clusterIDs[clus2]).second;
            clus2_isolated = (isoVal_2 >= isoMin && isoVal_2 < isoMax);
        }

        if (verbose)
        {
            std::cout << "clus1_isolated: " << clus1_isolated
                      << ", clus2_isolated: " << clus2_isolated << std::endl;
        }

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // 1) Fill the "WITHOUT ShowerShapeCuts" histogram
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (clus1_isolated || clus2_isolated)
        {
            if (verbose)
            {
                std::cout << "At least one cluster is isolated in this isoEt range." << std::endl;
            }

            std::string histName_noCuts =
                "isolatedPhotonCount_E" + formatFloatForFilename(minClusEnergy) +
                "_Chi" + formatFloatForFilename(maxChi2) +
                "_Asym" + formatFloatForFilename(maxAsym) +
                massWindowLabel +
                "_isoEt_" + formatFloatForFilename(isoMin) + "to" + formatFloatForFilename(isoMax) +
                "_pT_" + formatFloatForFilename(pT_min) + "to" + formatFloatForFilename(pT_max) +
                "_withoutShowerShapeCuts_" + triggerName;

            TH1F* hist_noCuts = dynamic_cast<TH1F*>(cutHistMap[pT_bin][histName_noCuts]);
            if (hist_noCuts)
            {
                // cluster #1 fill
                if (clus1_isolated)
                {
                    // Only fill if we haven't used clusterIDs[clus1] for this hist
                    if (!m_usedClustersThisEvent[histName_noCuts].count(clusterIDs[clus1]))
                    {
                        hist_noCuts->Fill(1);
                        m_usedClustersThisEvent[histName_noCuts].insert(clusterIDs[clus1]);
                        filledHistogram = true;
                        if (verbose)
                        {
                            std::cout << "[DEBUG] Filled "
                                      << histName_noCuts
                                      << " for clusterID=" << clusterIDs[clus1]
                                      << std::endl;
                        }
                    }
                }
                // cluster #2 fill
                if (clus2_isolated)
                {
                    if (!m_usedClustersThisEvent[histName_noCuts].count(clusterIDs[clus2]))
                    {
                        hist_noCuts->Fill(1);
                        m_usedClustersThisEvent[histName_noCuts].insert(clusterIDs[clus2]);
                        filledHistogram = true;
                        if (verbose)
                        {
                            std::cout << "[DEBUG] Filled "
                                      << histName_noCuts
                                      << " for clusterID=" << clusterIDs[clus2]
                                      << std::endl;
                        }
                    }
                }
            }

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // 2) Fill "WITH ShowerShapeCuts" histogram
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            bool clus1_passCuts = false;
            bool clus2_passCuts = false;
            if (clusterPassedShowerCuts.count(clusterIDs[clus1]))
                clus1_passCuts = clusterPassedShowerCuts.at(clusterIDs[clus1]);
            if (clusterPassedShowerCuts.count(clusterIDs[clus2]))
                clus2_passCuts = clusterPassedShowerCuts.at(clusterIDs[clus2]);

            if ((clus1_isolated && clus1_passCuts) || (clus2_isolated && clus2_passCuts))
            {
                std::string histName_withCuts =
                    "isolatedPhotonCount_E" + formatFloatForFilename(minClusEnergy) +
                    "_Chi" + formatFloatForFilename(maxChi2) +
                    "_Asym" + formatFloatForFilename(maxAsym) +
                    massWindowLabel +
                    "_isoEt_" + formatFloatForFilename(isoMin) + "to" + formatFloatForFilename(isoMax) +
                    "_pT_" + formatFloatForFilename(pT_min) + "to" + formatFloatForFilename(pT_max) +
                    "_withShowerShapeCuts_" + triggerName;

                TH1F* hist_withCuts = dynamic_cast<TH1F*>(cutHistMap[pT_bin][histName_withCuts]);
                if (hist_withCuts)
                {
                    // cluster #1
                    if (clus1_isolated && clus1_passCuts)
                    {
                        if (!m_usedClustersThisEvent[histName_withCuts].count(clusterIDs[clus1]))
                        {
                            hist_withCuts->Fill(1);
                            m_usedClustersThisEvent[histName_withCuts].insert(clusterIDs[clus1]);
                            filledHistogram = true;
                            if (verbose)
                            {
                                std::cout << "[DEBUG] Filled "
                                          << histName_withCuts
                                          << " (with shape cuts) for clusterID="
                                          << clusterIDs[clus1] << std::endl;
                            }
                        }
                    }
                    // cluster #2
                    if (clus2_isolated && clus2_passCuts)
                    {
                        if (!m_usedClustersThisEvent[histName_withCuts].count(clusterIDs[clus2]))
                        {
                            hist_withCuts->Fill(1);
                            m_usedClustersThisEvent[histName_withCuts].insert(clusterIDs[clus2]);
                            filledHistogram = true;
                            if (verbose)
                            {
                                std::cout << "[DEBUG] Filled "
                                          << histName_withCuts
                                          << " (with shape cuts) for clusterID="
                                          << clusterIDs[clus2] << std::endl;
                            }
                        }
                    }
                }
            }
        } // end if (clus1_isolated || clus2_isolated)
    } // end isoRange loop
}


void caloTreeGen::fillHistogramsForTriggers(
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
    const std::map<int, std::pair<float, float>>& clusterEtIsoMap_unsubtracted,
    const std::vector<std::string>& activeTriggerNames,
    bool& filledHistogram,
    const std::unordered_map<int,bool> &clusterPassedShowerCuts) {
    
    /*
     temporary -- next should do the invariant mass analysis in first passs go back through using proper mass window for each pT bin
     */
    const float defaultPionMass = 0.15;
    const float defaultPionMassWindow = 0.02;
    const float defaultEtaMass = 0.59;
    const float defaultEtaMassWindow = 0.05;
    
    for (const auto& firedShortName : activeTriggerNames) {
        if (verbose) {
            std::cout << "Processing trigger: " << firedShortName << std::endl;
        }
        /*
         Filling Histograms outside of pT binning
         */
        auto& noPtBinHistMap = massAndIsolationHistogramsNoPtBins[firedShortName];
        
        std::string invMassHistName_noBinsOfPt_name = "invMass_noPtBins_E" + formatFloatForFilename(minClusEnergy) +
                               "_Chi" + formatFloatForFilename(maxChi2) +
                               "_Asym" + formatFloatForFilename(maxAsym) +
                               "_" + firedShortName;
        
        if (verbose) {
            std::cout << "Attempting to fill histogram: " << invMassHistName_noBinsOfPt_name << std::endl;
        }

        // Attempt to find and fill the histogram
        TH1F* invMassHistName_noBinsOfPt = dynamic_cast<TH1F*>(noPtBinHistMap[invMassHistName_noBinsOfPt_name]);

        if (!invMassHistName_noBinsOfPt) {
            std::cerr << "Error: Histogram " << invMassHistName_noBinsOfPt_name
                      << " not found when trying to fill for Trigger: "
                      << firedShortName << std::endl;
            continue;
        }
        
        invMassHistName_noBinsOfPt->Fill(mesonMass);
        filledHistogramCount++;
        
        if (verbose) {
            std::cout << "Filled histogram " << invMassHistName_noBinsOfPt_name << " with meson mass " << mesonMass << std::endl;
        }
        
        auto& cutHistMap = massAndIsolationHistograms[firedShortName][std::make_tuple(maxAsym, maxChi2, minClusEnergy)];

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
                std::string invMassHistName = "invMass_E" + formatFloatForFilename(minClusEnergy) +
                                              "_Chi" + formatFloatForFilename(maxChi2) +
                                              "_Asym" + formatFloatForFilename(maxAsym) +
                                              "_pT_" + formatFloatForFilename(pT_min) + "to" + formatFloatForFilename(pT_max) +
                                              "_" + firedShortName;

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

                auto massWindowKey = std::make_tuple(firedShortName, minClusEnergy, maxChi2, maxAsym, pT_min, pT_max);
                if (!mesonMassWindowsMap.empty() && mesonMassWindowsMap.count(massWindowKey)) {
                    const MesonMassWindow& massWindow = mesonMassWindowsMap[massWindowKey];

                    if (verbose) {
                        std::cout << "Mass window found in map for trigger " << firedShortName << ":" << std::endl;
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
                         minClusEnergy,
                         maxChi2,
                         maxAsym,
                         massWindowLabel,
                         pT_min,
                         pT_max,
                         firedShortName,
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
                    minClusEnergy,
                    maxChi2,
                    maxAsym,
                    massWindowLabel,
                    pT_min,
                    pT_max,
                    firedShortName,
                    clusterEtIsoMap_unsubtracted,
                    cutHistMap,
                    filledHistogram,
                    verbose,
                    pT_bin,
                    clusterPassedShowerCuts
                );
                // Fill all photons histogram
                std::string allPhotonHistName = "allPhotonCount_E" + formatFloatForFilename(minClusEnergy) +
                                                "_Chi" + formatFloatForFilename(maxChi2) +
                                                "_Asym" + formatFloatForFilename(maxAsym) +
                                                massWindowLabel +
                                                "_pT_" + formatFloatForFilename(pT_min) + "to" + formatFloatForFilename(pT_max) +
                                                "_" + firedShortName;

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
                std::string ptPhotonHistName = "ptPhoton_E" + formatFloatForFilename(minClusEnergy) +
                                               "_Chi" + formatFloatForFilename(maxChi2) +
                                               "_Asym" + formatFloatForFilename(maxAsym) +
                                               massWindowLabel +
                                               "_pT_" + formatFloatForFilename(pT_min) + "to" + formatFloatForFilename(pT_max) +
                                               "_" + firedShortName;

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
    const std::vector<std::string>& activeTriggerNames,
    const std::unordered_map<int,bool> &clusterPassedShowerCuts)
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
    size_t skippedPairsDueToEnergy = 0;
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
                    for (float minClusEnergy : clus_Energy_values) {
                        // Apply selection criteria and count specific cut skips, ensuring each pair is counted once
                        if (!pairSkipped && asym >= maxAsym) {
                            skippedPairsDueToAsymmetry++;
                            pairSkipped = true;
                        } else if (!pairSkipped && (chi1 >= maxChi2 || chi2 >= maxChi2)) {
                            skippedPairsDueToChi2++;
                            pairSkipped = true;
                        } else if (!pairSkipped && (E1 < minClusEnergy || E2 < minClusEnergy)) {
                            skippedPairsDueToEnergy++;
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
                            minClusEnergy,
                            maxChi2,
                            maxAsym,
                            filledHistogramCount,
                            clusterIDs,
                            clusterEtIsoMap_unsubtracted,
                            activeTriggerNames,
                            filledHistogram,
                            clusterPassedShowerCuts
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
        std::cout << "Pairs skipped due to Cluster Energy cuts: " << skippedPairsDueToEnergy
                  << " (" << (100.0 * skippedPairsDueToEnergy / totalPairs) << "%)" << std::endl;
        std::cout << "Pairs with no histograms filled: " << noHistogramsFilledCount
                  << " (" << (100.0 * noHistogramsFilledCount / totalPairs) << "%)" << std::endl;
        std::cout << "Histograms filled: " << filledHistogramCount << std::endl;
        std::cout << "Cluster pairs that passed all cuts with corresponding meson masses and cuts:" << std::endl;

    }
}



bool caloTreeGen::applyShowerShapeCuts(const ShowerShapeVars& s, float clusterEt)
{
    //------------------------------------------------------------------
    // 1) ratio = E_{3x7}/E_{7x7} < 0.9
    //------------------------------------------------------------------
    float ratio_3x7over7x7 = 0.f;
    if (s.e7x7 > 1e-6f) ratio_3x7over7x7 = s.e3x7 / s.e7x7;
    if (ratio_3x7over7x7 >= 0.9f) return false;

    //------------------------------------------------------------------
    // 2) E_T^{HCal, 3x3} / ( E_T^{HCal, 3x3} + E_T^{cluster} ) < 0.1
    //------------------------------------------------------------------
    float hcalEt33 = s.ihcal_et33 + s.ohcal_et33;
    float denom = hcalEt33 + clusterEt;
    if (denom > 1e-6f)
    {
        float frac = hcalEt33 / denom;
        if (frac >= 0.1f) return false;
    }

    //------------------------------------------------------------------
    // 3) E_{1x1}/E_{7x7} < 0.98
    //------------------------------------------------------------------
    float ratio_1x1 = 0.f;
    if (s.e7x7 > 1e-6f) ratio_1x1 = s.e1x1 / s.e7x7;
    if (ratio_1x1 >= 0.98f) return false;

    //------------------------------------------------------------------
    // 4) w72 > 0.75
    //------------------------------------------------------------------
    if (s.w72 <= 0.75f) return false;

    // If we survive all checks => pass
    return true;
}


//____________________________________________________________________________..
int caloTreeGen::process_event(PHCompositeNode *topNode)
{
    // 1) If user wants to process DATA events:
    if (wantData)
    {
      // Everything you do for data
      process_event_Data(topNode);
    }
    
    // 2) If user wants to process SIMULATION events:
    if (wantSim)
    {
      // Everything you do for simulation
      process_event_Sim(topNode);
    }

  return Fun4AllReturnCodes::EVENT_OK;
}


int caloTreeGen::process_event_Data(PHCompositeNode *topNode) {
    event_count++;

    std::cout << "\n========== Processing CALOTREEGEN DATA MODE -- Event " << event_count << " ==========\n";

    // 1) Decode triggers
    if (trigAna)
    {
        if (verbose)
        {
            std::cout << "[DEBUG] About to call trigAna->decodeTriggers()"
                      << " on topNode=" << topNode << std::endl;
        }
        trigAna->decodeTriggers(topNode);
    }
    else
    {
        std::cerr << "[ERROR] No TriggerAnalyzer pointer!\n";
        return Fun4AllReturnCodes::ABORTEVENT;
    }

    // 2) Build a vector of fired trigger *short* names using the DB-names in triggerNameMap.
    std::vector<std::string> activeTriggerNames;
    activeTriggerNames.reserve(triggerNameMap.size());

    for (const auto &kv : triggerNameMap)
    {
        const std::string &dbTriggerName   = kv.first;   // how it's known in the DB
        const std::string &histFriendlyStr = kv.second;  // short name for histograms

        // 3) Check if this DB trigger fired
        if (trigAna->didTriggerFire(dbTriggerName))
        {
            // 4) Save the *short* name for future reference
            activeTriggerNames.push_back(histFriendlyStr);

            if (verbose)
            {
                std::cout << "Trigger fired: \"" << dbTriggerName
                          << "\" => short name \"" << histFriendlyStr << "\"" << std::endl;
            }
        }
    }

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
     switch to non retowered emcal -- TOWERINFO_CALIB_CEMC OR TOWERINFO_CALIB_CEMC_RETOWER
     */
    TowerInfoContainer* emcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC"); // RETOWER MADE THE SHOWER PROFILE CLACULATION NOT WORK
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
                      &m_ihcTowE, &m_ihciTowEta, &m_ihciTowPhi, &m_emcal_good, &m_ohc_good, &m_ihc_good, activeTriggerNames);
    
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
    int skippedEnergyCount = 0; // Counter for clusters skipped due to Energy cut
    float max_energy_clus = 0.0;
    float max_isoEt = std::numeric_limits<float>::lowest();
    float min_isoEt = std::numeric_limits<float>::max();
    std::map<int, std::pair<float, float>> clusterEtIsoMap_unsubtracted; //to store cluster ID and corresponding et iso value
    std::map<int, std::pair<float, float>> clusterEtIsoMap_subtracted;
    std::vector<int> m_clusterIds; // Store cluster IDs
    
    //  We'll store shape results in a map, keyed by cluster->get_id() -- DECLARED LOCALLY DO NOT NEED TO CLEAR AFTER EVENT -- AUTO RESETS as is not a member
    std::map<int, ShowerShapeVars> shapeVarsMap;
    
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
        
        //        CLHEP::Hep3Vector Ecore_vec_cluster = RawClusterUtility::GetECoreVec(*cluster, vertex);  NOT USING ECore in pp
        CLHEP::Hep3Vector clusterEnergy_Full = RawClusterUtility::GetEVec(*cluster, vertex);
        
        //        float clusEcore = Ecore_vec_cluster.mag(); NOT USING ECore in pp
        
        float clusEnergy = clusterEnergy_Full.mag();
        
        if (clusEnergy < 1.0) { // cut on cluster energy
            skippedEnergyCount++;
            continue;
        }
        
        float clus_eta = clusterEnergy_Full.pseudoRapidity();
        float clus_phi = clusterEnergy_Full.phi();
        
        ShowerShapeVars ssv = computeShowerShapesForCluster(
            cluster,
            emcTowerContainer, geomEM,
            ihcTowerContainer, geomIH,
            ohcTowerContainer, geomOH,
            m_vz,              // the z-vertex
            clus_eta,          // pass previously computed
            clus_phi           // pass previously computed
        );
        
        float clus_eT = clusEnergy / std::cosh(clus_eta);
        float clus_pt = clusterEnergy_Full.perp();
        float clus_chi = cluster -> get_chi2();
        float maxTowerEnergy = getMaxTowerE(cluster,emcTowerContainer);
        int clusID = cluster->get_id();

        shapeVarsMap[clusID] = ssv;
        /*
         Save map of clus ID's that pass shower shape cuts for this event
         */
        bool pass = applyShowerShapeCuts(ssv, clus_eT);
        clusterPassedShowerCuts[clusID] = pass;
        
        m_clusterIds.push_back(clusID);
        m_clusterEnergy.push_back(clusEnergy);
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
        
        
        if (clusEnergy > max_energy_clus) {
            max_energy_clus = clusEnergy; // Update the maximum cluster energy
        }
    }
    
    if (verbose) {
        std::cout << "\n--- Cluster Processing Summary ---" << std::endl;
        std::cout << "Clusters processed: " << m_clusterIds.size() << std::endl;
        std::cout << "Clusters skipped due to Cluster Energy < 1: " << skippedEnergyCount << std::endl;

        std::cout << "\nVector Sizes:" << std::endl;
        std::cout << "  m_clusterIds size: " << m_clusterIds.size() << std::endl;
        std::cout << "  m_clusterEnergy size: " << m_clusterEnergy.size() << std::endl;
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
        std::cout << std::setw(12) << "Cluster ID" << std::setw(20) << "Cluster Energy" << std::setw(20) << "Isolation Energy (et_iso)" << std::endl;
        std::cout << "------------------------------------------------" << std::endl;
        for (const auto& entry : clusterEtIsoMap_unsubtracted) {
            std::cout << std::setw(12) << entry.first << std::setw(20) << entry.second.first << std::setw(20) << entry.second.second << std::endl;
        }
        std::cout << "------------------------------------------------\n";
    }
    processClusterInvariantMass(m_clusterEnergy, m_clusterPt, m_clusterChi, m_clusterEta, m_clusterPhi, m_clusterIds, clusterEtIsoMap_unsubtracted, activeTriggerNames, clusterPassedShowerCuts);

    
    /*
     Loop over active trigger names
     */
    for (const std::string &firedShortName : activeTriggerNames) {
        
        auto& qaHistograms = qaHistogramsByTrigger[firedShortName];
        std::map<std::pair<float, float>, std::map<std::string, TObject*>>& qaIsolationHistograms = qaIsolationHistogramsByTriggerAndPt[firedShortName];

        // Check if the histogram exists before filling it
        std::string histName = "hVtxZ_" + firedShortName;
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
        TH1F* hTriggerCount = (TH1F*)qaHistogramsByTrigger[firedShortName]["hTriggerCount_" + firedShortName];
        if (hTriggerCount) {
            hTriggerCount->Fill(0.5); // Fill the single bin to count the trigger occurrence
        }
        
        if (verbose) {
            std::cout << "Processing Trigger: " << firedShortName << std::endl;
        }
        // ------------------------------------------------
        // Fill nClustersInEvent_pT_... histogram for each pT bin
        // ------------------------------------------------
        for (const auto &pT_bin : pT_bins)
        {
          float pT_min = pT_bin.first;
          float pT_max = pT_bin.second;
          std::pair<float,float> pT_range = {pT_min, pT_max};

          // Count how many clusters in this pT bin
          int nClustersInBin = 0;
          for (size_t iclus = 0; iclus < m_clusterPt.size(); ++iclus)
          {
            float thisPt = m_clusterPt[iclus];
            if (thisPt >= pT_min && thisPt < pT_max)
            {
              nClustersInBin++;
            }
          }

          // Build the histogram name
          std::string nClusHistName = "nClustersInEvent_pT_" +
              formatFloatForFilename(pT_min) + "to" +
              formatFloatForFilename(pT_max) + "_" + firedShortName;

          // Retrieve from qaIsolationHistograms (the same place you created it)
          TH1F* nClusHist = dynamic_cast<TH1F*>(
              qaIsolationHistograms[pT_range][nClusHistName]
          );

          if (nClusHist)
          {
            nClusHist->Fill(nClustersInBin);
          }
          else
          {
            if (verbose)
            {
              std::cerr << "[WARNING] Could not find histogram '"
                        << nClusHistName << "' for trigger '"
                        << firedShortName << "'" << std::endl;
            }
          }
        }
      

        //----------------------------------------------------------------------
        // Define pointers to all relevant histograms for this trigger ONCE for shower shape variable filling QA
        //----------------------------------------------------------------------
        TH1F* h_noCut_ratio         = (TH1F*) qaHistograms["E3by7_over_E7by7_NoShowerShapeCuts_" + firedShortName];
        TH1F* h_withCut_ratio       = (TH1F*) qaHistograms["E3by7_over_E7by7_withShowerShapeCuts_" + firedShortName];

        TH1F* h_noCut_w72           = (TH1F*) qaHistograms["w72_NoShowerShapeCuts_" + firedShortName];
        TH1F* h_withCut_w72         = (TH1F*) qaHistograms["w72_withShowerShapeCuts_" + firedShortName];

        TH1F* h_noCut_1x1overE      = (TH1F*) qaHistograms["E1x1_over_ClusterE_NoShowerShapeCuts_" + firedShortName];
        TH1F* h_withCut_1x1overE    = (TH1F*) qaHistograms["E1x1_over_ClusterE_withShowerShapeCuts_" + firedShortName];

        TH1F* h_noCut_1x1over3x3    = (TH1F*) qaHistograms["E1x1_over_E3x3_NoShowerShapeCuts_" + firedShortName];
        TH1F* h_withCut_1x1over3x3  = (TH1F*) qaHistograms["E1x1_over_E3x3_withShowerShapeCuts_" + firedShortName];

        TH1F* h_noCut_3x2over3x5    = (TH1F*) qaHistograms["E3x2_over_E3x5_NoShowerShapeCuts_" + firedShortName];
        TH1F* h_withCut_3x2over3x5  = (TH1F*) qaHistograms["E3x2_over_E3x5_withShowerShapeCuts_" + firedShortName];

        TH1F* h_noCut_1by7over7by7  = (TH1F*) qaHistograms["E1by7_over_E7by7_NoShowerShapeCuts_" + firedShortName];
        TH1F* h_withCut_1by7over7by7= (TH1F*) qaHistograms["E1by7_over_E7by7_withShowerShapeCuts_" + firedShortName];
        
        TH1F* h_noCut_3x3overE   = (TH1F*) qaHistograms["E3x3_over_ClusterE_NoShowerShapeCuts_" + firedShortName];
        TH1F* h_withCut_3x3overE = (TH1F*) qaHistograms["E3x3_over_ClusterE_withShowerShapeCuts_" + firedShortName];

        TH1F* h_noCut_3x3over3x7 = (TH1F*) qaHistograms["E3x3_over_E3x7_NoShowerShapeCuts_" + firedShortName];
        TH1F* h_withCut_3x3over3x7 = (TH1F*) qaHistograms["E3x3_over_E3x7_withShowerShapeCuts_" + firedShortName];

        TH1F* h_noCut_weta   = (TH1F*) qaHistograms["weta_NoShowerShapeCuts_" + firedShortName];
        TH1F* h_withCut_weta = (TH1F*) qaHistograms["weta_withShowerShapeCuts_" + firedShortName];

        TH1F* h_noCut_wphi   = (TH1F*) qaHistograms["wphi_NoShowerShapeCuts_" + firedShortName];
        TH1F* h_withCut_wphi = (TH1F*) qaHistograms["wphi_withShowerShapeCuts_" + firedShortName];


        //----------------------------------------------------------------------
        // Now loop over all clusters and fill SHOWER SHAPE histograms
        //----------------------------------------------------------------------
        for (size_t ic = 0; ic < m_clusterIds.size(); ++ic)
        {
            int   clusID   = m_clusterIds[ic];
            float clusE    = m_clusterEnergy[ic];

            // Retrieve the shape variables
            const auto &ssv   = shapeVarsMap[clusID];
            bool  passCuts    = clusterPassedShowerCuts[clusID]; // Did it pass your shape cuts?

            float weta_val = ssv.weta;
            float wphi_val = ssv.wphi;
            //----------------------------
            // Compute each ratio
            //----------------------------
            // E3x3
            float e3x3 = ssv.e3x3;          // from computeShowerShapesForCluster
            
            // (A) E3x7 / E7x7
            float ratio_3x7_over_7x7 = 0.f;
            if (ssv.e7x7 > 1e-6f) ratio_3x7_over_7x7 = ssv.e3x7 / ssv.e7x7;

            // (B) w72 => ssv.w72

            // (C) E1x1 / E_cluster
            float ratio_1x1_over_E = 0.f;
            if (clusE > 1e-6f) ratio_1x1_over_E = ssv.e1x1 / clusE;

            // (D) E1x1 / E3x3
            float ratio_1x1_over_3x3 = 0.f;
            if (ssv.e3x3 > 1e-6f) ratio_1x1_over_3x3 = ssv.e1x1 / ssv.e3x3;

            // (E) E1x7 / E7x7
            float ratio_1x7_over_7x7 = 0.f;
            if (ssv.e7x7 > 1e-6f) ratio_1x7_over_7x7 = ssv.e17 / ssv.e7x7;

            // (F) E3x2 / E3x5
            float ratio_3x2_over_3x5 = 0.f;
            if (ssv.e35 > 1e-6f) ratio_3x2_over_3x5 = ssv.e3x2 / ssv.e35;

            //(G) E3x3/Ecluster
            float ratio_3x3_over_E = 0.f;
            if (clusE > 1e-6f)
            {
              ratio_3x3_over_E = e3x3 / clusE;
            }

            //(H) E3x3/E3x7
            float ratio_3x3_over_3x7 = 0.f;
            if (ssv.e3x7 > 1e-6f)
            {
              ratio_3x3_over_3x7 = e3x3 / ssv.e3x7;
            }
            
            //---------------------------------------------------
            // Fill "No Cuts" histograms for every cluster
            //---------------------------------------------------
            if (h_noCut_ratio)               h_noCut_ratio->Fill(ratio_3x7_over_7x7);
            if (h_noCut_w72)                 h_noCut_w72->Fill(ssv.w72);
            if (h_noCut_1x1overE)            h_noCut_1x1overE->Fill(ratio_1x1_over_E);
            if (h_noCut_1x1over3x3)          h_noCut_1x1over3x3->Fill(ratio_1x1_over_3x3);
            if (h_noCut_1by7over7by7)        h_noCut_1by7over7by7->Fill(ratio_1x7_over_7x7);
            if (h_noCut_3x2over3x5)          h_noCut_3x2over3x5->Fill(ratio_3x2_over_3x5);
            if (h_noCut_3x3overE)            h_noCut_3x3overE->Fill(ratio_3x3_over_E);
            if (h_noCut_3x3over3x7)          h_noCut_3x3over3x7->Fill(ratio_3x3_over_3x7);
            if (h_noCut_weta)                h_noCut_weta->Fill(weta_val);
            if (h_noCut_wphi)                h_noCut_wphi->Fill(wphi_val);

            //---------------------------------------------------
            // If cluster passes shape cuts => fill "With Cuts"
            //---------------------------------------------------
            if (passCuts)
            {
                if (h_withCut_ratio)         h_withCut_ratio->Fill(ratio_3x7_over_7x7);
                if (h_withCut_w72)           h_withCut_w72->Fill(ssv.w72);
                if (h_withCut_1x1overE)      h_withCut_1x1overE->Fill(ratio_1x1_over_E);
                if (h_withCut_1x1over3x3)    h_withCut_1x1over3x3->Fill(ratio_1x1_over_3x3);
                if (h_withCut_1by7over7by7)  h_withCut_1by7over7by7->Fill(ratio_1x7_over_7x7);
                if (h_withCut_3x2over3x5)    h_withCut_3x2over3x5->Fill(ratio_3x2_over_3x5);
                if (h_withCut_3x3overE)      h_withCut_3x3overE->Fill(ratio_3x3_over_E);
                if (h_withCut_3x3over3x7)    h_withCut_3x3over3x7->Fill(ratio_3x3_over_3x7);
                if (h_withCut_weta)          h_withCut_weta->Fill(weta_val);
                if (h_withCut_wphi)          h_withCut_wphi->Fill(wphi_val);
            }
        } // end cluster loop

        // Process towers and fill histograms
        if (emcTowerContainer) {
            for (size_t i = 0; i < m_emcTowE.size(); ++i) {
                ((TH2F*)qaHistograms["h2_EMCal_TowerEtaPhi_Weighted_2D_" + firedShortName])->Fill(m_emciEta[i], m_emciPhi[i], m_emcTowE[i]);
            }
            // Check if the total EMCal energy is positive or negative and fill the appropriate histogram
            if (totalCaloEEMCal >= 0) {
                ((TH1F*)qaHistograms["hTotalCaloEEMCal_" + firedShortName])->Fill(totalCaloEEMCal);
            } else {
                ((TH1F*)qaHistograms["hTotalCalo_Negative_EEMCal_" + firedShortName])->Fill(totalCaloEEMCal);
            }
            //maximum chi2 of all the towers fill the histograms for EMCAL chi2 for this event -- similar for other calos
            if (!m_emcChi2.empty()) {
                float max_emcalChi2 = *std::max_element(m_emcChi2.begin(), m_emcChi2.end());
                ((TH1F*)qaHistograms["h_emcalChi2_" + firedShortName])->Fill(max_emcalChi2);
            }

        }
        if (ihcTowerContainer) {
            for (size_t i = 0; i < m_ihcTowE.size(); ++i) {
                ((TH2F*)qaHistograms["h2_IHCal_TowerEnergy_" + firedShortName])->Fill(m_ihciTowEta[i], m_ihciTowPhi[i], m_ihcTowE[i]);
            }
            ((TH1F*)qaHistograms["hTotalCaloEIHCal_" + firedShortName])->Fill(totalCaloEIHCal);
            if (!m_ihcChi2.empty()) {
                float max_ihcalChi2 = *std::max_element(m_ihcChi2.begin(), m_ihcChi2.end());
                ((TH1F*)qaHistograms["h_ihcalChi2_" + firedShortName])->Fill(max_ihcalChi2);
            }
        }
        
        if (ohcTowerContainer) {
            for (size_t i = 0; i < m_ohcTowE.size(); ++i) {
                ((TH2F*)qaHistograms["h2_OHCal_TowerEnergy_" + firedShortName])->Fill(m_ohciTowEta[i], m_ohciTowPhi[i], m_ohcTowE[i]);
            }
            ((TH1F*)qaHistograms["hTotalCaloEOHCal_" + firedShortName])->Fill(totalCaloEOHCal);
            if (!m_ohcChi2.empty()) {
                float max_ohcalChi2 = *std::max_element(m_ohcChi2.begin(), m_ohcChi2.end());
                ((TH1F*)qaHistograms["h_ohcalChi2_" + firedShortName])->Fill(max_ohcalChi2);
            }
        }
        
        ((TH1F*)qaHistograms["h8by8TowerEnergySum_" + firedShortName])->Fill(energyMaps.max_8by8energy_emcal);
        ((TH1F*)qaHistograms["h_hcal_energy_" + firedShortName])->Fill(energyMaps.max_energy_hcal);
        ((TH1F*)qaHistograms["h_jet_energy_" + firedShortName])->Fill(energyMaps.max_energy_jet);

        ((TH1F*)qaHistograms["h_jet_emcal_energy_" + firedShortName])->Fill(energyMaps.energymap_jet_emcal[energyMaps.jet_ebin][energyMaps.jet_pbin]);
        ((TH1F*)qaHistograms["h_jet_hcalin_energy_" + firedShortName])->Fill(energyMaps.energymap_jet_hcalin[energyMaps.jet_ebin][energyMaps.jet_pbin]);
        ((TH1F*)qaHistograms["h_jet_hcalout_energy_" + firedShortName])->Fill(energyMaps.energymap_jet_hcalout[energyMaps.jet_ebin][energyMaps.jet_pbin]);
        

        
        for (size_t i = 0; i < m_clusterIds.size(); ++i) {
            // Fill the eT histogram
            TH1F* h_ET = (TH1F*)qaHistograms["h_ET_" + firedShortName];
            if (!h_ET) {
                std::cerr << "Error: h_ET_" << firedShortName << " histogram is null!" << std::endl;
            } else {
                h_ET->Fill(m_clusterEt[i]);
            }
            
            // Check if histograms exist and fill them
            TH1F* h_maxTowE = (TH1F*)qaHistograms["h_maxTowerEnergy_" + firedShortName];
            if (!h_maxTowE) {
                std::cerr << "Error: h_maxTowerEnergy_" << firedShortName << " histogram is null!" << std::endl;
            } else {
                h_maxTowE->Fill(m_maxTowEnergy[i]);
            }
            
            // Access histograms from the map
            TH1F* hPt = (TH1F*)qaHistograms["hClusterPt_" + firedShortName];
            TH1F* hChi2 = (TH1F*)qaHistograms["hClusterChi2_" + firedShortName];

            // Check if histograms exist and fill them
            if (!hPt) {
                std::cerr << "Error: hClusterPt_" << firedShortName << " histogram is null!" << std::endl;
            } else {
                hPt->Fill(m_clusterPt[i]);
            }

            if (!hChi2) {
                std::cerr << "Error: hClusterChi2_" << firedShortName << " histogram is null!" << std::endl;
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
                std::string hist2DName = "h2_cluster_iso_Et_unsubtracted_" + firedShortName;
                std::string hist1DName = "h1_isoEt_unsubtracted_" + firedShortName;

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
                                formatFloatForFilename(pT_max) + "_" + firedShortName;
                            std::string hist1DName = "h1_isoEt_unsubtracted_pT_" +
                                formatFloatForFilename(pT_min) + "to" +
                                formatFloatForFilename(pT_max) + "_" + firedShortName;

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
                                          << firedShortName << std::endl;
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
                std::string hist2DName = "h2_cluster_iso_Et_subtracted_" + firedShortName;
                std::string hist1DName = "h1_isoEt_subtracted_" + firedShortName;

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
                                formatFloatForFilename(pT_max) + "_" + firedShortName;
                            std::string hist1DName = "h1_isoEt_subtracted_pT_" +
                                formatFloatForFilename(pT_min) + "to" +
                                formatFloatForFilename(pT_max) + "_" + firedShortName;

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
                                          << firedShortName << std::endl;
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
        TH1F* h_maxEnergyClus = (TH1F*)qaHistograms["h_maxEnergyClus_" + firedShortName];

        if (!h_maxEnergyClus) {
            std::cerr << "Error: Histogram h_maxEnergyClus_" << firedShortName << " is null and cannot be filled." << std::endl;
        } else {
            h_maxEnergyClus->Fill(max_energy_clus);
            
            if (verbose) {
                std::cout << "Filled histogram h_maxEnergyClus_" << firedShortName
                          << " with value: " << max_energy_clus << std::endl;
            }
        }
    }
    return Fun4AllReturnCodes::EVENT_OK;
}


int caloTreeGen::process_event_Sim(PHCompositeNode *topNode)
{
  //----------------------------------------------------------------------
  // [0] Optional debug print
  // If 'verbose' is ON, log that we have entered this function
  //----------------------------------------------------------------------
  if (verbose)
  {
    std::cout << "\n[DEBUG] Entering caloTreeGen::process_event_Sim()"
              << " => event_count=" << event_count
              << ", topNode="      << topNode << std::endl;
  }

  //----------------------------------------------------------------------
  // [1] Update event counter and check user limit
  //----------------------------------------------------------------------
  event_count++;

  std::cout << "\n========== Processing CALOTREEGEN SIMN MODE -- Event " << event_count << " ==========\n";
  //----------------------------------------------------------------------
  // [2] Retrieve next entry from TTree "slimTree"
  //     We keep a static 'iEntry' to remember our position
  //----------------------------------------------------------------------
  static Long64_t iEntry = 0;
  if (!slimTree || iEntry >= slimTree->GetEntries())
  {
    // If TTree pointer is null or no more entries => stop
    std::cerr << "[WARNING] 'slimTree' is null or no more entries to read.\n";
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Actually load the data arrays for this event
  slimTree->GetEntry(iEntry++);
  
  //----------------------------------------------------------------------
  // [2A] If verbose, let’s print info about the TTree structure once,
  //      e.g., its branches and leaves
  //----------------------------------------------------------------------
  static bool printedTreeInfo = false;
  if (verbose && !printedTreeInfo)
  {
    printedTreeInfo = true; // so we only do this once

    std::cout << "[VERBOSE] Printing available branches in 'slimTree':\n";
    TObjArray* branchList = slimTree->GetListOfBranches();
    if (branchList)
    {
      for (int ib = 0; ib < branchList->GetEntries(); ib++)
      {
        TBranch* br = (TBranch*) branchList->At(ib);
        if (!br) continue;
        std::cout << "  Branch name: " << br->GetName()
                  << ", Title: "      << br->GetTitle() << std::endl;
      }
    }
    else
    {
      std::cout << "  [WARNING] No branches found in slimTree?\n";
    }
  }

  //----------------------------------------------------------------------
  // [3] Fill global isolation-based histograms: hIsoFromPi0Eta, hIsoNotPi0Eta
  //----------------------------------------------------------------------
  for (int ic = 0; ic < ncluster_CEMC_SIM; ++ic)
  {
    float isoVal = cluster_iso_04_CEMC_SIM[ic];
    float cEt    = cluster_Et_CEMC_SIM[ic];
    int   pid    = cluster_pid_CEMC_SIM[ic]; // e.g. 22=photon,111=pi0,221=eta,...

    // We require isoVal <= 6 to be "isolated" in this definition
    if (isoVal <= 6.f)
    {
      // If cluster is from pi0 or eta => fill hIsoFromPi0Eta
      if (pid == 111 || pid == 221)
      {
        if (hIsoFromPi0Eta) hIsoFromPi0Eta->Fill(cEt);
      }
      // else fill hIsoNotPi0Eta
      else
      {
        if (hIsoNotPi0Eta) hIsoNotPi0Eta->Fill(cEt);
      }
    }
  }

  //----------------------------------------------------------------------
  // [4] We define "prompt cluster" => (isoVal<6, cEt>5).
  //     We'll measure how many such clusters are truly prompt (photonclass=1).
  //----------------------------------------------------------------------
  // Create counters for each pT bin
  std::map<std::pair<float, float>, int> nTaggedPrompt_byPt;
  std::map<std::pair<float, float>, int> nCorrectlyPrompt_byPt;
  for (auto &bin : pT_bins)
  {
    nTaggedPrompt_byPt[bin]       = 0;
    nCorrectlyPrompt_byPt[bin]    = 0;
  }

  //----------------------------------------------------------------------
  // [5] Additional counters: total clusters that are prompt, frag, decay, other
  //----------------------------------------------------------------------
  int totalPrompt = 0;
  int totalFrag   = 0;
  int totalDecay  = 0;
  int totalOther  = 0;

  //----------------------------------------------------------------------
  // [6] Main loop over clusters => classification
  //----------------------------------------------------------------------
  for (int ic = 0; ic < ncluster_CEMC_SIM; ++ic)
  {
    float cEt     = cluster_Et_CEMC_SIM[ic];
    float cEta    = cluster_Eta_CEMC_SIM[ic];
    float cPhi    = cluster_Phi_CEMC_SIM[ic];
    float isoVal  = cluster_iso_04_CEMC_SIM[ic];

    // [6.1] Find which pT bin we belong to
    bool inRange = false;
    std::pair<float,float> matchedPtBin;
    for (auto &bin : pT_bins)
    {
      if (cEt >= bin.first && cEt < bin.second)
      {
        matchedPtBin = bin;
        inRange = true;
        break;
      }
    }
    if (!inRange) continue; // skip cluster if no pT bin match

    // [6.2] Find best truth photon => minimal deltaR<0.05
    int   bestMatchIndex = -1;
    float bestDR         = 9999.f;
    for (int ip = 0; ip < nparticles_SIM; ++ip)
    {
      // must be pid=22 => real photon
      if (particle_pid_SIM[ip] != 22) continue;

      float dR = deltaR(cEta, particle_Eta_SIM[ip], cPhi, particle_Phi_SIM[ip]);
      if (dR < 0.05 && dR < bestDR)
      {
        bestDR = dR;
        bestMatchIndex = ip;
      }
    }

    // [6.3] Derive truthClass => 1=prompt,2=frag,3=decay,4=other
    int truthClass = 4;
    if (bestMatchIndex >= 0)
    {
      int photClass = particle_photonclass_SIM[bestMatchIndex];
      switch (photClass)
      {
        case 1: truthClass = 1; break; // prompt
        case 2: truthClass = 2; break; // frag
        case 3: truthClass = 3; break; // decay
        default: truthClass = 4; break;
      }
    }

    // [6.4] Increment counters for the entire event
    if      (truthClass == 1) totalPrompt++;
    else if (truthClass == 2) totalFrag++;
    else if (truthClass == 3) totalDecay++;
    else                      totalOther++;

    // [6.5] Fill “classification vs. pT bin” hist, e.g. "hClusterTruthClass_pT_2.0to3.0"
    if (hClusterTruthClass_pTbin.count(matchedPtBin))
    {
      hClusterTruthClass_pTbin[matchedPtBin]->Fill(truthClass);
    }

    // [6.6] Check if cluster is "prompt candidate" => isoVal<6, cEt>5
    bool clusterIsPromptCandidate = (isoVal < 6.f && cEt > 5.f);
    if (clusterIsPromptCandidate)
    {
      nTaggedPrompt_byPt[matchedPtBin]++;
      if (truthClass == 1) // truly prompt
      {
        nCorrectlyPrompt_byPt[matchedPtBin]++;
      }
    }
  } // end cluster loop

  //----------------------------------------------------------------------
  // [7] Now compute the “prompt purity” => (# truly prompt) / (# tagged)
  //----------------------------------------------------------------------
  for (auto &bin : pT_bins)
  {
    int nTagged  = nTaggedPrompt_byPt[bin];
    int nCorrect = nCorrectlyPrompt_byPt[bin];
    if (nTagged > 0)
    {
      float ratio = float(nCorrect) / float(nTagged);
      // e.g. "hPromptPurity_pT_2.0to3.0"
      if (hPromptPurity_pTbin.count(bin))
      {
        hPromptPurity_pTbin[bin]->Fill(ratio);
      }
    }
  }

  //----------------------------------------------------------------------
  // [8] Also fill distribution of *all* prompt photons in pT
  //     => look for (pid=22, photonclass=1)
  //----------------------------------------------------------------------
  int nPromptPhotons = 0;
  for (int ip = 0; ip < nparticles_SIM; ++ip)
  {
    if (particle_pid_SIM[ip] == 22 && particle_photonclass_SIM[ip] == 1)
    {
      float truthPt = particle_Pt_SIM[ip];
      if (hPhotonPtPrompt) hPhotonPtPrompt->Fill(truthPt);
      nPromptPhotons++;
    }
  }

  //----------------------------------------------------------------------
  // [9] If 'verbose', print a summary for this event
  //----------------------------------------------------------------------
  if (verbose)
  {
    int totalClusters   = ncluster_CEMC_SIM;
    int totalBackground = totalFrag + totalDecay + totalOther;

    std::cout << "[SIM] Event " << event_count
              << " (iEntry=" << iEntry << ") summary:\n"
              << "   => #clusters=" << totalClusters
              << " => #prompt=" << totalPrompt
              << ", #frag="   << totalFrag
              << ", #decay="  << totalDecay
              << ", #other="  << totalOther
              << "  (Total BG=" << totalBackground << ")\n"
              << "   => # prompt truth photons in event=" << nPromptPhotons
              << "\n   => Done with process_event_Sim.\n";
  }

  //----------------------------------------------------------------------
  // All done => success
  //----------------------------------------------------------------------
  return Fun4AllReturnCodes::EVENT_OK;
}




int caloTreeGen::ResetEvent(PHCompositeNode *topNode)
{
    if (verbose)
    {
        std::cout << ANSI_COLOR_BLUE_BOLD
                  << "ResetEvent called: wantData=" << wantData
                  << ", wantSim=" << wantSim << ANSI_COLOR_RESET << std::endl;
    }

    // If user wants data
    if (wantData)
    {
        resetEvent_Data(topNode);
    }
    
    // If user wants simulation
    if (wantSim)
    {
        resetEvent_Sim(topNode);
    }

    // If both wantData & wantSim are set true,
    // we clear both sets of vectors in parallel.

    // Return
    return Fun4AllReturnCodes::EVENT_OK;
}


//____________________________________________________________________________..
int caloTreeGen::resetEvent_Data(PHCompositeNode *topNode) {
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
    m_clusterEnergy.clear();
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
    
    clusterPassedShowerCuts.clear();
    m_usedClustersThisEvent.clear();
    
    return Fun4AllReturnCodes::EVENT_OK;

    if (verbose) {
        std::cout << ANSI_COLOR_GREEN_BOLD << "Event reset complete." << ANSI_COLOR_RESET << std::endl;
    }
}

int caloTreeGen::resetEvent_Sim(PHCompositeNode *topNode)
{
    if (verbose)
    {
        std::cout << ANSI_COLOR_BLUE_BOLD
                  << "[ResetEvent-Sim] Clearing simulation-related arrays..."
                  << ANSI_COLOR_RESET << std::endl;
    }


    if (verbose)
    {
        std::cout << ANSI_COLOR_GREEN_BOLD
                  << "[ResetEvent-Sim] Done clearing sim vectors."
                  << ANSI_COLOR_RESET << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
}


int caloTreeGen::End(PHCompositeNode *topNode)
{
    if (verbose)
    {
        std::cout << ANSI_COLOR_BLUE_BOLD
                  << "caloTreeGen::End() => All events done. Starting final steps..."
                  << ANSI_COLOR_RESET << std::endl;
    }

    // 1) If user wants data, do the data end-of-job
    if (wantData)
    {
        endData(topNode);
    }

    // 2) If user wants sim, do the sim end-of-job
    if (wantSim)
    {
        endSim(topNode);
    }

    // If both are true, it executes both. If only one is true,
    // it executes just that one.
    // Then we can return success code:
    return Fun4AllReturnCodes::EVENT_OK;
}


int caloTreeGen::endData(PHCompositeNode *topNode)
{
    // Always announce we’re beginning final steps
    std::cout << ANSI_COLOR_BLUE_BOLD
              << "caloTreeGen::End(...) => All events have been processed. Beginning final analysis..."
              << ANSI_COLOR_RESET << std::endl;

    // 1) Ensure the output file is open
    if (!out || !out->IsOpen())
    {
        std::cerr << ANSI_COLOR_RED_BOLD
                  << "Error: Output file is not open. Exiting End() without writing histograms."
                  << ANSI_COLOR_RESET << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
    }
    out->cd();  // go to the root dir in the output file
    gDirectory = out; // explicitly set the global ROOT directory

    // 2) Data structures to keep track of duplicates + final hist list to write
    std::unordered_set<std::string> usedNames; // tracks histogram names we’ve encountered
    bool foundDuplicate = false;

    // We’ll store all histograms we intend to write in a vector so we can do a progress print
    // Each entry is: (triggerDirName, histName, pointerToHist)
    struct HistToWrite
    {
      std::string triggerDir;
      std::string histName;
      TObject*    histPtr;
    };
    std::vector<HistToWrite> finalHistList; // everything we actually want to write

    // 3) Helper to see if a TH1 is zero‐entry => skip
    auto isZeroEntry = [&](TObject* obj)
    {
      TH1* checkHist = dynamic_cast<TH1*>(obj);
      if (checkHist)
      {
        return (checkHist->GetEntries() == 0);
      }
      // For non‐TH1 types, you can choose to treat them as never empty or also skip.
      return false;
    };

    // 4) "Dry run" collecting histograms from each trigger:
    //    - We first check if `hTriggerCount_<trigger>` is zero.
    //    - If so => skip the entire trigger folder quickly.
    //    - If not zero => we gather all *non-empty* hists, skip duplicates, etc.

    // This function collects the main QA histograms for a single trigger.
    // A lambda that collects QA histograms from the "main" trigger map
    //  but can also handle map<std::string, TH1F*> or TH2F*, etc.
    auto collectTriggerDir = [&](const std::string& trigName, auto & histMap)
    {
      // Deduce the map's mapped_type (could be TObject*, TH1F*, etc.)
      using MMap    = std::remove_reference_t<decltype(histMap)>;
      using ValuePtr = typename MMap::mapped_type;  // e.g. TH1F*

      // 1) Check if hTriggerCount_<trigName> exists and is >0
      std::string trigCountName = "hTriggerCount_" + trigName;
      ValuePtr trigCountObj = nullptr;
      if (auto it = histMap.find(trigCountName); it != histMap.end())
      {
        trigCountObj = it->second;  // e.g. TH1F*
      }

      // If no trigger count histogram or zero entries => skip entire folder
      bool skipEntireTrigger = true;
      if (trigCountObj)
      {
        // cast from TH1F* (or T*) to TH1*
        //  (TH1 inherits from TObject, so dynamic_cast is safe)
        TH1* hTrigCount = dynamic_cast<TH1*>(static_cast<TObject*>(trigCountObj));
        if (hTrigCount && hTrigCount->GetEntries() > 0)
        {
          skipEntireTrigger = false;
        }
      }
      if (skipEntireTrigger)
      {
        if (verbose)
        {
          std::cout << "[INFO] Skipping entire trigger directory \""
                    << trigName << "\" => hTriggerCount is zero or missing.\n";
        }
        return; // do not collect anything for this trigger
      }

      // 2) Otherwise, gather non-empty histograms, skip duplicates
      int countNonEmpty = 0;
      for (auto &kv : histMap)
      {
        const std::string& hName = kv.first;
        ValuePtr hObj            = kv.second; // e.g. TH1F* pointer

        // If pointer is null => skip
        if (!hObj) continue;

        // cast to TH1* for the "get entries" check
        TH1* checkHist = dynamic_cast<TH1*>(static_cast<TObject*>(hObj));
        if (checkHist && checkHist->GetEntries() == 0)
        {
          // zero‐entry => skip
          if (verbose)
          {
            std::cout << "[INFO] Trigger=" << trigName
                      << ": skipping histogram \"" << hName
                      << "\" => zero entries.\n";
          }
          continue;
        }

        // Check duplicates by name
        if (usedNames.count(hName) > 0)
        {
          std::cerr << ANSI_COLOR_RED_BOLD
                    << "[ERROR] Duplicate histogram name: \""
                    << hName << "\" => aborting."
                    << ANSI_COLOR_RESET << std::endl;
          foundDuplicate = true;
          return; // immediate stop
        }
        usedNames.insert(hName);

        countNonEmpty++;

        // We'll store it in finalHistList as a TObject*
        TObject* castObj = static_cast<TObject*>(hObj);
        finalHistList.push_back({ trigName, hName, castObj });
      }

      // If countNonEmpty==0 => everything was zero => skip entire dir
      if (countNonEmpty == 0 && verbose)
      {
        std::cout << "[INFO] Trigger \"" << trigName
                  << "\": all histograms zero => skipping directory.\n";
      }
    };

    // For pT maps or cut combos => gather them
    auto collectTriggerPtMap = [&](const std::string &trigName, auto &ptMap)
    {
      // outer loop => pT bins or cut combos
      // we do *not* skip the entire trigger if zero, because it may partially fill some sub-bins
      for (auto &outerKV : ptMap)
      {
        // outerKV.second => a map<string, TObject*>
        auto & innerHistMap = outerKV.second;
        for (auto & kv : innerHistMap)
        {
          const std::string& hName = kv.first;
          TObject* hObj           = kv.second;

          if (!hObj || isZeroEntry(hObj))
          {
            if (verbose)
            {
              std::cout << "[INFO] Trigger=" << trigName
                        << ": skipping hist \"" << hName
                        << "\" => null or zero entries.\n";
            }
            continue;
          }
          if (usedNames.count(hName) > 0)
          {
            std::cerr << ANSI_COLOR_RED_BOLD
                      << "[ERROR] Duplicate histogram name: \""
                      << hName << "\" => aborting."
                      << ANSI_COLOR_RESET << std::endl;
            foundDuplicate = true;
            return;
          }
          usedNames.insert(hName);
          finalHistList.push_back( { trigName, hName, hObj } );
        }
        if (foundDuplicate) return;
      }
    };

    // 5) Now gather from the main "qaHistogramsByTrigger" maps
    for (auto & trigKV : qaHistogramsByTrigger)
    {
      if (foundDuplicate) break; // stop if we found a conflict
      const std::string & shortTriggerName = trigKV.first;
      auto & histMap = trigKV.second;
      collectTriggerDir(shortTriggerName, histMap);
    }

    // Then gather from the isolation hist
    for (auto & trigKV : qaIsolationHistogramsByTriggerAndPt)
    {
      if (foundDuplicate) break;
      const std::string & shortTriggerName = trigKV.first;
      auto & ptMap = trigKV.second;
      collectTriggerPtMap(shortTriggerName, ptMap);
    }

    // Then from massAndIsolationHistograms
    for (auto & trigKV : massAndIsolationHistograms)
    {
      if (foundDuplicate) break;
      const std::string & shortTriggerName = trigKV.first;

      // cut combos => pT bins => hist
      auto & outerMap = trigKV.second;
      for (auto & cutKV : outerMap)
      {
        auto & pTHistMap = cutKV.second;
        for (auto & ptBinKV : pTHistMap)
        {
          auto & histMap2 = ptBinKV.second;
          for (auto & kv : histMap2)
          {
            const std::string & hName = kv.first;
            TObject* hObj            = kv.second;
            if (!hObj || isZeroEntry(hObj))
            {
              if (verbose)
              {
                std::cout << "[INFO] Trigger=" << shortTriggerName
                          << ": skipping \"" << hName
                          << "\" => null or zero entries.\n";
              }
              continue;
            }
            if (usedNames.count(hName) > 0)
            {
              std::cerr << ANSI_COLOR_RED_BOLD
                        << "[ERROR] Duplicate name \"" << hName
                        << "\" => aborting.\n"
                        << ANSI_COLOR_RESET;
              foundDuplicate = true;
              break;
            }
            usedNames.insert(hName);
            finalHistList.push_back( { shortTriggerName, hName, hObj } );
          }
          if (foundDuplicate) break;
        }
        if (foundDuplicate) break;
      }
    }

    // Then massAndIsolationHistogramsNoPtBins
    for (auto & trigKV : massAndIsolationHistogramsNoPtBins)
    {
      if (foundDuplicate) break;
      const std::string & shortTriggerName = trigKV.first;
      auto & histMapNoPt = trigKV.second;
      collectTriggerDir(shortTriggerName, histMapNoPt);
    }

    // 6) If found a duplicate => abort now
    if (foundDuplicate)
    {
      std::cerr << ANSI_COLOR_RED_BOLD
                << "Aborting => duplicate histogram found. Closing file..."
                << ANSI_COLOR_RESET << std::endl;

      out->Close();
      delete out;
      out = nullptr;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    // 7) We now have finalHistList of all histograms to write => do a progress tracker
    size_t totalHists = finalHistList.size();
    if (verbose)
    {
      std::cout << "[INFO] Ready to write " << totalHists
                << " non-empty histograms in total.\n";
    }
    size_t hCounter = 0;

    // 8) Write them
    for (auto & item : finalHistList)
    {
      hCounter++;
      // Print progress every 500 if verbose
      if (verbose && (hCounter % 500 == 0))
      {
        std::cout << "[INFO] Wrote " << hCounter
                  << " / " << totalHists
                  << " histograms so far...\n";
      }

      // cd to the directory if it exists
      TDirectory* trigDir = out->GetDirectory(item.triggerDir.c_str());
      if (!trigDir)
      {
        if (verbose)
        {
          std::cout << "[INFO] Directory \"" << item.triggerDir
                    << "\" not found => skipping histogram \""
                    << item.histName << "\".\n";
        }
        continue;
      }
      trigDir->cd();

      // Attempt to write
      if (item.histPtr->Write() == 0)
      {
        std::cerr << ANSI_COLOR_RED_BOLD
                  << "[ERROR] Failed to write hist \""
                  << item.histName << "\".\n"
                  << ANSI_COLOR_RESET;
      }
      else if (verbose)
      {
        std::cout << "[INFO] Wrote histogram \"" << item.histName
                  << "\" in directory \"" << item.triggerDir << "\".\n";
      }
    }

    // final progress message
    if (verbose)
    {
      std::cout << "[INFO] Done writing all " << hCounter
                << " histograms.\n";
    }

    // Return to root dir
    out->cd();

    // 9) Now close the file
    std::cout << ANSI_COLOR_BLUE_BOLD
              << "Closing output file and cleaning up..."
              << ANSI_COLOR_RESET << std::endl;
    if (out)
    {
        out->Close();
        delete out;
        out = nullptr;
        if (verbose)
        {
            std::cout << ANSI_COLOR_GREEN_BOLD
                      << "Output file successfully closed and deleted."
                      << ANSI_COLOR_RESET << std::endl;
        }
    }
    else
    {
        std::cerr << ANSI_COLOR_RED_BOLD
                  << "Warning: Output file was already null, possibly already closed or deleted."
                  << ANSI_COLOR_RESET << std::endl;
    }

    std::cout << ANSI_COLOR_GREEN_BOLD
              << "End of caloTreeGen::End(PHCompositeNode *topNode). Exiting smoothly."
              << ANSI_COLOR_RESET << std::endl;

    return Fun4AllReturnCodes::EVENT_OK;
}


int caloTreeGen::endSim(PHCompositeNode *topNode)
{
  // [1] Print a message that we are writing out simulation histograms
  std::cout << ANSI_COLOR_BLUE_BOLD
            << "[SIM End] Writing simulation histograms..."
            << ANSI_COLOR_RESET << std::endl;

  // [2] Check that outSim is valid
  if (!outSim || !outSim->IsOpen())
  {
    std::cerr << "[ERROR] outSim is null or not open. No SIM histos can be written.\n";
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // [3] Make sure to go to outSim root directory
  outSim->cd();

  // [4] A helper lambda to write a histogram if it exists
  auto writeHistIfExists = [&](TH1F* histPtr, const std::string& histName)
  {
    if (!histPtr) return; // skip if pointer is null
    if (histPtr->Write() == 0)
    {
      std::cerr << "[ERROR] Failed to write " << histName << ".\n";
    }
    else if (verbose)
    {
      std::cout << "[INFO] Successfully wrote " << histName
                << " to outSim." << std::endl;
    }
  };

  // [5] Write the “global” histograms that are not pT-binned
  writeHistIfExists(hIsoFromPi0Eta,   "hIsoFromPi0Eta");
  writeHistIfExists(hIsoNotPi0Eta,    "hIsoNotPi0Eta");
  writeHistIfExists(hPhotonPtPrompt,  "hPhotonPtPrompt");

  // [6] Write out the pT-binned histograms: classification & purity
  for (auto &bin : pT_bins)
  {
    // classification
    if (hClusterTruthClass_pTbin.count(bin))
    {
      TH1F* hClass = hClusterTruthClass_pTbin[bin];
      if (hClass)
      {
        std::string histName = "hClusterTruthClass_pT_" +
          formatFloatForFilename(bin.first) + "to" + formatFloatForFilename(bin.second);
        writeHistIfExists(hClass, histName);
      }
    }
    // purity
    if (hPromptPurity_pTbin.count(bin))
    {
      TH1F* hPurity = hPromptPurity_pTbin[bin];
      if (hPurity)
      {
        std::string histName = "hPromptPurity_pT_" +
          formatFloatForFilename(bin.first) + "to" + formatFloatForFilename(bin.second);
        writeHistIfExists(hPurity, histName);
      }
    }
  }

  // [7] Close the simulation file
  outSim->Close();
  delete outSim;
  outSim = nullptr;

  if (verbose)
  {
    std::cout << "[INFO] endSim: Finished writing SIM histos & closed outSim.\n";
  }
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

double caloTreeGen::getTowerEta(RawTowerGeom *tower_geom, double vx, double vy, double vz)
{
  float r;
  if (vx == 0 && vy == 0 && vz == 0)
  {
    r = tower_geom->get_eta();
  }
  else
  {
    double radius = sqrt((tower_geom->get_center_x() - vx) * (tower_geom->get_center_x() - vx) + (tower_geom->get_center_y() - vy) * (tower_geom->get_center_y() - vy));
    double theta = atan2(radius, tower_geom->get_center_z() - vz);
    r = -log(tan(theta / 2.));
  }
  return r;
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
