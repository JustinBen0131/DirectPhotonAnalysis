#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllUtils.h>
#include <ffamodules/CDBInterface.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>

#include <caloreco/CaloTowerStatus.h>

#include <phool/recoConsts.h>
#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <jetbase/FastJetAlgo.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <g4jets/TruthJetInput.h>

#include <jetbackground/CopyAndSubtractJets.h>
#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/DetermineTowerRho.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/SubtractTowers.h>
#include <jetbackground/SubtractTowersCS.h>
#include <jetbackground/TowerRho.h>
#include <calotrigger/TriggerRunInfoReco.h>
#include "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/clusterIsoCopy_src/ClusterIso.h"

#include <calotreegen/caloTreeGen.h>

#include <Calo_Calib.C>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffarawobjects.so)
R__LOAD_LIBRARY(libcaloTreeGen.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbase.so)
R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(/sphenix/user/patsfan753/install/lib/libclusteriso.so)

namespace Enable
{
  bool HIJETS = true;
  int HIJETS_VERBOSITY = 0;
  bool HIJETS_MC = false;
  bool HIJETS_TRUTH = false;
}  // namespace Enable

namespace HIJETS
{
  bool do_flow = false; // should be set to true once the EPD event plane correction is implemented
  bool do_CS = false;
  bool is_pp = false;  // turn off functionality only relevant for nucleon collisions
  std::string tower_prefix = "TOWERINFO_CALIB";
}  // namespace HIJETS

#endif



static const bool WANT_VERBOSE = true;
static const bool WANT_MANUAL_EVENT_LIMIT = true;  // or false if you want unlimited
static const int MANUAL_EVENT_LIMIT       = 10000;  // e.g. 5k

/**
 * @brief Macro to run either Data pipeline or Simulation pipeline for caloTreeGen
 *
 * @param nEvents   # of events to process (0 => process all)
 * @param listFile  A text file containing your data file list (data) or first line for run extraction
 * @param inName    Output name for data (if runData) or user-chosen label
 * @param runSim    If true => runs the SIM pipeline only
 * @param runData   If true => runs the DATA pipeline
 */
void Fun4All_CaloTreeGen(const int nEvents = 0,
                         const char *listFile = "input_files.list",
                         const char *inName = "commissioning.root",
                         bool runSim = false,
                         bool runData = false)
{
    
    if (WANT_VERBOSE) {
        std::cout << "\n[DEBUG] Entering Fun4All_CaloTreeGen function...\n";
        std::cout << "    nEvents = " << nEvents << std::endl;
        std::cout << "    listFile = " << listFile << std::endl;
        std::cout << "    inName   = " << inName << std::endl;
        std::cout << "    runSim   = " << (runSim ? "true" : "false") << std::endl;
        std::cout << "    runData  = " << (runData ? "true" : "false") << std::endl;
        
        std::cout << "[INFO] Starting Fun4All_CaloTreeGen..." << std::endl;
    }

    
    Fun4AllServer *se = Fun4AllServer::instance();
    if (WANT_VERBOSE) {
        std::cout << "[DEBUG] Fun4AllServer instance acquired: " << se << std::endl;
    }
    
    
    // If user wants MC truth jets
    if (Enable::HIJETS_MC && Enable::HIJETS_TRUTH)
    {
        std::cout << "[DEBUG] HIJETS_MC && HIJETS_TRUTH both true. Creating truthjetreco...\n";
        JetReco *truthjetreco = new JetReco();
        TruthJetInput *tji = new TruthJetInput(Jet::PARTICLE);
        tji->add_embedding_flag(0);  // changes depending on signal vs. embedded
        truthjetreco->add_input(tji);
        truthjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.2), "AntiKt_Truth_r02");
        truthjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.3), "AntiKt_Truth_r03");
        truthjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.4), "AntiKt_Truth_r04");
        truthjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.5), "AntiKt_Truth_r05");
        truthjetreco->set_algo_node("ANTIKT");
        truthjetreco->set_input_node("TRUTH");
        truthjetreco->Verbosity(0);
        std::cout << "[INFO] Registering truthjetreco subsystem..." << std::endl;
        se->registerSubsystem(truthjetreco);
    }
    gSystem->Load("libg4dst");
    
    // Basic run config
    recoConsts *rc = recoConsts::instance();
    if (WANT_VERBOSE) {
        std::cout << "[DEBUG] Setting CDB_GLOBALTAG to 'ProdA_2024'..." << std::endl;
    }
    rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
    
    // Read the first filename to extract the run number
    if (WANT_VERBOSE) {
        std::cout << "[DEBUG] Attempting to open listFile: " << listFile << std::endl;
    }
    std::ifstream infile(listFile);
    std::string firstFilename;
    if (!infile.is_open())
    {
        std::cerr << "[ERROR] Could not open input file list: " << listFile << std::endl;
        return;
    }
    if (!std::getline(infile, firstFilename))
    {
        std::cerr << "[ERROR] Input file list is empty: " << listFile << std::endl;
        return;
    }
    std::cout << "[DEBUG] First filename read: " << firstFilename << std::endl;
    
    int runnumber = -1;  // Declare at an outer scope
    if (runData)
    {
        // Extract run number from the first filename
        std::pair<int, int> runseg = Fun4AllUtils::GetRunSegment(firstFilename);
        runnumber = runseg.first;
        int segnumber = runseg.second;
        if (WANT_VERBOSE) {
            std::cout << "[DEBUG] Extracted run: " << runnumber
            << " segment: " << segnumber << std::endl;
        }
        
        if (runnumber <= 0)
        {
            std::cerr << "[ERROR] Invalid run number extracted from first file: "
            << runnumber << ". Exiting..." << std::endl;
            return;
        }
        rc->set_uint64Flag("TIMESTAMP", runnumber);
        
        if (WANT_VERBOSE) {
            std::cout << "[DEBUG] Setting up CaloTowerStatus modules...\n";
        }
        CaloTowerStatus *statusEMC = new CaloTowerStatus("CEMCSTATUS");
        statusEMC->set_detector_type(CaloTowerDefs::CEMC);
        statusEMC->set_time_cut(1);
        statusEMC->set_inputNodePrefix("TOWERINFO_CALIB_");
        statusEMC->Verbosity(0);
        if (WANT_VERBOSE) {
            std::cout << "[INFO] Registering statusEMC subsystem..." << std::endl;
        }
        se->registerSubsystem(statusEMC);
        
        CaloTowerStatus *statusHCalIn = new CaloTowerStatus("HCALINSTATUS");
        statusHCalIn->set_detector_type(CaloTowerDefs::HCALIN);
        statusHCalIn->set_time_cut(2);
        statusHCalIn->set_inputNodePrefix("TOWERINFO_CALIB_");
        statusHCalIn->Verbosity(0);
        if (WANT_VERBOSE) {
            std::cout << "[INFO] Registering towerjetreco subsystem (statusHCalIn)..." << std::endl;
        }
        se->registerSubsystem(statusHCalIn);
        
        CaloTowerStatus *statusHCALOUT = new CaloTowerStatus("HCALOUTSTATUS");
        statusHCALOUT->set_detector_type(CaloTowerDefs::HCALOUT);
        statusHCALOUT->set_time_cut(2);
        statusHCALOUT->set_inputNodePrefix("TOWERINFO_CALIB_");
        statusHCALOUT->Verbosity(0);
        if (WANT_VERBOSE) {
            std::cout << "[INFO] Registering statusHCALOUT subsystem..." << std::endl;
        }
        se->registerSubsystem(statusHCALOUT);
        
        RetowerCEMC *rcemc = new RetowerCEMC();
        rcemc->Verbosity(0);
        rcemc->set_towerinfo(true);
        rcemc->set_frac_cut(0.5); // fraction of retower that must be masked
        rcemc->set_towerNodePrefix(HIJETS::tower_prefix);
        if (WANT_VERBOSE) {
            std::cout << "[INFO] Registering RetowerCEMC subsystem..." << std::endl;
        }
        se->registerSubsystem(rcemc);
        
        if (WANT_VERBOSE) {
            std::cout << "[DEBUG] Creating EmcRawClusterBuilderTemplate (ClusterBuilder)..." << std::endl;
        }
        RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
        ClusterBuilder->Detector("CEMC");
        ClusterBuilder->set_threshold_energy(0.030); // 30 MeV threshold
        std::string emc_prof = getenv("CALIBRATIONROOT");
        emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
        ClusterBuilder->LoadProfile(emc_prof);
        ClusterBuilder->set_UseTowerInfo(1);  // using TowerInfo objects
        std::cout << "[INFO] Registering ClusterBuilder subsystem..." << std::endl;
        se->registerSubsystem(ClusterBuilder);
        
        if (WANT_VERBOSE) {
            std::cout << "[DEBUG] Creating towerjetreco (JetReco) with Tower inputs...\n";
        }
        JetReco *towerjetreco = new JetReco();
        towerjetreco->add_input(new TowerJetInput(Jet::CEMC_TOWERINFO_RETOWER, HIJETS::tower_prefix));
        towerjetreco->add_input(new TowerJetInput(Jet::HCALIN_TOWERINFO, HIJETS::tower_prefix));
        towerjetreco->add_input(new TowerJetInput(Jet::HCALOUT_TOWERINFO, HIJETS::tower_prefix));
        towerjetreco->add_algo(new FastJetAlgoSub(Jet::ANTIKT, 0.2), "AntiKt_TowerInfo_HIRecoSeedsRaw_r02");
        towerjetreco->set_algo_node("ANTIKT");
        towerjetreco->set_input_node("TOWER");
        towerjetreco->Verbosity(0);
        if (WANT_VERBOSE) {
            std::cout << "[INFO] Registering towerjetreco subsystem..." << std::endl;
        }
        se->registerSubsystem(towerjetreco);
        
        if (WANT_VERBOSE) {
            std::cout << "[DEBUG] Creating DetermineTowerBackground (dtb)..." << std::endl;
        }
        
        DetermineTowerBackground *dtb = new DetermineTowerBackground();
        dtb->SetBackgroundOutputName("TowerInfoBackground_Sub1");
        dtb->SetFlow(HIJETS::do_flow);
        dtb->SetSeedType(0);
        dtb->SetSeedJetD(3);
        dtb->set_towerinfo(true);
        dtb->Verbosity(0);
        dtb->set_towerNodePrefix(HIJETS::tower_prefix);
        if (WANT_VERBOSE) {
            std::cout << "[INFO] Registering dtb subsystem..." << std::endl;
        }

        se->registerSubsystem(dtb);
        
        CopyAndSubtractJets *casj = new CopyAndSubtractJets();
        casj->SetFlowModulation(HIJETS::do_flow);
        casj->Verbosity(0);
        casj->set_towerinfo(true);
        casj->set_towerNodePrefix(HIJETS::tower_prefix);
        if (WANT_VERBOSE) {
            std::cout << "[INFO] Registering CopyAndSubtractJets subsystem (casj)..." << std::endl;
        }

        se->registerSubsystem(casj);
        
        DetermineTowerBackground *dtb2 = new DetermineTowerBackground();
        dtb2->SetBackgroundOutputName("TowerInfoBackground_Sub2");
        dtb2->SetFlow(HIJETS::do_flow);
        dtb2->SetSeedType(1);
        dtb2->SetSeedJetPt(7);
        dtb2->Verbosity(0);
        dtb2->set_towerinfo(true);
        dtb2->set_towerNodePrefix(HIJETS::tower_prefix);
        if (WANT_VERBOSE) {
            std::cout << "[INFO] Registering dtb2 subsystem..." << std::endl;
        }
        se->registerSubsystem(dtb2);
        
        SubtractTowers *st = new SubtractTowers();
        st->SetFlowModulation(HIJETS::do_flow);
        st->Verbosity(0);
        st->set_towerinfo(true);
        st->set_towerNodePrefix(HIJETS::tower_prefix);
        if (WANT_VERBOSE) {
            std::cout << "[INFO] Registering SubtractTowers subsystem (st)..." << std::endl;
        }

        se->registerSubsystem(st);
        
        
        //  ClusterIso(const std::string&, float eTCut, int coneSize, bool do_subtracted, bool do_unsubtracted);
        ClusterIso *makeClusterEt = new ClusterIso("CaloTreeGen", 0, 3, 1, 1);
        makeClusterEt->Verbosity(0);
        
        if (WANT_VERBOSE) {
            std::cout << "[INFO] ClusterIso subsystem created with name 'CaloTreeGen'"
            << " (only for data mode)..." << std::endl;
        }
        se->registerSubsystem(makeClusterEt);
    } //END IF RUNDATA
    
    
    // Finally, register caloTreeGen
    if (WANT_VERBOSE) {
        std::cout << "[DEBUG] Creating caloTreeGen with Outfile name: " << inName << std::endl;
    }

    caloTreeGen *eval = new caloTreeGen(inName);
    
    if (WANT_VERBOSE) {
        std::cout << "[INFO] Registering caloTreeGen subsystem...\n";
    }
    eval->setVerbose(WANT_VERBOSE);
    
    /*
     SET RUNNING OF DATA AND SIM from inputted condor submit arguments
     */
    eval->setWantSim(runSim);
    eval->setWantData(runData);
    if (runData)
    {
        if (WANT_VERBOSE) {
            std::cout << "[DEBUG] Still using run: " << runnumber << std::endl;
        }
        eval->setRunNumber(runnumber);
    }
    if (runSim)
    {
        if (WANT_VERBOSE) {
            std::cout << "[SIM] Using single sim file => " << firstFilename << std::endl;
        }

        eval->setSimInputFileName(firstFilename);
    }
    se->registerSubsystem(eval);
    
    
    // If wantData == true
    if (runData)
    {
        if (WANT_VERBOSE) {
            std::cout << "[DEBUG] Creating TriggerRunInfoReco subsystem...\n";
        }

        TriggerRunInfoReco *triggerruninforeco = new TriggerRunInfoReco();
        se->registerSubsystem(triggerruninforeco);
        
        if (WANT_VERBOSE) {
            std::cout << "[DEBUG] Creating Fun4AllDstInputManager (DSTcalo)...\n";
        }

        Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTcalo");
        
        if (WANT_VERBOSE) {
            // Reset the file stream to read all filenames from the beginning
            std::cout << "[DEBUG] Rewinding input file list to line 0...\n";
        }

        infile.clear();
        infile.seekg(0, std::ios::beg);
        
        // Read all filenames and add them to the input manager
        std::string filename;
        while (std::getline(infile, filename))
        {
            if (filename.empty()) continue;
            in->AddFile(filename.c_str());
            if (WANT_VERBOSE) {
                std::cout << "[INFO] Added input file: " << filename << std::endl;
            }

        }
        infile.close();
        se->registerInputManager(in);
    }
    
    if (!WANT_MANUAL_EVENT_LIMIT)
    {
      if (WANT_VERBOSE)
      {
        std::cout << "[DEBUG] calling se->run(" << nEvents << ")...\n";
      }
        if (runSim)
        {
          // If user said nEvents>0 => process that many
          // else => keep calling run(1) until ABORTRUN
          long long eventCount = 0;
          if (nEvents > 0)
          {
            while (eventCount < nEvents)
            {
              int retval = se->run(1);
              if (retval != 0) break;
              eventCount++;
            }
          }
          else
          {
            // nEvents=0 => process all TTree
            while (true)
            {
              int retval = se->run(1);
              if (retval != 0) break; // ABORTRUN => TTree done
            }
          }
        }
        else // runData
        {
          se->run(nEvents);
        }
        se->End();
    }
    else
    {
      // The manual event loop approach:
      int eventCount = 0;
      while (true)
      {
        int retval = se->run(1);
        if (retval != 0)
        {
          // e.g. ABORTRUN => break
          if (WANT_VERBOSE)
          {
            std::cout << "[INFO] run(1) returned " << retval
                      << " => stopping.\n";
          }
          break;
        }
        eventCount++;
        if (eventCount >= MANUAL_EVENT_LIMIT)
        {
          if (WANT_VERBOSE) {
            std::cout << "[INFO] Reached manual limit => stopping.\n";
          }
          break;
        }
      }
        se->End();
    }

    // Once we exit that block, we do the final End()
    if (WANT_VERBOSE)
    {
      std::cout << "[DEBUG] event loop finished => calling se->End()...\n";
    }


    // Now done => exit
    if (WANT_VERBOSE)
    {
      std::cout << "[DEBUG] Done => exiting.\n";
    }
    gSystem->Exit(0);
}
