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
#include "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/clusterIsoCopy_src/ClusterIso.h"

#include <calotreegen/caloTreeGen.h>

#include <Calo_Calib.C>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffarawobjects.so)
R__LOAD_LIBRARY(libcaloTreeGen.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libJetValidation.so)
R__LOAD_LIBRARY(libjetbase.so)
R__LOAD_LIBRARY(/sphenix/user/patsfan753/install/lib/libclusteriso.so)


//R__LOAD_LIBRARY(libclusteriso.so)

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

void Fun4All_CaloTreeGen(const int nEvents = 0, const char *listFile = "DST_CALO_run2pp_ana430_2024p007-00049070-00057.root", const char *inName = "commissioning.root") {
    std::cout << "[INFO] Starting Fun4All_CaloTreeGen..." << std::endl;

    Fun4AllServer *se = Fun4AllServer::instance();
    if (Enable::HIJETS_MC && Enable::HIJETS_TRUTH)
      {
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
    
    recoConsts *rc = recoConsts::instance();

    rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
    
    pair<int, int> runseg = Fun4AllUtils::GetRunSegment(listFile);
    int runnumber = runseg.first;
    if (runnumber <= 0) {
        std::cerr << "[ERROR] Invalid run number extracted: " << runnumber << ". Exiting..." << std::endl;
        return;
    }
    rc->set_uint64Flag("TIMESTAMP", runnumber);
    
    std::cout << "status setters" << std::endl;
    CaloTowerStatus *statusEMC = new CaloTowerStatus("CEMCSTATUS");
    statusEMC->set_detector_type(CaloTowerDefs::CEMC);
    statusEMC->set_time_cut(1);
    statusEMC->set_inputNodePrefix("TOWERINFO_CALIB_");
    statusEMC->Verbosity(Fun4AllBase::VERBOSITY_MORE);
    std::cout << "[INFO] Registering statusEMC subsystem..." << std::endl;
    se->registerSubsystem(statusEMC);
    
    CaloTowerStatus *statusHCalIn = new CaloTowerStatus("HCALINSTATUS");
    statusHCalIn->set_detector_type(CaloTowerDefs::HCALIN);
    statusHCalIn->set_time_cut(2);
    statusHCalIn->set_inputNodePrefix("TOWERINFO_CALIB_");
    statusHCalIn->Verbosity(Fun4AllBase::VERBOSITY_MORE);
    std::cout << "[INFO] Registering towerjetreco subsystem..." << std::endl;
    se->registerSubsystem(statusHCalIn);

    CaloTowerStatus *statusHCALOUT = new CaloTowerStatus("HCALOUTSTATUS");
    statusHCALOUT->set_detector_type(CaloTowerDefs::HCALOUT);
    statusHCALOUT->set_time_cut(2);
    statusHCALOUT->set_inputNodePrefix("TOWERINFO_CALIB_");
    statusHCALOUT->Verbosity(Fun4AllBase::VERBOSITY_MORE);
    se->registerSubsystem(statusHCALOUT);
    
    RetowerCEMC *rcemc = new RetowerCEMC();
    rcemc->Verbosity(0);
    rcemc->set_towerinfo(true);
    rcemc->set_frac_cut(0.5); //fraction of retower that must be masked to mask the full retower
    rcemc->set_towerNodePrefix(HIJETS::tower_prefix);
    std::cout << "[INFO] Registering RetowerCEMC subsystem..." << std::endl;
    se->registerSubsystem(rcemc);
    /*
     Relevent code from Calo_Calib.C since production p007 do not contain the TOWERS_CEMCnode that the Process_Calo_Calib() call needs
     */
    RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
    ClusterBuilder->Detector("CEMC");
    ClusterBuilder->set_threshold_energy(0.030);  // for when using basic calibration
    std::string emc_prof = getenv("CALIBRATIONROOT");
    emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
    ClusterBuilder->LoadProfile(emc_prof);
    ClusterBuilder->set_UseTowerInfo(1);  // to use towerinfo objects rather than old RawTower
    std::cout << "[INFO] Registering ClusterBuilder subsystem..." << std::endl;
    se->registerSubsystem(ClusterBuilder);

    JetReco *towerjetreco = new JetReco();
    towerjetreco->add_input(new TowerJetInput(Jet::CEMC_TOWERINFO_RETOWER,HIJETS::tower_prefix));
    towerjetreco->add_input(new TowerJetInput(Jet::HCALIN_TOWERINFO,HIJETS::tower_prefix));
    towerjetreco->add_input(new TowerJetInput(Jet::HCALOUT_TOWERINFO,HIJETS::tower_prefix));
    towerjetreco->add_algo(new FastJetAlgoSub(Jet::ANTIKT, 0.2), "AntiKt_TowerInfo_HIRecoSeedsRaw_r02");
    towerjetreco->set_algo_node("ANTIKT");
    towerjetreco->set_input_node("TOWER");
    towerjetreco->Verbosity(0);
    se->registerSubsystem(towerjetreco);

    DetermineTowerBackground *dtb = new DetermineTowerBackground();
    dtb->SetBackgroundOutputName("TowerInfoBackground_Sub1");
    dtb->SetFlow(HIJETS::do_flow);
    dtb->SetSeedType(0);
    dtb->SetSeedJetD(3);
    dtb->set_towerinfo(true);
    dtb->Verbosity(0);
    dtb->set_towerNodePrefix(HIJETS::tower_prefix);
    se->registerSubsystem(dtb);

    CopyAndSubtractJets *casj = new CopyAndSubtractJets();
    casj->SetFlowModulation(HIJETS::do_flow);
    casj->Verbosity(0);
    casj->set_towerinfo(true);
    casj->set_towerNodePrefix(HIJETS::tower_prefix);
    se->registerSubsystem(casj);

    DetermineTowerBackground *dtb2 = new DetermineTowerBackground();
    dtb2->SetBackgroundOutputName("TowerInfoBackground_Sub2");
    dtb2->SetFlow(HIJETS::do_flow);
    dtb2->SetSeedType(1);
    dtb2->SetSeedJetPt(7);
    dtb2->Verbosity(0);
    dtb2->set_towerinfo(true);
    dtb2->set_towerNodePrefix(HIJETS::tower_prefix);
    se->registerSubsystem(dtb2);

    SubtractTowers *st = new SubtractTowers();
    st->SetFlowModulation(HIJETS::do_flow);
    st->Verbosity(0);
    st->set_towerinfo(true);
    st->set_towerNodePrefix(HIJETS::tower_prefix);
    se->registerSubsystem(st);
    
    //  ClusterIso(const std::string&, float eTCut, int coneSize, bool do_subtracted, bool do_unsubtracted);
    ClusterIso *makeClusterEt = new ClusterIso("CaloTreeGen", 0, 3, 1, 1);
    makeClusterEt->Verbosity(0);
    se->registerSubsystem(makeClusterEt);
    std::cout << "[INFO] ClusterIso subsystem created and registered successfully." << std::endl;
    
    caloTreeGen *eval = new caloTreeGen(inName);
    eval->setTriggerNameMapForRun(runnumber);
    se -> registerSubsystem(eval);
    
    Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTcalo");
    in->AddFile(listFile);    // condor
    se->registerInputManager(in);
    
    se->run(nEvents);
    se->End();
    std::cout << "All done!" << std::endl;

    gSystem->Exit(0);
}
