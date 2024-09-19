#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <calotreegen/caloTreeGen.h>
#include <ClusterIso.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffarawobjects.so)
R__LOAD_LIBRARY(libcaloTreeGen.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libclusteriso.so)


#endif

void Fun4All_CaloTreeGen(const int nEvents = 0, const char *listFile = "/sphenix/lustre01/sphnxpro/physics/slurp/caloy2calib/ana430_2024p007/run_00047100_00047200/DST_CALO_run2pp_ana430_2024p007-00047114-00456.root", const char *inName = "commissioning.root")
{
    std::cout << "working ... " << std::endl;
  Fun4AllServer *se = Fun4AllServer::instance();
  recoConsts *rc = recoConsts::instance();

  gSystem->Load("libg4dst");
    
  std::cout << "Creating caloTreeGen module..." << std::endl;

  caloTreeGen *calo = new caloTreeGen(inName);

  //What subsystems do you want?
  calo -> doEMCal(1,"TOWERINFO_CALIB_CEMC");
  //Store HCal information?
  calo -> doHCals(1,"TOWERINFO_CALIB_HCALOUT","TOWERINFO_CALIB_HCALIN");
  //Store EMCal clusters?
  calo -> doClusters(1,"CLUSTERINFO_CEMC");
  //Store tower information for each EMCal cluster?
  calo -> doClusterDetails(1);
  //Store ZDC information?
  calo -> doZDC(1,"TOWERINFO_CALIB_ZDC");
    
  
  se->registerSubsystem(calo);


  Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTcalo");

  
  in->AddFile(listFile);    // condor

  se->registerInputManager(in);

  se->run(nEvents);
  se->End();
  se->PrintTimer();
  std::cout << "All done!" << std::endl;

  gSystem->Exit(0);
}
