#include <calotrigger/TriggerDefs.h>
#include <calobase/TowerInfoDefs.h>
#include "dlUtility.h"
#include <TMath.h>
#include <sstream>
#include <fstream>
#include <string>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h> 

R__LOAD_LIBRARY(libcalotrigger_io.so)
R__LOAD_LIBRARY(libcalo_io.so)

void get_scaledowns(int runnumber, int scaledowns[])
{

  TSQLServer *db = TSQLServer::Connect("pgsql://sphnxdaqdbreplica:5432/daq","phnxro","");
 
  if (db)
    {
      printf("Server info: %s\n", db->ServerInfo());
    }
  else
    {
      printf("bad\n");
    }


  TSQLRow *row;
  TSQLResult *res;
  TString cmd = "";
  char sql[1000];
  
  
  for (int is = 0; is < 64; is++)
    {
      sprintf(sql, "select scaledown%02d from gl1_scaledown where runnumber = %d;", is, runnumber);
      printf("%s \n" , sql);
 
      res = db->Query(sql);
      
      int nrows = res->GetRowCount();

      int nfields = res->GetFieldCount();
      for (int i = 0; i < nrows; i++) {
	row = res->Next();
	for (int j = 0; j < nfields; j++) {
	  scaledowns[is] = stoi(row->GetField(j));
	}
	delete row;
      }
 
      delete res;
    }
  delete db;
}
// Declare two pointers to TF1 objects
TF1 *f1, *f2;

// Define a function 'finter' that calculates the absolute difference between two functions at a given point 'x'
double finter(double *x, double*par) {
  // Evaluate the functions f1 and f2 at the point 'x' with parameters 'par'
  // Calculate the absolute difference between f1 and f2 at 'x'
  return TMath::Abs(f1->EvalPar(x, par) - f2->EvalPar(x, par));
}

// Main function 'fint' to create, draw, and find the intersection of two functions
void fint() {
  // Create the first function f1 with the expression "1 + 2*x + 0.2*x*x" in the range [0, 10]
  f1 = new TF1("f1", "1 + 2*x + 0.2*x*x", 0, 10);
  
  // Create the second function f2 with the expression "6 + 3*x - 0.3*x*x" in the range [0, 10]
  f2 = new TF1("f2", "6 + 3*x - 0.3*x*x", 0, 10);
  
  // Draw the first function f1
  f1->Draw();
  
  // Draw the second function f2 on the same canvas
  f2->Draw("same");
  
  // Create a new function fint that represents the absolute difference between f1 and f2
  // The function 'finter' is defined above and used here
  TF1 *fint = new TF1("fint", finter, 0, 10, 0);
  
  // Find the x-coordinate where the absolute difference fint is minimum
  double xint = fint->GetMinimumX();
  
  // Draw the function fint with the option "lsame" to overlay it on the existing canvas
  fint->Draw("lsame");
  
  // Create a marker at the intersection point with x-coordinate 'xint' and y-coordinate 'f1->Eval(xint)'
  TMarker *m = new TMarker(xint, f1->Eval(xint), 24);
  
  // Set the color of the marker to red
  m->SetMarkerColor(kRed);
  
  // Set the size of the marker to 3
  m->SetMarkerSize(3);
  
  // Draw the marker on the canvas
  m->Draw();
  
  // Print the x-coordinate of the intersection point
  printf("xint=%g\n", xint);
}

// Function to draw events based on the given run number, trigger bit, number of events, and segment
void Draw_Events(const int runnumber, const int triggerbit, const int nevents = 1, const int segment = 0) {
    
    
  std::string savename = "event_disp_run_" + to_string(runnumber) + "_" + to_string(segment) + "_gl1_" + to_string(triggerbit) + ".pdf";
    
    
  std::string savename2 = "event_disp_ll1_run_" + to_string(runnumber) + "_" + to_string(segment) + "_gl1_" + to_string(triggerbit) + ".pdf";

  std::string filename = Form("/sphenix/tg/tg01/commissioning/CaloCalibWG/dlis/TREE_CALO_TRIGGER_EMU-%08d-%04d.root", runnumber, segment);  
  TFile *f = new TFile(filename.c_str(),"r");
  
  TTree *t = (TTree*) f->Get("ttree");
  if (!t) {
      std::cout << "no tree" << std::endl;
      return;
    }
    
    
  // Initialize a variable to hold the type of trigger
  int triggertype = 0;
    
  // Determine the trigger type based on the trigger bit
  if (triggerbit >= 16 && triggerbit < 24)
    {
      triggertype = 1;
    }
  else if (triggerbit >= 24 && triggerbit < 32)
    {
      triggertype = 2;
    }
    
  // Initialize an array to hold scaledown values for 64 trigger bits, with all values initially set to -1
  int scaledowns[64]={-1};
    
  // Variables to hold scaledown data read from the file
  int a, b;
  std::string scaledownfile = "scaledown_" + to_string(runnumber) +".txt";
    
  // Open the scaledown file for reading
  std::ifstream infile(scaledownfile.c_str());
    
  // Read the scaledown values from the file
  while (infile >> a >> b) {
      // Check if the scaledown index is valid
      if (a < 0 || a >= 64) {
	  std::cout << "bad scaledowns " << std::endl;
	  return;
	}
      // Store the scaledown value in the array
      scaledowns[a] = b;      // process pair (a,b)
    }

  ULong64_t b_gl1_scaledvec; // Variable to hold the GL1 scaled vector
  ULong64_t b_gl1_clock;  // Variable to hold the GL1 clock
 
  // Define pointers to vectors to hold EMCal data
  std::vector<short>* b_emcal_good= nullptr;
  std::vector<float>* b_emcal_energy= nullptr;
  std::vector<float>* b_emcal_time= nullptr;
  std::vector<float>* b_emcal_etabin= nullptr;
  std::vector<float>* b_emcal_phibin= nullptr;

    
  // Define pointers to vectors to hold HCal inner data
  std::vector<short>* b_hcalin_good= nullptr;
  std::vector<float>* b_hcalin_energy= nullptr;
  std::vector<float>* b_hcalin_time= nullptr;
  std::vector<float>* b_hcalin_etabin= nullptr;
  std::vector<float>* b_hcalin_phibin= nullptr;

  // Define pointers to vectors to hold HCal outer data
  std::vector<short>* b_hcalout_good= nullptr;
  std::vector<float>* b_hcalout_energy= nullptr;
  std::vector<float>* b_hcalout_time= nullptr;
  std::vector<float>* b_hcalout_etabin= nullptr;
  std::vector<float>* b_hcalout_phibin= nullptr;

  // Define variables to hold cluster data
  int b_cluster_n;
  std::vector<float>* b_cluster_prob = nullptr;
  std::vector<float>* b_cluster_chi2 = nullptr;
  std::vector<float>* b_cluster_ecore = nullptr;
  std::vector<float>* b_cluster_pt = nullptr;
  std::vector<float>* b_cluster_phi = nullptr;
  std::vector<float>* b_cluster_eta = nullptr;
  std::vector<float>* b_cluster_iso = nullptr;

  // Define arrays to hold trigger sums for EMCal LL1 and EMCal
  unsigned int b_trigger_sum_emcal_ll1[384];
  unsigned int b_trigger_sumkey_emcal_ll1[384];
  unsigned int b_trigger_sum_emcal[6144];
  unsigned int b_trigger_sumkey_emcal[6144];

  // Define pointers to vectors to hold triggered sums for photon and jet
  std::vector<unsigned int>* b_triggered_sums_photon = nullptr;
  std::vector<unsigned int>* b_triggered_sums_jet= nullptr;

  // Set the branch addresses for the TTree branches to the corresponding variables
  t->SetBranchAddress("trigger_sum_emcal_ll1",&b_trigger_sum_emcal_ll1);
  t->SetBranchAddress("trigger_sumkey_emcal_ll1",&b_trigger_sumkey_emcal_ll1);
  t->SetBranchAddress("trigger_sum_emcal",&b_trigger_sum_emcal);
  t->SetBranchAddress("trigger_sumkey_emcal",&b_trigger_sumkey_emcal);

  t->SetBranchAddress("triggered_sums_photon",&b_triggered_sums_photon);
  t->SetBranchAddress("triggered_sums_jet",&b_triggered_sums_jet);
   
  t->SetBranchAddress("gl1_clock",&b_gl1_clock);//gl1_clock/l");
  t->SetBranchAddress("gl1_scaledvec",&b_gl1_scaledvec);//gl1_scaledvec/l");

  t->SetBranchAddress("emcal_good",&b_emcal_good);
  t->SetBranchAddress("emcal_energy",&b_emcal_energy);
  t->SetBranchAddress("emcal_time",&b_emcal_time);
  t->SetBranchAddress("emcal_phibin",&b_emcal_phibin);
  t->SetBranchAddress("emcal_etabin",&b_emcal_etabin);
  t->SetBranchAddress("hcalin_good",&b_hcalin_good);
  t->SetBranchAddress("hcalin_energy",&b_hcalin_energy);
  t->SetBranchAddress("hcalin_time",&b_hcalin_time);
  t->SetBranchAddress("hcalin_phibin",&b_hcalin_phibin);
  t->SetBranchAddress("hcalin_etabin",&b_hcalin_etabin);
  t->SetBranchAddress("hcalout_good",&b_hcalout_good);
  t->SetBranchAddress("hcalout_energy",&b_hcalout_energy);
  t->SetBranchAddress("hcalout_time",&b_hcalout_time);
  t->SetBranchAddress("hcalout_phibin",&b_hcalout_phibin);
  t->SetBranchAddress("hcalout_etabin",&b_hcalout_etabin);

  // Define vectors of 2D histograms for storing event data for various detector components
  std::vector<TH2D> h_events_emcal;           // Vector to store histograms for EMCal events
  std::vector<TH2D> h_events_emcal_retwr2;    // Vector to store histograms for EMCal retower 2x2 events
  std::vector<TH2D> h_events_emcal_retwr4;    // Vector to store histograms for EMCal retower 4x4 events
  std::vector<TH2D> h_events_emcal_retwr8;    // Vector to store histograms for EMCal retower 8x8 events
  std::vector<TH2D> h_events_emcal_ll1;       // Vector to store histograms for EMCal LL1 events
  std::vector<TH2D> h_events_emcal_2x2;       // Vector to store histograms for EMCal 2x2 events
  std::vector<TH2D> h_events_hcalin;          // Vector to store histograms for inner HCal events
  std::vector<TH2D> h_events_hcalout;         // Vector to store histograms for outer HCal events
  std::vector<std::vector<unsigned int>> h_events_ll1;  // Vector of vectors to store LL1 events data
  std::vector<int> event_numbers;             // Vector to store event numbers

  // Get the number of entries (events) in the TTree
  int entries = t->GetEntries();
    
  std::cout << " running through " << entries << " events" << std::endl;
    
  // Loop through each event in the TTree
  for (int i = 0; i < entries; i++)
    {
      t->GetEntry(i);
      bool skip = true; // Initialize a flag to determine whether to skip the current event
        
      // Vector to store trigger bits
      std::vector<int> trig_bits{};
        
      // Loop through all possible trigger bits (0 to 63)
      for (unsigned int bit = 0; bit < 64; bit++) {
          /*
          This line checks if the bit at position 'bit' in 'b_gl1_scaledvec' is set (i.e., is 1).
          1. 'b_gl1_scaledvec >> bit' shifts the bits of 'b_gl1_scaledvec' to the right by 'bit' positions.
             This moves the bit at the 'bit' position to the least significant bit (rightmost position).
          2. The result is then bitwise ANDed with '0x1U' (which is 0001 in binary).
             This masks out all bits except the least significant bit.
          3. The result of the AND operation is compared to '0x1U'.
             If it equals '0x1U', it means the original 'bit'-th bit in 'b_gl1_scaledvec' was 1 (set).
          */
          if (((b_gl1_scaledvec >> bit ) & 0x1U ) == 0x1U) {
              
          // If the current bit matches the trigger bit, do not skip this event
	      if (bit == triggerbit) skip = false;
	    }
	}
      // If the event does not contain the desired trigger bit, skip it
      if (skip) continue;
        
      // Create 2D histograms for the current event
      TH2D h_emcal_2x2(Form("h_emcal_2x2_%d", i), ";#eta EMCAL;#phi EMCAL", 48, -0.5, 47.5, 128, -0.5, 127.5);
      TH2D h_emcal_ll1(Form("h_emcal_ll1_%d", i), ";#eta EMCAL;#phi EMCAL", 12, -0.5, 11.5, 32, -0.5, 31.5);
      TH2D h_emcal(Form("h_emcal_%d", i), ";#eta EMCAL;#phi EMCAL", 96, -0.5, 95.5, 256, -0.5, 255.5);
      TH2D h_emcal_retwr2(Form("h_emcal_retwr2_%d", i), ";#eta EMCAL;#phi EMCAL", 48, -0.5, 95.5, 128, -0.5, 255.5);
      TH2D h_emcal_retwr4(Form("h_emcal_retwr4_%d", i), ";#eta EMCAL;#phi EMCAL", 24, -0.5, 95.5, 64, -0.5, 255.5);
      TH2D h_emcal_retwr8(Form("h_emcal_retwr8_%d", i), ";#eta EMCAL;#phi EMCAL", 12, -0.5, 95.5, 32, -0.5, 255.5);
      TH2D h_hcalin(Form("h_hcalin_%d", i), ";#eta HCALIN;#phi HCALIN", 24, -0.5, 23.5, 64, -0.5, 63.5);
      TH2D h_hcalout(Form("h_hcalout_%d", i), ";#eta HCALOUT;#phi HCALOUT", 24, -0.5, 23.5, 64, -0.5, 63.5);

      // Fill the EMCal LL1 histogram with data
      for (int ie = 0; ie < 384; ie++) {
          unsigned int eta = TriggerDefs::getSumEtaId(b_trigger_sumkey_emcal_ll1[ie]);
          unsigned int phi = TriggerDefs::getSumPhiId(b_trigger_sumkey_emcal_ll1[ie]) + 2*  TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(b_trigger_sumkey_emcal_ll1[ie]);
          h_emcal_ll1.Fill(eta, phi, b_trigger_sum_emcal_ll1[ie]);
	  }
      for (int ie = 0; ie < 6144; ie++) {
          unsigned int eta = TriggerDefs::getSumEtaId(b_trigger_sumkey_emcal[ie]) + 4*TriggerDefs::getPrimitiveEtaId_from_TriggerSumKey(b_trigger_sumkey_emcal[ie]);
          unsigned int phi = TriggerDefs::getSumPhiId(b_trigger_sumkey_emcal[ie]) + 4*TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(b_trigger_sumkey_emcal[ie]);
          h_emcal_2x2.Fill(eta, phi, b_trigger_sum_emcal[ie]);
      }

      for (int ie = 0; ie < b_emcal_energy->size(); ie++) {
          if (b_emcal_energy->at(ie) < 0.01) continue;
          //	  if (!b_emcal_good->at(ie)) continue;
          h_emcal.Fill(b_emcal_etabin->at(ie), b_emcal_phibin->at(ie), b_emcal_energy->at(ie));
          h_emcal_retwr2.Fill(b_emcal_etabin->at(ie), b_emcal_phibin->at(ie), b_emcal_energy->at(ie));
          h_emcal_retwr4.Fill(b_emcal_etabin->at(ie), b_emcal_phibin->at(ie), b_emcal_energy->at(ie));
          h_emcal_retwr8.Fill(b_emcal_etabin->at(ie), b_emcal_phibin->at(ie), b_emcal_energy->at(ie));
	}
      for (int ie = 0; ie < b_hcalin_energy->size(); ie++) {
	  if (!b_hcalin_good->at(ie)) continue;
	  if (b_hcalin_energy->at(ie) < 0.01) continue;
	  h_hcalin.SetBinContent(b_hcalin_etabin->at(ie), b_hcalin_phibin->at(ie), b_hcalin_energy->at(ie));
	}
      for (int ie = 0; ie < b_hcalout_energy->size(); ie++)
	{
	  if (!b_hcalout_good->at(ie)) continue;
	  if (b_hcalout_energy->at(ie) < 0.01) continue;
	  h_hcalout.SetBinContent(b_hcalout_etabin->at(ie), b_hcalout_phibin->at(ie), b_hcalout_energy->at(ie));
	}

      std::cout << "FOUND ONE" << std::endl;

        

      if (triggertype == 1) 
	{ 
	  for (auto &j : *b_triggered_sums_jet)
	    v.push_back(j);
	}
      else if (triggertype == 2) 
	{ 
	  for (auto &j : *b_triggered_sums_photon)
	    {
	      std::cout <<j << std::endl;
	      v.push_back(j);
	    }
	}

      h_events_ll1.push_back(v);
      if (h_events_emcal.size() == nevents) break;      
    }

  TString rundir = Form("/sphenix/user/dlis/Projects/CaloTriggerEmulator/plotted/run%d/", runnumber);
  gSystem->mkdir(rundir);

  TCanvas *canvas2 = new TCanvas("c2", "c2", 1400, 700);
  canvas2->cd();
  TPad *pll1 = new TPad("pll1","pll1", 0, 0, 0.25, 1.);
  pll1->SetLeftMargin(0.1);
  pll1->SetRightMargin(0.2);
  pll1->Draw();

  canvas2->cd();
  TPad *p2x2 = new TPad("p2x2","p2x2", 0.25,  0, 0.5, 1.);
  p2x2->SetLeftMargin(0.1);
  p2x2->SetRightMargin(0.2);
  p2x2->Draw();
  canvas2->cd();

  TPad *prll1 = new TPad("prll1","prll1", 0.5, 0, 0.75, 1.);
  prll1->SetLeftMargin(0.1);
  prll1->SetRightMargin(0.2);
  prll1->Draw();

  canvas2->cd();
  TPad *pr2x2 = new TPad("pr2x2","pr2x2", 0.75,  0, 1.0, 1.);
  pr2x2->SetLeftMargin(0.1);
  pr2x2->SetRightMargin(0.2);
  pr2x2->Draw();
  canvas2->cd();

  gStyle->SetOptStat(0);
  gStyle->SetPalette(kCool);
  gPad->SetTicks(1,1);
  for (int i = 0; i < h_events_emcal.size(); i++)
    {
      pll1->cd();
      h_events_emcal_ll1.at(i).Draw("colz");
      p2x2->cd();
      h_events_emcal_2x2.at(i).Draw("colz");
      prll1->cd();
      h_events_emcal_retwr8.at(i).Draw("colz");
      pr2x2->cd();
      h_events_emcal_retwr2.at(i).Draw("colz");

      if (h_events_emcal.size() == 1)
	{
	  canvas2->Print(Form("%s", savename2.c_str()),"pdf");
	}
      else if (i == 0)
	{
	  canvas2->Print(Form("%s(", savename2.c_str()),"pdf");
	}
      else if (i == h_events_emcal.size() - 1)
	{
	  canvas2->Print(Form("%s)", savename2.c_str()),"pdf");
	}
      else
	{
	  canvas2->Print(savename2.c_str(),"pdf");
	}

    }

  TCanvas *canvas = new TCanvas("c", "c", 1000, 700);
  canvas->cd();
  TPad *pe = new TPad("pe","pe", 0, 0, 0.33, 1.);
  pe->SetLeftMargin(0.1);
  pe->SetRightMargin(0.2);
  pe->Draw();

  canvas->cd();
  TPad *phi = new TPad("phi","phi", 0.33, 0, 0.67, 1.);
  phi->SetLeftMargin(0.1);
  phi->SetRightMargin(0.2);
  phi->Draw();
  canvas->cd();
  TPad *pho = new TPad("pho","pho", 0.67, 0, 1., 1.);
  pho->SetLeftMargin(0.1);
  pho->SetRightMargin(0.2);
  pho->Draw();
  canvas->cd();

  gStyle->SetOptStat(0);
  gStyle->SetPalette(kCool);
  gPad->SetTicks(1,1);


  for (int i = 0; i < h_events_emcal.size(); i++)
    {
      pe->cd();

      if (triggerbit >= 24)
	{
	  h_events_emcal_retwr4.at(i).GetZaxis()->SetRangeUser(0, 6);
	}
      else
	{
	  h_events_emcal_retwr4.at(i).GetZaxis()->SetRangeUser(0, 10);
	}

      h_events_emcal_retwr4.at(i).Draw("colz");

      if (triggertype == 1)
	{
	  for (auto &j : h_events_ll1.at(i))
	    {
	      unsigned int phi = j & 0xffffU;
	      unsigned int eta = (j>>16U) & 0xffffU;
	      unsigned int topphi = phi + 4;
	      if (phi > 28)
		{
		  topphi = 32;
		}
	      DrawBox(eta*8 - 0.5, phi*8 - 0.5, (eta + 4)*8 - 0.5, topphi*8 - 0.5);
	      if (phi > 28)
		{
		  DrawBox(eta*8 - 0.5, 0.5, (eta + 4)*8 - 0.5, (phi - 28)*8 - 0.5);
		}	      
	    }
	}
      else
	{
	  for (auto &j : h_events_ll1.at(i))
	    {

	      
	      unsigned int phi = TriggerDefs::getSumPhiId(j) + TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(j)*2;
	      unsigned int eta = TriggerDefs::getSumEtaId(j);
	      std::cout << eta <<" "<<phi<< " -- " << eta*8 << " " << phi*8 << std::endl;
	      DrawBox(eta*8 - 0.5, phi*8 - 0.5, (eta + 1)*8 - 0.5, (phi + 1)*8 - 0.5);
	    }
	}
	  	  
      drawText(Form("Run %d - Event %d", runnumber, event_numbers.at(i)), 0.1, 0.95);
      phi->cd();
      h_events_hcalin.at(i).GetZaxis()->SetRangeUser(0, 10);
      h_events_hcalin.at(i).Draw("colz");
      if (triggertype == 1)
	{
	  for (auto &j : h_events_ll1.at(i))
	    {
	      unsigned int phi = j & 0xffffU;
	      unsigned int eta = (j>>16U) & 0xffffU;
	      unsigned int topphi = phi + 4;
	      if (phi > 28)
		{
		  topphi = 32;
		}
	      TBox *hbox = new TBox(eta*2 - 0.5, phi*2 - 0.5, (eta + 4)*2 - 0.5, topphi*2 - 0.5);

	      hbox->SetFillStyle(4000);
	      hbox->SetLineWidth(2);
	      hbox->SetLineColor(2);
	      hbox->Draw("same");

	      if (phi > 28)
		{
		  TBox *hbox2 = new TBox(eta*2 - 0.5,  -0.5, (eta+4)*2  - 0.5, (phi - 28)*2 - 0.5);

		  hbox2->SetFillStyle(4000);
		  hbox2->SetLineWidth(2);
		  hbox2->SetLineColor(2);
		  hbox2->Draw("same");
	      
		}	      
	    }
	}
	  
      pho->cd();
      h_events_hcalout.at(i).GetZaxis()->SetRangeUser(0, 10);
      h_events_hcalout.at(i).Draw("colz");
      if (triggertype == 1)
	{
	  for (auto &j : h_events_ll1.at(i))
	    {
	      unsigned int phi = j & 0xffffU;
	      unsigned int eta = (j>>16U) & 0xffffU;
	      unsigned int topphi = phi+4;
	      if (phi > 28)
		{
		  topphi = 32;
		}
	      TBox *hbox = new TBox(eta*2 - 0.5, phi*2 - 0.5, (eta + 4)*2 - 0.5, topphi*2 - 0.5);

	      hbox->SetFillStyle(4000);
	      hbox->SetLineWidth(2);
	      hbox->SetLineColor(2);
	      hbox->Draw("same");

	      if (phi > 28)
		{
		  TBox *hbox2 = new TBox(eta*2 - 0.5,  -0.5, (eta+4)*2  - 0.5, (phi - 28)*2 - 0.5);

		  hbox2->SetFillStyle(4000);
		  hbox2->SetLineWidth(2);
		  hbox2->SetLineColor(2);
		  hbox2->Draw("same");
	      
		}	      
	    }
	}
	  

      if (h_events_emcal.size() == 1)
	{
	  canvas->Print(Form("%s", savename.c_str()),"pdf");
	}
      else if (i == 0)
	{
	  canvas->Print(Form("%s(", savename.c_str()),"pdf");
	}
      else if (i == h_events_emcal.size() - 1)
	{
	  canvas->Print(Form("%s)", savename.c_str()),"pdf");
	}
      else
	{
	  canvas->Print(savename.c_str(),"pdf");
	}
    }

  return;
}
void Draw_Trigger_TurnOn(const int runnumber, std::string trigger)
{
  int threshold = 3 ;  
  const std::string filename = "/sphenix/tg/tg01/commissioning/CaloCalibWG/dlis/HIST_" + trigger + "_TRIGGER_QA-000" + to_string(runnumber) + ".root";
  TFile *f = new TFile(filename.c_str(),"r");

  TH2D *h_lut_energy = (TH2D*) f->Get("h_emcal_8x8_energy_lutsum");
  h_lut_energy->GetYaxis()->SetRangeUser(0, 16);
  // We want to get the slices of this histogram 
  const int nbinsy = h_lut_energy->GetYaxis()->GetNbins();
  const int nbinsx = h_lut_energy->GetXaxis()->GetNbins();
  
  TF1 *divider = new TF1("divider","[0] + [1]*x");
  divider->FixParameter(0, 3.);
  divider->FixParameter(1, (7.-3.)/6.);
  
  TH1D *h_lut_energy_proj = (TH1D*) h_lut_energy->ProjectionX();
  TH2D *h_lut_energy_normalized = (TH2D*) h_lut_energy->Clone();
  h_lut_energy_normalized->SetName("h_lut_energy_normalized");
  h_lut_energy_normalized->Reset();

  for (int i = 0; i < nbinsx; i++)
    {
      if (i < threshold)
	{
	  for (int j = 0; j < nbinsy; j++)
	    {
	      int bin = h_lut_energy_normalized->GetBin(i+1, j+1);
	      h_lut_energy_normalized->SetBinContent(bin, 0);
	    }
	  continue;
	}
      
      for (int j = 0; j < nbinsy; j++)
	{
	  int bin = h_lut_energy_normalized->GetBin(i+1, j+1);
	  // float x = h_lut_energy_normalized->GetXaxis()->GetBinCenter(j+1);
	  // float y = h_lut_energy_normalized->GetYaxis()->GetBinCenter(i+1);
	  // if (divider->Eval(x) < y) continue; 
	  h_lut_energy_normalized->SetBinContent(bin, h_lut_energy->GetBinContent(bin));
	  h_lut_energy_normalized->SetBinError(bin, h_lut_energy->GetBinError(bin));
	}
    }
  
  TH1D *hproj = (TH1D*)  h_lut_energy_normalized->ProjectionX();
  for (int i = 0; i < nbinsx; i++)
    {

      if (i == 0)
	{
	  for (int j = 0; j < nbinsy; j++)
	    {
	      int bin = h_lut_energy_normalized->GetBin(i+1, j+1);
	      h_lut_energy_normalized->SetBinContent(bin, 0);
	    }
	  continue;
	}

      float scale = hproj->GetBinContent(i+1);

      if (scale == 0) continue;
      for (int j = 0; j < nbinsy; j++)
	{
	  int bin = h_lut_energy_normalized->GetBin(i+1, j+1);
	  h_lut_energy_normalized->SetBinContent(bin, h_lut_energy_normalized->GetBinContent(bin)/scale);
	  h_lut_energy_normalized->SetBinError(bin, h_lut_energy_normalized->GetBinError(bin)/scale);
	}

    }

  //Get The Step Function;

  TH1D *h_turn_on = (TH1D*) h_lut_energy_normalized->ProjectionY();
  h_turn_on->SetName("h_turn_on");

  TH1D* h_turn_on_norm = (TH1D*) h_turn_on->Clone();
  h_turn_on_norm->SetName("h_turn_on_norm");
  
  TF1 *flatline = new TF1("flatline", "[0]", 5, 8);
  h_turn_on->Fit("flatline","R");
  float level = flatline->GetParameter(0);
  h_turn_on_norm->Scale(1./level);
  // TF1 *trigeffcurve = new TF1("trigeffcurve","1.0-TMath::Exp(-pow((x/[0]), [1]))",0.0,10.);
  // trigeffcurve->SetParameters(0.732,0.532); // just initial value guesses
  // //  trigeffcurve->SetParLimits(0,0.2,20.0);    
  // //trigeffcurve->SetParLimits(1,0.2,5.0);
  // h_turn_on_norm->Fit(trigeffcurve,"R","",1.0,7.0);


  float plateau = 0;
  for (int i = 0; i < h_turn_on_norm->GetNbinsX(); i++)
    {
      if (h_turn_on_norm->GetBinContent(i+1) > 0.95)
	{
	  plateau = h_turn_on_norm->GetBinCenter(i+1);
	  break;
	}
    }
  gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("c","c", 500, 500);
  c->SetTicks(1,1);
  h_turn_on_norm->SetTitle("; Energy [GeV]; #epsilon(E)");
  h_turn_on_norm->SetLineColor(kBlue);
  h_turn_on_norm->SetLineWidth(2);

  h_turn_on_norm->SetMaximum(1.5);
  h_turn_on_norm->GetXaxis()->SetRangeUser(0, 10);
  h_turn_on_norm->Draw();

  DrawSPHENIXraw(0.16, 0.83);
  drawText("pp - #sqrt{s} = 200 GeV", 0.16, 0.78, 0, kBlack, 0.03); 
  drawText(Form("Run %d", runnumber), 0.16, 0.74, 0, kBlack, 0.03); 
  drawText(Form("%s %d", trigger.c_str(), threshold), 0.16, 0.70, 0, kBlack, 0.03); 
  drawText(Form("Plateau ~ %2.1f GeV", plateau), 0.16, 0.66, 0, kBlack, 0.03); 
  TLine *l = new TLine(0, 1, 10, 1);
  l->SetLineWidth(1);
  l->SetLineColor(kBlack);
  l->SetLineStyle(2);
  l->Draw("same");
  TCanvas *c2 = new TCanvas("c2","c2", 500, 500);
  h_turn_on->Draw("hist");

  TCanvas *c3 = new TCanvas("c3","c3", 500, 500);
  c3->SetLogz();
  h_lut_energy_normalized->GetXaxis()->SetRangeUser(0, 20);
  h_lut_energy_normalized->Draw("colz");

}
void Draw_Trigger_Basic(const std::string filename)
{

  ofstream outf;
  outf.open ("example.txt");

  TFile *fin = new TFile("/sphenix/user/dlis/Projects/macros/CDBTest/optmask_051624.root","r");

  TNtuple *tn = (TNtuple*) fin->Get("tn_optmask");

  std::vector<uint32_t> v_optmask{};
  float mask = 0;
  tn->SetBranchAddress("primkey", &mask);
  for (int i = 0 ; i < tn->GetEntries(); i++)
    {
      tn->GetEntry(i);
      //      v_optmask.push_back(static_cast<uint32_t>((int)mask));
    }
  
  TFile *f = new TFile(filename.c_str(),"r");
  
  TTree *t = (TTree*) f->Get("ttree");


  unsigned int b_trigger_sum_smpl_emcal[6144];
  unsigned int b_trigger_sumkey_emcal[6144];
  unsigned int b_trigger_sum_emcal[6144];  
  unsigned int b_trigger_sum_smpl_hcalin[384];
  unsigned int b_trigger_sumkey_hcalin[384];
  unsigned int b_trigger_sum_hcalin[384];
  unsigned int b_trigger_sum_smpl_hcalout[384];
  unsigned int b_trigger_sumkey_hcalout[384];
  unsigned int b_trigger_sum_hcalout[384];

  std::vector<uint32_t> b_triggered_sums;
  ULong64_t b_gl1_rawvec;
  ULong64_t b_gl1_livevec;
  ULong64_t b_gl1_scaledvec;
  ULong64_t b_gl1_clock;
 
  unsigned int b_trigger_sum_smpl_emcal_ll1[384];
  unsigned int b_trigger_sumkey_emcal_ll1[384];
  unsigned int b_trigger_sum_emcal_ll1[384];
  unsigned int b_trigger_sum_smpl_hcal_ll1[384];
  unsigned int b_trigger_sumkey_hcal_ll1[384];
  unsigned int b_trigger_sum_hcal_ll1[384];

  unsigned int b_trigger_sum_smpl_jet[288];
  unsigned int b_trigger_sumkey_jet[288];
  unsigned int b_trigger_sum_jet[288];

  unsigned int b_trigger_sum_smpl_jet_input[384];
  unsigned int b_trigger_sumkey_jet_input[384];
  unsigned int b_trigger_sum_jet_input[384];

  std::vector<unsigned int>* b_trigger_bits = nullptr;;
  std::vector<unsigned int>* b_trigger_raw_bits= nullptr;

  std::vector<float>* b_emcal_energy= nullptr;
  std::vector<float>* b_emcal_time= nullptr;
  std::vector<float>* b_emcal_etabin= nullptr;
  std::vector<float>* b_emcal_phibin= nullptr;

  std::vector<float>* b_hcalin_energy= nullptr;
  std::vector<float>* b_hcalin_time= nullptr;
  std::vector<float>* b_hcalin_etabin= nullptr;
  std::vector<float>* b_hcalin_phibin= nullptr;

  std::vector<float>* b_hcalout_energy= nullptr;
  std::vector<float>* b_hcalout_time= nullptr;
  std::vector<float>* b_hcalout_etabin= nullptr;
  std::vector<float>* b_hcalout_phibin= nullptr;

  t->SetBranchAddress("trigger_sum_emcal",b_trigger_sum_emcal);
  t->SetBranchAddress("trigger_sumkey_emcal",b_trigger_sumkey_emcal);
  t->SetBranchAddress("trigger_sum_smpl_emcal",b_trigger_sum_smpl_emcal);
  t->SetBranchAddress("trigger_sum_hcalin",b_trigger_sum_hcalin);
  t->SetBranchAddress("trigger_sumkey_hcalin",b_trigger_sumkey_hcalin);//trigger_sumkey_hcalin[384]/i");
  t->SetBranchAddress("trigger_sum_smpl_hcalin",b_trigger_sum_smpl_hcalin);//trigger_sum_smpl_hcalin[384]/i");
  t->SetBranchAddress("trigger_sum_hcalout",b_trigger_sum_hcalout);//trigger_sum_hcalout[384]/i");
  t->SetBranchAddress("trigger_sumkey_hcalout",b_trigger_sumkey_hcalout);//trigger_sumkey_hcalout[384]/i");
  t->SetBranchAddress("trigger_sum_smpl_hcalout",b_trigger_sum_smpl_hcalout);//trigger_sum_smpl_hcalout[384]/i");
  

  t->SetBranchAddress("trigger_sum_emcal_ll1",b_trigger_sum_emcal_ll1);//trigger_sum_emcal_ll1[384]/i");
  t->SetBranchAddress("trigger_sumkey_emcal_ll1",b_trigger_sumkey_emcal_ll1);//trigger_sumkey_emcal_ll1[384]/i");
  t->SetBranchAddress("trigger_sum_smpl_emcal_ll1",b_trigger_sum_smpl_emcal_ll1);//trigger_sum_smpl_emcal_ll1[384]/i");
  t->SetBranchAddress("trigger_sum_hcal_ll1",b_trigger_sum_hcal_ll1);//trigger_sum_hcal_ll1[384]/i");
  t->SetBranchAddress("trigger_sumkey_hcal_ll1",b_trigger_sumkey_hcal_ll1);//trigger_sumkey_hcal_ll1[384]/i");
  t->SetBranchAddress("trigger_sum_smpl_hcal_ll1",b_trigger_sum_smpl_hcal_ll1);//trigger_sum_smpl_hcal_ll1[384]/i");

  t->SetBranchAddress("triggered_sums",&b_triggered_sums);//triggered_sums/i");

  t->SetBranchAddress("gl1_clock",&b_gl1_clock);//gl1_clock/l");
  t->SetBranchAddress("gl1_rawvec",&b_gl1_rawvec);//gl1_triggervec/l");
  t->SetBranchAddress("gl1_livevec",&b_gl1_livevec);//gl1_livevec/l");
  t->SetBranchAddress("gl1_scaledvec",&b_gl1_scaledvec);//gl1_scaledvec/l");

  t->SetBranchAddress("trigger_bits",&b_trigger_bits);
  t->SetBranchAddress("trigger_raw_bits",&b_trigger_raw_bits);

  t->SetBranchAddress("trigger_sum_jet",b_trigger_sum_jet);//trigger_sum_jet[288]/i");
  t->SetBranchAddress("trigger_sum_smpl_jet",b_trigger_sum_smpl_jet);//trigger_sum_smpl_jet[288]/i");
  t->SetBranchAddress("trigger_sumkey_jet",b_trigger_sumkey_jet);//trigger_sumkey_jet[288]/i");

  t->SetBranchAddress("trigger_sum_smpl_jet_input",b_trigger_sum_smpl_jet_input);//trigger_sum_smpl_jet_input[384]/i");
  t->SetBranchAddress("trigger_sum_jet_input",b_trigger_sum_jet_input);//trigger_sum_jet_input[384]/i");
  t->SetBranchAddress("trigger_sumkey_jet_input",b_trigger_sumkey_jet_input);//trigger_sumkey_jet_input[384]/i");

  t->SetBranchAddress("emcal_energy",&b_emcal_energy);
  t->SetBranchAddress("emcal_time",&b_emcal_time);
  t->SetBranchAddress("emcal_phibin",&b_emcal_phibin);
  t->SetBranchAddress("emcal_etabin",&b_emcal_etabin);
  t->SetBranchAddress("hcalin_energy",&b_hcalin_energy);
  t->SetBranchAddress("hcalin_time",&b_hcalin_time);
  t->SetBranchAddress("hcalin_phibin",&b_hcalin_phibin);
  t->SetBranchAddress("hcalin_etabin",&b_hcalin_etabin);
  t->SetBranchAddress("hcalout_energy",&b_hcalout_energy);
  t->SetBranchAddress("hcalout_time",&b_hcalout_time);
  t->SetBranchAddress("hcalout_phibin",&b_hcalout_phibin);
  t->SetBranchAddress("hcalout_etabin",&b_hcalout_etabin);

  // Historgram Time
  float energymap[12][32] = {0};
  int clock_diff_bins = 1000;
  float clock_diff_max = 100000;
  int n_primitives = 384;
  int n_sums = 6144;
  int energy_bins = 100;
  float energy_max = 20.0;
  int turnon_bins = 40;
  float turnon_max = 10.0;
  int lut_bins = 256;


  // clock v sum

  TH2D *h_clock_sum = new TH2D("h_clock_sum", ";Clock Diff;Sum", clock_diff_bins, 0, clock_diff_max, n_sums, -0.5, n_sums - 0.5);
  // clock v primitive
  TH2D *h_clock_primitive = new TH2D("h_clock_primitive", ";Clock Diff; Primitive",clock_diff_bins, 0, clock_diff_max, n_primitives, -0.5, n_primitives - 0.5);

  //channel by channel scatter plots;
  TH2D *h_prim_energy[384];
  for (int i = 0; i < 384; i++)
    {
      h_prim_energy[i] = new TH2D(Form("h_prim_energy_%d", i), "; Energy; Lut", energy_bins, 0, energy_max, lut_bins, -0.5, lut_bins - 0.5);
    }

  TH1D *h_sum_energy_slopes = new TH1D("h_sum_energy_slopes", "; Sum/Energy Slopes", energy_bins, 0, energy_max);
  TH1D *h_emcal_lut = new TH1D("h_emcal_lut", "; Lut;", lut_bins, 0, lut_bins - 0.5);
  TH1D *h_emcal_energy = new TH1D("h_emcal_energy", "; Energy;", 200, -0.5,  19.5);
  TH1D *h_emcal_energy_patch = new TH1D("h_emcal_energy_patch", "; Energy;", 200, -0.5,  19.5);
  TH1D *h_emcal_time = new TH1D("h_emcal_time", "; Time;", 12, -0.5,  11.5);

  TH1D *h_emcal_lut_gl1[64];// = new TH1D("h_emcal_lut", "; Lut;", lut_bins, 0, lut_bins - 0.5);
  TH1D *h_emcal_energy_gl1[64];// = new TH1D("h_emcal_energy", "; Energy;", 200, -0.5,  19.5);
  TH1D *h_emcal_energy_raw_gl1[64];// = new TH1D("h_emcal_energy", "; Energy;", 200, -0.5,  19.5);
  TH1D *h_emcal_energy_patch_gl1[64];// = new TH1D("h_emcal_energy_patch", "; Energy;", 200, -0.5,  19.5);
  TH1D *h_emcal_time_gl1[64];// = new TH1D("h_emcal_time", "; Time;", 12, -0.5,  11.5);

  for (int i = 0; i < 64; i++)
    {
      h_emcal_lut_gl1[i] = new TH1D(Form("h_emcal_lut_%d",i), "; Lut;", lut_bins, 0, lut_bins - 0.5);
      h_emcal_energy_gl1[i] = new TH1D(Form("h_emcal_energy_%d",i), "; Energy;", 200, -0.5,  19.5);
      h_emcal_energy_raw_gl1[i] = new TH1D(Form("h_emcal_energy_raw_%d",i), "; Energy;", 200, -0.5,  19.5);
      h_emcal_energy_patch_gl1[i] = new TH1D(Form("h_emcal_energy_patch_%d",i), "; Energy;", 200, -0.5,  19.5);
      h_emcal_time_gl1[i] = new TH1D(Form("h_emcal_time_%d",i), "; Time;", 12, -0.5,  11.5);
    }

  TH2D *h_emcal_lut_energy = new TH2D("h_emcal_lut_energy", "; Energy;Lut;", 200, -0.5, 19.5, lut_bins, 0, lut_bins - 0.5);
  //Turn on curve
  TEfficiency *he_turn_on_curves[20];
  for (int i = 0; i < 20; i++)
    {
      he_turn_on_curves[i] = new TEfficiency(Form("he_turn_on_curve_%d", i), "; Energy [GeV]; Efficiency", turnon_bins, 0, turnon_max); 
    }
  TGraph * g_clock_primitive[64];
  TGraph * g_clock_primitive_above_10 [64];
  for (int i = 0 ; i < 64; i++)
    {
      g_clock_primitive[i] = new TGraph();
      g_clock_primitive[i]->SetName(Form("g_clock_primitive_%d", i));
      g_clock_primitive_above_10[i] = new TGraph();
      g_clock_primitive_above_10[i]->SetName(Form("g_clock_primitive_above_10%d", i));

    }

  
  float e_cut = 0.03;
  float e_cut_patch = 0.1;
  int entries = t->GetEntries();
  std::cout << " running through " << entries << " events" << std::endl;
  for (int i = 0; i < entries; i++)
    {
      t->GetEntry(i);
      if (entries < 10)
	{
	  std::cout << "event " << i << std::endl;
	}
      if (i % 100 == 0)
	{
	  std::cout << "\rEvent "<<std::dec << i << " ...   "<< std::hex << b_gl1_scaledvec <<std::flush;
	}
      std::vector<int> trig_bits{};
      for (uint16_t bit = 0; bit < 64; bit++)
	{
	  if (((b_gl1_scaledvec >> bit) & 0x1U) == 0x1U)
	    {
	      trig_bits.push_back(bit);
	    }
	}

      for (int j = 0; j < 32; j++)
	{
	  for (int k =0 ; k < 12; k++)
	    {
	      energymap[k][j] = 0.0;
	    }
	}      
      for (int ie = 0; ie < b_emcal_energy->size(); ie++)
	{
	  int ebin = b_emcal_etabin->at(ie)/8;
	  int pbin = b_emcal_phibin->at(ie)/8;
	  if (b_emcal_time->at(ie) < 0) continue;	  
	  energymap[ebin][pbin] += b_emcal_energy->at(ie);
	  h_emcal_energy->Fill(b_emcal_energy->at(ie));
	  if (b_emcal_energy->at(ie) > e_cut)
	    h_emcal_time->Fill(b_emcal_time->at(ie));
	  for (auto &b : trig_bits)
	    {
	      h_emcal_energy_gl1[b]->Fill(b_emcal_energy->at(ie));
	      if (b_emcal_energy->at(ie) > e_cut)
		h_emcal_time_gl1[b]->Fill(b_emcal_time->at(ie));
	    }
	}

      for (int j = 0; j < 32; j++)
	{
	  if (entries <=10) std::cout << j << " \t |";

	  for (int k =0 ; k < 12; k++)
	    {

	      if (entries <= 10) std::cout << energymap[k][j] << "\t";
	      if (energymap[k][j] > e_cut_patch)
		h_emcal_energy_patch->Fill(energymap[k][j]);
	      for (auto &b : trig_bits)
		{
		  h_emcal_energy_patch_gl1[b]->Fill(energymap[k][j]);
		}
	    }

	  if (entries <=10) std::cout << " | " <<std::endl;;

	}      

      // look through sums
      for (int j = 0; j < 384; j++)
	{

	  // Now look and see if sumkey matches any prim mask;
	  uint32_t key = (b_trigger_sumkey_emcal_ll1[j] & 0xffffffffU);
	  if (std::find(v_optmask.begin(), v_optmask.end(), key) != v_optmask.end()) continue;

	  int phibin = TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(key) * 2 + TriggerDefs::getSumPhiId(key);
	  int etabin = TriggerDefs::getSumEtaId(key);

	  for (int it = 0; it < 20; it++)
	    {
	      he_turn_on_curves[it]->Fill(b_trigger_sum_emcal_ll1[j] > it, energymap[etabin][phibin]);
	    }

	  if (b_trigger_sum_emcal_ll1[j] == 0 && energymap[etabin][phibin] > 4)
	    {
	      outf << " Event " << i << " | primitive " << j << " | phi/eta " << phibin << "/" << etabin << " | Energy " << energymap[etabin][phibin] << " GeV " <<std::endl;
	    }
	  for (auto &b : trig_bits)
	    {
	      h_emcal_lut_gl1[b]->Fill(b_trigger_sum_emcal_ll1[j]);
	      if (b_trigger_sum_emcal_ll1[j] > 10)
		g_clock_primitive_above_10[b]->SetPoint(g_clock_primitive_above_10[b]->GetN(), (float)(b_gl1_clock & 0xffffffffff), j);
	      else if (b_trigger_sum_emcal_ll1[j] > 0)
		g_clock_primitive[b]->SetPoint(g_clock_primitive[b]->GetN(), (float)(b_gl1_clock & 0xffffffffff), j);
	    }
	  for (auto &b : trig_bits)
	    {
	      h_emcal_lut_gl1[b]->Fill(b_trigger_sum_emcal_ll1[j]);
	    }

	  h_emcal_lut_energy->Fill( energymap[etabin][phibin], b_trigger_sum_emcal_ll1[j]);
	  h_prim_energy[j]->Fill(energymap[etabin][phibin], b_trigger_sum_emcal_ll1[j]);
	  h_emcal_lut->Fill( b_trigger_sum_emcal_ll1[j]);
	  // Now see if the channel is masked.
	 	  
	}
      
    }

  
  TFile *fout = new TFile("outhist.root","recreate");
  for (int i = 0; i < 64; i++)
    {

      h_emcal_lut_gl1[i]->Write();
      h_emcal_energy_gl1[i]->Write();
      h_emcal_time_gl1[i]->Write();
      h_emcal_energy_patch_gl1[i]->Write();
      g_clock_primitive[i]->Write();
      g_clock_primitive_above_10[i]->Write();
    }
  for (int i = 0; i < 20; i++)
    {
      he_turn_on_curves[i]->Write();
    }
  h_emcal_lut->Write();
  h_emcal_energy->Write();
  h_emcal_time->Write();
  h_emcal_energy_patch->Write();
  fout->Close();
  outf.close();
  std::cout <<" Done." << std::endl;
  return;
}
void Draw_Trigger_Sim_Efficiency(int runnumber = -1, const std::string filename = "nothing", std::string trigger = "SIM")
{
  
  gStyle->SetOptStat(0);
  std::string outfilename = Form("/sphenix/user/dlis/Projects/CaloTriggerEmulator/plots/HIST_%s-%08d.root",trigger.c_str(), runnumber);
  std::string outfilenamenew = Form("/sphenix/user/dlis/Projects/CaloTriggerEmulator/plots/HIST_%sNEW-%08d.root",trigger.c_str(), runnumber);
  if (runnumber < 0)
    {
      outfilename = filename;
    }

  TFile *f = new TFile(outfilename.c_str(), "r");
  if (!f)
    {
      std::cout << "nope" << std::endl;
      return;
    }
  TH1D *h_pt_thresh[100];
  for (int i = 0; i < 100; i++)
    {
      h_pt_thresh[i] = (TH1D*) f->Get(Form("h_pt_thresh_%d", i));
    }
  TH1D *h_ntower_thresh[100][5];
  for (int i = 0; i < 100; i++)
    {
      for (int j = 0; j < 5; j++)
	{
	  h_ntower_thresh[i][j] = (TH1D*) f->Get(Form("h_ntower_thresh_%d_%d", i, j));
	}
    }

  TH1D *h_pt_thresh_ratio[100];
  for (int i = 0; i < 100; i++)
    {
      h_pt_thresh_ratio[i] = (TH1D*) h_pt_thresh[i]->Clone();
      //      h_pt_thresh_ratio[i]->SetName(Form("ratio_%d"), i);
      h_pt_thresh_ratio[i]->Divide(h_pt_thresh[0]);
    }

  TH1D *h_ntower_thresh_ratio[100][5];
  for (int i = 0; i < 100; i++)
    {
      for (int j = 0; j < 5; j++)
	{
	  h_ntower_thresh_ratio[i][j] = (TH1D*) h_ntower_thresh[i][j]->Clone();
	  //      h_pt_thresh_ratio[i]->SetName(Form("ratio_%d"), i);
	  h_ntower_thresh_ratio[i][j]->Divide(h_ntower_thresh[0][j]);
	}
    }


  TFile *fn = new TFile(outfilenamenew.c_str(), "r");
  if (!fn)
    {
      std::cout << "nope" << std::endl;
      return;
    }
  TH1D *h_pt_thresh_n[100];
  for (int i = 0; i < 100; i++)
    {
      h_pt_thresh_n[i] = (TH1D*) fn->Get(Form("h_pt_thresh_%d", i));
    }
  TH1D *h_ntower_thresh_n[100][5];
  for (int i = 0; i < 100; i++)
    {
      for (int j = 0; j < 5; j++)
	{
	  h_ntower_thresh_n[i][j] = (TH1D*) fn->Get(Form("h_ntower_thresh_%d_%d", i, j));
	}
    }

  TH1D *h_pt_thresh_ratio_n[100];
  for (int i = 0; i < 100; i++)
    {
      h_pt_thresh_ratio_n[i] = (TH1D*) h_pt_thresh_n[i]->Clone();
      //      h_pt_thresh_ratio[i]->SetName(Form("ratio_%d"), i);
      h_pt_thresh_ratio_n[i]->Divide(h_pt_thresh_n[0]);
    }
  TH1D *h_ntower_thresh_ratio_n[100][5];
  for (int i = 0; i < 100; i++)
    {
      for (int j = 0; j < 5; j++)
	{
	  h_ntower_thresh_ratio_n[i][j] = (TH1D*) h_ntower_thresh_n[i][j]->Clone();
	  //      h_pt_thresh_ratio[i]->SetName(Form("ratio_%d"), i);
	  h_ntower_thresh_ratio_n[i][j]->Divide(h_ntower_thresh_n[0][j]);
	}
    }

  
  int colors[5] = {kBlack, kRed, kBlue, kGreen, kOrange};

  for (int i = 0; i < 5; i++)
    {
      //      SetLineAtt(h_pt_thresh_ratio[i*2+6], colors[i], 2, 1);
      SetMarkerAtt(h_pt_thresh_ratio[i*2+4], colors[i], 1, 4);
      SetMarkerAtt(h_pt_thresh_ratio_n[i*4+6], colors[i], 1, 8);
    }
  
  SetyjPadStyle();
  TCanvas *c1 = new TCanvas("c1","c1", 500, 500);
  h_pt_thresh_ratio[4]->SetMaximum(1.4);
  h_pt_thresh_ratio[4]->SetMinimum(0);
  h_pt_thresh_ratio[4]->SetTitle(";#gamma p_{T} [ GeV];Efficieny");
  h_pt_thresh_ratio[4]->Draw("p");
  h_pt_thresh_ratio_n[6]->Draw("same p");
  for (int i = 1 ; i < 5 ; i++)
    {
      h_pt_thresh_ratio[i*2 + 4]->Draw("same p");
      h_pt_thresh_ratio_n[i*4 + 6]->Draw("same p");
    } 
  drawText("#gamma Sim 5 - 25 GeV", 0.19, 0.85);

  TLegend *leg = new TLegend(0.75, 0.75, 0.9, 0.9);
  for (int i = 0 ; i < 5 ; i++)
    {
      leg->AddEntry(h_pt_thresh_ratio[i*2 + 4], Form("Old Threshold - %d", i*2+4));
    }
  leg->SetBorderSize(0);
  leg->SetFillStyle(4000);

  leg->Draw("same");

  leg = new TLegend(0.6, 0.75, 0.75, 0.9);
  for (int i = 0 ; i < 5 ; i++)
    {
      leg->AddEntry(h_pt_thresh_ratio_n[i*4 + 6], Form("New.2 Threshold - %d", i*4+6));
    }
  leg->SetBorderSize(0);
  leg->SetFillStyle(4000);

  leg->Draw("same");

  TLine *l = new TLine(0, 1, 25, 1);
  SetLineAtt(l, kBlack, 2, 2);

  
  SetyjPadStyle();
  TCanvas *c2 = new TCanvas("c2","c2", 500, 500);
  int i = 2;
  for (int j = 0; j < 5; j++)
    {
      //      SetLineAtt(h_ntower_thresh_ratio[i*2+6], colors[i], 2, 1);
	SetMarkerAtt(h_ntower_thresh_ratio[j*2+4][i-1], colors[j], 1, 4);
	SetMarkerAtt(h_ntower_thresh_ratio_n[j*4+6][i-1], colors[j], 1, 8);
    }
  h_ntower_thresh_ratio[4][i-1]->GetXaxis()->SetRangeUser(-0.5, 25.5);
  h_ntower_thresh_ratio[4][i-1]->SetMaximum(1.4);
  h_ntower_thresh_ratio[4][i-1]->SetMinimum(0);
  h_ntower_thresh_ratio[4][i-1]->SetTitle(";nTower > 100 MeV;Efficieny");
  h_ntower_thresh_ratio[4][i-1]->Draw("p");
  //h_ntower_thresh_ratio_n[6][i-1]->Draw("same p");
  for (int j = 1 ; j < 5 ; j++)
    {
      h_ntower_thresh_ratio[j*2 + 4][i-1]->Draw("same p");
      //h_ntower_thresh_ratio_n[j*4 + 6][i-1]->Draw("same p");
    } 
  drawText("#gamma Sim 5 - 25 GeV", 0.19, 0.85);
  drawText(Form("%d GeV <= pt < %d GeV", i*5, (i+1)*5), 0.19, 0.8);
  leg = new TLegend(0.65, 0.65, 0.9, 0.9);
  for (int j = 0 ; j < 5 ; j++)
    {
      leg->AddEntry(h_ntower_thresh_ratio[j*2 + 4][i-1], Form("Old Threshold - %d", j*2+4));
    }
  leg->SetBorderSize(0);
  leg->SetFillStyle(4000);
  
  leg->Draw("same");
}


void Draw_Trigger_Efficieny(int runnumber = -1, const std::string filename = "nothing")
{
  std::string trigger = "PHOTON";
  gStyle->SetOptStat(0);
  std::string outfilename = Form("/sphenix/user/dlis/Projects/CaloTriggerEmulator/plots/HIST_%s-%08d.root",trigger.c_str(), runnumber);
  if (runnumber < 0)
    {
      outfilename = filename;
    }

  TFile *f = new TFile(outfilename.c_str(), "r");
  if (!f)
    {
      std::cout << "nope" << std::endl;
      return;
    }
  TH1D *h_emcal_energy_gl1[64];
  TH1D *h_jet_energy_gl1[64];

  std::string histtrigname;
  for (int i = 0; i < 64; i++)
    {
      histtrigname = "h_emcal_energy_raw_gl1_" + to_string(i);
      h_emcal_energy_gl1[i] = (TH1D*) f->Get(histtrigname.c_str());
      histtrigname = "h_jet_energy_raw_gl1_" + to_string(i);
      h_jet_energy_gl1[i] = (TH1D*) f->Get(histtrigname.c_str());
          
    }
  TH1D *h_emcal_energy_ref = (TH1D*) f->Get("h_emcal_energy_ref");
  TH1D *h_jet_energy_ref = (TH1D*) f->Get("h_jet_energy_ref");

  int jet_bit = 18;
  int emcal_bit = 22;
  int mbd_bit = 10;
  TH1D *h_turn_on_emcal[4];
  TH1D *h_turn_on_jet[4];
  for (int i = 0; i < 4; i++)
    {
      h_turn_on_emcal[i] = (TH1D*) h_emcal_energy_gl1[emcal_bit + i]->Clone();
      h_turn_on_jet[i] = (TH1D*) h_jet_energy_gl1[jet_bit + i]->Clone();

      h_turn_on_emcal[i]->Divide(h_emcal_energy_ref);
      h_turn_on_jet[i]->Divide(h_jet_energy_ref);
    }

  TCanvas *c1 = new TCanvas("c2","c2", 500, 700);
  
  TPad *p12 = new TPad("p12","p12", 0.0, 0.4, 1.0, 1.0);
  p12->SetTicks(1,1);
  p12->Draw();	
  p12->cd();
  p12->SetLogy();  
  h_emcal_energy_ref->SetTitle(";Maximum 8x8 Energy Sum (EMCAL) [GeV]; Prescale Counts");
  SetLineAtt(h_emcal_energy_ref, kBlack, 2, 1);
  SetLineAtt(h_emcal_energy_gl1[emcal_bit + 0], kBlue, 2, 1);
  SetLineAtt(h_emcal_energy_gl1[emcal_bit + 1], kRed, 2, 1);
  SetLineAtt(h_emcal_energy_gl1[emcal_bit + 2], kSpring+2, 2, 1);
  SetLineAtt(h_emcal_energy_gl1[emcal_bit + 3], kBlue, 2, 1);
  h_emcal_energy_ref->SetMinimum(1);
  h_emcal_energy_ref->Draw("hist");
  //  h_emcal_energy_gl1[emcal_bit + 0]->Draw("hist same");
  h_emcal_energy_gl1[emcal_bit + 1]->Draw("hist same");
  h_emcal_energy_gl1[emcal_bit + 2]->Draw("hist same");
  h_emcal_energy_gl1[emcal_bit + 3]->Draw("hist same");
  TLegend *l = new TLegend(0.65, 0.4, 0.85, 0.7);
  l->SetBorderSize(0);  
  l->AddEntry(h_emcal_energy_ref, "MBD >= 1");
  //  l->AddEntry(h_emcal_energy_gl1[emcal_bit + 0], "Photon 1");
  l->AddEntry(h_emcal_energy_gl1[emcal_bit + 1], "Photon 2 GeV");
  l->AddEntry(h_emcal_energy_gl1[emcal_bit + 2], "Photon 3 GeV");
  l->AddEntry(h_emcal_energy_gl1[emcal_bit + 3], "Photon 4 GeV");
  l->Draw("same");

  float x = 0.6;
  float y = 0.8;
  float dy = 0.05;
  DrawSPHENIXraw(x, y);
  drawText("pp - #sqrt{s} = 200 GeV", x, y - dy, 0, kBlack, 0.03); 
  drawText(Form("Run %d", runnumber), x , y - 2*dy, 0, kBlack, 0.03); 

  c1->cd();
  TPad *p2 = new TPad("p2","p2", 0.0, 0.0, 1.0,0.4);
  p2->SetTicks(1,1);
  p2->Draw();
  p2->cd();
  std::string triggername = "PHOTON";
  h_turn_on_emcal[0]->SetTitle(Form(";Maximum 8x8 Energy Sum (EMCAL) [GeV]; MBD Coincidence/%s" , triggername.c_str()));
  SetLineAtt(h_turn_on_emcal[0], kBlue, 2, 1);
  SetLineAtt(h_turn_on_emcal[1], kRed, 2, 1);
  SetLineAtt(h_turn_on_emcal[2], kSpring + 2, 2, 1);
  SetLineAtt(h_turn_on_emcal[3], kBlue, 2, 1);
  h_turn_on_emcal[1]->GetXaxis()->SetRangeUser(0, 10);
  h_turn_on_emcal[1]->SetMaximum(2.0);
  h_turn_on_emcal[1]->SetMinimum(0);


  TF1 *f1 = new TF1("f1","[0]*(1 + TMath::Erf((x-[1])/[2]))",0, 7);
  f1->SetParameter(0, 0.3);
  f1->SetParameter(1, 4);
  f1->SetParameter(2, 1);
  //  h_turn_on_emcal[1]->Fit(f1, "RN");
  //  h_turn_on_emcal[0]->Draw();
  h_turn_on_emcal[1]->Draw("");
  h_turn_on_emcal[2]->Draw("same");
  h_turn_on_emcal[3]->Draw("same");
  // f1->SetLineColor(kBlack);
  // f1->Draw("same");

  // double plateau = f1->GetParameter(0)*2;
  // double half = 0;
  // double y1 = 0;
  // while (y1 < plateau/2.)
  //   {
  //     half += 0.01;
  //     y1 = f1->Eval(half);
  //     if (half == 20) break;
  //   }

  // drawText("Photon 1 :", 0.17, 0.8, 0, kBlack, 0.06);
  // drawText(Form("   plateau     : %1.2f", plateau), 0.17, 0.73, 0, kBlack, 0.06);
  // drawText(Form("   Half height : %1.2f GeV", half), 0.17, 0.66, 0, kBlack, 0.06);

  // Jet

  TCanvas *c = new TCanvas("c1","c1", 500, 700);

  TPad *p1 = new TPad("p1","p1", 0.0, 0.4, 1.0, 1.0);
  p1->SetTicks(1,1);
  p1->Draw();	
  p1->cd();
  p1->SetLogy();  
  h_jet_energy_ref->SetTitle(";Maximum 0.8x0.8 Energy Sum (EMCAL + HCAL) [GeV]; Prescale Counts");
  SetLineAtt(h_jet_energy_ref, kBlack, 2, 1);
  SetLineAtt(h_jet_energy_gl1[jet_bit + 0], kBlue, 2, 1);
  SetLineAtt(h_jet_energy_gl1[jet_bit + 1], kRed, 2, 1);
  SetLineAtt(h_jet_energy_gl1[jet_bit + 2], kSpring+2, 2, 1);
  SetLineAtt(h_jet_energy_gl1[jet_bit + 3], kBlue, 2, 1);
  h_jet_energy_ref->SetMinimum(1);
  h_jet_energy_ref->Draw("hist");
  drawText("06/12/2024", 0.75, 0.92);
  //  h_jet_energy_gl1[jet_bit + 0]->Draw("hist same");
  h_jet_energy_gl1[jet_bit + 1]->Draw("hist same");
  h_jet_energy_gl1[jet_bit + 2]->Draw("hist same");
  h_jet_energy_gl1[jet_bit + 3]->Draw("hist same");
  l = new TLegend(0.65, 0.4, 0.85, 0.7);
  l->SetBorderSize(0);  
  l->AddEntry(h_jet_energy_ref, "MBD >= 1");
  //  l->AddEntry(h_jet_energy_gl1[jet_bit + 0], "Jet 6 GeV");
  l->AddEntry(h_jet_energy_gl1[jet_bit + 1], "Jet 6 GeV");
  l->AddEntry(h_jet_energy_gl1[jet_bit + 2], "Jet 8 GeV");
  l->AddEntry(h_jet_energy_gl1[jet_bit + 3], "Jet 10 GeV");
  l->Draw("same");

  DrawSPHENIXraw(x, y);
  drawText("pp - #sqrt{s} = 200 GeV", x, y - dy, 0, kBlack, 0.03); 
  drawText(Form("Run %d", runnumber), x , y - 2*dy, 0, kBlack, 0.03); 

  c->cd();
  TPad *p22 = new TPad("p22","p22", 0.0, 0.0, 1.0,0.4);
  p22->SetTicks(1,1);
  p22->Draw();
  p22->cd();
  h_turn_on_jet[0]->SetTitle(";Maximum 0.8x0.8 Energy Sum (EMCAL + HCAL) [GeV]; MBD Coincidence/JET");
  SetLineAtt(h_turn_on_jet[0], kBlue, 2, 1);
  SetLineAtt(h_turn_on_jet[1], kRed, 2, 1);
  SetLineAtt(h_turn_on_jet[2], kSpring + 2, 2, 1);
  SetLineAtt(h_turn_on_jet[3], kBlue, 2, 1);
  
  h_turn_on_jet[1]->SetMaximum(1.4);
  h_turn_on_jet[1]->SetMinimum(0);

  TF1 *f2 = new TF1("f2","[0]*(1 + TMath::Erf((x-[1])/[2]))",0, 20);
  f2->SetParameter(0, 0.4);
  f2->SetParameter(1, 5);
  f2->SetParameter(2, 1);
  h_turn_on_jet[1]->GetXaxis()->SetRangeUser(0, 20);
  //  h_turn_on_jet[1]->Fit(f2, "RN");
  h_turn_on_jet[1]->Draw();
  //h_turn_on_jet[1]->Draw("same");
  h_turn_on_jet[2]->Draw("same");
  h_turn_on_jet[3]->Draw("same");
  // f2->SetLineColor(kBlack);
  // f2->Draw("same");

  // plateau = f2->GetParameter(0)*2;
  // half = 0;
  // y1 = 0;
  // while (y1 < plateau/2.)
  //   {
  //     half += 0.01;
  //     y1 = f2->Eval(half);
  //     if (half == 10) break;
  //   }

  // drawText("Jet 6 GeV :", 0.17, 0.8, 0, kBlack, 0.06);
  // drawText(Form("   plateau     : %1.2f", plateau), 0.17, 0.73, 0, kBlack, 0.06);
  // drawText(Form("   Half height : %1.2f GeV", half), 0.17, 0.66, 0, kBlack, 0.06);

  
}

void Draw_Trigger_Turn_On(int runnumber = -1, const std::string trigger = "PHOTON", int triggerbit = 24, const std::string filename = "nothing")
{

  gStyle->SetOptStat(0);
  std::string outfilename = Form("/sphenix/user/dlis/Projects/CaloTriggerEmulator/plots/HIST_%s-%08d.root",trigger.c_str(), runnumber);
  if (runnumber < 0)
    {
      outfilename = filename;
    }
  TFile *f = new TFile(outfilename.c_str(), "r");
  if (!f)
    {
      std::cout << "nope" << std::endl;
      return;
    }

  int scaledowns[64]={-1};
  int a, b;
  std::string scaledownfile = "scaledown_" + to_string(runnumber) +".txt";
  std::ifstream infile(scaledownfile.c_str());
  while (infile >> a >> b)
    {
      if (a < 0 || a >= 64)
	{
	  std::cout << "bad scaledowns " << std::endl;
	  return;
	}
      scaledowns[a] = b;      // process pair (a,b)
    }
  TH1D *h_emcal_energy_gl1[64];
  TH1D *h_jet_energy_gl1[64];
  std::string histtrigname;
  for (int i = 0; i < 64; i++)
    {
      histtrigname = "h_emcal_energy_scaled_gl1_" + to_string(i);
      h_emcal_energy_gl1[i] = (TH1D*) f->Get(histtrigname.c_str());
      histtrigname = "h_jet_energy_scaled_gl1_" + to_string(i);
      h_jet_energy_gl1[i] = (TH1D*) f->Get(histtrigname.c_str());
          
    }
  int jet_bit = 16;
  int emcal_bit = 24;
  int mbd_bit = 10;
  TH1D *h_turn_on_emcal[4];
  TH1D *h_turn_on_jet[4];
  for (int i = 0; i < 4; i++)
    {
      h_turn_on_emcal[i] = (TH1D*) h_emcal_energy_gl1[emcal_bit + i]->Clone();
      h_turn_on_jet[i] = (TH1D*) h_jet_energy_gl1[jet_bit + i]->Clone();

      h_turn_on_emcal[i]->Divide(h_emcal_energy_gl1[mbd_bit]);
      h_turn_on_jet[i]->Divide(h_jet_energy_gl1[mbd_bit]);
    }

  TCanvas *c1 = new TCanvas("c2","c2", 500, 700);
  
  TPad *p12 = new TPad("p12","p12", 0.0, 0.4, 1.0, 1.0);
  p12->SetTicks(1,1);
  p12->Draw();	
  p12->cd();
  p12->SetLogy();  
  h_emcal_energy_gl1[mbd_bit]->SetTitle(";Maximum 0.8x0.8 Energy Sum (EMCAL) [GeV]; Prescale Counts");
  SetLineAtt(h_emcal_energy_gl1[mbd_bit], kBlack, 2, 1);
  SetLineAtt(h_emcal_energy_gl1[emcal_bit + 0], kBlue, 2, 1);
  SetLineAtt(h_emcal_energy_gl1[emcal_bit + 1], kRed, 2, 1);
  SetLineAtt(h_emcal_energy_gl1[emcal_bit + 2], kSpring+2, 2, 1);
  SetLineAtt(h_emcal_energy_gl1[emcal_bit + 3], kBlue, 2, 1);
  h_emcal_energy_gl1[mbd_bit]->SetMinimum(1);
  h_emcal_energy_gl1[mbd_bit]->Draw("hist");
  //  h_emcal_energy_gl1[emcal_bit + 0]->Draw("hist same");
  h_emcal_energy_gl1[emcal_bit + 1]->Draw("hist same");
  h_emcal_energy_gl1[emcal_bit + 2]->Draw("hist same");
  h_emcal_energy_gl1[emcal_bit + 3]->Draw("hist same");
  TLegend *l = new TLegend(0.65, 0.4, 0.85, 0.7);
  l->SetBorderSize(0);  
  l->AddEntry(h_emcal_energy_gl1[mbd_bit], "MBD >= 1");
  //  l->AddEntry(h_emcal_energy_gl1[emcal_bit + 0], "Photon 1");
  l->AddEntry(h_emcal_energy_gl1[emcal_bit + 1], "Photon 2 GeV");
  l->AddEntry(h_emcal_energy_gl1[emcal_bit + 2], "Photon 3 GeV");
  l->AddEntry(h_emcal_energy_gl1[emcal_bit + 3], "Photon 4 GeV");
  l->Draw("same");

  float x = 0.6;
  float y = 0.8;
  float dy = 0.05;
  DrawSPHENIXraw(x, y);
  drawText("pp - #sqrt{s} = 200 GeV", x, y - dy, 0, kBlack, 0.03); 
  drawText(Form("Run %d", runnumber), x , y - 2*dy, 0, kBlack, 0.03); 

  c1->cd();
  TPad *p2 = new TPad("p2","p2", 0.0, 0.0, 1.0,0.4);
  p2->SetTicks(1,1);
  p2->Draw();
  p2->cd();
  std::string triggername = "PHOTON";
  h_turn_on_emcal[0]->SetTitle(Form(";Maximum 8x8 Energy Sum (EMCAL) [GeV]; MBD Coincidence/%s" , triggername.c_str()));
  SetLineAtt(h_turn_on_emcal[0], kBlue, 2, 1);
  SetLineAtt(h_turn_on_emcal[1], kRed, 2, 1);
  SetLineAtt(h_turn_on_emcal[2], kSpring + 2, 2, 1);
  SetLineAtt(h_turn_on_emcal[3], kBlue, 2, 1);
  h_turn_on_emcal[1]->GetXaxis()->SetRangeUser(0, 10);
  h_turn_on_emcal[1]->SetMaximum(2.0);
  h_turn_on_emcal[1]->SetMinimum(0);


  TF1 *f1 = new TF1("f1","[0]*(1 + TMath::Erf((x-[1])/[2]))",0, 7);
  f1->SetParameter(0, 0.3);
  f1->SetParameter(1, 4);
  f1->SetParameter(2, 1);
  h_turn_on_emcal[1]->Fit(f1, "RN");
  //  h_turn_on_emcal[0]->Draw();
  h_turn_on_emcal[1]->Draw("");
  h_turn_on_emcal[2]->Draw("same");
  h_turn_on_emcal[3]->Draw("same");
  f1->SetLineColor(kBlack);
  f1->Draw("same");

  double plateau = f1->GetParameter(0)*2;
  double half = 0;
  double y1 = 0;
  while (y1 < plateau/2.)
    {
      half += 0.01;
      y1 = f1->Eval(half);
      if (half == 20) break;
    }

  drawText("Photon 1 :", 0.17, 0.8, 0, kBlack, 0.06);
  drawText(Form("   plateau     : %1.2f", plateau), 0.17, 0.73, 0, kBlack, 0.06);
  drawText(Form("   Half height : %1.2f GeV", half), 0.17, 0.66, 0, kBlack, 0.06);

  // Jet

  TCanvas *c = new TCanvas("c1","c1", 500, 700);
  
  TPad *p1 = new TPad("p1","p1", 0.0, 0.4, 1.0, 1.0);
  p1->SetTicks(1,1);
  p1->Draw();	
  p1->cd();
  p1->SetLogy();  
  h_jet_energy_gl1[mbd_bit]->SetTitle(";Maximum 0.8x0.8 Energy Sum (EMCAL + HCAL) [GeV]; Prescale Counts");
  SetLineAtt(h_jet_energy_gl1[mbd_bit], kBlack, 2, 1);
  SetLineAtt(h_jet_energy_gl1[jet_bit + 0], kBlue, 2, 1);
  SetLineAtt(h_jet_energy_gl1[jet_bit + 1], kRed, 2, 1);
  SetLineAtt(h_jet_energy_gl1[jet_bit + 2], kSpring+2, 2, 1);
  SetLineAtt(h_jet_energy_gl1[jet_bit + 3], kBlue, 2, 1);
  h_jet_energy_gl1[mbd_bit]->SetMinimum(1);
  h_jet_energy_gl1[mbd_bit]->Draw("hist");
  //  h_jet_energy_gl1[jet_bit + 0]->Draw("hist same");
  h_jet_energy_gl1[jet_bit + 1]->Draw("hist same");
  h_jet_energy_gl1[jet_bit + 2]->Draw("hist same");
  h_jet_energy_gl1[jet_bit + 3]->Draw("hist same");
  l = new TLegend(0.65, 0.4, 0.85, 0.7);
  l->SetBorderSize(0);  
  l->AddEntry(h_jet_energy_gl1[mbd_bit], "MBD >= 1");
  //  l->AddEntry(h_jet_energy_gl1[jet_bit + 0], "Jet 6 GeV");
  l->AddEntry(h_jet_energy_gl1[jet_bit + 1], "Jet 6 GeV");
  l->AddEntry(h_jet_energy_gl1[jet_bit + 2], "Jet 8 GeV");
  l->AddEntry(h_jet_energy_gl1[jet_bit + 3], "Jet 10 GeV");
  l->Draw("same");

  DrawSPHENIXraw(x, y);
  drawText("pp - #sqrt{s} = 200 GeV", x, y - dy, 0, kBlack, 0.03); 
  drawText(Form("Run %d", runnumber), x , y - 2*dy, 0, kBlack, 0.03); 

  c->cd();
  TPad *p22 = new TPad("p22","p22", 0.0, 0.0, 1.0,0.4);
  p22->SetTicks(1,1);
  p22->Draw();
  p22->cd();
  h_turn_on_jet[0]->SetTitle(";Maximum 8x8 Energy Sum (EMCAL + HCAL) [GeV]; MBD Coincidence/JET");
  SetLineAtt(h_turn_on_jet[0], kBlue, 2, 1);
  SetLineAtt(h_turn_on_jet[1], kRed, 2, 1);
  SetLineAtt(h_turn_on_jet[2], kSpring + 2, 2, 1);
  SetLineAtt(h_turn_on_jet[3], kBlue, 2, 1);
  
  h_turn_on_jet[1]->SetMaximum(1.4);
  h_turn_on_jet[1]->SetMinimum(0);

  TF1 *f2 = new TF1("f2","[0]*(1 + TMath::Erf((x-[1])/[2]))",0, 20);
  f2->SetParameter(0, 0.4);
  f2->SetParameter(1, 5);
  f2->SetParameter(2, 1);
  h_turn_on_jet[1]->GetXaxis()->SetRangeUser(0, 20);
  h_turn_on_jet[1]->Fit(f2, "RN");
  h_turn_on_jet[1]->Draw();
  //h_turn_on_jet[1]->Draw("same");
  h_turn_on_jet[2]->Draw("same");
  h_turn_on_jet[3]->Draw("same");
  f2->SetLineColor(kBlack);
  f2->Draw("same");

  plateau = f2->GetParameter(0)*2;
  half = 0;
  y1 = 0;
  while (y1 < plateau/2.)
    {
      half += 0.01;
      y1 = f2->Eval(half);
      if (half == 10) break;
    }

  drawText("Jet 6 GeV :", 0.17, 0.8, 0, kBlack, 0.06);
  drawText(Form("   plateau     : %1.2f", plateau), 0.17, 0.73, 0, kBlack, 0.06);
  drawText(Form("   Half height : %1.2f GeV", half), 0.17, 0.66, 0, kBlack, 0.06);

}

void Draw_Trigger_Turn_On(const std::string outfilename, const std::string trigger = "PHOTON")
{
  std::string triggername = "PHOTON";
  gStyle->SetOptStat(0);
  TFile *f = new TFile(outfilename.c_str(), "r");
  if (!f)
    {
      std::cout << "nope" << std::endl;
      return;
    }

  TH1D *h_emcal_energy_gl1[64];
  TH1D *h_cluster_energy_gl1[64];
  TH1D *h_jet_energy_gl1[64];
  std::string histtrigname;
  for (int i = 0; i < 64; i++)
    {
      histtrigname = "h_emcal_energy_scaled_gl1_" + to_string(i);
      h_emcal_energy_gl1[i] = (TH1D*) f->Get(histtrigname.c_str());
      histtrigname = "h_cluster_energy_scaled_gl1_" + to_string(i);
      h_cluster_energy_gl1[i] = (TH1D*) f->Get(histtrigname.c_str());

      histtrigname = "h_jet_energy_scaled_gl1_" + to_string(i);
      h_jet_energy_gl1[i] = (TH1D*) f->Get(histtrigname.c_str());
          
    }
  int jet_bit = 16;
  int emcal_bit = 24;
  int mbd_bit = 10;
  TH1D *h_turn_on_emcal[4];
  TH1D *h_turn_on_jet[4];
  TH1D *h_turn_on_cluster[4];
  for (int i = 0; i < 4; i++)
    {
      h_turn_on_emcal[i] = (TH1D*) h_emcal_energy_gl1[emcal_bit + i]->Clone();
      h_turn_on_cluster[i] = (TH1D*) h_cluster_energy_gl1[emcal_bit + i]->Clone();
      h_turn_on_jet[i] = (TH1D*) h_jet_energy_gl1[jet_bit + i]->Clone();

      h_turn_on_cluster[i]->Divide(h_cluster_energy_gl1[mbd_bit]);
      h_turn_on_emcal[i]->Divide(h_emcal_energy_gl1[mbd_bit]);
      h_turn_on_jet[i]->Divide(h_jet_energy_gl1[mbd_bit]);
    }

  TCanvas *c1 = new TCanvas("c2","c2", 500, 700);
  
  TPad *p12 = new TPad("p12","p12", 0.0, 0.4, 1.0, 1.0);
  p12->SetTicks(1,1);
  p12->SetLeftMargin(0.2);
  p12->Draw();	
  p12->cd();
  p12->SetLogy();  
  h_emcal_energy_gl1[mbd_bit]->SetTitle(";Maximum 8x8 Energy Sum (EMCAL) [GeV]; Prescale Counts");
  SetLineAtt(h_emcal_energy_gl1[mbd_bit], kBlack, 2, 1);
  SetLineAtt(h_emcal_energy_gl1[emcal_bit + 0], kBlue, 2, 1);
  SetLineAtt(h_emcal_energy_gl1[emcal_bit + 1], kRed, 2, 1);
  SetLineAtt(h_emcal_energy_gl1[emcal_bit + 2], kSpring+2, 2, 1);
  SetLineAtt(h_emcal_energy_gl1[emcal_bit + 3], kBlue, 2, 1);
  h_emcal_energy_gl1[mbd_bit]->SetMinimum(1);
  h_emcal_energy_gl1[mbd_bit]->Draw("hist");
  drawText("06/12/2024", 0.75, 0.92);
  //  h_emcal_energy_gl1[emcal_bit + 0]->Draw("hist same");
  h_emcal_energy_gl1[emcal_bit + 1]->Draw("hist same");
  h_emcal_energy_gl1[emcal_bit + 2]->Draw("hist same");
  h_emcal_energy_gl1[emcal_bit + 3]->Draw("hist same");
  TLegend *l = new TLegend(0.6, 0.4, 0.85, 0.7);
  l->SetBorderSize(0);
  l->SetTextSize(0.033);
  l->SetTextFont(42);
 
  l->AddEntry(h_emcal_energy_gl1[mbd_bit], "MBD >= 1");
  //  l->AddEntry(h_emcal_energy_gl1[emcal_bit + 0], "Photon 1");
  l->AddEntry(h_emcal_energy_gl1[emcal_bit + 1], "Photon 2 GeV");
  l->AddEntry(h_emcal_energy_gl1[emcal_bit + 2], "Photon 3 GeV");
  l->AddEntry(h_emcal_energy_gl1[emcal_bit + 3], "Photon 4 GeV");
  l->Draw("same");

  float x = 0.54;
  float y = 0.8;
  float dy = 0.05;
  DrawSPHENIXprelim(x, y);
  drawText("p+p - #sqrt{s} = 200 GeV", x, y - dy, 0, kBlack, 0.03); 
  //  drawText(Form("Run %d", runnumber), x , y - 2*dy, 0, kBlack, 0.03); 

  c1->cd();
  TPad *p2 = new TPad("p2","p2", 0.0, 0.0, 1.0,0.4);
  p2->SetTicks(1,1);
  p2->SetLeftMargin(0.2);
  p2->SetBottomMargin(0.2);
  p2->Draw();
  p2->cd();

  h_turn_on_emcal[1]->SetTitle(Form(";Maximum 8x8 Energy Sum (EMCAL) [GeV]; MBD Coincidence/%s" , triggername.c_str()));
  h_turn_on_emcal[1]->SetLabelSize(0.05, "XY");
  h_turn_on_emcal[1]->SetTitleSize(0.05, "XY");

  SetLineAtt(h_turn_on_emcal[0], kBlue, 2, 1);
  SetLineAtt(h_turn_on_emcal[1], kRed, 2, 1);
  SetLineAtt(h_turn_on_emcal[2], kSpring + 2, 2, 1);
  SetLineAtt(h_turn_on_emcal[3], kBlue, 2, 1);
  h_turn_on_emcal[1]->GetXaxis()->SetRangeUser(0, 10);
  h_turn_on_emcal[1]->SetMaximum(2.0);
  h_turn_on_emcal[1]->SetMinimum(0);
  

  // TF1 *f1 = new TF1("f1","[0]*(1 + TMath::Erf((x-[1])/[2]))",0, 10);
  // f1->SetParameter(0, 0.3);
  // f1->SetParameter(1, 4);
  // f1->SetParameter(2, 1);
  // h_turn_on_emcal[1]->Fit(f1, "RN");
  // //  h_turn_on_emcal[0]->Draw();
  h_turn_on_emcal[1]->Draw("");
  TLine *dotted31 = new TLine(0, 1, 10, 1);
  SetLineAtt(dotted31, kBlack, 1, 4);
  dotted31->Draw("same");
  h_turn_on_emcal[1]->Draw("same");
  h_turn_on_emcal[2]->Draw("same");
  h_turn_on_emcal[3]->Draw("same");


  // f1->SetLineColor(kBlack);
  // f1->Draw("same");

  // double plateau = f1->GetParameter(0)*2;
  // double half = 0;
  // double y1 = 0;
  // while (y1 < plateau/2.)
  //   {
  //     half += 0.01;
  //     y1 = f1->Eval(half);
  //     if (half == 20) break;
  //   }

  // drawText("Photon 2 GeV :", 0.17, 0.8, 0, kBlack, 0.06);
  // drawText(Form("   plateau     : %1.2f", plateau), 0.17, 0.73, 0, kBlack, 0.06);
  // drawText(Form("   Half height : %1.2f GeV", half), 0.17, 0.66, 0, kBlack, 0.06);

  c1->Print(Form("../plotted/photon_turn_on.pdf"));
  c1->Print(Form("../plotted/photon_turn_on.png"));
  // Jet
  TCanvas *cc1 = new TCanvas("cc2","cc2", 500, 700);
  
  TPad *pc1 = new TPad("pc1","pc1", 0.0, 0.4, 1.0, 1.0);
  pc1->SetTicks(1,1);
  pc1->SetLeftMargin(0.2);
  pc1->Draw();	
  pc1->cd();
  pc1->SetLogy();  
  h_cluster_energy_gl1[mbd_bit]->SetTitle(";Maximum Cluster Energy Sum (EMCAL) [GeV]; Prescale Counts");
  SetLineAtt(h_cluster_energy_gl1[mbd_bit], kBlack, 2, 1);
  SetLineAtt(h_cluster_energy_gl1[emcal_bit + 0], kBlue, 2, 1);
  SetLineAtt(h_cluster_energy_gl1[emcal_bit + 1], kRed, 2, 1);
  SetLineAtt(h_cluster_energy_gl1[emcal_bit + 2], kSpring+2, 2, 1);
  SetLineAtt(h_cluster_energy_gl1[emcal_bit + 3], kBlue, 2, 1);
  h_cluster_energy_gl1[mbd_bit]->SetMinimum(1);
  h_cluster_energy_gl1[mbd_bit]->Draw("hist");
  //drawText("06/12/2024", 0.75, 0.92);
  //  h_cluster_energy_gl1[cluster_bit + 0]->Draw("hist same");
  h_cluster_energy_gl1[emcal_bit + 1]->Draw("hist same");
  h_cluster_energy_gl1[emcal_bit + 2]->Draw("hist same");
  h_cluster_energy_gl1[emcal_bit + 3]->Draw("hist same");
  l = new TLegend(0.6, 0.4, 0.85, 0.7);
  l->SetBorderSize(0);
  l->SetTextSize(0.033);
  l->SetTextFont(42);
 
  l->AddEntry(h_cluster_energy_gl1[mbd_bit], "MBD >= 1");
  //  l->AddEntry(h_cluster_energy_gl1[emcal_bit + 0], "Photon 1");
  l->AddEntry(h_cluster_energy_gl1[emcal_bit + 1], "Photon 2 GeV");
  l->AddEntry(h_cluster_energy_gl1[emcal_bit + 2], "Photon 3 GeV");
  l->AddEntry(h_cluster_energy_gl1[emcal_bit + 3], "Photon 4 GeV");
  l->Draw("same");

  DrawSPHENIXraw(x, y);
  drawText("p+p - #sqrt{s} = 200 GeV", x, y - dy, 0, kBlack, 0.03); 
  //  drawText(Form("Run %d", runnumber), x , y - 2*dy, 0, kBlack, 0.03); 

  cc1->cd();
  TPad *pc2 = new TPad("pc2","pc2", 0.0, 0.0, 1.0,0.4);
  pc2->SetTicks(1,1);
  pc2->SetLeftMargin(0.2);
  pc2->SetBottomMargin(0.2);
  pc2->Draw();
  pc2->cd();

  h_turn_on_cluster[1]->SetTitle(Form(";Maximum Cluster Energy Sum [GeV]; MBD Coincidence/%s" , triggername.c_str()));
  h_turn_on_cluster[1]->SetLabelSize(0.05, "XY");
  h_turn_on_cluster[1]->SetTitleSize(0.05, "XY");

  SetLineAtt(h_turn_on_cluster[0], kBlue, 2, 1);
  SetLineAtt(h_turn_on_cluster[1], kRed, 2, 1);
  SetLineAtt(h_turn_on_cluster[2], kSpring + 2, 2, 1);
  SetLineAtt(h_turn_on_cluster[3], kBlue, 2, 1);
  h_turn_on_cluster[1]->GetXaxis()->SetRangeUser(0, 10);
  h_turn_on_cluster[1]->SetMaximum(2.0);
  h_turn_on_cluster[1]->SetMinimum(0);
  

  // TF1 *f1 = new TF1("f1","[0]*(1 + TMath::Erf((x-[1])/[2]))",0, 10);
  // f1->SetParameter(0, 0.3);
  // f1->SetParameter(1, 4);
  // f1->SetParameter(2, 1);
  // h_turn_on_emcal[1]->Fit(f1, "RN");
  // //  h_turn_on_emcal[0]->Draw();
  h_turn_on_cluster[1]->Draw("");
  TLine *dotted1 = new TLine(0, 1, 10, 1);
  SetLineAtt(dotted1, kBlack, 1, 4);
  dotted1->Draw("same");
  h_turn_on_cluster[1]->Draw("same");
  h_turn_on_cluster[2]->Draw("same");
  h_turn_on_cluster[3]->Draw("same");


  // f1->SetLineColor(kBlack);
  // f1->Draw("same");

  // double plateau = f1->GetParameter(0)*2;
  // double half = 0;
  // double y1 = 0;
  // while (y1 < plateau/2.)
  //   {
  //     half += 0.01;
  //     y1 = f1->Eval(half);
  //     if (half == 20) break;
  //   }

  // drawText("Photon 2 GeV :", 0.17, 0.8, 0, kBlack, 0.06);
  // drawText(Form("   plateau     : %1.2f", plateau), 0.17, 0.73, 0, kBlack, 0.06);
  // drawText(Form("   Half height : %1.2f GeV", half), 0.17, 0.66, 0, kBlack, 0.06);

  c1->Print(Form("../plotted/cluster_turn_on.pdf"));
  c1->Print(Form("../plotted/cluster_turn_on.png"));
  // Jet

  TCanvas *c = new TCanvas("c1","c1", 500, 700);
  
  TPad *p1 = new TPad("p1","p1", 0.0, 0.4, 1.0, 1.0);
  p1->SetTicks(1,1);
  p1->SetLeftMargin(0.2);
  p1->Draw();	
  p1->cd();
  p1->SetLogy();  
  h_jet_energy_gl1[mbd_bit]->SetTitle(";Maximum 0.8x0.8 Energy Sum (EMCAL + HCAL) [GeV]; Prescale Counts");
  SetLineAtt(h_jet_energy_gl1[mbd_bit], kBlack, 2, 1);
  SetLineAtt(h_jet_energy_gl1[jet_bit + 0], kBlue, 2, 1);
  SetLineAtt(h_jet_energy_gl1[jet_bit + 1], kRed, 2, 1);
  SetLineAtt(h_jet_energy_gl1[jet_bit + 2], kSpring+2, 2, 1);
  SetLineAtt(h_jet_energy_gl1[jet_bit + 3], kBlue, 2, 1);
  h_jet_energy_gl1[mbd_bit]->SetMinimum(1);
  h_jet_energy_gl1[mbd_bit]->GetXaxis()->SetRangeUser(0, 35);
  h_jet_energy_gl1[mbd_bit]->Draw("hist");
  drawText("06/12/2024", 0.75, 0.92);
  //  h_jet_energy_gl1[jet_bit + 0]->Draw("hist same");
  h_jet_energy_gl1[jet_bit + 1]->Draw("hist same");
  h_jet_energy_gl1[jet_bit + 2]->Draw("hist same");
  h_jet_energy_gl1[jet_bit + 3]->Draw("hist same");
  l = new TLegend(0.6, 0.4, 0.85, 0.7);
  l->SetBorderSize(0);
  l->SetTextSize(0.033);
  l->SetTextFont(42);  
  l->AddEntry(h_jet_energy_gl1[mbd_bit], "MBD >= 1");
  //  l->AddEntry(h_jet_energy_gl1[jet_bit + 0], "Jet 1");
  l->AddEntry(h_jet_energy_gl1[jet_bit + 1], "Jet 6 GeV");
  l->AddEntry(h_jet_energy_gl1[jet_bit + 2], "Jet 8 GeV");
  l->AddEntry(h_jet_energy_gl1[jet_bit + 3], "Jet 10 GeV");
  l->Draw("same");

  DrawSPHENIXprelim(x, y);
  drawText("p+p - #sqrt{s} = 200 GeV", x, y - dy, 0, kBlack, 0.03); 
  //  drawText(Form("Run %d", runnumber), x , y - 2*dy, 0, kBlack, 0.03); 

  c->cd();
  TPad *p22 = new TPad("p22","p22", 0.0, 0.0, 1.0,0.4);
  p22->SetTicks(1,1);
  p22->Draw();
  p22->SetLeftMargin(0.2);
  p22->SetBottomMargin(0.2);
  p22->cd();
  h_turn_on_jet[1]->SetTitle(";Maximum 8x8 Energy Sum (EMCAL + HCAL) [GeV]; MBD Coincidence/JET");
  h_turn_on_jet[1]->SetLabelSize(0.05, "XY");
  h_turn_on_jet[1]->SetTitleSize(0.05, "XY");
  SetLineAtt(h_turn_on_jet[0], kBlue, 2, 1);
  SetLineAtt(h_turn_on_jet[1], kRed, 2, 1);
  SetLineAtt(h_turn_on_jet[2], kSpring + 2, 2, 1);
  SetLineAtt(h_turn_on_jet[3], kBlue, 2, 1);
  h_turn_on_jet[1]->GetXaxis()->SetRangeUser(0, 20);
  h_turn_on_jet[1]->SetMaximum(2);
  h_turn_on_jet[1]->SetMinimum(0);

  // TF1 *f2 = new TF1("f2","[0]*(1 + TMath::Erf((x-[1])/[2]))",0, 20);
  // f2->SetParameter(0, 0.4);
  // f2->SetParameter(1, 5);
  // f2->SetParameter(2, 1);
  // h_turn_on_jet[1]->GetXaxis()->SetRangeUser(0, 30);
  // h_turn_on_jet[1]->Fit(f2, "RN");
  h_turn_on_jet[1]->Draw();
  //h_turn_on_jet[1]->Draw("same");
  TLine *dotted2 = new TLine(0, 1, 30, 1);
  SetLineAtt(dotted2, kBlack, 1, 4);
  dotted2->Draw("same");
  h_turn_on_jet[1]->Draw("same");
  h_turn_on_jet[2]->Draw("same");
  h_turn_on_jet[3]->Draw("same");

  // f2->SetLineColor(kBlack);
  // f2->Draw("same");

  // plateau = f2->GetParameter(0)*2;
  // half = 0;
  // y1 = 0;
  // while (y1 < plateau/2.)
  //   {
  //     half += 0.01;
  //     y1 = f2->Eval(half);
  //     if (half == 10) break;
  //   }

  // drawText("Jet 6 GeV :", 0.17, 0.8, 0, kBlack, 0.06);
  // drawText(Form("   plateau     : %1.2f", plateau), 0.17, 0.73, 0, kBlack, 0.06);
  // drawText(Form("   Half height : %1.2f GeV", half), 0.17, 0.66, 0, kBlack, 0.06);

  c->Print(Form("../plotted/jet_turn_on.pdf"));
  c->Print(Form("../plotted/jet_turn_on.png"));

}

void Draw_Trigger_ThresholdScan()
{
  const int nruns = 4;
  int runnumbers[4] = {44042, 44043, 44045, 44047};
  float scaledowns[4] = {50., 60., 70., 90.};
  float thresholds[4][4] ={{2,3,4,5},{3,5,6,8}, {4,6,8,10}, {6,10,14,18}};
  float elut[4] = {2.0, 1.5, 1.0, 0.5};

  gStyle->SetOptStat(0);

  float peaks_hcal[4][4] = {0};
  float peaks_emcal[4][4] = {0};
  for (int i = 0; i< nruns; i++)
    {
      int runnumber = runnumbers[i];
      std::string outfilename = Form("/sphenix/user/dlis/Projects/CaloTriggerEmulator/plots/HIST_PHOTON-%08d.root", runnumber);

      TFile *f = new TFile(outfilename.c_str(), "r");
      if (!f)
	{
	  std::cout << "nope" << std::endl;
	  return;
	}
      float scaledown = scaledowns[i] + 1.;
      int triggerbit = 24;
      std::string histtrigname = "h_emcal_energy_gl1_" + to_string(triggerbit);
      std::string triggername = "PHOTON";
      TH1D *h_emcal_energy_mbd = (TH1D*) f->Get("h_emcal_energy_gl1_10");
      TH1D *h_jet_energy_mbd = (TH1D*) f->Get("h_jet_energy_gl1_10");
      TH1D *h_hcal_energy_mbd = (TH1D*) f->Get("h_hcal_energy_gl1_10");
      h_emcal_energy_mbd->Scale(scaledown);
      h_jet_energy_mbd->Scale(scaledown);
      h_hcal_energy_mbd->Scale(scaledown);
      
      //TH1D *h_emcal_energy_ref = (TH1D*) f->Get("h_emcal_energy_ref");
      
      //  h_emcal_energy_mbd->Add(h_emcal_energy_ref);
      TH1D *h_emcal_energy_trig0 = (TH1D*) f->Get(histtrigname.c_str());
      histtrigname = "h_emcal_energy_gl1_" + to_string(triggerbit + 1);
      TH1D *h_emcal_energy_trig1 = (TH1D*) f->Get(histtrigname.c_str());
      histtrigname = "h_emcal_energy_gl1_" + to_string(triggerbit + 2);
      TH1D *h_emcal_energy_trig2 = (TH1D*) f->Get(histtrigname.c_str());
      histtrigname = "h_emcal_energy_gl1_" + to_string(triggerbit + 3);
      TH1D *h_emcal_energy_trig3 = (TH1D*) f->Get(histtrigname.c_str());
      

      int jetbit = 16;
      histtrigname = "h_jet_energy_gl1_" + to_string(jetbit);
      TH1D *h_jet_energy_trig0 = (TH1D*) f->Get(histtrigname.c_str());
      histtrigname = "h_jet_energy_gl1_" + to_string(jetbit + 1);
      TH1D *h_jet_energy_trig1 = (TH1D*) f->Get(histtrigname.c_str());
      histtrigname = "h_jet_energy_gl1_" + to_string(jetbit + 2);
      TH1D *h_jet_energy_trig2 = (TH1D*) f->Get(histtrigname.c_str());
      histtrigname = "h_jet_energy_gl1_" + to_string(jetbit + 3);
      TH1D *h_jet_energy_trig3 = (TH1D*) f->Get(histtrigname.c_str());
      //  h_emcal_energy_mbd->Add(h_emcal_energy_ref);
      
      histtrigname = "h_hcal_energy_gl1_" + to_string(jetbit);
      TH1D *h_hcal_energy_trig0 = (TH1D*) f->Get(histtrigname.c_str());
      histtrigname = "h_hcal_energy_gl1_" + to_string(jetbit + 1);
      TH1D *h_hcal_energy_trig1 = (TH1D*) f->Get(histtrigname.c_str());
      histtrigname = "h_hcal_energy_gl1_" + to_string(jetbit + 2);
      TH1D *h_hcal_energy_trig2 = (TH1D*) f->Get(histtrigname.c_str());
      histtrigname = "h_hcal_energy_gl1_" + to_string(jetbit + 3);
      TH1D *h_hcal_energy_trig3 = (TH1D*) f->Get(histtrigname.c_str());
      
      
      TH1D *h_turn_on = (TH1D*) h_emcal_energy_trig0->Clone();
      
      h_turn_on->Divide(h_emcal_energy_mbd);
      
      TH1D *h_turn_on2 = (TH1D*) h_emcal_energy_trig1->Clone();
      
      h_turn_on2->Divide(h_emcal_energy_mbd);
      
      TH1D *h_turn_on3 = (TH1D*) h_emcal_energy_trig2->Clone();
      
      h_turn_on3->Divide(h_emcal_energy_mbd);
      
      TH1D *h_turn_on4 = (TH1D*) h_emcal_energy_trig3->Clone();
      
      h_turn_on4->Divide(h_emcal_energy_mbd);
      
      
      TH1D *h_turn_on_j = (TH1D*) h_jet_energy_trig0->Clone();

      h_turn_on_j->Divide(h_jet_energy_mbd);

      TH1D *h_turn_on_j2 = (TH1D*) h_jet_energy_trig1->Clone();

      h_turn_on_j2->Divide(h_jet_energy_mbd);

      TH1D *h_turn_on_j3 = (TH1D*) h_jet_energy_trig2->Clone();

      h_turn_on_j3->Divide(h_jet_energy_mbd);

      TH1D *h_turn_on_j4 = (TH1D*) h_jet_energy_trig3->Clone();

      h_turn_on_j4->Divide(h_jet_energy_mbd);


      TH1D *h_turn_on_h = (TH1D*) h_hcal_energy_trig0->Clone();

      h_turn_on_h->Divide(h_hcal_energy_mbd);

      TH1D *h_turn_on_h2 = (TH1D*) h_hcal_energy_trig1->Clone();

      h_turn_on_h2->Divide(h_hcal_energy_mbd);

      TH1D *h_turn_on_h3 = (TH1D*) h_hcal_energy_trig2->Clone();

      h_turn_on_h3->Divide(h_hcal_energy_mbd);

      TH1D *h_turn_on_h4 = (TH1D*) h_hcal_energy_trig3->Clone();

      h_turn_on_h4->Divide(h_hcal_energy_mbd);
      h_emcal_energy_trig0->GetXaxis()->SetRangeUser(1, 10);
      h_emcal_energy_trig1->GetXaxis()->SetRangeUser(1, 10);
      h_emcal_energy_trig2->GetXaxis()->SetRangeUser(1, 10);
      h_emcal_energy_trig3->GetXaxis()->SetRangeUser(1, 10);
      h_hcal_energy_trig0->GetXaxis()->SetRangeUser(1, 10);
      h_hcal_energy_trig1->GetXaxis()->SetRangeUser(1, 10);
      h_hcal_energy_trig2->GetXaxis()->SetRangeUser(1, 10);
      h_hcal_energy_trig3->GetXaxis()->SetRangeUser(1, 10);
      peaks_emcal[i][0] = h_emcal_energy_trig0->GetBinCenter(h_emcal_energy_trig0->GetMaximumBin());
      peaks_emcal[i][1] = h_emcal_energy_trig1->GetBinCenter(h_emcal_energy_trig1->GetMaximumBin());
      peaks_emcal[i][2] = h_emcal_energy_trig2->GetBinCenter(h_emcal_energy_trig2->GetMaximumBin());
      peaks_emcal[i][3] = h_emcal_energy_trig3->GetBinCenter(h_emcal_energy_trig3->GetMaximumBin());
      peaks_hcal[i][0] = h_hcal_energy_trig0->GetBinCenter(h_hcal_energy_trig0->GetMaximumBin());
      peaks_hcal[i][1] = h_hcal_energy_trig1->GetBinCenter(h_hcal_energy_trig1->GetMaximumBin());
      peaks_hcal[i][2] = h_hcal_energy_trig2->GetBinCenter(h_hcal_energy_trig2->GetMaximumBin());
      peaks_hcal[i][3] = h_hcal_energy_trig3->GetBinCenter(h_hcal_energy_trig3->GetMaximumBin());

    }

  int color[4] = {kViolet, kBlue, kRed, kOrange};
  TGraph *gEmcal[4];
  TGraph *gHcal[4];
  for (int i = 0; i < 4; i++)
    {
      std::cout << i << " : "<<std::endl;
      for (int j = 0; j < 4; j++)
	{
	  std::cout << "    " << thresholds[i][j] << " -- " << peaks_emcal[i][j] <<" - " << peaks_hcal[i][j] << std::endl;
	}

      gEmcal[i] = new TGraph(4, thresholds[i], peaks_emcal[i]);
      SetMarkerAtt(gEmcal[i], color[i], 2, 29);
      SetLineAtt(gEmcal[i], color[i], 1, 1);
      gHcal[i] = new TGraph(4, thresholds[i], peaks_hcal[i]);
      SetMarkerAtt(gHcal[i], color[i], 2, 30);
      SetLineAtt(gHcal[i], color[i], 1, 2);


    }


  float slopes_emcal[4] = {0};
  float slopes_hcal[4] = {0};

  float ints_emcal[4] = {0};
  float ints_hcal[4] = {0};

  TF1 *flines[8];
  for (int i = 0; i<4;i++)
    {
      flines[i]= new TF1(Form("fline_%d", i),"[0] + x*[1]", 0, 20);
      flines[i]->SetParameter(0, 0);
      flines[i]->SetParameter(1, 1);
      gEmcal[i]->Fit(Form("fline_%d", i), "NDORQ");
      slopes_emcal[i] = flines[i]->GetParameter(1);
      ints_emcal[i] = flines[i]->GetParameter(0);
      flines[i]->SetLineColor(color[i]);
      flines[i]->SetLineStyle(1);
      flines[i+4]= new TF1(Form("fline_%d", i+4),"[0] + x*[1]", 0, 20);
      flines[i+4]->SetParameter(0, 0);
      flines[i+4]->SetParameter(1, 1);
      gHcal[i]->Fit(Form("fline_%d", i+4), "NDORQ");
      slopes_hcal[i] = flines[i+4]->GetParameter(1);
      ints_hcal[i] = flines[i+4]->GetParameter(0);
      flines[i+4]->SetLineColor(color[i]);
      flines[i+4]->SetLineStyle(2);
    }


  TGraph *g_slopes_emcal = new TGraph(4, elut, slopes_emcal);
  TGraph *g_slopes_hcal = new TGraph(4, elut, slopes_hcal);
  TGraph *g_ints_emcal = new TGraph(4, elut, ints_emcal);
  TGraph *g_ints_hcal = new TGraph(4, elut, ints_hcal);

  TF1 *f_slope_emcal = new TF1("f_slope_emcal", "[0] + x*[1]", 0, 20);
  TF1 *f_int_emcal = new TF1("f_int_emcal", "[0] + x*[1]", 0, 20);
  TF1 *f_slope_hcal = new TF1("f_slope_hcal", "[0] + x*[1]", 0, 20);
  TF1 *f_int_hcal = new TF1("f_int_hcal", "[0] + x*[1]", 0, 20);

  f_slope_emcal->SetParameters(1, 1);
  f_int_emcal->SetParameters(1, 1);
  f_slope_hcal->SetParameters(1, 1);
  f_int_hcal->SetParameters(1, 1);

  g_slopes_emcal->Fit("f_slope_emcal", "NDORQ");
  g_slopes_hcal->Fit("f_slope_hcal", "NDORQ");
  g_ints_emcal->Fit("f_int_emcal", "NDORQ");
  g_ints_hcal->Fit("f_int_hcal", "NDORQ");

  SetyjPadStyle();
  TCanvas *c1 = new TCanvas("c1","c1", 700, 500);
  TPad *pad = new TPad("p","p", 0, 0, 1, 1);
  pad->SetRightMargin(0.3);
  pad->Draw();
  pad->cd();
  TH1D *h = new TH1D("h","; Digital Threshold; Offline 8x8 Tower Energy [GeV]", 1, 0, 20);
  h->SetMaximum(10);
  h->SetMinimum(0);
  h->Draw();
  for (int i = 0; i < 4; i++)
    {
      TLine *l = new TLine(0, i+1, 20, i+1);
      SetLineAtt(l, kBlack, 1, 2);
      l->Draw("same");
    }
  for (int i = 0; i < 4; i++)
    {
      gEmcal[i]->Draw("P");
      gHcal[i]->Draw("P");
      flines[i+4]->Draw("same");
      flines[i]->Draw("same");
    }

  // Find 4 Gev;
  f1 = new TF1("f1", "[0] + [1]*x", 0, 20);
  f2 = new TF1("f2", "[0]", 0, 20);
  TGraph *gthresh[4];
  TGraph *gthresh_hcal[4];
  for (int i =0; i < 4; i++)
    {
      gthresh[i] = new TGraph(15);
      gthresh_hcal[i] = new TGraph(15);
    }
  SetLineAtt(f1, kBlack, 1, 4);
  for (int i = 0; i < 15; i++)
    {
      float elutp = elut[0] - 0.1*i;
      float a1 = f_slope_emcal->Eval(elutp);
      float a2 = f_slope_hcal->Eval(elutp);
      float b1 = f_int_emcal->Eval(elutp);
      float b2 = f_int_hcal->Eval(elutp);
      TF1 *flin = new TF1("flin", "[0] + [1]*x", 0, 20);
      flin->SetParameters(b1, a1);
      f1->SetParameters(b1, a1);

      TF1 *fint = new TF1("fint", finter, 0, 10, 0);
      SetLineAtt(flin, kBlack, 1, 4);
      for (int j = 0; j < 4; j++)
	{
	  f2->SetParameter(0, j+1);

	  double found = fint->GetMinimumX();

	  gthresh[j]->SetPoint(i,  elutp, found);
	  SetLineAtt(gthresh[j], kBlack, 2, 1);
	}

      f1->SetParameters(b2, a2);

      for (int j = 0; j < 4; j++)
	{
	  f2->SetParameter(0, j+1);

	  double found = fint->GetMinimumX();

	  gthresh_hcal[j]->SetPoint(i,  elutp, found);
	  SetLineAtt(gthresh_hcal[j], kBlack, 2, 3);	  
	}

      flin->Draw("same");
    }

  TLegend *l = new TLegend(0.75, 0.2, 1.0, 0.6);
  for (int i = 0; i < 4; i++)
    {
      l->AddEntry(gEmcal[i], Form("EMCAL - %1.1f GeV", elut[i]));
      l->AddEntry(gHcal[i], Form("HCAL - %1.1f GeV", elut[i]));
    }
  l->SetHeader("Detector - Est. Energy Step");
  l->SetBorderSize(0); 
  l->Draw("same");

  DrawSPHENIXpp(0.75, 0.9);
  drawText("Threshold Scan", 0.75, 0.78);
  drawText("8x8 Non-overlap", 0.75, 0.72);

  TCanvas *c2 = new TCanvas("c2","c2", 500, 500);
  TPad *padb = new TPad("p","p", 0, 0.0, 1, 1);
  padb->Draw();
  padb->cd();

  TMultiGraph *mg = new TMultiGraph();
  for (int i = 0; i < 4; i++)
    {
      SetMarkerAtt(gthresh[i], color[i], 1, 34);
      mg->Add(gthresh[i]);
      SetMarkerAtt(gthresh_hcal[i], color[i], 1, 28);
      mg->Add(gthresh_hcal[i]);

    }
  mg->GetYaxis()->SetTitle("Digital threshold");
  mg->GetXaxis()->SetTitle("Calibration Correction");

  mg->Draw("AP");
  mg->GetYaxis()->SetRangeUser(0, 12);
  c2->Update();
  DrawSPHENIXpp(0.6, 0.85);
  drawText("Threshold Scan", 0.6, 0.74);
  drawText("8x8 Non-overlap", 0.6, 0.68);
  TLegend *l2 = new TLegend(0.65, 0.45, 0.85, 0.62);
  l2->AddEntry(gthresh[0], "1 GeV","p");
  l2->AddEntry(gthresh[1], "2 GeV","p");
  l2->AddEntry(gthresh[2], "3 GeV","p");
  l2->AddEntry(gthresh[3], "4 GeV","p");
  l2->SetBorderSize(0);
  l2->Draw("same");
  TLegend *l3 = new TLegend(0.4, 0.75, 0.6, 0.85);
  l3->AddEntry(gthresh[0], "EMCAL","lp");
  l3->AddEntry(gthresh_hcal[0], "HCAL","lp");  
  l3->SetBorderSize(0);
  l3->Draw("same");

  TLine *thresh_emcal = new TLine(0.8, 0, 0.8, 12);
  SetLineAtt(thresh_emcal, kBlack, 2, 1);
  thresh_emcal->Draw("same");
  TLine *thresh_hcal = new TLine(0.9, 0, 0.9, 12);
  SetLineAtt(thresh_hcal, kBlack, 2, 3);
  thresh_hcal->Draw("same");
  mg->Draw("P");
}
// void Draw_Trigger_ThresholdScan()
// {
//   const int nruns = 4;
//   int runnumbers[4] = {44042, 44043, 44045, 44047};
//   float scaledowns[4] = {50., 60., 70., 90.}
//   int thresholds[4][4] ={{2,3,4,5},{3,5,6,8}, {4,6,8,10}, {6,10,14,18}};

//   gStyle->SetOptStat(0);

//   float energy_thresholds[4][4] = {0};
//   float efficiencies[4][4] = {0};
//   for (int irun = 0; irun < nruns ; irun++)
//     {
//       int runnumber = runnumbers[irun];
//       float scaledown = scaledowns[irun] + 1.;
//       std::string outfilename = Form("/sphenix/user/dlis/Projects/CaloTriggerEmulator/plots/HIST_%s-%08d.root",trigger.c_str(), runnumber);

//       TFile *f = new TFile(outfilename.c_str(), "r");
//       if (!f)
// 	{
// 	  std::cout << "nope" << std::endl;
// 	  return;
// 	}

//       std::string histtrigname = "h_emcal_energy_gl1_" + to_string(triggerbit);
//       std::string triggername = "PHOTON";
//       TH1D *h_emcal_energy_mbd = (TH1D*) f->Get("h_emcal_energy_gl1_10");
//       TH1D *h_hcal_energy_mbd = (TH1D*) f->Get("h_hcal_energy_gl1_10");
//       h_emcal_energy_mbd->Scale(scaledown);
//       h_hcal_energy_mbd->Scale(scaledown);
  
//       //TH1D *h_emcal_energy_ref = (TH1D*) f->Get("h_emcal_energy_ref");
      
//       //  h_emcal_energy_mbd->Add(h_emcal_energy_ref);
//       TH1D *h_emcal_energy_trig0 = (TH1D*) f->Get(histtrigname.c_str());
//       histtrigname = "h_emcal_energy_gl1_" + to_string(triggerbit + 1);
//       TH1D *h_emcal_energy_trig1 = (TH1D*) f->Get(histtrigname.c_str());
//       histtrigname = "h_emcal_energy_gl1_" + to_string(triggerbit + 2);
//       TH1D *h_emcal_energy_trig2 = (TH1D*) f->Get(histtrigname.c_str());
//       histtrigname = "h_emcal_energy_gl1_" + to_string(triggerbit + 3);
//       TH1D *h_emcal_energy_trig3 = (TH1D*) f->Get(histtrigname.c_str());

//       int jetbit = 16;
//       histtrigname = "h_jet_energy_gl1_" + to_string(jetbit);
//       TH1D *h_jet_energy_trig0 = (TH1D*) f->Get(histtrigname.c_str());
//       histtrigname = "h_jet_energy_gl1_" + to_string(jetbit + 1);
//       TH1D *h_jet_energy_trig1 = (TH1D*) f->Get(histtrigname.c_str());
//       histtrigname = "h_jet_energy_gl1_" + to_string(jetbit + 2);
//       TH1D *h_jet_energy_trig2 = (TH1D*) f->Get(histtrigname.c_str());
//       histtrigname = "h_jet_energy_gl1_" + to_string(jetbit + 3);
//       TH1D *h_jet_energy_trig3 = (TH1D*) f->Get(histtrigname.c_str());

//       int hcalbit = 16;
//       histtrigname = "h_hcal_energy_gl1_" + to_string(hcalbit);
//       TH1D *h_hcal_energy_trig0 = (TH1D*) f->Get(histtrigname.c_str());
//       histtrigname = "h_hcal_energy_gl1_" + to_string(hcalbit + 1);
//       TH1D *h_hcal_energy_trig1 = (TH1D*) f->Get(histtrigname.c_str());
//       histtrigname = "h_hcal_energy_gl1_" + to_string(hcalbit + 2);
//       TH1D *h_hcal_energy_trig2 = (TH1D*) f->Get(histtrigname.c_str());
//       histtrigname = "h_hcal_energy_gl1_" + to_string(hcalbit + 3);
//       TH1D *h_hcal_energy_trig3 = (TH1D*) f->Get(histtrigname.c_str());


//       TH1D *h_turn_on = (TH1D*) h_emcal_energy_trig0->Clone();

//       h_turn_on->Divide(h_emcal_energy_mbd);

//       TH1D *h_turn_on2 = (TH1D*) h_emcal_energy_trig1->Clone();

//       h_turn_on2->Divide(h_emcal_energy_mbd);

//       TH1D *h_turn_on3 = (TH1D*) h_emcal_energy_trig2->Clone();

//       h_turn_on3->Divide(h_emcal_energy_mbd);

//       TH1D *h_turn_on4 = (TH1D*) h_emcal_energy_trig3->Clone();

//       h_turn_on4->Divide(h_emcal_energy_mbd);


//       TH1D *h_turn_on_j = (TH1D*) h_jet_energy_trig0->Clone();

//       h_turn_on_j->Divide(h_jet_energy_mbd);

//       TH1D *h_turn_on_j2 = (TH1D*) h_jet_energy_trig1->Clone();

//       h_turn_on_j2->Divide(h_jet_energy_mbd);

//       TH1D *h_turn_on_j3 = (TH1D*) h_jet_energy_trig2->Clone();

//       h_turn_on_j3->Divide(h_jet_energy_mbd);

//       TH1D *h_turn_on_j4 = (TH1D*) h_jet_energy_trig3->Clone();

//       h_turn_on_j4->Divide(h_jet_energy_mbd);



//       TH1D *h_turn_on_h = (TH1D*) h_hcal_energy_trig0->Clone();

//       h_turn_on_h->Divide(h_hcal_energy_mbd);

//       TH1D *h_turn_on_h2 = (TH1D*) h_hcal_energy_trig1->Clone();

//       h_turn_on_h2->Divide(h_hcal_energy_mbd);

//       TH1D *h_turn_on_h3 = (TH1D*) h_hcal_energy_trig2->Clone();

//       h_turn_on_h3->Divide(h_hcal_energy_mbd);

//       TH1D *h_turn_on_h4 = (TH1D*) h_hcal_energy_trig3->Clone();

//       h_turn_on_h4->Divide(h_hcal_energy_mbd);


//       TCanvas *c1 = new TCanvas("c2","c2", 500, 700);
  
//       TPad *p12 = new TPad("p12","p12", 0.0, 0.4, 1.0, 1.0);
//       p12->SetTicks(1,1);
//       p12->Draw();	
//       p12->cd();
//       p12->SetLogy();  
//       h_emcal_energy_mbd->SetTitle(";Maximum 8x8 Energy Sum (EMCAL) [GeV]; Prescale Counts");
//       SetLineAtt(h_emcal_energy_mbd, kBlack, 2, 1);
//       SetLineAtt(h_emcal_energy_trig0, kBlue, 2, 1);
//       SetLineAtt(h_emcal_energy_trig1, kRed, 2, 1);
//       SetLineAtt(h_emcal_energy_trig2, kSpring+2, 2, 1);
//       SetLineAtt(h_emcal_energy_trig3, kOrange, 2, 1);
//       h_emcal_energy_mbd->SetMinimum(1);
//       h_emcal_energy_mbd->Draw();
//       h_emcal_energy_trig0->Draw("same");
//       h_emcal_energy_trig1->Draw("same");
//       h_emcal_energy_trig2->Draw("same");
//       h_emcal_energy_trig3->Draw("same");
//       TLegend *l = new TLegend(0.65, 0.4, 0.85, 0.7);
//       l->SetBorderSize(0);  
//       l->AddEntry(h_emcal_energy_mbd, "MBD >= 1");
//       l->AddEntry(h_emcal_energy_trig0, "Photon 1");
//       l->AddEntry(h_emcal_energy_trig1, "Photon 2");
//       l->AddEntry(h_emcal_energy_trig2, "Photon 3");
//       l->AddEntry(h_emcal_energy_trig3, "Photon 4");
//       l->Draw("same");

//       float x = 0.6;
//       float y = 0.8;
//       float dy = 0.05;
//       DrawSPHENIXraw(x, y);
//       drawText("pp - #sqrt{s} = 200 GeV", x, y - dy, 0, kBlack, 0.03); 
//       drawText(Form("Run %d", runnumber), x , y - 2*dy, 0, kBlack, 0.03); 

//       c1->cd();
//       TPad *p2 = new TPad("p2","p2", 0.0, 0.0, 1.0,0.4);
//       p2->SetTicks(1,1);
//       p2->Draw();
//       p2->cd();
//       h_turn_on->SetTitle(Form(";Maximum 8x8 Energy Sum (EMCAL) [GeV]; MBD Coincidence/%s" , triggername.c_str()));
//       SetLineAtt(h_turn_on, kBlue, 2, 1);
//       SetLineAtt(h_turn_on2, kRed, 2, 1);
//       SetLineAtt(h_turn_on3, kSpring + 2, 2, 1);
//       SetLineAtt(h_turn_on4, kOrange, 2, 1);
//       h_turn_on->GetXaxis()->SetRangeUser(0, 10);
//       h_turn_on->SetMaximum(1.4);
//       h_turn_on->SetMinimum(0);


//       TF1 *f1 = new TF1("f1","[0]*(1 + TMath::Erf((x-[1])/[2]))",2, 7);
//       f1->SetParameter(0, 0.3);
//       f1->SetParameter(1, 4);
//       f1->SetParameter(2, 1);
//       h_turn_on->Fit(f1, "RN");
//       h_turn_on->Draw();
//       h_turn_on2->Draw("same");
//       h_turn_on3->Draw("same");
//       h_turn_on4->Draw("same");
//       f1->SetLineColor(kBlack);
//       f1->Draw("same");

//       double plateau = f1->GetParameter(0)*2;
//       double half = 0;
//       double y1 = 0;
//       while (y1 < plateau/2.)
// 	{
// 	  half += 0.01;
// 	  y1 = f1->Eval(half);
// 	  if (half == 20) break;
// 	}

//       drawText("Photon 1 :", 0.17, 0.8, 0, kBlack, 0.06);
//       drawText(Form("   plateau     : %1.2f", plateau), 0.17, 0.73, 0, kBlack, 0.06);
//       drawText(Form("   Half height : %1.2f GeV", half), 0.17, 0.66, 0, kBlack, 0.06);

//       // Hcal

//       TCanvas *c = new TCanvas("c1","c1", 500, 700);
  
//       TPad *p1 = new TPad("p1","p1", 0.0, 0.4, 1.0, 1.0);
//       p1->SetTicks(1,1);
//       p1->Draw();	
//       p1->cd();
//       p1->SetLogy();  
//       h_hcal_energy_mbd->SetTitle(";Maximum 8x8 Energy Sum (HCAL) [GeV]; Prescale Counts");
//       SetLineAtt(h_hcal_energy_mbd, kBlack, 2, 1);
//       SetLineAtt(h_hcal_energy_trig0, kBlue, 2, 1);
//       SetLineAtt(h_hcal_energy_trig1, kRed, 2, 1);
//       SetLineAtt(h_hcal_energy_trig2, kSpring+2, 2, 1);
//       SetLineAtt(h_hcal_energy_trig3, kOrange, 2, 1);
//       h_hcal_energy_mbd->SetMinimum(1);
//       h_hcal_energy_mbd->Draw();
//       h_hcal_energy_trig0->Draw("same");
//       h_hcal_energy_trig1->Draw("same");
//       h_hcal_energy_trig2->Draw("same");
//       h_hcal_energy_trig3->Draw("same");
//       l = new TLegend(0.65, 0.4, 0.85, 0.7);
//       l->SetBorderSize(0);  
//       l->AddEntry(h_hcal_energy_mbd, "MBD >= 1");
//       l->AddEntry(h_hcal_energy_trig0, "Hcal 1");
//       l->AddEntry(h_hcal_energy_trig1, "Hcal 2");
//       l->AddEntry(h_hcal_energy_trig2, "Hcal 3");
//       l->AddEntry(h_hcal_energy_trig3, "Hcal 4");
//       l->Draw("same");

//       DrawSPHENIXraw(x, y);
//       drawText("pp - #sqrt{s} = 200 GeV", x, y - dy, 0, kBlack, 0.03); 
//       drawText(Form("Run %d", runnumber), x , y - 2*dy, 0, kBlack, 0.03); 

//       c->cd();
//       TPad *p23 = new TPad("p23","p23", 0.0, 0.0, 1.0,0.4);
//       p23->SetTicks(1,1);
//       p23->Draw();
//       p23->cd();
//       htu_rn_on_j->SetTitle(";Maximum 8x8 Energy Sum (EMCAL + HCAL) [GeV]; MBD Coincidence/HCAL");
//       SetLineAtt(h_turn_on_h, kBlue, 2, 1);
//       SetLineAtt(h_turn_on_h2, kRed, 2, 1);
//       SetLineAtt(h_turn_on_h3, kSpring + 2, 2, 1);
//       SetLineAtt(h_turn_on_h4, kOrange, 2, 1);
//       //  h_turn_on_h->GetXaxis()->SetRangeUser(0, 20);
//       h_turn_on_h->SetMaximum(1.4);
//       h_turn_on_h->SetMinimum(0);

//       TF1 *f3 = new TF1("f3","[0]*(1 + TMath::Erf((x-[1])/[2]))",0, 15);
//       f3->SetParameter(0, 0.4);
//       f3->SetParameter(1, 5);
//       f3->SetParameter(2, 1);
//       h_turn_on_h->GetXaxis()->SetRangeUser(0, 15);
//       h_turn_on_h->Fit(f3, "RN");
//       h_turn_on_h->Draw();
//       h_turn_on_h2->Draw("same");
//       h_turn_on_h3->Draw("same");
//       h_turn_on_h4->Draw("same");
//       f3->SetLineColor(kBlack);
//       f3->Draw("same");

//       plateau = f3->GetParameter(0)*2;
//       half = 0;
//       y1 = 0;
//       while (y1 < plateau/2.)
// 	{
// 	  half += 0.01;
// 	  y1 = f3->Eval(half);
// 	  if (half == 10) break;
// 	}

//       drawText("Hcal 1 :", 0.17, 0.8, 0, kBlack, 0.06);
//       drawText(Form("   plateau     : %1.2f", plateau), 0.17, 0.73, 0, kBlack, 0.06);
//       drawText(Form("   Half height : %1.2f GeV", half), 0.17, 0.66, 0, kBlack, 0.06);
//     }

// }

// void GetScaledowns(int runnumber, float scaledowns[])
// {
//   TString testcmd = "";
//   TString cmd = "";

//   for (int i = 0; i < 64; i++)
//     {
//       cmd = "psql -h sphnxdaqdbreplica daq -c 'select scaledown" + to_string(i) + " from gl1_scaledown where runnumber = " + to_string(runnumber) + ";'";
//       gSystem->Exec( cmd.Data() );
//     }
  
// }
void Draw_Trigger_Hists_Sim(int runnumber, const std::string file, int segment, const std::string trigger)
{
  std::string filename = file;
  TFile *f = new TFile(filename.c_str(),"r");
  
  TTree *t = (TTree*) f->Get("ttree");

  unsigned int b_trigger_sum_emcal_ll1[384];
  std::vector<float> *b_truth_particle_pt = nullptr;
  std::vector<float> *b_truth_particle_eta = nullptr;
  std::vector<float> *b_truth_particle_phi = nullptr;
  std::vector<float> *b_emcal_energy = nullptr;

  t->SetBranchAddress("emcal_energy", &b_emcal_energy);
  t->SetBranchAddress("trigger_sum_emcal_ll1", b_trigger_sum_emcal_ll1);
  t->SetBranchAddress("_truth_particle_pt", &b_truth_particle_pt);
  t->SetBranchAddress("_truth_particle_phi", &b_truth_particle_phi);
  t->SetBranchAddress("_truth_particle_eta", &b_truth_particle_eta);
  float cut = 3;
  TH1D *h_ntower_thresh[100][5];
  for (int i = 0; i < 100; i++)
    {
      for (int j = 0; j < 5; j++)
	{
	  h_ntower_thresh[i][j] = new TH1D(Form("h_ntower_thresh_%d_%d", i, j), "; ntower; counts", 50, -0.5, 49.5);
	}
    }

  TH1D *h_pt_thresh[100];
  for (int i = 0; i < 100; i++)
    {
      h_pt_thresh[i] = new TH1D(Form("h_pt_thresh_%d", i), "; pT [GeV]; counts", 50, 0, 25);
    }

  for (int i = 0; i < t->GetEntries(); i++)
    {
      t->GetEntry(i);
      unsigned int max = 0;

      auto maxpt = max_element(b_truth_particle_pt->begin(), b_truth_particle_pt->end());
      int index = std::distance(b_truth_particle_pt->begin(), maxpt);

      if (fabs(b_truth_particle_eta->at(index)) > 1.1) continue;
      if (fabs(b_truth_particle_pt->at(index)) < cut) continue;
      int count = 0;
      for (int j = 0; j < b_emcal_energy->size(); j++)
	{
	  if (b_emcal_energy->at(j) > 0.100)
	    count++;
	}
      for (int j = 0; j < 384; j++)
	{
	  if (b_trigger_sum_emcal_ll1[j] > max)
	    {
	      max = b_trigger_sum_emcal_ll1[j];
	    }
	}

      int hindex = floor(b_truth_particle_pt->at(index)/5.);
      if (hindex >= 5) hindex = 4;
      for (int j = 0; j < 100; j++)
	{
	  if (max >= j)
	    {
	      h_pt_thresh[j]->Fill(b_truth_particle_pt->at(index));
	      h_ntower_thresh[j][hindex]->Fill(count);
	      h_ntower_thresh[j][0]->Fill(count);
	    }
	}
    }
  
  std::string outfilename = Form("/sphenix/user/dlis/Projects/CaloTriggerEmulator/plots/HIST_%s-%08d-%04d.root",trigger.c_str(), runnumber, segment);  

  TFile *fout = new TFile(outfilename.c_str(),"recreate");
  for (int i = 0; i < 100; i++)
    {
      h_pt_thresh[i]->Write();
      for (int j = 0; j < 5; j++)
	h_ntower_thresh[i][j]->Write();
    }
  fout->Close();

}

void Draw_Trigger_Hists(int runnumber, const std::string file, int segment, const std::string trigger)
{
  std::string filename = file;
  TFile *f = new TFile(filename.c_str(),"r");
  
  TTree *t = (TTree*) f->Get("ttree");
  if (!t)
    {
      std::cout << "no tree" << std::endl;
      return;
    }
  
  int scaledowns[64]={-1};
  get_scaledowns(runnumber, scaledowns);
  
  ULong64_t b_gl1_scaledvec;
  ULong64_t b_gl1_rawvec;
  ULong64_t b_gl1_clock;
 
  std::vector<short>* b_emcal_good= nullptr;
  std::vector<float>* b_emcal_energy= nullptr;
  std::vector<float>* b_emcal_time= nullptr;
  std::vector<float>* b_emcal_etabin= nullptr;
  std::vector<float>* b_emcal_phibin= nullptr;

  std::vector<short>* b_hcalin_good= nullptr;
  std::vector<float>* b_hcalin_energy= nullptr;
  std::vector<float>* b_hcalin_time= nullptr;
  std::vector<float>* b_hcalin_etabin= nullptr;
  std::vector<float>* b_hcalin_phibin= nullptr;

  std::vector<short>* b_hcalout_good= nullptr;
  std::vector<float>* b_hcalout_energy= nullptr;
  std::vector<float>* b_hcalout_time= nullptr;
  std::vector<float>* b_hcalout_etabin= nullptr;
  std::vector<float>* b_hcalout_phibin= nullptr;

  int b_cluster_n;
  std::vector<float>* b_cluster_prob = nullptr;
  std::vector<float>* b_cluster_chi2 = nullptr;
  std::vector<float>* b_cluster_ecore = nullptr;
  std::vector<float>* b_cluster_pt = nullptr;
  std::vector<float>* b_cluster_phi = nullptr;
  std::vector<float>* b_cluster_eta = nullptr;
  std::vector<float>* b_cluster_iso = nullptr;

  t->SetBranchAddress("gl1_clock",&b_gl1_clock);//gl1_clock/l");
  t->SetBranchAddress("gl1_scaledvec",&b_gl1_scaledvec);//gl1_scaledvec/l");  
  t->SetBranchAddress("gl1_rawvec",&b_gl1_rawvec);//gl1_scaledvec/l");

  t->SetBranchAddress("emcal_good",&b_emcal_good);
  t->SetBranchAddress("emcal_energy",&b_emcal_energy);
  t->SetBranchAddress("emcal_time",&b_emcal_time);
  t->SetBranchAddress("emcal_phibin",&b_emcal_phibin);
  t->SetBranchAddress("emcal_etabin",&b_emcal_etabin);
  t->SetBranchAddress("hcalin_good",&b_hcalin_good);
  t->SetBranchAddress("hcalin_energy",&b_hcalin_energy);
  t->SetBranchAddress("hcalin_time",&b_hcalin_time);
  t->SetBranchAddress("hcalin_phibin",&b_hcalin_phibin);
  t->SetBranchAddress("hcalin_etabin",&b_hcalin_etabin);
  t->SetBranchAddress("hcalout_good",&b_hcalout_good);
  t->SetBranchAddress("hcalout_energy",&b_hcalout_energy);
  t->SetBranchAddress("hcalout_time",&b_hcalout_time);
  t->SetBranchAddress("hcalout_phibin",&b_hcalout_phibin);
  t->SetBranchAddress("hcalout_etabin",&b_hcalout_etabin);
  t->SetBranchAddress("cluster_n",&b_cluster_n);
  t->SetBranchAddress("cluster_ecore",&b_cluster_ecore);
  t->SetBranchAddress("cluster_pt",&b_cluster_pt);
  t->SetBranchAddress("cluster_prob",&b_cluster_prob);
  t->SetBranchAddress("cluster_chi2",&b_cluster_chi2);
  t->SetBranchAddress("cluster_phi",&b_cluster_phi);
  t->SetBranchAddress("cluster_eta",&b_cluster_eta);
  t->SetBranchAddress("cluster_iso",&b_cluster_iso);

  // Historgram Time
  float energymap[12][32] = {0};
  float energymap_hcalin[12][35] = {0};
  float energymap_hcalout[12][35] = {0};
  float energymap_emcal[12][35] = {0};
  float energymap_jet[9][32] = {0};
  float energymap_jet_emcal[9][32] = {0};
  float energymap_jet_hcalin[9][32] = {0};
  float energymap_jet_hcalout[9][32] = {0};
  float energymap_extend[12][35] = {0};

  int goodmap[12][32] = {0};
  int goodmap_hcalin[12][35] = {0};
  int goodmap_hcalout[12][35] = {0};
  int goodmap_emcal[12][35] = {0};
  int goodmap_jet[9][32] = {0};
  int goodmap_jet_emcal[9][32] = {0};
  int goodmap_jet_hcalin[9][32] = {0};
  int goodmap_jet_hcalout[9][32] = {0};
  int goodmap_extend[12][35] = {0};

  int n_primitives = 384;
  int n_sums = 6144;
  int energy_bins = 100;
  float energy_max = 30.0;
  int turnon_bins = 40;
  float turnon_max = 10.0;
  int lut_bins = 256;


  // clock v sum
  TH1D *h_triggers = new TH1D("h_triggers", "", 64, -0.5, 63.5);
  TH1D *h_mbd_triggers = new TH1D("h_mbd_triggers", "", 64, -0.5, 63.5);
  TH1D *h_emcal_energy = new TH1D("h_emcal_energy", "; Energy;", 40, 0,  20);
  TH1D *h_cluster_energy = new TH1D("h_cluster_energy", "; Energy;", 40, 0,  20);
  TH1D *h_hcal_energy = new TH1D("h_hcal_energy", "; Energy;", 40, 0,  20);
  TH1D *h_jet_energy = new TH1D("h_hmcal_energy", "; Energy;", 50, 0,  50);
  TH1D *h_emcal_energy_ref = new TH1D("h_emcal_energy_ref", "; Energy;", 40, 0,20);
  TH1D *h_jet_energy_ref = new TH1D("h_jet_energy_ref", "; Energy;", 50, 0,  50);

  TH1D *h_cluster_phi = new TH1D("h_cluster_phi", "; #phi;", 64, -1*TMath::Pi(), TMath::Pi());
  TH1D *h_cluster_eta = new TH1D("h_cluster_eta", "; #eta;", 64, -1.2, 1.2);
  TH1D *h_cluster_iso = new TH1D("h_cluster_iso", "; Energy;", 100, -50, 50);

  TH1D *h_cluster_phi_gl1[64];
  TH1D *h_cluster_eta_gl1[64];
  TH1D *h_cluster_iso_gl1[64];

  TH1D *h_cluster_phi_chi2 = new TH1D("h_cluster_phi_chi2", "; #phi;", 64, -1*TMath::Pi(), TMath::Pi());
  TH1D *h_cluster_eta_chi2 = new TH1D("h_cluster_eta_chi2", "; #eta;", 64, -1.2, 1.2);
  TH1D *h_cluster_iso_chi2 = new TH1D("h_cluster_iso_chi2", "; Energy;", 100, -50, 50);

  TH1D *h_cluster_phi_gl1_chi2[64];
  TH1D *h_cluster_eta_gl1_chi2[64];
  TH1D *h_cluster_iso_gl1_chi2[64];

  TH1D *h_emcal_energy_gl1[64];
  TH1D *h_emcal_energy_raw_gl1[64];
  TH1D *h_cluster_energy_gl1[64];
  TH1D *h_cluster_energy_raw_gl1[64];
  TH1D *h_emcal_energy_scaled_gl1[64];
  TH1D *h_cluster_energy_scaled_gl1[64];
  TH1D *h_hcal_energy_gl1[64];  
  TH1D *h_jet_emcal_energy_gl1[64];
  TH1D *h_jet_hcalin_energy_gl1[64];
  TH1D *h_jet_hcalout_energy_gl1[64];
  TProfile *h_jet_emcal_fraction_gl1[64];  
  TProfile *h_jet_hcalin_fraction_gl1[64];  
  TProfile *h_jet_hcalout_fraction_gl1[64];  
  TH2F *h2_good = new TH2F("h2_good","; Good towers; energy", 65, -0.5, 64.5, 200, 0, 20);
  TH2F *h2_good_jet = new TH2F("h2_good_jet","; Good towers; energy", 64*16, -0.5, 64*16 + .5, 200, 0, 20);
  TH2F *h2_good_gl1[64];
  TH2F *h2_good_jet_gl1[64];
  TH1D *h_jet_energy_gl1[64];
  TH1D *h_jet_energy_scaled_gl1[64];
  TH1D *h_jet_energy_raw_gl1[64];
  TH1D *h_hcalout_imbalance[64];
  TH1D *h_emcal_imbalance[64];
  TProfile *h_emcal_spread = new TProfile("h","h", 96*256, -0.5, 96*256 - 0.5);
  for (int i = 0 ; i < 64; i++)
    {

      h_cluster_phi_gl1[i] = new TH1D(Form("h_cluster_phi_gl1_%d", i), "; #phi;", 64, -1*TMath::Pi(), TMath::Pi());
      h_cluster_eta_gl1[i] = new TH1D(Form("h_cluster_eta_gl1_%d", i), "; #eta;", 64, -1.2, 1.2);
      h_cluster_iso_gl1[i] = new TH1D(Form("h_cluster_iso_gl1_%d", i), "; Energy;", 100, -50, 50);
      
      h2_good_gl1[i] = new TH2F(Form("h2_good_gl1_%d", i),"; Good towers; energy", 65, -0.5, 64.5, 200, 0, 20);
      h2_good_jet_gl1[i] = new TH2F(Form("h2_good_jet_gl1_%d", i),"; Good towers; energy", 64*16, -0.5, 64*16 + .5, 40, 0, 20);
      h_jet_emcal_fraction_gl1[i] = new TProfile(Form("h_jet_emcal_fraction_gl1_%d", i), "; Energy; Fraction", 40, 0, 20);
      h_jet_hcalin_fraction_gl1[i] = new TProfile(Form("h_jet_hcalin_fraction_gl1_%d", i), "; Energy; Fraction", 40, 0, 20);
      h_jet_hcalout_fraction_gl1[i] = new TProfile(Form("h_jet_hcalout_fraction_gl1_%d", i), "; Energy; Fraction", 40, 0, 20);
      h_emcal_energy_gl1[i] = new TH1D(Form("h_emcal_energy_gl1_%d", i), "; Energy;", 40, 0,  20);
      h_emcal_energy_raw_gl1[i] = new TH1D(Form("h_emcal_energy_raw_gl1_%d", i), "; Energy;", 40, 0,  20);
      h_cluster_energy_gl1[i] = new TH1D(Form("h_cluster_energy_gl1_%d", i), "; Energy;", 40, 0,  20);
      h_cluster_energy_raw_gl1[i] = new TH1D(Form("h_cluster_energy_raw_gl1_%d", i), "; Energy;", 40, 0,  20);
      h_hcal_energy_gl1[i] = new TH1D(Form("h_hcal_energy_gl1_%d", i), "; Energy;", 40, 0,  20);

      h_hcalout_imbalance[i] = new TH1D(Form("h_hcalout_imbalance_%d", i), "; Energy;", 40, 0,  20);
      h_emcal_imbalance[i] = new TH1D(Form("h_emcal_imbalance_%d", i), "; Energy;", 40, 0,  20);
      h_jet_emcal_energy_gl1[i] = new TH1D(Form("h_jet_emcal_energy_gl1_%d", i), "; Energy;", 40, 0,  20);
      h_jet_hcalin_energy_gl1[i] = new TH1D(Form("h_jet_hcalin_energy_gl1_%d", i), "; Energy;", 40, 0,  20);
      h_jet_hcalout_energy_gl1[i] = new TH1D(Form("h_jet_hcalout_energy_gl1_%d", i), "; Energy;", 40, 0,  20);
      h_jet_energy_gl1[i] = new TH1D(Form("h_jet_energy_gl1_%d", i), "; Energy;", 50, 0,  50);
      h_jet_energy_raw_gl1[i] = new TH1D(Form("h_jet_energy_raw_gl1_%d", i), "; Energy;", 50, 0,  50);
    }

  int entries = t->GetEntries();
  std::cout << " running through " << entries << " events" << std::endl;
  for (int i = 0; i < entries; i++)
    {
      t->GetEntry(i);

      //      std::cout << "Event "<<i << " ...   "<<std::endl;

      std::vector<int> trig_bits{};
      std::vector<int> trig_bits_raw{};
      int raw_mbd = 0;
      //      ULong64_t triggervec  = 0;      
      for (unsigned int bit = 0; bit < 64; bit++)
	{
	  if (((b_gl1_scaledvec >> bit ) & 0x1U ) == 0x1U)
	    {	      
	      trig_bits.push_back(bit);
	    }
	  if (((b_gl1_rawvec >> bit ) & 0x1U ) == 0x1U)
	    {
	      if (bit == 8 || bit == 9) raw_mbd++;
	      h_triggers->Fill(bit); 
	      trig_bits_raw.push_back(bit);
	    }

	}
      // clusters here
      //      std::cout << std::hex << triggervec << std::endl;
      for (int j = 0; j < 35; j++)
	{
	  for (int k =0 ; k < 12; k++)
	    {
	      if (j < 32)
		{
		  energymap[k][j] = 0.0;
		  goodmap[k][j] = 0;
		}

	      energymap_hcalin[k][j] = 0.0;
	      energymap_hcalout[k][j] = 0.0;
	      energymap_emcal[k][j] = 0.0;
	      energymap_extend[k][j] = 0.0;
	      goodmap_extend[k][j] = 0;
	    }
	}      
      for (int ie = 0; ie < b_emcal_energy->size(); ie++)
	{
	  int ebin = b_emcal_etabin->at(ie)/8;
	  int pbin = b_emcal_phibin->at(ie)/8;
	  //	  if (ebin >= 6 && pbin == 27) continue;
	  goodmap[ebin][pbin] += b_emcal_good->at(ie);
	  goodmap_extend[ebin][pbin] += b_emcal_good->at(ie);
	  if (pbin < 3)
	    {
	      goodmap_extend[ebin][pbin+32] += b_emcal_good->at(ie);
	    }
	  
	  if (!b_emcal_good->at(ie)) continue;	  
	  energymap[ebin][pbin] += b_emcal_energy->at(ie);
	  energymap_emcal[ebin][pbin] += b_emcal_energy->at(ie);
	  energymap_extend[ebin][pbin] += b_emcal_energy->at(ie);
	  h_emcal_spread->Fill(ie, b_emcal_energy->at(ie));

	  if (pbin < 3)
	    {
	      energymap_emcal[ebin][pbin+32] += b_emcal_energy->at(ie);
	      energymap_extend[ebin][pbin+32] += b_emcal_energy->at(ie);
	    }
	}
      for (int ie = 0; ie < b_hcalin_energy->size(); ie++)
	{
	  int ebin = b_hcalin_etabin->at(ie)/2;
	  int pbin = b_hcalin_phibin->at(ie)/2;
	  if (pbin < 3)
	    {
	      goodmap_extend[ebin][pbin+32] += b_hcalin_good->at(ie);
	    }
	  goodmap_extend[ebin][pbin] += b_hcalin_good->at(ie);
	  if (!b_hcalin_good->at(ie)) continue;	  
	  energymap_extend[ebin][pbin] += b_hcalin_energy->at(ie);
	  energymap_hcalin[ebin][pbin] += b_hcalin_energy->at(ie);
	  if (pbin < 3)
	    {
	      energymap_extend[ebin][pbin+32] += b_hcalin_energy->at(ie);
	      energymap_hcalin[ebin][pbin+32] += b_hcalin_energy->at(ie);
	    }
	}
      for (int ie = 0; ie < b_hcalout_energy->size(); ie++)
	{
	  int ebin = b_hcalout_etabin->at(ie)/2;
	  int pbin = b_hcalout_phibin->at(ie)/2;
	  if (pbin < 3)
	    {
	      goodmap_extend[ebin][pbin+32] += b_hcalout_good->at(ie);
	    }
	  goodmap_extend[ebin][pbin] += b_hcalout_good->at(ie);
	  if (!b_hcalout_good->at(ie)) continue;	  
	  energymap_extend[ebin][pbin] += b_hcalout_energy->at(ie);
	  energymap_hcalout[ebin][pbin] += b_hcalout_energy->at(ie);
	  if (pbin < 3)
	    {
	      energymap_hcalout[ebin][pbin+32] += b_hcalout_energy->at(ie);
	      energymap_extend[ebin][pbin+32] += b_hcalout_energy->at(ie);
	    }
	}

      for (int ie = 0; ie< 9; ie++)
	{
	  for (int ip = 0 ; ip < 32; ip++)
	    {
	      energymap_jet[ie][ip] = 0.0;
	      energymap_jet_emcal[ie][ip] = 0.0;
	      energymap_jet_hcalin[ie][ip] = 0.0;
	      energymap_jet_hcalout[ie][ip] = 0.0;
	      goodmap_jet[ie][ip] = 0;
	      for (int is = 0; is < 16; is++)
		{
		  goodmap_jet[ie][ip] += goodmap_extend[ie + is%4][ip + is/4];
		  energymap_jet[ie][ip] += energymap_extend[ie + is%4][ip + is/4];
		  energymap_jet_emcal[ie][ip] += energymap_emcal[ie + is%4][ip + is/4];
		  energymap_jet_hcalin[ie][ip] += energymap_hcalin[ie + is%4][ip + is/4];
		  energymap_jet_hcalout[ie][ip] += energymap_hcalout[ie + is%4][ip + is/4];
		}
	    }
	}

      float max_energy = 0.0;
      float max_energy_clus = 0.0;
      float max_energy_hcal = 0.0;
      float max_energy_jet = 0.0;
      int jet_ebin = 0;
      int jet_pbin = 0;
      int ebin = 0;
      int pbin = 0;
      int hcal_ebin = 0;
      int hcal_pbin = 0;

      for (int j  = 0 ; j < b_cluster_n; j++)
      	{
      	  if (b_cluster_ecore->at(j) > max_energy_clus)
      	    {
      	      max_energy_clus = b_cluster_ecore->at(j);
      	    }
      	}

      //      std::cout << "max cluster " << max_energy_clus << " at " << b_cluster_n << std::endl;
      for (int j = 0; j < 32; j++)
	{
	  if (entries <=10) std::cout << j << " \t |";

	  for (int k =0 ; k < 12; k++)
	    {
	      if (k < 9)
		{
		  if (energymap_jet[k][j] > max_energy_jet)
		    {
		      max_energy_jet = energymap_jet[k][j];
		      jet_ebin = k;
		      jet_pbin = j;
		    }
		}
	      if (entries <= 10) std::cout << energymap[k][j] << "\t";
	      if (energymap[k][j] > max_energy)
		{
		  max_energy = energymap[k][j];
		  ebin = k;
		  pbin = j;
		}
	      if (energymap_hcalout[k][j] + energymap_hcalin[k][j] > max_energy_hcal)
		{
		  max_energy_hcal = energymap_hcalin[k][j] + energymap_hcalout[k][j];
		  hcal_ebin = k;
		  hcal_pbin = j;
		}

	    }

	  if (entries <=10) std::cout << " | " <<std::endl;;

	}      

      h_emcal_energy->Fill(max_energy);
      h_cluster_energy->Fill(max_energy_clus);
      h_hcal_energy->Fill(max_energy_hcal);
      h_jet_energy->Fill(max_energy_jet);
      h2_good->Fill(goodmap[ebin][pbin], max_energy);
      h2_good_jet->Fill(goodmap_jet[jet_ebin][jet_pbin], max_energy_jet);
      if (raw_mbd == 2)
	{

	  h_emcal_energy_ref->Fill(max_energy);
	  h_jet_energy_ref->Fill(max_energy_jet);
	  for (auto &b : trig_bits_raw)
	    {
	      h_mbd_triggers->Fill(b); 	      
	      h_emcal_energy_raw_gl1[b]->Fill(max_energy);
	      h_cluster_energy_raw_gl1[b]->Fill(max_energy_clus);
	      h_jet_energy_raw_gl1[b]->Fill(max_energy_jet);
	    }
	}

      for (auto &b : trig_bits)
	{
	  //	  std::cout << "    " << std::dec <<b << " " << max_energy << std::endl;
	  h2_good_gl1[b]->Fill(goodmap[ebin][pbin], max_energy);
	  h2_good_jet_gl1[b]->Fill(goodmap_jet[jet_ebin][jet_pbin], max_energy_jet);
	  //if ( goodmap[ebin][pbin] == 64)
	  h_cluster_energy_gl1[b]->Fill(max_energy_clus);
	  h_emcal_energy_gl1[b]->Fill(max_energy);
	  h_hcal_energy_gl1[b]->Fill(max_energy_hcal);
	  h_jet_energy_gl1[b]->Fill(max_energy_jet);
	  h_jet_emcal_energy_gl1[b]->Fill(energymap_jet_emcal[jet_ebin][jet_pbin]);
	  h_jet_hcalin_energy_gl1[b]->Fill(energymap_jet_hcalin[jet_ebin][jet_pbin]);
	  h_jet_hcalout_energy_gl1[b]->Fill(energymap_jet_hcalout[jet_ebin][jet_pbin]);
	  h_jet_emcal_fraction_gl1[b]->Fill(max_energy_jet, energymap_jet_emcal[jet_ebin][jet_pbin]/max_energy_jet);
	  h_jet_hcalin_fraction_gl1[b]->Fill(max_energy_jet, energymap_jet_hcalin[jet_ebin][jet_pbin]/max_energy_jet);
	  h_jet_hcalout_fraction_gl1[b]->Fill(max_energy_jet, energymap_jet_hcalout[jet_ebin][jet_pbin]/max_energy_jet);
	  if (energymap_jet_emcal[jet_ebin][jet_pbin] < 1)
	    {
	      h_hcalout_imbalance[b]->Fill(energymap_jet_hcalout[jet_ebin][jet_pbin]);
	    }
	  if (energymap_jet_hcalout[jet_ebin][jet_pbin] < 1)
	    {
	      h_emcal_imbalance[b]->Fill(energymap_jet_emcal[jet_ebin][jet_pbin]);
	    }
	}

    }

  //  h_emcal_energy_gl1[10]->Scale(scaledown);  
  //h_jet_energy_gl1[10]->Scale(scaledown);  
  //h_jet_emcal_energy_gl1[10]->Scale(scaledown);  
  //h_jet_hcalin_energy_gl1[10]->Scale(scaledown);  
  //h_jet_hcalout_energy_gl1[10]->Scale(scaledown);  
  std::string outfilename = Form("/sphenix/user/dlis/Projects/CaloTriggerEmulator/plots/HIST_%s-%08d-%04d.root",trigger.c_str(), runnumber, segment);  

  TFile *fout = new TFile(outfilename.c_str(),"recreate");
  h_mbd_triggers->Write();
  h_triggers->Write();
  std::string name;
  h_emcal_spread->Write();

  for (int i = 0; i < 64; i++)
    {
      h2_good_gl1[i]->Write();
      h2_good_jet_gl1[i]->Write();
      h_emcal_energy_gl1[i]->Write();
      h_emcal_energy_raw_gl1[i]->Write();
      h_emcal_energy_scaled_gl1[i] = (TH1D*) h_emcal_energy_gl1[i]->Clone();
      name = "h_emcal_energy_scaled_gl1_" + to_string(i);
      h_emcal_energy_scaled_gl1[i]->SetName(name.c_str());
      h_emcal_energy_scaled_gl1[i]->Scale(scaledowns[i] + 1.);
      h_emcal_energy_scaled_gl1[i]->Write();

      h_cluster_energy_gl1[i]->Write();
      h_cluster_energy_raw_gl1[i]->Write();
      h_cluster_energy_scaled_gl1[i] = (TH1D*) h_cluster_energy_gl1[i]->Clone();
      name = "h_cluster_energy_scaled_gl1_" + to_string(i);
      h_cluster_energy_scaled_gl1[i]->SetName(name.c_str());
      h_cluster_energy_scaled_gl1[i]->Scale(scaledowns[i] + 1.);
      h_cluster_energy_scaled_gl1[i]->Write();

      h_hcal_energy_gl1[i]->Write();
      h_jet_energy_gl1[i]->Write();
      h_jet_energy_raw_gl1[i]->Write();
      h_jet_energy_scaled_gl1[i] = (TH1D*) h_jet_energy_gl1[i]->Clone();
      name = "h_jet_energy_scaled_gl1_" + to_string(i);
      h_jet_energy_scaled_gl1[i]->SetName(name.c_str());
      h_jet_energy_scaled_gl1[i]->Scale(scaledowns[i] + 1.);
      h_jet_energy_scaled_gl1[i]->Write();

      h_jet_emcal_energy_gl1[i]->Write();
      h_jet_hcalin_energy_gl1[i]->Write();
      h_jet_hcalout_energy_gl1[i]->Write();
      h_jet_emcal_fraction_gl1[i]->Write();
      h_jet_hcalin_fraction_gl1[i]->Write();
      h_jet_hcalout_fraction_gl1[i]->Write();
      h_hcalout_imbalance[i]->Write();
      h_emcal_imbalance[i]->Write();
    }
  h2_good->Write();
  h2_good_jet->Write();
  h_emcal_energy_ref->Write();
  h_jet_energy_ref->Write();
  h_emcal_energy->Write();
  h_cluster_energy->Write();
  h_jet_energy->Write();
  h_hcal_energy->Write();
  fout->Close();
  
  std::cout <<" Done." << std::endl;
  return;
}

void Draw_TriggerQA(int runnumber, const std::string file ,int segment, const std::string trigger)
{

  Draw_Trigger_Hists_Sim(runnumber, file, segment, trigger);
  std::cout << "DONE" <<std::endl;
  return;
}
