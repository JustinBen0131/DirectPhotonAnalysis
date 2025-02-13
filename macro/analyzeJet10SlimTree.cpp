#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TLatex.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <sys/stat.h>
#include <sys/types.h>
#include <map>
#include <fstream>

static const int kMaxClusters   = 2000000;
static const int kMaxParticles  = 2000000;

struct ClusterRecord {
    float cPt;
    float isoAll;
    bool  isTruePrompt;  // cat=0 or not
    int   ptBinIndex;    // bPur => which pT bin in [5..60]
};

static std::vector<ClusterRecord> gAllClusters; // global to store them

// A debug toggle to control how much we print:
static const bool kDebug = false;   // set to 'true' to see extra debug lines

// -------------------------------------------------------------------
// Existing arrays/variables for the iso-efficiency study
// -------------------------------------------------------------------
int   ncluster;
int   cluster_pid[kMaxClusters];
float cluster_Et[kMaxClusters];
float cluster_E[kMaxClusters];
float cluster_Eta[kMaxClusters];
float cluster_Phi[kMaxClusters];
float cluster_iso_04[kMaxClusters];
float cluster_iso_04_emcal[kMaxClusters];
static float cluster_prob[kMaxClusters];
static float cluster_MergedProb[kMaxClusters];
int   particle_photon_mother_pid[kMaxParticles];
static const int kMaxBinsPur = 50;
static const int kNumExpandedCats = 10;
static long long catCountExpanded[kNumExpandedCats][kMaxBinsPur];
static long long catCount[6][kMaxBinsPur];
static int particle_converted[kMaxParticles];
float cluster_iso_04_hcalin[kMaxClusters];
float cluster_iso_04_hcalout[kMaxClusters];

int   nparticles;
int   ndaughters;
static float vertexz_truth;
int   particle_pid[kMaxParticles];
int   daughter_pid[kMaxParticles];
int   particle_photonclass[kMaxParticles];
float particle_Pt[kMaxParticles];
float particle_Eta_[kMaxParticles];
float particle_Phi_[kMaxParticles];
float particle_truth_iso_04[kMaxParticles];
float particle_E[kMaxParticles];      // Energy for each particle

// For storing kinematic information for daughter particles:
float daughter_E[kMaxParticles];        // Energy for each daughter
float daughter_Pt[kMaxParticles];       // Transverse momentum for each daughter
float daughter_Phi_[kMaxParticles];     // Phi for each daughter

// If daughter_vtx_z is a single value per event, declare it as a float:
float daughter_vtx_z;

static std::map<int, long long> otherPidCounts;
// For storing the isoAll distribution of truly prompt clusters in each pT bin
static TH1F* hIsoPromptInBin[kMaxBinsPur] = {nullptr};

// -------------------------------------------------------------------
static long long gCutFlow_totalPairs = 0;      // total pairs considered
static long long gCutFlow_skipVertexZ = 0;     // entire event had |z|>30 => skip all pairs
static long long gCutFlow_skipProb1 = 0;       // cluster1 prob < 0.05
static long long gCutFlow_skipProb2 = 0;       // cluster2 prob < 0.05
static long long gCutFlow_skipPt1neg = 0;      // cluster1 had negative pt/E
static long long gCutFlow_skipPt2neg = 0;      // cluster2 had negative pt/E
static long long gCutFlow_skipMesonPt = 0;     // meson pT < 1.0
static long long gCutFlow_skipAsym = 0;        // asym>0.7
static long long gCutFlow_skipDR = 0;          // dR not in (1/(10pT),1/pT)
static long long gCutFlow_skipLeadPt = 0;      // leading pT outside [5..40]
static long long gCutFlow_passAll = 0;         // final accepted pairs

// -------------------------------------------------------------------
// pT bin edges and iso-thresholds for the iso-efficiency
// -------------------------------------------------------------------
static std::vector<float> pTedges = {5, 6, 8, 10, 12, 15, 20, 25, 30, 40};
static int nBins = 0;
static std::vector<double> isoCutsGeV = {
    0.0, 1.0, 2.0, 3.0, 4.0,
    6.0, 8.0, 10.0, 12.0, 15.0, 20.0
};
static int nIsoCuts = 0;

// isoCountByClass[classIdx][pTbin][isoCutIndex]
static long long isoCountByClass[2][10][50];
static long long nPhotonsInBin[2][10];

// -------------------------------------------------------------------
// pT bin edges for the prompt-photon purity measurement
// -------------------------------------------------------------------
static std::vector<float> pTedges_purEff = {5, 6, 8, 10, 12, 15, 20, 25, 30, 40};
static int nBins_pur = 0;
static long long nFakeCount[6][kMaxBinsPur]; // number of clusters that pass "prompt" cuts but are cat=3,4,5
static long long nTagAll[kMaxBinsPur];
static long long nTagCorrect[kMaxBinsPur];
// -- FRAGMENTATION PHOTON PURITY ARRAYS:
static long long nTagAllFrag[kMaxBinsPur];
static long long nTagCorrectFrag[kMaxBinsPur];
static long long nTagAllAlt[3][kMaxBinsPur];
static long long nTagCorrAlt[3][kMaxBinsPur];

static long long nTagAll_varIso[kMaxBinsPur];
static long long nTagCorr_varIso[kMaxBinsPur];


static bool passIsoCombo(float isoAll, float isoEM, int iCombo)
{
    // iCombo=0 => isoAll<6 only
    // iCombo=1 => isoAll<4, isoEM<2
    // iCombo=2 => isoAll<6, isoEM<2
    switch(iCombo){
      case 0: // isoAll<6 only
        if(isoAll>6.f) return false;
        return true; // ignoring isoEM
      case 1: // isoAll<4, isoEM<2
        if(isoAll>4.f) return false;
        if(isoEM>2.f)  return false;
        return true;
      case 2: // isoAll<6, isoEM<2
        if(isoAll>6.f) return false;
        if(isoEM>2.f)  return false;
        return true;
    }
    return false;
}

// -------------------------------------------------------------------
int photonClassIndex(int photClass)
{
    if(photClass == 1) return 0;  // prompt => index 0
    if(photClass == 2) return 1;  // frag   => index 1
    return -1;                    // skip decay/other
}

// -------------------------------------------------------------------
void ensureOutputDirectory(const std::string & outDir)
{
    struct stat info;
    if(stat(outDir.c_str(), &info)!=0){
        std::cout << "[INFO] Creating directory: " << outDir << std::endl;
        std::string cmd = "mkdir -p " + outDir;
        system(cmd.c_str());
    } else {
        std::cout << "[INFO] Output directory already exists: " << outDir << std::endl;
    }
}

// -------------------------------------------------------------------
int findPtBin(float pt)
{
    for(int i=0; i<(int)pTedges.size()-1; i++){
        if(pt >= pTedges[i] && pt < pTedges[i+1]) return i;
    }
    return -1;
}

// -------------------------------------------------------------------
float deltaR_func(float eta1, float eta2, float phi1, float phi2)
{
    float dEta = eta1 - eta2;
    float dPhi = phi1 - phi2;
    while(dPhi >  TMath::Pi()) dPhi -= 2.f*TMath::Pi();
    while(dPhi < -TMath::Pi()) dPhi += 2.f*TMath::Pi();
    return std::sqrt(dEta*dEta + dPhi*dPhi);
}

// -------------------------------------------------------------------
bool passClusterPromptCuts(float isoAll, float isoEM)
{
    // The normal "prompt candidate" cuts used in purity
    if(isoAll > 6.0) return false;
    if(isoEM  > 4.0) return false;
    return true;
}

// -------------------------------------------------------------------
int findPurityBin(float cPt)
{
    for(int ib=0; ib<nBins_pur; ib++){
        if(cPt >= pTedges_purEff[ib] && cPt < pTedges_purEff[ib+1]){
            return ib;
        }
    }
    return -1;
}

int findBestMatchAnyParticle(float cEta, float cPhi, float cE)
{
    const float maxDR = 0.05;
    float bestDR = 9999.f;
    int   bestIdx= -1;

    for(int ip=0; ip<nparticles; ip++){
        // We can match to ANY pid.  If you only want stable hadrons or something,
        // you could skip e.g. short-lived stuff here. We'll allow all for now.

        float dr = deltaR_func(cEta, particle_Eta_[ip], cPhi, particle_Phi_[ip]);
        if(dr > maxDR) continue; // must be within 0.05

        float eTruth = particle_Pt[ip] * std::cosh(particle_Eta_[ip]);
        if(eTruth < 1e-3) continue;

        float ratio = cE / eTruth;
        if(ratio < 0.5f) continue; // want cluster to carry at least half the truth

        if(dr < bestDR){
            bestDR  = dr;
            bestIdx = ip;
        }
    }

    return bestIdx;
}



static void makeSlicesAndOverlays_3cat(
    TH2F &h2Prompt,
    TH2F &h2Frag,
    TH2F &h2Meson,            // NEW
    const std::string & outFolderLabel
)
{
    // We can keep the same slice boundaries or refine them
    std::vector<std::pair<double,double>> truthIsoSlices = {
       {0.0,0.5}, {0.5,1.0}, {1.0,2.0}, {2.0,3.0}, {3.0,4.0}, {4.0,5.0}, {5.0,8.0}, {8.0,12.0}
    };

    int sliceIndex = 0;
    for(const auto &slice : truthIsoSlices)
    {
        double xLo = slice.first;
        double xHi = slice.second;

        int binLo = h2Prompt.GetXaxis()->FindBin(xLo + 1e-6);
        int binHi = h2Prompt.GetXaxis()->FindBin(xHi - 1e-6);

        // Project Y for prompt
        TH1D *hPromptProj = h2Prompt.ProjectionY(
            Form("hPromptProj3_%d", sliceIndex),
            binLo, binHi
        );
        // Project Y for frag
        TH1D *hFragProj = h2Frag.ProjectionY(
            Form("hFragProj3_%d", sliceIndex),
            binLo, binHi
        );
        // Project Y for meson
        TH1D *hMesonProj = h2Meson.ProjectionY(
            Form("hMesonProj3_%d", sliceIndex),
            binLo, binHi
        );

        // Normalize each
        double intP = hPromptProj->Integral();
        double intF = hFragProj->Integral();
        double intM = hMesonProj->Integral();
        if(intP>0) hPromptProj->Scale(1./intP);
        if(intF>0) hFragProj->Scale(1./intF);
        if(intM>0) hMesonProj->Scale(1./intM);

        // Canvas
        TCanvas cSlice(
          Form("cSlice3_%s_%d", outFolderLabel.c_str(), sliceIndex),
          Form("Reco-iso slice (truth iso=%.1f..%.1f)", xLo, xHi),
          800,600
        );
        cSlice.SetGrid();

        // line colors
        hPromptProj->SetLineColor(kRed);
        hPromptProj->SetLineWidth(2);

        hFragProj->SetLineColor(kBlue);
        hFragProj->SetLineWidth(2);

        hMesonProj->SetLineColor(kGreen+2);
        hMesonProj->SetLineWidth(2);

        // Title
        hPromptProj->SetTitle(
            Form("Reco iso slices (Prompt/Frag/Meson); iso_{reco} [GeV]; Normalized entries (truth iso=%.1f..%.1f)",
                 xLo, xHi)
        );

        // find max
        double maxY = std::max( {hPromptProj->GetMaximum(),
                                 hFragProj->GetMaximum(),
                                 hMesonProj->GetMaximum()} );
        hPromptProj->SetMaximum(1.2 * maxY);

        // Draw
        hPromptProj->Draw("HIST");
        hFragProj->Draw("HIST SAME");
        hMesonProj->Draw("HIST SAME");

        // Legend
        TLegend leg(0.65,0.70,0.9,0.9);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.AddEntry(hPromptProj, "Prompt (#gamma_{class}=1)", "l");
        leg.AddEntry(hFragProj,   "Frag (#gamma_{class}=2)",   "l");
        leg.AddEntry(hMesonProj,  "Meson (#pi^{0}/#eta)",       "l");
        leg.Draw();

        // A TLatex label for the slice
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.SetTextFont(42);
        latex.DrawLatex(0.2, 0.85,
             Form("#bf{Truth iso range: %.1f - %.1f GeV}", xLo, xHi)
        );

        // Save
        std::string outName = Form("%s/RecoIso_1Doverlay_3cat_slice%d_%.1fto%.1f.png",
                                   outFolderLabel.c_str(),
                                   sliceIndex, xLo, xHi);
        cSlice.SaveAs(outName.c_str());
        std::cout << "[INFO] Wrote => " << outName << "\n";

        ++sliceIndex;
    }
}


// -------------------------------------------------------------------
// The main function that saves 2D + 1D slices to subfolders
// -------------------------------------------------------------------
void buildAndSaveTruthVsRecoIsoPlots(
    const std::string & outDir,
    TH2F & h2TruthVsReco_prompt,
    TH2F & h2TruthVsReco_frag,
    TH2F & h2TruthVsReco_meson)
{
    gStyle->SetOptStat(0);

    //
    // A) "totalIsolation" subfolder
    //
    {
        // 1) Make sure the subdirectory exists
        std::string totalIsoDir = outDir + "/totalIsolation";
        gSystem->mkdir(totalIsoDir.c_str(), /*mode*/true);  // in case it doesn't exist

        // 2) Draw the 2D with profile for prompt
        {
            TCanvas cPromptProfile("cPromptProfile","Prompt: 2D + ProfileX",800,600);
            cPromptProfile.SetGrid();
            cPromptProfile.SetLogz();
            
            h2TruthVsReco_prompt.Draw("COLZ");
            TProfile *profPrompt = h2TruthVsReco_prompt.ProfileX("profPrompt", 1, -1, "");
            profPrompt->SetLineColor(kMagenta+1);
            profPrompt->SetLineWidth(2);
            profPrompt->Draw("SAME");
            std::string outNameProfPrompt = totalIsoDir + "/TruthVsRecoIso_Prompt_WithProfile.png";
            cPromptProfile.SaveAs(outNameProfPrompt.c_str());
            std::cout << "[INFO] Wrote " << outNameProfPrompt << std::endl;
        }
        // 3) Draw the 2D with profile for frag
        {
            TCanvas cFragProfile("cFragProfile","Frag: 2D + ProfileX",800,600);
            cFragProfile.SetGrid();
            cFragProfile.SetLogz();
            h2TruthVsReco_frag.Draw("COLZ");
            TProfile *profFrag = h2TruthVsReco_frag.ProfileX("profFrag", 1, -1, "");
            profFrag->SetLineColor(kMagenta+1);
            profFrag->SetLineWidth(2);
            profFrag->Draw("SAME");
            std::string outNameProfFrag = totalIsoDir + "/TruthVsRecoIso_Frag_WithProfile.png";
            cFragProfile.SaveAs(outNameProfFrag.c_str());
            std::cout << "[INFO] Wrote " << outNameProfFrag << std::endl;
        }
        
        {
            TCanvas cMesonProfile("cMesonProfile","Meson decay: 2D + ProfileX",800,600);
            cMesonProfile.SetGrid();
            cMesonProfile.SetLogz();
            h2TruthVsReco_meson.Draw("COLZ");
            TProfile* profM = h2TruthVsReco_meson.ProfileX("profMeson",1,-1,"");
            profM->SetLineColor(kMagenta+1);
            profM->SetLineWidth(2);
            profM->Draw("SAME");
            std::string outName = totalIsoDir + "/TruthVsRecoIso_Meson_WithProfile.png";
            cMesonProfile.SaveAs(outName.c_str());
        }
        
        // 4) Now do the slice‐and‐overlay in the same folder
        //    i.e. pass "totalIsoDir" to makeSlicesAndOverlays
        gSystem->cd(totalIsoDir.c_str()); // so the 1D overlays also go here
        makeSlicesAndOverlays_3cat(h2TruthVsReco_prompt, h2TruthVsReco_frag, h2TruthVsReco_meson, totalIsoDir);
    }

    std::cout << "[INFO] Done building & saving truth-vs-reco iso plots + 1D slices.\n";
}


// -------------------------------------------------------------------
// findMatchedPhotonClass => for a cluster with leading PID=22, check
// all MC photons with PID=22. If ∆R<0.05 and E_ratio>0.5, we pick the
// best match => return that photon’s class in [1..4]. If no match => -1
// -------------------------------------------------------------------
int findMatchedPhotonClass(float cEta, float cPhi, float cE)
{
    float bestDR= 9999.f;
    int   bestIdx=-1;
    const float maxDR= 0.05;
    for(int ip=0; ip<nparticles; ip++){
        if(particle_pid[ip] != 22) continue; // must be gamma
        float dr = deltaR_func(cEta, particle_Eta_[ip], cPhi, particle_Phi_[ip]);
        if(dr> maxDR) continue;

        float eTruth= particle_Pt[ip]* std::cosh(particle_Eta_[ip]);
        if(eTruth<1e-3) continue;

        float ratio= cE/eTruth;
        if(ratio< 0.5f) continue;

        if(dr< bestDR){
            bestDR=dr;
            bestIdx= ip;
        }
    }
    if(bestIdx<0) return -1; // no match
    return particle_photonclass[bestIdx]; // [1..4]
}

// -------------------------------------------------------------------
int findMatchedPromptPhoton(float cEta, float cPhi, float cE)
{
    // used for purity => must be class=1
    int phClass= findMatchedPhotonClass(cEta,cPhi,cE);
    if(phClass==1) return 12345; // matched => positive sentinel
    return -1; // no match or not class=1
}

int findMatchedFragPhoton(float cEta, float cPhi, float cE)
{
    // Identify a truth photon with class=2 (frag)
    int phClass = findMatchedPhotonClass(cEta, cPhi, cE);
    if(phClass == 2) return 12345; // matched => sentinel
    return -1; // no match or not class=2
}

int mapOtherPidToSubcat(int pidMatch)
{
    // If we had "no match" => stored sentinel 999999
    if(pidMatch == 999999 || pidMatch < 0){
        // negative means no matched index, etc.
        return 5; // "No-match"
    }

    // absolute value for easier grouping
    int ab = std::abs(pidMatch);

    // pions
    if(ab == 211) return 6;

    // kaons: ±321, ±311, maybe also 310, 130, etc.
    //        you can choose to expand or contract as needed
    if(ab == 321 || ab == 311 || ab == 310 || ab == 130) return 7;

    // nucleons: ±2212 (protons), ±2112 (neutrons)
    if(ab == 2212 || ab == 2112) return 8;

    // otherwise catch‐all
    return 9;
}

void buildAndSaveIsoEfficiencyPlot_2Classes(const std::string & outDir)
{
    // nBins is computed as pTedges.size()-1.
    // With: static std::vector<float> pTedges = {5, 7, 10, 15, 20, 30, 40};
    // we have 6 pT bins.
    std::cout << "[INFO] Building iso-efficiency per pT bin overlaying prompt and frag curves (6 plots)...\n";

    // Fixed colors: prompt always red, frag always blue.
    Color_t promptColor = kRed+1;
    Color_t fragColor   = kBlue+1;

    // Loop over each pT bin.
    for (int b = 0; b < nBins; b++) {
        // Check if there is any data for either prompt or frag in this bin.
        bool promptExists = (nPhotonsInBin[0][b] > 0);
        bool fragExists   = (nPhotonsInBin[1][b] > 0);
        if (!promptExists && !fragExists) continue; // Skip if no data.

        // Create a canvas for this pT bin.
        TCanvas cBin("cBin", Form("Iso Efficiency: %.0f < p_{T} < %.0f GeV", pTedges[b], pTedges[b+1]), 800,600);
        cBin.SetGrid();
        // Change here: set y-axis lower bound to 0.0 (instead of 0.93)
        TH1F* frameBin = cBin.DrawFrame(0., 0.0, 20., 1.02);
        frameBin->SetXTitle("Truth isoE_{T} cut [GeV]");
        frameBin->SetYTitle("Efficiency");
        frameBin->SetTitle(Form("%.0f < p_{T} < %.0f GeV", pTedges[b], pTedges[b+1]));

        // --- Prompt curve (class index 0) ---
        if (promptExists) {
            std::vector<double> vx(nIsoCuts, 0.);
            std::vector<double> vy(nIsoCuts, 0.);
            std::vector<double> vex(nIsoCuts, 0.);
            std::vector<double> vey(nIsoCuts, 0.);
            long long denom = nPhotonsInBin[0][b];
            for (int ic = 0; ic < nIsoCuts; ic++) {
                vx[ic] = isoCutsGeV[ic];
                long long num = isoCountByClass[0][b][ic];
                double eff = 0., err = 0.;
                if (denom > 0) {
                    eff = double(num) / double(denom);
                    err = std::sqrt(eff * (1.-eff) / double(denom));
                }
                vy[ic]  = eff;
                vey[ic] = err;
            }
            TGraphErrors* grPrompt = new TGraphErrors(nIsoCuts, &vx[0], &vy[0],
                                                      &vex[0], &vey[0]);
            grPrompt->SetLineColor(promptColor);
            grPrompt->SetMarkerColor(promptColor);
            grPrompt->SetMarkerStyle(20); // circle marker for prompt
            grPrompt->SetMarkerSize(1.3);
            grPrompt->Draw("P SAME");
        }

        // --- Frag curve (class index 1) ---
        if (fragExists) {
            std::vector<double> vx(nIsoCuts, 0.);
            std::vector<double> vy(nIsoCuts, 0.);
            std::vector<double> vex(nIsoCuts, 0.);
            std::vector<double> vey(nIsoCuts, 0.);
            long long denom = nPhotonsInBin[1][b];
            for (int ic = 0; ic < nIsoCuts; ic++) {
                vx[ic] = isoCutsGeV[ic];
                long long num = isoCountByClass[1][b][ic];
                double eff = 0., err = 0.;
                if (denom > 0) {
                    eff = double(num) / double(denom);
                    err = std::sqrt(eff * (1.-eff) / double(denom));
                }
                vy[ic]  = eff;
                vey[ic] = err;
            }
            TGraphErrors* grFrag = new TGraphErrors(nIsoCuts, &vx[0], &vy[0],
                                                    &vex[0], &vey[0]);
            grFrag->SetLineColor(fragColor);
            grFrag->SetMarkerColor(fragColor);
            grFrag->SetMarkerStyle(21); // square marker for frag
            grFrag->SetMarkerSize(1.3);
            grFrag->Draw("P SAME");
        }

        // --- Legend ---
        TLegend leg(0.6, 0.15, 0.85, 0.50);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        TGraph* dummyPrompt = new TGraph();
        dummyPrompt->SetMarkerStyle(20);
        dummyPrompt->SetMarkerColor(promptColor);
        leg.AddEntry(dummyPrompt, "Prompt (#gamma_{class}=1)", "p");
        TGraph* dummyFrag = new TGraph();
        dummyFrag->SetMarkerStyle(21);
        dummyFrag->SetMarkerColor(fragColor);
        leg.AddEntry(dummyFrag, "Frag (#gamma_{class}=2)", "p");
        leg.Draw();

        // Save the canvas.
        std::string fileName = outDir + Form("/IsoEfficiency_pTBin_%d_%.0f-%.0f.png",
                                             b, pTedges[b], pTedges[b+1]);
        cBin.SaveAs(fileName.c_str());
        std::cout << "[INFO] Wrote: " << fileName << std::endl;
    }
}

void buildExpandedCategoryPlots(const std::string & outDir)
{
    // 1) Summation across all pT bins => total pT‐independent
    long long catTotalExp[kNumExpandedCats];
    std::memset(catTotalExp, 0, sizeof(catTotalExp));
    long long totalClustersExp = 0;

    for(int cat=0; cat<kNumExpandedCats; cat++){
        for(int ib=0; ib<nBins_pur; ib++){
            catTotalExp[cat] += catCountExpanded[cat][ib];
        }
        totalClustersExp += catTotalExp[cat];
    }

    if(kDebug){
        std::cout << "[DEBUG] Summation of catCountExpanded across pT bins:\n";
        for(int cat=0; cat<kNumExpandedCats; cat++){
            std::cout << "   cat=" << cat
                      << " => total=" << catTotalExp[cat] << "\n";
        }
        std::cout << "   totalClustersExp=" << totalClustersExp << "\n";
    }

    // 2) Define user-friendly labels & colors for the 10 categories
    const char* catLabelsExp[kNumExpandedCats] = {
      "Prompt #gamma",        // cat=0
      "Frag #gamma",          // cat=1
      "Decay #gamma",         // cat=2
      "#pi^{0}-led",          // cat=3
      "#eta-led",             // cat=4
      "nucleon-led",      // cat=5
      "#pm pion-led",      // cat=6
      "kaon-led",      // cat=7
      "no-match",   // cat=8
      "remainder"      // cat=9
    };

    Color_t catColorsExp[kNumExpandedCats] = {
      kGreen+2,   // prompt
      kBlue+1,    // frag
      kMagenta,   // decay
      kOrange+1,  // pi0
      kAzure+1,   // eta
      kGray+1,    // no-match
      kRed+1,     // pions
      kYellow+2,  // kaons
      kSpring+5,  // nucleons
      kCyan+2     // leftover
    };

    // 3) Make the single bar chart (pT‐integrated)
    TH1F* hBarExp = new TH1F("hBarExp",
        "Cluster Categories (Expanded) pT-indep",
        kNumExpandedCats, 0.5, kNumExpandedCats+0.5 );

    hBarExp->SetYTitle("Fraction of total clusters");
    hBarExp->SetStats(0);
    hBarExp->GetXaxis()->SetLabelSize(0.05);
    hBarExp->GetXaxis()->SetTitleSize(0.05);

    // Fill the bin contents for each category
    double maxFracExp=0.;
    for(int cat=0; cat<kNumExpandedCats; cat++){
        double frac=0.;
        if(totalClustersExp>0){
            frac = double(catTotalExp[cat]) / double(totalClustersExp);
        }
        hBarExp->SetBinContent(cat+1, frac);
        if(frac>maxFracExp) maxFracExp= frac;
        hBarExp->GetXaxis()->SetBinLabel(cat+1, catLabelsExp[cat]);
    }

    // 4) Draw them as colored TBoxes
    TCanvas cCatExp("cCatExp","Expanded Categories pT-indep",800,600);
    cCatExp.SetGrid();

    double yMaxExp = std::min(1.0, 1.1*maxFracExp);
    if(yMaxExp < 0.2) yMaxExp=0.2;
    hBarExp->SetMaximum(yMaxExp);

    // Just draw the axis
    hBarExp->SetFillStyle(0);
    hBarExp->SetLineColor(kBlack);
    hBarExp->Draw("hist");

    // Over‐draw boxes for color fill
    double barWidth=1.0;
    for(int cat=0; cat<kNumExpandedCats; cat++){
        double frac = hBarExp->GetBinContent(cat+1);
        double xCenter = hBarExp->GetBinCenter(cat+1);
        double xLow  = xCenter - 0.5*barWidth;
        double xHigh = xCenter + 0.5*barWidth;

        TBox* box = new TBox(xLow, 0., xHigh, frac);
        box->SetFillColor( catColorsExp[cat] );
        box->SetLineColor( catColorsExp[cat] );
        box->SetLineWidth(0);
        box->Draw("same");
    }

    // 5) Save
    std::string outCatExp = outDir + "/ClusterCategories_Expanded_SingleBar.png";
    cCatExp.SaveAs(outCatExp.c_str());
    std::cout << "[INFO] Wrote single-bar (expanded) category plot => " << outCatExp << std::endl;

    // 6) Print numeric results
    std::cout << "\n===== pT-Independent Expanded Category Breakdown =====\n";
    std::cout << "  total clusters across all bins="<< totalClustersExp <<"\n";
    for(int cat=0; cat<kNumExpandedCats; cat++){
        double frac = (totalClustersExp > 0
                       ? double(catTotalExp[cat])/double(totalClustersExp)
                       : 0.);
        std::cout << "    cat=" << cat <<" ("<< catLabelsExp[cat] <<") => "
                  << catTotalExp[cat] <<" => fraction="<< frac <<"\n";
    }
    std::cout << "===============================================\n\n";

    // 7) Produce separate bar charts for each pT bin (just as before)
    for(int ib=0; ib<nBins_pur; ib++){
        // Sum how many clusters in this bin (any cat)
        long long binTotalExp = 0;
        for(int cat=0; cat<kNumExpandedCats; cat++){
            binTotalExp += catCountExpanded[cat][ib];
        }

        // Book a TH1F with 10 bins
        TH1F* hBarBinExp = new TH1F(Form("hBarExp_bin%d", ib),
            Form("Exp. Categories: %.0f #leq p_{T} < %.0f GeV",
                 pTedges_purEff[ib], pTedges_purEff[ib+1]),
            kNumExpandedCats, 0.5, kNumExpandedCats+0.5);

        hBarBinExp->SetYTitle("Fraction of clusters in category");
        hBarBinExp->SetStats(0);
        hBarBinExp->GetXaxis()->SetLabelSize(0.05);
        hBarBinExp->GetXaxis()->SetTitleSize(0.05);

        double maxFracBinExp=0.;
        for(int cat=0; cat<kNumExpandedCats; cat++){
            double frac=0.;
            if(binTotalExp>0){
                frac= double(catCountExpanded[cat][ib]) / double(binTotalExp);
            }
            hBarBinExp->SetBinContent(cat+1, frac);
            if(frac > maxFracBinExp) maxFracBinExp= frac;
            hBarBinExp->GetXaxis()->SetBinLabel(cat+1, catLabelsExp[cat]);
        }

        TCanvas* cCatBinExp = new TCanvas(Form("cCatExp_bin%d", ib),
            Form("Cluster Breakdown (expanded) in pT bin %d", ib),
            800,600);
        cCatBinExp->SetGrid();

        double yMaxBinExp = std::min(1.0, 1.1*maxFracBinExp);
        if(yMaxBinExp<0.2) yMaxBinExp=0.2;
        hBarBinExp->SetMaximum(yMaxBinExp);

        // Draw axes only
        hBarBinExp->SetFillStyle(0);
        hBarBinExp->SetLineColor(kBlack);
        hBarBinExp->Draw("hist");

        // Over‐draw colored boxes
        for(int cat=0; cat<kNumExpandedCats; cat++){
            double frac = hBarBinExp->GetBinContent(cat+1);
            double xCenter = hBarBinExp->GetBinCenter(cat+1);
            double xLow  = xCenter - 0.5*barWidth;
            double xHigh = xCenter + 0.5*barWidth;

            TBox* box = new TBox(xLow, 0., xHigh, frac);
            box->SetFillColor( catColorsExp[cat] );
            box->SetLineColor( catColorsExp[cat] );
            box->SetLineWidth(0);
            box->Draw("same");
        }

        double binLo = pTedges_purEff[ib];
        double binHi = pTedges_purEff[ib+1];

        // 8) Save
        std::string outBinNameExp =
            outDir + Form("/ClusterCatExp_Bin%d_%.0fto%.0f.png", ib, binLo, binHi);
        cCatBinExp->SaveAs(outBinNameExp.c_str());
        std::cout << "[INFO] Wrote bin="<< ib
                  <<" expanded cat bar chart => "<< outBinNameExp << "\n";
    }

    std::cout << "[INFO] Done building expanded category plots.\n";
}




static void buildPromptPurityPlot(const std::string & outDir,
                                  const std::vector<double>& isoCutBin)
{
    // (A) Recompute the purity vs pT as before

    // We assume you already computed:
    //   isoCutBin[b]   (the 95% iso cut in each bin)
    //   nTagAll_varIso[b], nTagCorr_varIso[b]
    //   hIsoPromptInBin[b]  => your iso distribution for true prompts in bin b

    std::vector<double> vBinC(nBins_pur, 0.);
    std::vector<double> vPur(nBins_pur,  0.);
    std::vector<double> vErr(nBins_pur,  0.);
    std::vector<double> vIsoCut(nBins_pur, 0.);

    double maxPur = 0.0;

    // 1) Recompute purity in each pT bin
    for(int ib = 0; ib < nBins_pur; ib++)
    {
        double cCenter = 0.5*( pTedges_purEff[ib] + pTedges_purEff[ib+1] );
        vBinC[ib] = cCenter;

        long long denom = nTagAll_varIso[ib];
        long long num   = nTagCorr_varIso[ib];

        double pur  = 0.;
        double ePur = 0.;
        if(denom > 0) {
            pur  = double(num) / double(denom);
            ePur = std::sqrt( pur*(1.0 - pur) / double(denom) );
            if(pur > maxPur) maxPur = pur;
        }
        vPur[ib] = pur;
        vErr[ib] = ePur;

        // also store the isoCut we used for that bin
        vIsoCut[ib] = isoCutBin[ib];
    }

    // 2) Fix the y‐axis range => [0..1.5]
    //    (since user wants from 0 to 1.5, we do that explicitly)
    double yMax = 1.5; // Hard-coded upper bound

    // 3) Make the purity TCanvas
    TCanvas cPur("cPur","Prompt Photon Purity (Var Iso)",800,600);
    gStyle->SetOptStat(0);
    cPur.SetGrid();

    double xMin = pTedges_purEff.front() - 0.5;
    double xMax = pTedges_purEff.back()   + 0.5;

    // We draw a frame from 0..1.5 in purity
    TH1F* hFr = cPur.DrawFrame(xMin, 0.0, xMax, yMax);
    hFr->SetXTitle("Cluster p_{T} [GeV]");
    hFr->SetYTitle("Purity = N_{correct}/N_{tagged}");
    hFr->SetTitle("Prompt Photon Purity vs. p_{T} (bin-by-bin iso @ 95% eff)");

    // 4) TGraphErrors for purity
    TGraphErrors* grPur = new TGraphErrors(nBins_pur,
                                           &vBinC[0],
                                           &vPur[0],
                                           nullptr,
                                           &vErr[0]);
    grPur->SetMarkerStyle(20);
    grPur->SetMarkerColor(kAzure+2);
    grPur->SetLineColor(kAzure+2);
    grPur->SetMarkerSize(1.3);
    grPur->Draw("P SAME");

    // 5) Some text annotation
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.03);
    latex.SetTextFont(42);
    latex.DrawLatex(0.2, 0.25,
        "#it{Variable iso cut at 95% prompt efficiency per bin (Purity up to 1.5)}");

    // 6) Print summary to console
    std::cout << "\n=== Prompt Purity Summary (bin-by-bin isoAll @ 95% eff) ===\n";
    for(int ib=0; ib<nBins_pur; ib++){
        long long denom = nTagAll_varIso[ib];
        long long num   = nTagCorr_varIso[ib];
        double pur = (denom > 0 ? double(num)/double(denom) : 0.);
        std::cout << " pT bin [" << pTedges_purEff[ib] << ","
                  << pTedges_purEff[ib+1] << "]"
                  << " isoCut=" << isoCutBin[ib]
                  << " #tagged=" << denom
                  << " #correct=" << num
                  << " => purity=" << pur << "\n";
    }
    std::cout << "===========================================================\n\n";

    // 7) Save the purity plot
    std::string outName = outDir + "/PromptPurity.png";
    cPur.SaveAs(outName.c_str());
    std::cout << "[INFO] Wrote variable-iso Purity plot => " << outName << "\n";


    // (B) For each pT bin, also produce a plot of hIsoPromptInBin[ib]
    // with a dashed line at isoCutBin[ib].

    for(int ib=0; ib<nBins_pur; ib++)
    {
        TH1F* hIso = hIsoPromptInBin[ib];
        if(!hIso) {
            // If somehow that histogram isn't available, skip
            std::cerr << "[WARN] hIsoPromptInBin["<<ib<<"] is null => skip.\n";
            continue;
        }

        double binLo = pTedges_purEff[ib];
        double binHi = pTedges_purEff[ib+1];
        double isoCut = isoCutBin[ib];

        // Make a canvas for the isolation distribution in this bin
        TCanvas cIso(Form("cIso_bin%d", ib),
                     Form("Prompt iso distribution in pT bin %d", ib),
                     700,600);
        cIso.SetGrid();

        // Optionally set some axis range. For example, do 20% above max:
        double maxY = 1.2 * hIso->GetMaximum();
        hIso->SetTitle(
            Form("Isolation for true prompt, %.0f < pT < %.0f", binLo, binHi)
        );
        hIso->GetXaxis()->SetTitle("Isolation energy [GeV]");
        hIso->GetYaxis()->SetTitle("Entries");
        hIso->SetLineColor(kBlue+1);
        hIso->SetLineWidth(2);
        hIso->SetStats(false);
        hIso->SetMaximum(maxY);

        hIso->Draw("HIST");

        // Now draw a dashed line at isoCut
        TLine cutLine(isoCut, 0.0, isoCut, maxY);
        cutLine.SetLineColor(kRed+1);
        cutLine.SetLineWidth(3);
        cutLine.SetLineStyle(2); // dashed
        cutLine.Draw("SAME");

        // Optionally label it
        TLatex latCut;
        latCut.SetTextFont(42);
        latCut.SetTextSize(0.03);
        latCut.SetTextColor(kRed+1);
        latCut.SetNDC(false); // draw in data-coord
        latCut.DrawLatex(isoCut+0.5, 0.7*maxY, Form("isoCut=%.2f", isoCut));

        // Save each bin's iso distribution plot in the *same* outDir
        // as the purity plot. Name them e.g. "IsoCut_bX.png"
        std::string isoPlotName = Form("%s/IsoCutBin_%d.png", outDir.c_str(), ib);
        cIso.SaveAs(isoPlotName.c_str());
        std::cout << "[INFO] Wrote iso distribution w/ cut => " << isoPlotName << "\n";
    }

    // Done
}



static void buildMesonFakeRatePlot(const std::string & outDir)
{
    // EXACT SAME CODE as your existing block:

    // "Fake Rate" = (N_{meson passing prompt cuts}) / (N_{meson total}),
    // where "meson" means pi0 + eta combined.
    std::vector<double> vx(nBins_pur, 0.);
    std::vector<double> vy(nBins_pur, 0.);
    std::vector<double> vex(nBins_pur, 0.);
    std::vector<double> vey(nBins_pur, 0.);

    double maxFake = 0.0; // track maximum y-value

    for(int ib = 0; ib < nBins_pur; ib++){
        double ptLo  = pTedges_purEff[ib];
        double ptHi  = pTedges_purEff[ib+1];
        double ptMid = 0.5*(ptLo + ptHi);

        // denominator = total (pi0 + eta)
        long long denom = catCount[3][ib] + catCount[4][ib];
        // numerator = total passing prompt cuts from pi0 + eta
        long long num   = nFakeCount[3][ib] + nFakeCount[4][ib];

        double frate = 0.;
        double err   = 0.;
        if(denom > 0){
            frate = double(num) / double(denom);
            err   = std::sqrt(frate*(1.-frate)/ double(denom));
            if(frate > maxFake) maxFake = frate;
        }

        vx[ib]  = ptMid;
        vy[ib]  = frate;
        vey[ib] = err;
    }

    TGraphErrors* grFake = new TGraphErrors(nBins_pur, &vx[0], &vy[0],
                                            &vex[0], &vey[0]);
    grFake->SetMarkerStyle(21);
    grFake->SetMarkerSize(1.3);
    grFake->SetMarkerColor(kOrange+1);
    grFake->SetLineColor(kOrange+1);

    TCanvas cFake("cFake","Combined Meson Fake Rate vs pT",800,600);
    cFake.SetGrid();
    gStyle->SetOptStat(0);

    double xMin = pTedges_purEff.front() - 0.5;
    double xMax = pTedges_purEff.back()  + 0.5;
    double yMax = (maxFake < 0.01 ? 0.1 : 1.2 * maxFake);
    if(yMax > 1.0) yMax = 1.0;

    TH1F* hFrFake = cFake.DrawFrame(xMin, 0., xMax, yMax);
    hFrFake->SetTitle("Meson Fake Rate vs. p_{T} (#pi^{0} + #eta combined)");
    hFrFake->SetXTitle("Cluster p_{T} [GeV]");
    hFrFake->SetYTitle("Fake Rate = #frac{N_{(#pi^{0}+#eta) passing}}{N_{(#pi^{0}+#eta) total}}");

    grFake->Draw("P SAME");

    TLatex latexFake;
    latexFake.SetNDC();
    latexFake.SetTextSize(0.03);
    latexFake.SetTextFont(42);
    latexFake.DrawLatex(0.2, 0.2,
        "Cuts: #DeltaR<0.05, E_{reco}/E_{truth}>0.5, iso_{All}<6, iso_{EM}<4");
    latexFake.DrawLatex(0.2, 0.28,
        "Fake Rate = fraction of (#pi^{0}+#eta) clusters passing 'prompt' ID");

    std::string outFake = outDir + "/FakeRateVsPt_MesonsCombined.png";
    cFake.SaveAs(outFake.c_str());
    std::cout<<"[INFO] Wrote combined meson fake-rate plot => "<< outFake <<std::endl;
}


// -------------------------------------------------------------------
// Global bin edges for leading cluster pT in the invariant-mass analysis
// -------------------------------------------------------------------
static const std::vector<float> pTleadingBins_invM = {5,7,10,15,20,30,40};
static const int nBins_invM = 6;  // since we have 7 edges => 6 intervals

// -------------------------------------------------------------------
// Histograms for each leading-pT bin, each covering inv. mass [0..1] GeV
// -------------------------------------------------------------------
static TH1F* hInvMassBins[nBins_invM] = {nullptr};


void initInvariantMassHistograms()
{
    for(int i = 0; i < nBins_invM; i++){
        float binLo = pTleadingBins_invM[i];
        float binHi = pTleadingBins_invM[i+1];
        
        // Build a histogram name and title
        TString hName  = Form("hInvM_pt%.0f_to_%.0f", binLo, binHi);
        TString hTitle = Form("Invariant Mass for %.0f #leq p_{T}^{leading} < %.0f;M_{#gamma#gamma} [GeV];Entries",
                              binLo, binHi);
        
        hInvMassBins[i] = new TH1F(hName, hTitle, 100, 0.0, 1.0);
        hInvMassBins[i]->Sumw2();

        // <-- The crucial line to ensure ROOT won't delete them when the file closes:
        hInvMassBins[i]->SetDirectory(nullptr);  // or hInvMassBins[i]->SetDirectory(0);
    }
    std::cout << "[INFO] Initialized " << nBins_invM
              << " invariant-mass histograms (SetDirectory=0).\n";
}


// -------------------------------------------------------------------
// Utility to find the correct pT-bin index based on the leading pT
// returns -1 if out of range [5..40]
// -------------------------------------------------------------------
int findLeadingPtBin_invM(float ptLead)
{
    for(int i = 0; i < nBins_invM; i++){
        if(ptLead >= pTleadingBins_invM[i] && ptLead < pTleadingBins_invM[i+1]){
            return i;
        }
    }
    return -1;
}

void printInvMassCutFlowSummary()
{
    std::cout << "\n=== Invariant-Mass Cut-Flow Summary ===\n";
    std::cout << "  Total pairs encountered: " << gCutFlow_totalPairs << "\n";
    std::cout << "  *Skipped* because event |vertexZ|>30 => "
              << gCutFlow_skipVertexZ << "\n";
    std::cout << "  *Skipped* cluster1 prob<0.05 => "
              << gCutFlow_skipProb1 << "\n";
    std::cout << "  *Skipped* cluster2 prob<0.05 => "
              << gCutFlow_skipProb2 << "\n";
    std::cout << "  *Skipped* clus1 had negative pt/E => "
              << gCutFlow_skipPt1neg << "\n";
    std::cout << "  *Skipped* clus2 had negative pt/E => "
              << gCutFlow_skipPt2neg << "\n";
    std::cout << "  *Skipped* mesonPt<1.0 => "
              << gCutFlow_skipMesonPt << "\n";
    std::cout << "  *Skipped* asym>0.7 => "
              << gCutFlow_skipAsym << "\n";
    std::cout << "  *Skipped* dR not in (1/(10pT),1/pT) => "
              << gCutFlow_skipDR << "\n";
    std::cout << "  *Skipped* leading pT outside [5..40] => "
              << gCutFlow_skipLeadPt << "\n";
    std::cout << "  *PASSED* all cuts => "
              << gCutFlow_passAll << "\n";
    std::cout << "=========================================\n\n";
}

void fillInvariantMassDistributions(float vertexz_truth)
{
    // Quick boundary check on ncluster
    if(ncluster < 0) {
        std::cerr << "[ERROR] ncluster < 0 => skipping.\n";
        return;
    }
    if(ncluster > kMaxClusters) {
        std::cerr << "[ERROR] ncluster=" << ncluster
                  << " > kMaxClusters => skipping.\n";
        return;
    }

    // If the event's vertex Z is out of range => skip entire event
    // BUT let's also count how many pairs that *would* have existed
    // so we can credit them to skipVertexZ
    if(std::fabs(vertexz_truth) > 30.0) {
        // number of distinct pairs = N*(N-1)/2
        long long nPairsInEvent = (long long)ncluster * (long long)(ncluster-1)/2;
        gCutFlow_skipVertexZ += nPairsInEvent;
        return;
    }

    // If no clusters, no pairs
    if(ncluster == 0) {
        return; // no pairs to count
    }

    // Loop over all pairs of clusters
    for(int clus1 = 0; clus1 < ncluster; clus1++)
    {
        // cluster1
        float prob1 = cluster_MergedProb[clus1];
        float pt1   = cluster_Et[clus1];
        float E1    = cluster_E[clus1];

        for(int clus2 = clus1 + 1; clus2 < ncluster; clus2++)
        {
            // increment total pairs
            gCutFlow_totalPairs++;

            // cluster2
            float prob2 = cluster_MergedProb[clus2];
            float pt2   = cluster_Et[clus2];
            float E2    = cluster_E[clus2];

            // 1) skip if prob1 < 0.05
            if(prob1 < 0.001f) {
                gCutFlow_skipProb1++;
                continue;
            }

            // 2) skip if prob2 < 0.05
            if(prob2 < 0.001f) {
                gCutFlow_skipProb2++;
                continue;
            }

            // 3) skip if negative pt/E
            if(pt1 < 0.f || E1 < 0.f) {
                gCutFlow_skipPt1neg++;
                continue;
            }
            if(pt2 < 0.f || E2 < 0.f) {
                gCutFlow_skipPt2neg++;
                continue;
            }

            // Build TLorentzVector for each photon
            float eta1 = cluster_Eta[clus1];
            float phi1 = cluster_Phi[clus1];
            TLorentzVector photon1;
            photon1.SetPtEtaPhiE(pt1, eta1, phi1, E1);

            float eta2 = cluster_Eta[clus2];
            float phi2 = cluster_Phi[clus2];
            TLorentzVector photon2;
            photon2.SetPtEtaPhiE(pt2, eta2, phi2, E2);

            // build meson
            TLorentzVector meson = photon1 + photon2;
            float mesonPt   = meson.Pt();
            float mesonMass = meson.M();

            // 4) cut: meson pT > 1.0
            if(mesonPt < 1.0f) {
                gCutFlow_skipMesonPt++;
                continue;
            }

            // 5) asymmetry cut
            float denomE = E1 + E2;
            float asym   = (denomE > 0.f)
                         ? std::fabs(E1 - E2)/denomE
                         : 0.f;
            if(asym > 0.7f) {
                gCutFlow_skipAsym++;
                continue;
            }

            // 6) Delta-R cut
            float dR    = photon1.DeltaR(photon2);
            float dRmin = 1.0f / (10.f*mesonPt);
            float dRmax = 1.0f / mesonPt;
            if(dR <= dRmin || dR >= dRmax) {
                gCutFlow_skipDR++;
                continue;
            }

            // 7) leading pT in [5..40]
            float ptLead = (pt1 > pt2? pt1: pt2);
            int binIdx   = findLeadingPtBin_invM(ptLead);
            if(binIdx < 0) {
                gCutFlow_skipLeadPt++;
                continue;
            }

            // If we reach here => pass all cuts => fill histogram
            if(!hInvMassBins[binIdx]) {
                std::cerr << "[ERROR] hInvMassBins["<<binIdx<<"] is null.\n";
                continue;
            }
            hInvMassBins[binIdx]->Fill(mesonMass);

            // final pass
            gCutFlow_passAll++;
        } // clus2
    } // clus1
}

static double findIsoCutAtEff(TH1F* hist, double eff)
{
    if(!hist) return 0.0;

    double total = hist->Integral();
    if(total <= 1e-6) {
        // If no entries, return 0.0 or some fallback
        return 0.0;
    }

    double goal = eff * total;
    double sum  = 0.0;

    // scan the bins from low to high
    for(int ib = 1; ib <= hist->GetNbinsX(); ib++){
        sum += hist->GetBinContent(ib);
        if(sum >= goal){
            // return the upper edge of this bin
            return hist->GetXaxis()->GetBinUpEdge(ib);
        }
    }

    // if we never reached the goal, return the last bin’s upper edge
    return hist->GetXaxis()->GetBinUpEdge( hist->GetNbinsX() );
}



void saveInvariantMassPlots(const std::string & outDir)
{
    // Ensure the directory exists or can be created
    struct stat info;
    if(stat(outDir.c_str(), &info) != 0){
        std::cout << "[INFO] Creating directory: " << outDir << std::endl;
        std::string cmd = "mkdir -p " + outDir;
        int retCode = system(cmd.c_str());
        if(retCode != 0) {
            std::cerr << "[ERROR] Could not create output directory: "
                      << outDir << "\n";
            return; // or throw, or handle differently
        }
    }

    TCanvas cInvM("cInvM","Invariant Mass Distributions",800,600);
    for(int i = 0; i < nBins_invM; i++){
        cInvM.Clear();

        if(!hInvMassBins[i]){
            // If the histogram was never initialized or is null
            std::cerr << "[WARN] hInvMassBins[" << i
                      << "] is null. Skipping.\n";
            continue;
        }

        hInvMassBins[i]->SetLineColor(kBlue+1);
        hInvMassBins[i]->SetLineWidth(2);

        hInvMassBins[i]->Draw("HIST E");

        float binLo = pTleadingBins_invM[i];
        float binHi = pTleadingBins_invM[i+1];

        // e.g. "InvMass_pt5_7.png"
        TString outName = Form("%s/InvMass_pt%.0f_%.0f.png",
                               outDir.c_str(), binLo, binHi);

        cInvM.SaveAs(outName);
        std::cout << "[INFO] Saved inv-mass plot => " << outName << "\n";
    }
    std::cout << "[INFO] Done saving all invariant-mass plots.\n";
}

/**
 * plotTruthPhotonYields_6to30_Overlay
 *
 * Overlays the yields of "prompt" (class=1) vs "frag" (class=2) photons
 * in 12 bins: [6..8), [8..10), ... [28..30),
 * along with a ratio pad below (Prompt / Frag).
 *
 * In your main code, you should keep a pair of arrays/vectors:
 *   promptCount_2GeV[12], fragCount_2GeV[12]
 * and fill them only for photons with 6 <= pT < 30 (in increments of 2).
 *
 * Usage:
 *   // Suppose we have promptCount_2GeV[12], fragCount_2GeV[12].
 *   // Then, after filling them in the event loop:
 *   plotTruthPhotonYields_6to30_Overlay(promptCount_2GeV, fragCount_2GeV, outDir);
 */
static void plotTruthPhotonYields_6to30_Overlay(
    const std::vector<long long> & promptCount,
    const std::vector<long long> & fragCount,
    const std::string & outDir)
{
    // 1) We have 12 bins from pT=6..30 in steps of 2
    const int nLocalBins = 12;
    if ( (int)promptCount.size()!=nLocalBins || (int)fragCount.size()!=nLocalBins )
    {
        std::cerr << "[ERROR] plotTruthPhotonYields_6to30_Overlay: input vector sizes"
                  << " must be 12. Aborting.\n";
        return;
    }

    // 2) Build the bin edges for 6..8,8..10,...,28..30
    double ptEdgesLocal[nLocalBins+1];
    for(int i=0; i<=nLocalBins; i++){
        ptEdgesLocal[i] = 6.0 + 2.0*double(i);  // 6,8,10,...,30
    }

    // 3) Prepare histograms for the bar chart
    //    We'll do hBarPrompt (red) and hBarFrag (blue) with these 12 bins
    TH1F *hBarPrompt = new TH1F("hBarPrompt","Prompt yields (bar)",
                                nLocalBins, ptEdgesLocal);
    TH1F *hBarFrag   = new TH1F("hBarFrag",  "Frag yields (bar)",
                                nLocalBins, ptEdgesLocal);

    // 4) We'll also prepare arrays for the ratio points
    std::vector<double> vx(nLocalBins, 0.);
    std::vector<double> vyRatio(nLocalBins, 0.);
    std::vector<double> veRatio(nLocalBins, 0.);

    // Fill histograms from promptCount/fragCount, track maximum for top pad
    double maxY = 1.0;
    for(int i=0; i<nLocalBins; i++){
        long long Np = promptCount[i];
        long long Nf = fragCount[i];
        double dp = double(Np);
        double df = double(Nf);

        hBarPrompt->SetBinContent(i+1, dp);
        hBarFrag->SetBinContent(i+1,   df);

        if(dp>maxY) maxY=dp;
        if(df>maxY) maxY=df;

        // ratio = prompt/frag
        double rVal=0.0, rErr=0.0;
        if(df > 0.0) {
            rVal = dp / df;
            if(dp>0.0) {
                // error: ratio * sqrt(1/Np + 1/Nf)
                rErr = rVal * std::sqrt(1.0/dp + 1.0/df);
            }
        }
        vyRatio[i] = rVal;
        veRatio[i] = rErr;

        // x at the center of the bin
        vx[i] = 0.5*(ptEdgesLocal[i] + ptEdgesLocal[i+1]);
    }
    maxY *= 1.3; // pad 30% on top

    // 5) Create a canvas with top/bottom pads
    TCanvas *cCan = new TCanvas("cCan","Truth Photon Yields (Prompt vs Frag), 6 < p_{T} < 30",800,800);
    cCan->SetMargin(0,0,0,0);

    TPad *padTop = new TPad("padTop","padTop", 0.0,0.3, 1.0,1.0);
    padTop->SetBottomMargin(0.02);
    padTop->SetLeftMargin(0.13);
    padTop->Draw();

    TPad *padBot = new TPad("padBot","padBot", 0.0,0.0, 1.0,0.3);
    padBot->SetTopMargin(0.02);
    padBot->SetBottomMargin(0.25);
    padBot->SetLeftMargin(0.13);
    padBot->Draw();

    // === TOP PAD => Bar charts
    padTop->cd();
    padTop->SetGrid();
    gStyle->SetOptStat(0);

    // We'll explicitly set the range from 6..30
    TH1F *hFrameTop = padTop->DrawFrame(6.0, 0.0, 30.0, maxY);
    hFrameTop->SetTitle("Truth Photon Yields in Bins of p_{T} (6..30)");
    hFrameTop->GetXaxis()->SetTitle("");
    hFrameTop->GetXaxis()->SetLabelSize(0.0); // hide X labels on top
    hFrameTop->GetYaxis()->SetTitle("Photon count (#gamma_{class}=1 or 2)");
    hFrameTop->GetYaxis()->SetTitleSize(0.045);
    hFrameTop->GetYaxis()->SetTitleOffset(1.25);

    // Style for prompt bars
    hBarPrompt->SetFillColor(kRed+1);
    hBarPrompt->SetLineColor(kRed+1);
    hBarPrompt->SetBarWidth(0.4);
    hBarPrompt->SetBarOffset(0.1);

    // Style for frag bars
    hBarFrag->SetFillColor(kBlue+1);
    hBarFrag->SetLineColor(kBlue+1);
    hBarFrag->SetBarWidth(0.4);
    hBarFrag->SetBarOffset(0.5);

    // Draw them as bars
    hBarPrompt->Draw("SAME BAR");
    hBarFrag->Draw("SAME BAR");

    // A legend
    TLegend *leg = new TLegend(0.60,0.65,0.88,0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hBarPrompt,"Prompt (#gamma_{class}=1)","f");
    leg->AddEntry(hBarFrag,  "Frag (#gamma_{class}=2)","f");
    leg->Draw();

    // === BOTTOM PAD => ratio
    padBot->cd();
    padBot->SetGrid();

    // find ratio max
    double maxRatio = 1.0;
    for(int i=0; i<nLocalBins; i++){
        double up = vyRatio[i] + veRatio[i];
        if(up>maxRatio) maxRatio=up;
    }
    if(maxRatio<0.2) maxRatio=1.0;
    else maxRatio *= 1.2;

    TH1F *hFrameBot = padBot->DrawFrame(6.0, 0.0, 30.0, maxRatio);
    hFrameBot->SetTitle("");
    hFrameBot->GetXaxis()->SetTitle("p_{T} [GeV]");
    hFrameBot->GetYaxis()->SetTitle("Prompt / Frag");
    hFrameBot->GetYaxis()->SetTitleSize(0.08);
    hFrameBot->GetYaxis()->SetTitleOffset(0.55);
    hFrameBot->GetYaxis()->SetLabelSize(0.07);
    hFrameBot->GetXaxis()->SetLabelSize(0.08);
    hFrameBot->GetXaxis()->SetTitleSize(0.09);

    // TGraphErrors for ratio points
    TGraphErrors *grRatio = new TGraphErrors(nLocalBins, &vx[0], &vyRatio[0], nullptr, &veRatio[0]);
    grRatio->SetMarkerStyle(20);
    grRatio->SetMarkerColor(kGray+2);
    grRatio->SetLineColor(kGray+2);
    grRatio->SetMarkerSize(1.2);
    grRatio->Draw("P SAME");

    // 6) Optionally force x-ticks at 6,8,10,...,30
    hFrameBot->GetXaxis()->SetNdivisions(-12,  false);
    // Negative => "exactly that many" ticks, no minor divisions.
    // or you can do more advanced labeling with custom bin labels

    // 7) Save
    cCan->cd();
    cCan->Update();
    std::string outName = outDir + "/TruthPhotonYields_Overlay_6to30.png";
    cCan->SaveAs(outName.c_str());
    std::cout << "[INFO] Wrote bar-chart yields + ratio => " << outName << "\n";
    
    // If you don't want the canvas to persist in memory:
    // delete cCan;
}

static void plotTruthPhotonYields_6to30_ThreeBar(
    const std::vector<long long> & promptCount,      // red bars
    const std::vector<long long> & fragCount,        // blue bars
    const std::vector<long long> & mesonDecayCount,  // green bars
    const std::string & outDir)
{
    // 1) We have 12 bins from pT=6..30 in steps of 2
    const int nLocalBins = 12;
    // Check input sizes
    if (  promptCount.size() != size_t(nLocalBins)
       || fragCount.size()   != size_t(nLocalBins)
       || mesonDecayCount.size() != size_t(nLocalBins) )
    {
        std::cerr << "[ERROR] plotTruthPhotonYields_6to30_ThreeBar: input vector sizes"
                  << " must each be 12. Aborting.\n";
        return;
    }

    // 2) Build bin edges for [6..8), [8..10), ..., [28..30)
    double ptEdgesLocal[nLocalBins+1];
    for(int i=0; i<=nLocalBins; i++){
        ptEdgesLocal[i] = 6.0 + 2.0*double(i);  // 6,8,10,...,30
    }

    // 3) Prepare histograms for the bar chart
    //    We'll do hBarPrompt (red), hBarFrag (blue), hBarMeson (green).
    TH1F *hBarPrompt = new TH1F("hBarPrompt3","Prompt yields (bar)",
                                nLocalBins, ptEdgesLocal);
    TH1F *hBarFrag   = new TH1F("hBarFrag3",  "Frag yields (bar)",
                                nLocalBins, ptEdgesLocal);
    TH1F *hBarMeson  = new TH1F("hBarMeson3", "Meson-decay yields (bar)",
                                nLocalBins, ptEdgesLocal);

    // 4) Fill histograms, track maximum
    double maxY = 1.0;
    for(int i=0; i<nLocalBins; i++){
        double pVal = double( promptCount[i] );
        double fVal = double( fragCount[i] );
        double mVal = double( mesonDecayCount[i] );

        hBarPrompt->SetBinContent(i+1, pVal);
        hBarFrag->SetBinContent(i+1,   fVal);
        hBarMeson->SetBinContent(i+1,  mVal);

        if(pVal>maxY) maxY=pVal;
        if(fVal>maxY) maxY=fVal;
        if(mVal>maxY) maxY=mVal;
    }
    maxY *= 1.3; // pad 30% on top

    // 5) Create a canvas (no subpads)
    TCanvas *cCan = new TCanvas("cCan3","Truth Photon Yields (Prompt vs Frag vs Meson) 6..30",800,600);
    cCan->SetMargin(0.12,0.05,0.12,0.08); // left,right,bottom,top margins as you prefer
    cCan->SetGrid();
    gStyle->SetOptStat(0);

    // 6) Draw a frame from x=6..30
    TH1F *hFrame = (TH1F*)gPad->DrawFrame(6.0,0.0, 30.0, maxY);
    hFrame->SetTitle("Truth Photon Yields in Bins of p_{T} (6..30)");
    hFrame->GetXaxis()->SetTitle("p_{T} [GeV]");
    hFrame->GetYaxis()->SetTitle("Photon count");
    hFrame->GetYaxis()->SetTitleSize(0.045);
    hFrame->GetYaxis()->SetTitleOffset(1.1);

    // 7) Style & offsets for the three sets of bars
    hBarPrompt->SetFillColor(kRed+1);
    hBarPrompt->SetLineColor(kRed+1);
    hBarPrompt->SetBarWidth(0.25);
    hBarPrompt->SetBarOffset(0.05);

    hBarFrag->SetFillColor(kBlue+1);
    hBarFrag->SetLineColor(kBlue+1);
    hBarFrag->SetBarWidth(0.25);
    hBarFrag->SetBarOffset(0.35);

    hBarMeson->SetFillColor(kGreen+2);
    hBarMeson->SetLineColor(kGreen+2);
    hBarMeson->SetBarWidth(0.25);
    hBarMeson->SetBarOffset(0.65);

    // 8) Draw them as bars
    hBarPrompt->Draw("SAME BAR");
    hBarFrag->Draw("SAME BAR");
    hBarMeson->Draw("SAME BAR");

    // 9) Legend
    TLegend *leg = new TLegend(0.7,0.65,0.88,0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hBarPrompt, "Prompt (#gamma_{class}=1)","f");
    leg->AddEntry(hBarFrag,   "Frag (#gamma_{class}=2)",  "f");
    leg->AddEntry(hBarMeson,  "#pi^{0}/#eta decay",       "f");
    leg->Draw();

    // 10) Save
    cCan->Update();
    std::string outName = outDir + "/TruthPhotonYields_ThreeBar_6to30.png";
    cCan->SaveAs(outName.c_str());
    std::cout << "[INFO] Wrote 3-bar yields => " << outName << std::endl;

    // if you don’t want it in memory afterward:
    // delete cCan;
}


void analyzeJet10SlimTree()
{
    // 1) Input & Output
    std::string inputFile = "/Users/patsfan753/Desktop/caloana0130.root";
    std::string outDir    = "/Users/patsfan753/Desktop/SimOut";
    std::string catBarDir  = outDir + "/ClusterCatBarCharts";
    std::string qaDir      = outDir + "/QA";
    std::string invMdir    = outDir + "/InvariantMass";
    std::string purityDir  = outDir + "/Purity";
    ensureOutputDirectory(outDir);

    gAllClusters.clear();
    static std::vector<long long> gPromptCount_2GeV(12, 0);
    static std::vector<long long> gFragCount_2GeV(12, 0);
    static std::vector<long long> gMesonDecayCount_2GeV(12, 0);
    
    // 2) Check TTree
    {
        TFile* fTest = TFile::Open(inputFile.c_str(), "READ");
        if(!fTest || fTest->IsZombie()){
            std::cerr << "[ERROR] Cannot open file: " << inputFile << std::endl;
            return;
        }
        TTree* tTest = (TTree*) fTest->Get("slimtree");
        if(!tTest){
            std::cerr << "[ERROR] 'slimtree' not found.\n";
            fTest->Close();
            return;
        }
        std::cout << "\n=== Branches in 'slimtree' ===\n";
        TObjArray* branchList = tTest->GetListOfBranches();
        for(int iB = 0; iB < branchList->GetEntries(); iB++){
            TBranch* br = (TBranch*) branchList->At(iB);
            if(br) std::cout << "   " << br->GetName() << "\n";
        }
        std::cout << "================================\n\n";
        fTest->Close();
    }

    // 3) Initialize bin counts for iso-efficiency
    nBins   = (int)pTedges.size() - 1;   // => 4
    nIsoCuts= (int)isoCutsGeV.size();    // => 11
    std::memset(isoCountByClass, 0, sizeof(isoCountByClass));
    std::memset(nPhotonsInBin,   0, sizeof(nPhotonsInBin));

    // 4) Initialize bin counts for purity
    nBins_pur = (int)pTedges_purEff.size() - 1; // => 7
    std::memset(nTagAll,     0, sizeof(nTagAll));
    std::memset(nTagCorrect, 0, sizeof(nTagCorrect));
    
    //Create histograms to store isoAll for 'true prompt' clusters in each pT bin
    for(int ib = 0; ib < nBins_pur; ib++){
        float lo = pTedges_purEff[ib];
        float hi = pTedges_purEff[ib+1];
        TString hName = Form("hIsoPrompt_b%d", ib);
        TString hTit  = Form("isoAll for true prompt, %.0f < pT < %.0f", lo, hi);
        // Suppose isoAll can reach up to ~50. Adjust bins/range as needed:
        hIsoPromptInBin[ib] = new TH1F(hName, hTit, 200, 0.0, 50.0);
        hIsoPromptInBin[ib]->SetDirectory(nullptr); // keep it memory-resident
    }
    
    
    std::memset(nTagAllAlt,  0, sizeof(nTagAllAlt));
    std::memset(nTagCorrAlt, 0, sizeof(nTagCorrAlt));
    
    // Also zero out the fragmentation arrays:
    std::memset(nTagAllFrag,     0, sizeof(nTagAllFrag));
    std::memset(nTagCorrectFrag, 0, sizeof(nTagCorrectFrag));

    // 5) Initialize counters for the 6-category breakdown
    std::memset(catCount, 0, sizeof(catCount));
    std::memset(nFakeCount, 0, sizeof(nFakeCount));

    TFile* fIn = TFile::Open(inputFile.c_str(), "READ");
    if(!fIn || fIn->IsZombie()){
        std::cerr << "[ERROR] Could not open input file: " << inputFile << std::endl;
        return;
    }
    TTree* slim = (TTree*) fIn->Get("slimtree");
    if(!slim){
        std::cerr << "[ERROR] 'slimtree' not found.\n";
        fIn->Close();
        return;
    }

    static float vertexz = 9999.f;
    
    //identification -- truth level information
    slim->SetBranchAddress("particle_pid",       particle_pid);
    slim->SetBranchAddress("daughter_pid",       daughter_pid);
    slim->SetBranchAddress("particle_photonclass", particle_photonclass);
    slim->SetBranchAddress("particle_photon_mother_pid", particle_photon_mother_pid);
    slim->SetBranchAddress("particle_converted", particle_converted);
    
    //kinematic -- truth level information
    slim->SetBranchAddress("nparticles",         &nparticles);
    slim->SetBranchAddress("vertexz_truth", &vertexz_truth);
    slim->SetBranchAddress("particle_E",        particle_E);
    slim->SetBranchAddress("particle_Pt",        particle_Pt);
    slim->SetBranchAddress("particle_Eta",       particle_Eta_);
    slim->SetBranchAddress("particle_Phi",       particle_Phi_);
    
    slim->SetBranchAddress("daughter_vtx_z", &daughter_vtx_z);
    slim->SetBranchAddress("ndaughter",         &ndaughters);
    slim->SetBranchAddress("daughter_E",        daughter_E);
    slim->SetBranchAddress("daughter_Pt",       daughter_Pt);
    slim->SetBranchAddress("daughter_Eta",       particle_Eta_);
    slim->SetBranchAddress("daughter_Phi",       daughter_Phi_);
    
    
    //truth level isolation information
    slim->SetBranchAddress("particle_truth_iso_04", particle_truth_iso_04);

    
    /*
     reconstructed cluster MC informatioon
     */
    //identifying info
    slim->SetBranchAddress("cluster_pid_CLUSTERINFO_CEMC",          cluster_pid);
    slim->SetBranchAddress("cluster_prob_CLUSTERINFO_CEMC", cluster_prob);
    slim->SetBranchAddress("cluster_merged_prob_CLUSTERINFO_CEMC", cluster_MergedProb);
    
    //reco kinematic information
    slim->SetBranchAddress("ncluster_CLUSTERINFO_CEMC",    &ncluster);
    slim->SetBranchAddress("vertexz", &vertexz);
    slim->SetBranchAddress("cluster_Et_CLUSTERINFO_CEMC",  cluster_Et);
    slim->SetBranchAddress("cluster_Eta_CLUSTERINFO_CEMC", cluster_Eta);
    slim->SetBranchAddress("cluster_Phi_CLUSTERINFO_CEMC", cluster_Phi);
    slim->SetBranchAddress("cluster_E_CLUSTERINFO_CEMC",   cluster_E);
    
    //reco isolation energy information
    slim->SetBranchAddress("cluster_iso_04_CLUSTERINFO_CEMC",       cluster_iso_04);
    slim->SetBranchAddress("cluster_iso_04_emcal_CLUSTERINFO_CEMC", cluster_iso_04_emcal);
    slim->SetBranchAddress("cluster_iso_04_hcalin_CLUSTERINFO_CEMC",  cluster_iso_04_hcalin);
    slim->SetBranchAddress("cluster_iso_04_hcalout_CLUSTERINFO_CEMC", cluster_iso_04_hcalout);

    static TH2F h2TruthVsReco_prompt(
        "h2TruthVsReco_prompt",
        "Prompt photons; true iso [GeV]; reco iso [GeV]",
        60,0,12,
        60,0,12
    );

    static TH2F h2TruthVsReco_frag(
        "h2TruthVsReco_frag",
        "Frag photons; true iso [GeV]; reco iso [GeV]",
        60,0,12,
        60,0,12
    );

    static TH2F h2TruthVsReco_meson(
        "h2TruthVsReco_meson",
        "Meson decay (pi^{0} or #eta); true iso [GeV]; reco iso [GeV]",
        60,0,12,
        60,0,12
    );

    
    TH1F *hClusterProb = new TH1F("hClusterProb", "Cluster Probability Distribution; Cluster Probability; Entries", 100, 0.0, 1.0);
    hClusterProb->SetDirectory(nullptr);

    TH1F *hClusterMergedProb = new TH1F("hClusterMergedProb", "Cluster Merged Probability Distribution; Cluster Merged Probability; Entries", 100, 0.0, 1.0);
    hClusterMergedProb->SetDirectory(nullptr);

    initInvariantMassHistograms();
    
    TH1F *hIsoTotal    = new TH1F("hIsoTotal", "Total Isolation Distribution; Isolation [GeV]; Entries", 100, 0, 100);
    TH1F *hIsoEMCal    = new TH1F("hIsoEMCal", "EMCal Isolation Distribution; Isolation [GeV]; Entries", 100, 0, 100);
    TH1F *hIsoHCALin   = new TH1F("hIsoHCALin", "HCAL-in Isolation Distribution; Isolation [GeV]; Entries", 100, 0, 100);
    TH1F *hIsoHCALout  = new TH1F("hIsoHCALout", "HCAL-out Isolation Distribution; Isolation [GeV]; Entries", 100, 0, 100);
    // No need to set Directory (they are not saved in ROOT file)
    hIsoTotal->SetDirectory(nullptr);
    hIsoEMCal->SetDirectory(nullptr);
    hIsoHCALin->SetDirectory(nullptr);
    hIsoHCALout->SetDirectory(nullptr);

    
    // 8) Main loop => fill iso-efficiency + purity + categories
    Long64_t nEntries = slim->GetEntries();
    std::cout << "[INFO] nEntries in TTree: " << nEntries << "\n";

    for(Long64_t iEvt = 0; iEvt < nEntries; iEvt++){
        slim->GetEntry(iEvt);
        
        
        // -- Check cluster array size
        if (ncluster < 0 || ncluster > kMaxClusters) {
            std::cerr << "[ERROR] ncluster=" << ncluster
                      << " is out of bounds [0.."
                      << kMaxClusters << "]. Aborting.\n";
            return;  // or break; or throw, depending on your preference
        }

        // -- Check particle array size
        if (nparticles < 0 || nparticles > kMaxParticles) {
            std::cerr << "[ERROR] nparticles=" << nparticles
                      << " is out of bounds [0.."
                      << kMaxParticles << "]. Aborting.\n";
            return;  // or break; or throw
        }


        bool doPrint = (kDebug && iEvt < 10);

        if(doPrint){
            std::cout<<"\n----------------------\n";
            std::cout<<" [DEBUG] Event "<< iEvt
                     <<": nparticles="<< nparticles
                     <<", ncluster="<< ncluster <<"\n";
        }

        fillInvariantMassDistributions(vertexz_truth);
        
        // (A) iso-efficiency from particle info
        for(int ip = 0; ip < nparticles; ip++){
            if(particle_pid[ip] != 22) continue; // must be photon
            int cIdx = photonClassIndex(particle_photonclass[ip]); // 0=prompt,1=frag
            if(cIdx < 0) continue; // skip if decay or something
            float pt   = particle_Pt[ip];
            float isoV = particle_truth_iso_04[ip];
            int bIso   = findPtBin(pt); // [10..30], or -1 if out-of-range
            if(bIso < 0) continue;

            nPhotonsInBin[cIdx][bIso]++;
            for(int ic = 0; ic < nIsoCuts; ic++){
                if(isoV < isoCutsGeV[ic]){
                    isoCountByClass[cIdx][bIso][ic]++;
                }
            }
            if(pt >= 6.f && pt < 30.f){
                int binIdx = int( (pt - 6.f)/2.f );
                if(binIdx >=0 && binIdx<12) {
                    // cIdx=0 => prompt, cIdx=1 => frag
                    if(cIdx == 0) gPromptCount_2GeV[binIdx]++;
                    if(cIdx == 1) gFragCount_2GeV[binIdx]++;
                }

                int mom = std::abs( particle_photon_mother_pid[ip] );
                if(mom == 111 || mom == 221){
                    gMesonDecayCount_2GeV[binIdx]++;
                }
            }
        }

        // (B) cluster-based loop => do prompt-purity + new 6-category
        for(int ic = 0; ic < ncluster; ic++){
            
            hClusterProb->Fill(cluster_prob[ic]);
            hClusterMergedProb->Fill(cluster_MergedProb[ic]);
            
            if(std::fabs(vertexz) > 30.0) {
                // We skip the event entirely if the reconstructed vertex is out of range
                continue;
            }
            
            // NEW: we enforce a minimum cluster_Et to skip low-pT clusters
            float cPt    = cluster_Et[ic];
            if(cPt < 5.0f) {  // user-chosen threshold, e.g. 5 GeV
               continue;
            }

            float cEta   = cluster_Eta[ic];
            float cPhi   = cluster_Phi[ic];
            float cE     = cluster_E[ic];
            float isoAll = cluster_iso_04[ic];
            float isoEM  = cluster_iso_04_emcal[ic];
            int   pid    = cluster_pid[ic]; // leading PID

            hIsoTotal->Fill(cluster_iso_04[ic]);
            hIsoEMCal->Fill(cluster_iso_04_emcal[ic]);
            hIsoHCALin->Fill(cluster_iso_04_hcalin[ic]);
            hIsoHCALout->Fill(cluster_iso_04_hcalout[ic]);
            
            // find purity bin [5..60]
            int   bPur   = findPurityBin(cPt);
            int   cat    = -1;

            // If out of [5..60], skip
            if(bPur < 0) {
                if(doPrint){
                    std::cout<<"   [DEBUG-cluster "<<ic<<"] cPt="<<cPt
                             <<" => outside [5..60], skip cat count.\n";
                }
                continue;
            }

            // 1) If cat==3 => pi0-led, or cat==4 => eta-led, skip from prompt
            //    We'll fill cat later, so we can't check cat==3,4 yet:
            //    We do that after we figure out cat below.

            // A) Normal "prompt candidate" cuts => fill nTagAll + nTagCorrect
            bool passCuts = passClusterPromptCuts(isoAll, isoEM);
            if(passCuts) {
                nTagAll[bPur]++;
                int mm = findMatchedPromptPhoton(cEta, cPhi, cE);
                if(mm >= 0) {
                   nTagCorrect[bPur]++;
                }
            }

            // B) alt combos
            for (int iC = 0; iC < 3; iC++) {
                if (!passIsoCombo(isoAll, isoEM, iC)) continue;
                nTagAllAlt[iC][bPur]++;

                int mm2 = findMatchedPromptPhoton(cEta, cPhi, cE);
                if (mm2 >= 0) {
                    nTagCorrAlt[iC][bPur]++;
                }
            }

            // 2) frag-candidate clusters => same bin
            if(passCuts){ // same iso cuts
                nTagAllFrag[bPur]++;
                int matchFrag = findMatchedFragPhoton(cEta, cPhi, cE);
                if(matchFrag >= 0){
                    nTagCorrectFrag[bPur]++;
                }
            }

            int bestIdx = -1;
            float bestDR = 9999.f;
            
            // 3) Decide cat => 0..5 for the big 6-category classification
            if(pid == 22)
            {
                // If leading PID=22 => match => see if prompt,frag,decay, or other
                int phClass = findMatchedPhotonClass(cEta, cPhi, cE);

                {
                   for(int ip=0; ip<nparticles; ip++){
                       if(particle_pid[ip] != 22) continue;
                       float dr = deltaR_func(cEta, particle_Eta_[ip], cPhi, particle_Phi_[ip]);
                       if(dr>0.05) continue;
                       float eTruth= particle_E[ip];
                       if(eTruth<1e-3) continue;
                       float ratio = cE/eTruth;
                       if(ratio<0.5f) continue;
                       if(dr<bestDR){
                           bestDR=dr;
                           bestIdx=ip;
                       }
                   }
                }

                // Basic classification from photonclass
                if     (phClass == 1) cat = 0; // prompt
                else if(phClass == 2) cat = 1; // frag
                else if(phClass == 3) cat = 2; // decay
                else                  cat = 5; // other

                // If we found a “prompt” photon BUT mother is pi0/eta => override to decay
                if(cat == 0 && bestIdx >= 0)
                {
                   int momPid = particle_photon_mother_pid[bestIdx];
                   if(std::abs(momPid) == 111 || std::abs(momPid) == 221){
                      cat = 2; // re‐label as decay
                   }

                   // NEW: also skip conversions
                   if(particle_converted[bestIdx] == 1) {
                      cat = 2; // or cat=5 => whichever you prefer for "not a pure prompt"
                   }
                }
            }
            else if(pid == 111) cat = 3; // pi0-led
            else if(pid == 221) cat = 4; // eta-led
            else                cat = 5; // other

            
            ClusterRecord rec;
            rec.cPt         = cPt;
            rec.isoAll      = isoAll;
            rec.isTruePrompt= (cat == 0); // i.e. “cat==0 => truly prompt”
            rec.ptBinIndex  = bPur;

            gAllClusters.push_back(rec);

            
            // If cat=0 => "truly prompt" in your final logic, fill the isoAll distribution
            if(cat == 0) {
                // We already have bPur for the pT bin
                hIsoPromptInBin[bPur]->Fill( isoAll );
            }
            
            // Then see if cat=3 or cat=4 => those we skip from "prompt" logic up top
            // but for the big category histograms:
            // FILL THE 2D HISTOGRAMS FOR RECO vs TRUTH ISOLATION
            if(cat == 0) // cat==0 => Prompt
            {
                // Step 1) find actual matched photon index
                int bestMatch = -1;
                float bestDR  = 9999.0;
                for(int ip = 0; ip < nparticles; ip++){
                    if(particle_pid[ip] != 22) continue;
                    if(particle_photonclass[ip] != 1) continue; // must be "prompt"
                    float dr = deltaR_func(cEta, particle_Eta_[ip], cPhi, particle_Phi_[ip]);
                    if(dr > 0.05) continue;
                    float eTruth = particle_Pt[ip] * std::cosh(particle_Eta_[ip]);
                    if(eTruth < 1e-3) continue;
                    float ratio = cE / eTruth;
                    if(ratio < 0.5f) continue;
                    if(dr < bestDR){
                        bestDR = dr;
                        bestMatch = ip;
                    }
                }
                if(bestMatch >= 0){
                    float truthIso   = particle_truth_iso_04[bestMatch];
                    float recoIsoAll = cluster_iso_04[ic];
                    h2TruthVsReco_prompt.Fill(truthIso, recoIsoAll);
                }
            }
            else if(cat == 1) // cat==1 => Frag
            {
                int bestMatch = -1;
                float bestDR  = 9999.0;
                for(int ip = 0; ip < nparticles; ip++){
                    if(particle_pid[ip] != 22) continue;
                    if(particle_photonclass[ip] != 2) continue; // must be frag
                    float dr = deltaR_func(cEta, particle_Eta_[ip], cPhi, particle_Phi_[ip]);
                    if(dr > 0.05) continue;
                    float eTruth = particle_Pt[ip]*std::cosh(particle_Eta_[ip]);
                    if(eTruth < 1e-3) continue;
                    float ratio = cE/eTruth;
                    if(ratio<0.5f) continue;
                    if(dr < bestDR){
                        bestDR = dr;
                        bestMatch = ip;
                    }
                }
                if(bestMatch >= 0){
                    float truthIso   = particle_truth_iso_04[bestMatch];
                    float recoIsoAll = cluster_iso_04[ic];
                    h2TruthVsReco_frag.Fill(truthIso, recoIsoAll);
                }
            }
            
            if(cat == 2) // i.e. “decay” from pi0 or eta
            {
                // But we must verify that we actually matched a photon with mother=111 or 221.
                // We can either check the bestIdx again or we can simply fill if motherPid is ±111 or ±221.
                // Easiest: check bestIdx's mother PID:
                if(bestIdx >= 0) {
                   int momPid = particle_photon_mother_pid[bestIdx];
                   if(std::abs(momPid)==111 || std::abs(momPid)==221)
                   {
                      // fill the new 2D hist
                      float truthIso = particle_truth_iso_04[bestIdx];
                      float recoIso  = cluster_iso_04[ic];
                      h2TruthVsReco_meson.Fill(truthIso, recoIso);
                   }
                }
            }

            // Expand cat => subcategory
            int expandedCat = -1;
            if(cat < 5) {
                expandedCat = cat; // 0..4
            }
            else {
                // cat==5 => "other". Map to subcat => no-match, pions, kaons, etc.
                int idxMatch2 = findBestMatchAnyParticle(cEta, cPhi, cE);
                int subcat = mapOtherPidToSubcat( (idxMatch2<0? 999999 : particle_pid[idxMatch2]) );
                expandedCat = subcat;
            }

            catCount[cat][bPur]++;
            catCountExpanded[expandedCat][bPur]++;
        } // end cluster loop

        if(doPrint){
            std::cout<<" [DEBUG] End of event "<< iEvt <<" catCount so far (first 10 bins):\n";
            for(int cat = 0; cat < 6; cat++){
                long long sumCat = 0;
                for(int ib = 0; ib < std::min(nBins_pur,10); ib++){
                    sumCat += catCount[cat][ib];
                }
                std::cout<<"    cat="<<cat<<" => sumCat="<<sumCat<<"\n";
            }
        }


        if((iEvt > 0) && (iEvt % 100000 == 0)){
            std::cout << "[INFO] Processed " << iEvt << "/" << nEntries << " events...\n";
        }
    } // end event loop

    // Now close input
    fIn->Close();
    

    std::vector<double> isoCutBin(nBins_pur, 0.0);

    // e.g. target eff=95% => find the iso cut in each pT bin
    double targetEff = 0.95;
    for(int b = 0; b < nBins_pur; b++){
        // find the 95% iso cut from the hIsoPromptInBin histogram
        isoCutBin[b] = findIsoCutAtEff(hIsoPromptInBin[b], targetEff);
        std::cout << "[INFO] pT bin=" << b
                  << " => isoCut=" << isoCutBin[b] << std::endl;
    }

    // Zero out the global arrays for variable-iso counts
    std::memset(nTagAll_varIso,  0, sizeof(nTagAll_varIso));
    std::memset(nTagCorr_varIso, 0, sizeof(nTagCorr_varIso));

    // Now loop over the 'gAllClusters' in memory and apply per-bin iso cuts
    for(const auto & rec : gAllClusters)
    {
        int b = rec.ptBinIndex; // which pT bin
        if(b < 0 || b >= nBins_pur) continue; // safety check

        // "tag" as prompt-like if isoAll < isoCutBin[b]
        if(rec.isoAll < isoCutBin[b]) {
            nTagAll_varIso[b]++;
            if(rec.isTruePrompt) {
                nTagCorr_varIso[b]++;
            }
        }
    }

    // now call buildPromptPurityPlot(...) which uses
    // nTagAll_varIso[b], nTagCorr_varIso[b], and isoCutBin
    // to build the final Purity vs. pT plot
    buildPromptPurityPlot(purityDir, isoCutBin);


    // 9) Build iso-efficiency & purity plots
    buildAndSaveIsoEfficiencyPlot_2Classes(qaDir);
    buildAndSaveTruthVsRecoIsoPlots(qaDir,
        h2TruthVsReco_prompt,
        h2TruthVsReco_frag,
        h2TruthVsReco_meson
    );
    
    buildMesonFakeRatePlot(purityDir);
    
    buildExpandedCategoryPlots(catBarDir);
    
    printInvMassCutFlowSummary();
    saveInvariantMassPlots(invMdir);

    // Optionally save cluster_prob hist, cluster_MergedProb hist, etc.
    TCanvas cProb("cProb", "Cluster Probability", 800,600);
    hClusterProb->SetLineColor(kRed);
    hClusterProb->SetLineWidth(2);
    hClusterProb->Draw("HIST");
    TString outNameProb = Form("%s/ClusterProb.png", qaDir.c_str());
    cProb.SaveAs(outNameProb);
    std::cout << "[INFO] Saved cluster_prob plot => " << outNameProb << "\n";

    TCanvas cMergedProb("cMergedProb", "Cluster Merged Probability", 800,600);
    hClusterMergedProb->SetLineColor(kBlue);
    hClusterMergedProb->SetLineWidth(2);
    hClusterMergedProb->Draw("HIST");
    TString outNameMergedProb = Form("%s/ClusterMergedProb.png", qaDir.c_str());
    cMergedProb.SaveAs(outNameMergedProb);
    std::cout << "[INFO] Saved cluster_MergedProb plot => " << outNameMergedProb << "\n";


    // ********** New code: Build and save 1D isolation histograms **********

    // Define output directories for 1D isolation plots
    std::string totalIsoDir   = qaDir + "/totalIsolation";
    std::string emcalIsoDir   = qaDir + "/EMCalIsolation";
    std::string ihcalIsoDir   = qaDir + "/IHCalIsolation";
    std::string ohcalIsoDir   = qaDir + "/OHCalIsolation";
    ensureOutputDirectory(totalIsoDir);
    ensureOutputDirectory(emcalIsoDir);
    ensureOutputDirectory(ihcalIsoDir);
    ensureOutputDirectory(ohcalIsoDir);

    // Automatically adjust axes by scanning nonzero bins
    auto adjustAxis = [](TH1F* h) {
         double minVal = 1e9;
         double maxVal = -1e9;
         for (int i = 1; i <= h->GetNbinsX(); i++){
             double content = h->GetBinContent(i);
             if (content > 0){
                 double center = h->GetBinCenter(i);
                 if(center < minVal) minVal = center;
                 if(center > maxVal) maxVal = center;
             }
         }
         if(minVal < maxVal) {
             h->GetXaxis()->SetRangeUser(minVal, maxVal);
         }
         h->SetMaximum(h->GetMaximum() * 1.2);
    };

    // Draw and save each histogram in its corresponding directory
    TCanvas cIso1("cIso1", "Isolation Distribution", 800,600);
    cIso1.SetGrid();
    cIso1.SetLogy();

    adjustAxis(hIsoTotal);
    hIsoTotal->Draw("HIST");
    cIso1.SaveAs((totalIsoDir + "/TotalIsolation.png").c_str());
    cIso1.Clear();

    adjustAxis(hIsoEMCal);
    hIsoEMCal->Draw("HIST");
    cIso1.SaveAs((emcalIsoDir + "/EMCalIsolation.png").c_str());
    cIso1.Clear();

    adjustAxis(hIsoHCALin);
    hIsoHCALin->Draw("HIST");
    cIso1.SaveAs((ihcalIsoDir + "/HCALinIsolation.png").c_str());
    cIso1.Clear();

    adjustAxis(hIsoHCALout);
    hIsoHCALout->Draw("HIST");
    cIso1.SaveAs((ohcalIsoDir + "/HCALoutIsolation.png").c_str());
    cIso1.Clear();

    
    plotTruthPhotonYields_6to30_Overlay(gPromptCount_2GeV, gFragCount_2GeV, qaDir);
    
    plotTruthPhotonYields_6to30_ThreeBar(
        gPromptCount_2GeV,
        gFragCount_2GeV,
        gMesonDecayCount_2GeV,
        qaDir
    );
    
    std::cout << "[INFO] All done.\n";
}
