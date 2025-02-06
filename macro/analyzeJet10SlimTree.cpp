// analyze_Jet10SlimTree.cpp
// -------------------------------------------------------------------
// Single macro that preserves the exact functionality of the original
// isolation-efficiency code and the prompt-photon purity measurement,
// and ALSO categorizes each cluster into 6 categories (prompt, frag,
// decay, pi0-led, eta-led, other) all in ONE TTree loop.
//
// Outputs:
//   1) The iso-efficiency plot (unchanged),
//   2) The prompt-photon purity plot (unchanged),
//   3) A NEW single bar chart (pT-independent) that shows the fraction
//      of all clusters in each of the 6 categories.
//
// Run it as: root -b -q analyze_Jet10SlimTree.cpp
// -------------------------------------------------------------------

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

static const int kMaxClusters   = 200000;
static const int kMaxParticles  = 200000;

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

int   nparticles;
int   particle_pid[kMaxParticles];
int   particle_photonclass[kMaxParticles]; // 1=prompt,2=frag,3=decay,4=other
float particle_Pt[kMaxParticles];
float particle_Eta_[kMaxParticles];
float particle_Phi_[kMaxParticles];
float particle_truth_iso_04[kMaxParticles];

// -------------------------------------------------------------------
// pT bin edges and iso-thresholds for the iso-efficiency
// -------------------------------------------------------------------
static std::vector<float> pTedges = {10, 15, 20, 25, 30};
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
static std::vector<float> pTedges_purEff = {5, 7, 10, 15, 20, 30, 40, 60};
static int nBins_pur = 0;
static const int kMaxBinsPur = 50;
static long long nTagAll[kMaxBinsPur];
static long long nTagCorrect[kMaxBinsPur];
// -- FRAGMENTATION PHOTON PURITY ARRAYS:
static long long nTagAllFrag[kMaxBinsPur];
static long long nTagCorrectFrag[kMaxBinsPur];

// -------------------------------------------------------------------
// We'll define 6 categories for clusters:
//
//   0: Prompt photon cluster
//   1: Frag photon cluster
//   2: Decay photon cluster
//   3: pi0-led cluster (PID=111)
//   4: eta-led cluster (PID=221)
//   5: other
//
// catCount[cat][bin], where bin is the purity bin in [5..60] range
// -------------------------------------------------------------------
static long long catCount[6][kMaxBinsPur];

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

// -------------------------------------------------------------------
void buildAndSaveIsoEfficiencyPlot_2Classes(const std::string & outDir)
{
    std::cout << "[INFO] Building iso-efficiency for {prompt,frag} × 4 pT bins...\n";

    Color_t classColors[2] = { kRed+1, kBlue+1 };
    Style_t binMarkers[4]  = {20, 21, 24, 25};

    TCanvas cIso("cIso","Iso Efficiency: prompt vs frag in 4 pT bins",800,600);
    cIso.SetGrid();
    TH1F* hFrame = cIso.DrawFrame(0., 0.93, 20., 1.02);
    hFrame->SetXTitle("Truth isoE_{T} cut [GeV]");
    hFrame->SetYTitle("Efficiency");
    hFrame->SetTitle("");

    TGraphErrors* graphs[2][4];
    std::memset(graphs, 0, sizeof(graphs));

    for(int cIdx=0; cIdx<2; cIdx++){
        for(int b=0; b<nBins; b++){
            long long denom = nPhotonsInBin[cIdx][b];
            if(denom<1) continue;

            std::vector<double> vx(nIsoCuts,0.);
            std::vector<double> vy(nIsoCuts,0.);
            std::vector<double> vex(nIsoCuts,0.);
            std::vector<double> vey(nIsoCuts,0.);

            for(int ic=0; ic<nIsoCuts; ic++){
                vx[ic] = isoCutsGeV[ic];
                long long num = isoCountByClass[cIdx][b][ic];
                double eff=0., err=0.;
                if(denom>0){
                    eff = double(num)/ double(denom);
                    err = std::sqrt(eff*(1.-eff)/ double(denom));
                }
                vy[ic]  = eff;
                vey[ic] = err;
            }
            TGraphErrors* gr = new TGraphErrors(nIsoCuts, &vx[0], &vy[0],
                                                &vex[0], &vey[0]);
            gr->SetLineColor(classColors[cIdx]);
            gr->SetMarkerColor(classColors[cIdx]);
            gr->SetMarkerStyle(binMarkers[b]);
            gr->SetMarkerSize(1.3);
            gr->Draw("P SAME");
            graphs[cIdx][b] = gr;
        }
    }

    TLegend leg(0.6, 0.15, 0.85, 0.50);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);

    for(int b=0; b<nBins; b++){
        TGraph* dummyBin = new TGraph();
        dummyBin->SetMarkerStyle(binMarkers[b]);
        dummyBin->SetMarkerColor(kBlack);
        float ptLo = pTedges[b];
        float ptHi = pTedges[b+1];
        TString labelBin = Form("%.0f < p_{T} < %.0f GeV", ptLo, ptHi);
        leg.AddEntry(dummyBin, labelBin, "p");
    }
    {
        TGraph* d0 = new TGraph();
        d0->SetMarkerStyle(20);
        d0->SetMarkerColor(classColors[0]);
        leg.AddEntry(d0, "Prompt (#gamma_{class}=1)", "p");

        TGraph* d1 = new TGraph();
        d1->SetMarkerStyle(20);
        d1->SetMarkerColor(classColors[1]);
        leg.AddEntry(d1, "Frag (#gamma_{class}=2)", "p");
    }
    leg.Draw();

    std::string outName = outDir + "/IsoEfficiency_PromptVsFrag_4Bins.png";
    cIso.SaveAs(outName.c_str());
    std::cout << "[INFO] Wrote: " << outName << std::endl;
}

// -------------------------------------------------------------------
void analyzeJet10SlimTree()
{
    // 1) Input & Output
    std::string inputFile = "/Users/patsfan753/Desktop/caloana0130.root";
    std::string outDir    = "/Users/patsfan753/Desktop/SimOut";
    ensureOutputDirectory(outDir);

    // 2) Check TTree
    {
        TFile* fTest = TFile::Open(inputFile.c_str(),"READ");
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
        for(int iB=0; iB < branchList->GetEntries(); iB++){
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
    
    // Also zero out the fragmentation arrays:
    std::memset(nTagAllFrag,     0, sizeof(nTagAllFrag));
    std::memset(nTagCorrectFrag, 0, sizeof(nTagCorrectFrag));

    // 5) Initialize counters for the 6-category breakdown
    std::memset(catCount, 0, sizeof(catCount));

    // 6) Open input file & retrieve TTree
    TFile* fIn = TFile::Open(inputFile.c_str(),"READ");
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

    // 7) Set branch addresses
    //    Particles => iso-efficiency
    slim->SetBranchAddress("nparticles",         &nparticles);
    slim->SetBranchAddress("particle_pid",       particle_pid);
    slim->SetBranchAddress("particle_photonclass", particle_photonclass);
    slim->SetBranchAddress("particle_Pt",        particle_Pt);
    slim->SetBranchAddress("particle_Eta",       particle_Eta_);
    slim->SetBranchAddress("particle_Phi",       particle_Phi_);
    slim->SetBranchAddress("particle_truth_iso_04", particle_truth_iso_04);

    //    Clusters => purity & categories
    slim->SetBranchAddress("ncluster_CLUSTERINFO_CEMC",             &ncluster);
    slim->SetBranchAddress("cluster_Et_CLUSTERINFO_CEMC",           cluster_Et);
    slim->SetBranchAddress("cluster_Eta_CLUSTERINFO_CEMC",          cluster_Eta);
    slim->SetBranchAddress("cluster_Phi_CLUSTERINFO_CEMC",          cluster_Phi);
    slim->SetBranchAddress("cluster_E_CLUSTERINFO_CEMC",            cluster_E);
    slim->SetBranchAddress("cluster_iso_04_CLUSTERINFO_CEMC",       cluster_iso_04);
    slim->SetBranchAddress("cluster_iso_04_emcal_CLUSTERINFO_CEMC", cluster_iso_04_emcal);
    slim->SetBranchAddress("cluster_pid_CLUSTERINFO_CEMC", cluster_pid);

    // 8) Main loop => fill iso-efficiency + purity + categories
    Long64_t nEntries = slim->GetEntries();
    std::cout << "[INFO] nEntries in TTree: " << nEntries << "\n";

    for(Long64_t iEvt=0; iEvt<nEntries; iEvt++){
        slim->GetEntry(iEvt);

        // Optionally limit debug prints to first 10 events
        bool doPrint = (kDebug && iEvt<10);

        if(doPrint){
            std::cout<<"\n----------------------\n";
            std::cout<<" [DEBUG] Event "<< iEvt
                     <<": nparticles="<< nparticles
                     <<", ncluster="<< ncluster <<"\n";
        }

        // (A) iso-efficiency from particle info
        for(int ip=0; ip<nparticles; ip++){
            if(particle_pid[ip] != 22) continue; // must be photon
            int cIdx = photonClassIndex(particle_photonclass[ip]); // 0=prompt,1=frag
            if(cIdx<0) continue; // skip
            float pt = particle_Pt[ip];
            float isoV= particle_truth_iso_04[ip];
            int bIso= findPtBin(pt); // [10..30]
            if(bIso<0) continue;

            nPhotonsInBin[cIdx][bIso]++;
            for(int ic=0; ic<nIsoCuts; ic++){
                if(isoV< isoCutsGeV[ic]){
                    isoCountByClass[cIdx][bIso][ic]++;
                }
            }
        }

        // (B) cluster-based loop => do prompt-purity + new 6-category
        for(int ic=0; ic<ncluster; ic++){
            float cPt    = cluster_Et[ic];
            float cEta   = cluster_Eta[ic];
            float cPhi   = cluster_Phi[ic];
            float cE     = cluster_E[ic];
            float isoAll = cluster_iso_04[ic];
            float isoEM  = cluster_iso_04_emcal[ic];
            int   pid    = cluster_pid[ic]; // leading PID

            int bPur= findPurityBin(cPt);

            // 1) Prompt purity => same logic as before
            if(bPur>=0){
                // Tag if passes prompt-candidate cuts
                bool passCuts = passClusterPromptCuts(isoAll, isoEM);
                if(passCuts){
                    nTagAll[bPur]++;
                    // check if truly matched to prompt photon
                    int mm= findMatchedPromptPhoton(cEta,cPhi,cE);
                    if(mm>=0){
                        nTagCorrect[bPur]++;
                    }
                }
            }
            // NOW do the same for "frag-candidate" clusters in the same bin
            // We'll reuse the same cluster cuts (isoAll, isoEM) or define a new pass if you want.
            if(bPur >= 0){
                bool passFragCuts = passClusterPromptCuts(isoAll, isoEM);
                if(passFragCuts){
                    nTagAllFrag[bPur]++;
                    int matchFrag = findMatchedFragPhoton(cEta, cPhi, cE);
                    if(matchFrag >= 0){
                        nTagCorrectFrag[bPur]++;
                    }
                }
            }

            // 2) 6-category classification => also uses same [5..60] bin
            if(bPur<0) {
                if(doPrint){
                    std::cout<<"   [DEBUG-cluster "<<ic<<"] cPt="<<cPt
                             <<" => outside [5..60] => skipping cat count.\n";
                }
                continue;
            }

            // Decide cat => 0..5
            int cat= -1;
            if(pid==22){
                // If leading PID=22 => match it => see if prompt,frag,decay, or other
                int phClass= findMatchedPhotonClass(cEta, cPhi, cE);
                if(phClass==1)      cat=0; // prompt
                else if(phClass==2) cat=1; // frag
                else if(phClass==3) cat=2; // decay
                else                cat=5; // other
            }
            else if(pid==111) cat=3; // pi0-led
            else if(pid==221) cat=4; // eta-led
            else              cat=5; // other

            catCount[cat][bPur]++;

            if(doPrint){
                std::cout<<"   [DEBUG-cluster "<<ic<<"] cPt="<<cPt
                         <<", pid="<<pid
                         <<", isoAll="<<isoAll
                         <<", isoEM="<<isoEM
                         <<" => cat="<<cat
                         <<" (0=prompt,1=frag,2=decay,3=pi0,4=eta,5=other)\n";
            }
        } // end cluster loop

        // Optional summary per event
        if(doPrint){
            std::cout<<" [DEBUG] End of event "<< iEvt <<" catCount so far (first 10 bins):\n";
            for(int cat=0; cat<6; cat++){
                long long sumCat=0;
                for(int ib=0; ib< std::min(nBins_pur,10); ib++){
                    sumCat += catCount[cat][ib];
                }
                std::cout<<"    cat="<<cat<<" => sumCat="<<sumCat<<"\n";
            }
        }

        if((iEvt>0) && (iEvt%100000==0)){
            std::cout << "[INFO] Processed " << iEvt << "/" << nEntries << " events...\n";
        }
    } // end event loop
    fIn->Close();

    // 9) Build iso-efficiency plot
    buildAndSaveIsoEfficiencyPlot_2Classes(outDir);

    // 10) Build prompt-purity plot
    {
        std::vector<double> vBinC(nBins_pur, 0.);
        std::vector<double> vPur(nBins_pur, 0.);
        std::vector<double> vErr(nBins_pur, 0.);
        double maxPur= 0.;

        for(int ib=0; ib<nBins_pur; ib++){
            double cCenter= 0.5*( pTedges_purEff[ib]+ pTedges_purEff[ib+1] );
            vBinC[ib]= cCenter;
            double denom= double(nTagAll[ib]);
            double num  = double(nTagCorrect[ib]);
            if(denom>0.){
                double pur= num/ denom;
                double ePur= std::sqrt(pur*(1.-pur)/ denom);
                vPur[ib]= pur;
                vErr[ib]= ePur;
                if(pur> maxPur) maxPur= pur;
            }
        }

        double yMax= (maxPur<0.01 ? 0.1 : 1.1* maxPur);
        if(yMax>1.0) yMax=1.0;

        TCanvas cPur("cPur","Prompt Photon Purity",800,600);
        gStyle->SetOptStat(0);
        cPur.SetGrid();

        double xMin= pTedges_purEff.front()-0.5;
        double xMax= pTedges_purEff.back() +0.5;
        TH1F* hFr= cPur.DrawFrame(xMin,0., xMax,yMax);
        hFr->SetXTitle("Cluster p_{T} [GeV]");
        hFr->SetYTitle("Purity = N_{correct}/N_{tagged}");
        hFr->SetTitle("Prompt Photon Purity vs. p_{T}");

        TGraphErrors* grPur= new TGraphErrors(nBins_pur,&vBinC[0],&vPur[0],nullptr,&vErr[0]);
        grPur->SetMarkerStyle(20);
        grPur->SetMarkerColor(kAzure+2);
        grPur->SetLineColor(kAzure+2);
        grPur->SetMarkerSize(1.3);
        grPur->Draw("P SAME");

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.032);
        latex.SetTextFont(42);
        latex.DrawLatex(0.3,0.45,
                        "Cuts: #DeltaR<0.05, E_{reco}/E_{truth}>0.5, iso_{All}<6, iso_{EM}<4");

        std::string outName= outDir + "/PromptPurity.png";
        cPur.SaveAs(outName.c_str());
        std::cout<<"[INFO] Wrote Purity plot => "<< outName <<std::endl;
    }

    // === NOW BUILD FRAG PHOTON PURITY ===
    {
        std::vector<double> vBinC_frag(nBins_pur, 0.);
        std::vector<double> vPur_frag(nBins_pur, 0.);
        std::vector<double> vErr_frag(nBins_pur, 0.);
        double maxPurFrag = 0.;

        for(int ib=0; ib<nBins_pur; ib++){
            double cCenter = 0.5*( pTedges_purEff[ib] + pTedges_purEff[ib+1] );
            vBinC_frag[ib] = cCenter;
            double denom = double(nTagAllFrag[ib]);
            double num   = double(nTagCorrectFrag[ib]);
            if(denom > 0.){
                double pur = num / denom;
                double ePur= std::sqrt( pur*(1.-pur)/ denom );
                vPur_frag[ib]    = pur;
                vErr_frag[ib]    = ePur;
                if(pur> maxPurFrag) maxPurFrag= pur;
            }
        }

        double yMaxFrag= (maxPurFrag<0.01 ? 0.1 : 1.1*maxPurFrag);
        if(yMaxFrag>1.0) yMaxFrag=1.0;

        TCanvas cPurFrag("cPurFrag","Frag Photon Purity",800,600);
        cPurFrag.SetGrid();
        gStyle->SetOptStat(0);

        double xMin = pTedges_purEff.front()-0.5;
        double xMax = pTedges_purEff.back() +0.5;
        TH1F* hFrFrag= cPurFrag.DrawFrame(xMin,0., xMax,yMaxFrag);
        hFrFrag->SetXTitle("Cluster p_{T} [GeV]");
        hFrFrag->SetYTitle("Purity = N_{correct}/N_{tagged}");
        hFrFrag->SetTitle("Fragmentation Photon Purity vs. p_{T}");

        TGraphErrors* grPurFrag = new TGraphErrors(nBins_pur,
            &vBinC_frag[0], &vPur_frag[0], nullptr, &vErr_frag[0]);
        grPurFrag->SetMarkerStyle(21);
        grPurFrag->SetMarkerColor(kOrange+1);
        grPurFrag->SetLineColor(kOrange+1);
        grPurFrag->SetMarkerSize(1.3);
        grPurFrag->Draw("P SAME");

        TLatex latexFrag;
        latexFrag.SetNDC();
        latexFrag.SetTextSize(0.032);
        latexFrag.SetTextFont(42);
        latexFrag.DrawLatex(0.3,0.45,
            "Cuts: #DeltaR<0.05, E_{reco}/E_{truth}>0.5, iso_{All}<6, iso_{EM}<4");

        std::string outNameFrag= outDir + "/FragPurity.png";
        cPurFrag.SaveAs(outNameFrag.c_str());
        std::cout<<"[INFO] Wrote frag purity plot => "<< outNameFrag <<std::endl;
    }
    
    // 11) Build the single bar chart for 6 categories => pT-independent
    {
        // Summation
        long long catTotal[6];
        std::memset(catTotal, 0, sizeof(catTotal));
        long long totalClusters=0;

        for(int cat=0; cat<6; cat++){
            for(int ib=0; ib<nBins_pur; ib++){
                catTotal[cat]+= catCount[cat][ib];
            }
            totalClusters += catTotal[cat];
        }

        if(kDebug){
            std::cout<<"[DEBUG] Summation of catCount across pT bins:\n";
            for(int cat=0; cat<6; cat++){
                std::cout<<"   cat="<<cat<<" => total="<< catTotal[cat] <<"\n";
            }
            std::cout<<"   totalClusters="<< totalClusters <<"\n";
        }

        const char* catLabels[6]={
          "Prompt #gamma","Frag #gamma","Decay #gamma",
          "#pi^{0}-led","#eta-led","Other"
        };
        Color_t catColors[6]={
          kGreen+2, kBlue+1, kMagenta, kOrange+1, kAzure+1, kGray+1
        };

        TH1F* hBar = new TH1F("hBar","Cluster Categories (pT-indep)",6,0.5,6.5);
        hBar->SetYTitle("Total Clusters in Category/Total Clusters");
        hBar->SetStats(0);
        // Increase the x-axis label and title sizes:
        hBar->GetXaxis()->SetLabelSize(0.05);   // Set label size to 0.05 (adjust as desired)
        hBar->GetXaxis()->SetTitleSize(0.05);   // Set title size to 0.05 (optional)

        double maxFrac = 0.;
        for(int cat=0; cat<6; cat++){
            double frac = 0.;
            if(totalClusters>0) {
                frac = double(catTotal[cat]) / double(totalClusters);
            }
            hBar->SetBinContent(cat+1, frac);
            if(frac > maxFrac) maxFrac = frac;
            // Label each bin as before
            hBar->GetXaxis()->SetBinLabel(cat+1, catLabels[cat]);
        }

        // 2) Draw the histogram with no fill, just as an axis/frame
        TCanvas cCat("cCat","6-Category pT-independent Breakdown",800,600);
        cCat.SetGrid();

        double yMax= std::min(1.0, 1.1* maxFrac);
        if(yMax< 0.2) yMax=0.2;
        hBar->SetMaximum(yMax);

        // Draw the histogram axes/labels only (no fill):
        hBar->SetFillStyle(0);
        hBar->SetLineColor(kBlack);
        hBar->Draw("hist");

        double barWidth = 1.0; // Fill the entire bin (the bin width is exactly 1)
        for(int cat=0; cat<6; cat++){
            double frac = hBar->GetBinContent(cat+1);
            double xCenter = hBar->GetBinCenter(cat+1);

            // The bin extends from xCenter - 0.5 to xCenter + 0.5 if barWidth=1.0:
            double xLow  = xCenter - 0.5*barWidth;
            double xHigh = xCenter + 0.5*barWidth;

            // Create the box, fill it, and remove the black outline
            TBox* box = new TBox(xLow, 0., xHigh, frac);
            box->SetFillColor(catColors[cat]);
            box->SetLineColor(catColors[cat]); // match fill color or kWhite
            box->SetLineWidth(0);             // or set to 1 if you want a colored border
            box->Draw("same");
        }

        // 5) Save the final result
        std::string outCat = outDir + "/ClusterCategories_6Way_SingleBar.png";
        cCat.SaveAs(outCat.c_str());
        std::cout << "[INFO] Wrote single-bar 6-category breakdown => " << outCat << std::endl;

        // Also print numeric results
        std::cout<<"\n===== pT-Independent 6-Category Breakdown =====\n";
        std::cout<<"  total clusters across all bins="<<totalClusters<<"\n";
        for(int cat=0; cat<6; cat++){
            double frac= (totalClusters>0 ?
                          double(catTotal[cat])/double(totalClusters) : 0.);
            std::cout<<"    cat="<<cat<<" ("<< catLabels[cat] <<") => "
                     << catTotal[cat] <<" => fraction="<< frac <<"\n";
        }
        std::cout<<"===============================================\n\n";
        
        {
            // Loop over each purity bin => produce a separate bar chart
            for(int ib=0; ib<nBins_pur; ib++){
                // 1) Sum how many clusters in this bin
                long long binTotal = 0;
                for(int cat=0; cat<6; cat++){
                    binTotal += catCount[cat][ib];
                }

                TH1F* hBarBin = new TH1F(Form("hBar_bin%d", ib),
                                         Form("Cluster Categories: %.0f #leq p_{T} < %.0f GeV",
                                              pTedges_purEff[ib], pTedges_purEff[ib+1]),
                                         6, 0.5, 6.5);
                hBarBin->SetYTitle("Fraction");
                hBarBin->SetStats(0);
                // Increase the x-axis label and title sizes:
                hBarBin->GetXaxis()->SetLabelSize(0.05);   // Increase label size for the x-axis
                hBarBin->GetXaxis()->SetTitleSize(0.05);   // Increase title size for the x-axis (optional)

                // 3) Fill each category's fraction
                double maxFracBin=0.;
                for(int cat=0; cat<6; cat++){
                    double frac = 0.;
                    if(binTotal>0){
                        frac = double(catCount[cat][ib]) / double(binTotal);
                    }
                    hBarBin->SetBinContent(cat+1, frac);
                    if(frac>maxFracBin) maxFracBin= frac;
                    hBarBin->GetXaxis()->SetBinLabel(cat+1, catLabels[cat]);
                }

                // 4) Draw it with TBoxes for color
                TCanvas * cCatBin = new TCanvas(Form("cCat_bin%d",ib),
                    Form("Cluster Breakdown in pT bin %d", ib), 800,600);
                cCatBin->SetGrid();

                double yMaxBin = std::min(1.0, 1.1* maxFracBin);
                if(yMaxBin<0.2) yMaxBin=0.2;
                hBarBin->SetMaximum(yMaxBin);

                // Draw the axis/labels only
                hBarBin->SetFillStyle(0);
                hBarBin->SetLineColor(kBlack);
                hBarBin->Draw("hist");

                double barWidth=1.0;
                for(int cat=0; cat<6; cat++){
                    double frac= hBarBin->GetBinContent(cat+1);
                    double xCenter= hBarBin->GetBinCenter(cat+1);
                    double xLow= xCenter - 0.5* barWidth;
                    double xHigh= xCenter + 0.5* barWidth;

                    TBox* box = new TBox(xLow, 0., xHigh, frac);
                    box->SetFillColor(catColors[cat]);
                    box->SetLineColor(catColors[cat]);
                    box->SetLineWidth(0);
                    box->Draw("same");
                }

                double binLo= pTedges_purEff[ib];
                double binHi= pTedges_purEff[ib+1];

                // 6) Save
                std::string outBinName=
                  outDir + Form("/ClusterCat_Bin%d_%.0fto%.0f.png", ib, binLo, binHi);
                cCatBin->SaveAs(outBinName.c_str());
                std::cout<<"[INFO] Wrote bin="<<ib<<" bar chart => "<< outBinName <<"\n";
            }
        }
    }

    std::cout << "[INFO] All done.\n";
}
