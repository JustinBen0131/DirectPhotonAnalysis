#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TFitResultPtr.h>
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

//______________________________________________________________________________
// Maximum possible arrays
static const int kMaxClusters   = 200000;
static const int kMaxParticles  = 200000;

//______________________________________________________________________________
// Branch variables
int   ncluster;
int   cluster_pid[kMaxClusters];        // PDG code for cluster's "leading" particle
float cluster_Et[kMaxClusters];         // Reconstructed cluster pT
float cluster_E[kMaxClusters];          // Reconstructed cluster E
float cluster_Eta[kMaxClusters];        // Reconstructed cluster Eta
float cluster_Phi[kMaxClusters];        // Reconstructed cluster Phi
float cluster_iso_04[kMaxClusters];     // isolation sum (R=0.4)
float cluster_iso_04_emcal[kMaxClusters];

int   nparticles;
int   particle_pid[kMaxParticles];         // PDG code for truth particles
int   particle_photonclass[kMaxParticles]; // 1=prompt,2=frag,3=decay,4=other
float particle_Pt[kMaxParticles];
float particle_Eta_[kMaxParticles];
float particle_Phi_[kMaxParticles];

//______________________________________________________________________________
// pT bin edges
static std::vector<float> pTedges = {5,7,10,15,20,30,45};

// We'll store the fraction-of-pi0/eta counters
static std::vector<long long> nAllClustersInBin;    // total clusters in bin
static std::vector<long long> nPi0EtaInBin;         // clusters truly from pi0 or eta

// For the prompt-photon purity (original method)
static std::vector<long long> nPromptTagCorrect;    // truly prompt & tagged
static std::vector<long long> nPromptTagFake;       // not prompt but tagged

// For the *enhanced* prompt-photon purity
static std::vector<long long> nPromptTagCorrectEnhanced;
static std::vector<long long> nPromptTagFakeEnhanced;

// We'll also store how many clusters are flagged by *our invariant-mass approach*
// as presumably from meson
static std::vector<long long> nFlaggedAsMesonInBin;

//______________________________________________________________________________
// Utility: ensure output directory exists
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

//______________________________________________________________________________
// Utility: find pT bin index for a given pT
int findPtBin(float pt)
{
    for(int i=0; i<(int)pTedges.size()-1; i++){
        if(pt >= pTedges[i] && pt < pTedges[i+1]) return i;
    }
    return -1; // out of range
}

//______________________________________________________________________________
// DeltaR function for cluster-particle matching
float deltaR_func(float eta1, float eta2, float phi1, float phi2)
{
    float dEta = eta1 - eta2;
    float dPhi = phi1 - phi2;
    while(dPhi >  TMath::Pi()) dPhi -= 2.f*TMath::Pi();
    while(dPhi < -TMath::Pi()) dPhi += 2.f*TMath::Pi();
    return std::sqrt(dEta*dEta + dPhi*dPhi);
}

//______________________________________________________________________________
// Example cut-based approach to “prompt photon candidate”
bool passPromptCandidateCuts(float cEt, float isoAll, float isoEM)
{
    if(cEt < 5.f) {
        return false;  // minimal pT
    }
    if(isoAll > 6.f) {
        return false;  // total hadronic isolation
    }
    if(isoEM  > 4.f) {
        return false;  // EMCal isolation
    }
    // Could add shower-shape cuts, etc.
    return true;
}

//______________________________________________________________________________
// Utility to format a float with 3 sig figs
std::string formatToThreeSigFigs(double val)
{
    if(val == 0.0) return "0.00";
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(2) << val;
    return oss.str();
}

//______________________________________________________________________________
// Data structure to hold the final fit results (mass windows) for each pT bin
struct FitResultData {
    int    binIndex;
    double pTmin;
    double pTmax;

    double meanPi0;
    double sigmaPi0;
    double meanEta;
    double sigmaEta;
};

//______________________________________________________________________________
// Detailed function that does: gaus(pi0) + gaus(eta) + pol4
TFitResultPtr PerformFitting(
    TH1* hMass,
    TF1*& totalFit, TF1*& gaussPi0Fit, TF1*& gaussEtaFit, TF1*& polyFit,
    double& fitStart, double& fitEnd)
{
    std::cout << "[DEBUG] Performing fit with double Gaussians + pol4..." << std::endl;
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);

    fitStart = 0.04;
    fitEnd   = 0.90;

    // Approx pi0, eta ranges
    double pi0Min = 0.10, pi0Max = 0.18;
    double etaMin = 0.50, etaMax = 0.70;

    // 1) Identify approximate pi0 peak
    int binPi0Min = hMass->GetXaxis()->FindBin(pi0Min);
    int binPi0Max = hMass->GetXaxis()->FindBin(pi0Max);
    int maxBinPi0 = binPi0Min;
    double maxPi0Content = 0.0;
    for(int b=binPi0Min; b<=binPi0Max; ++b){
        double c = hMass->GetBinContent(b);
        if(c>maxPi0Content){
            maxPi0Content = c;
            maxBinPi0 = b;
        }
    }
    double meanPi0Estimate = hMass->GetXaxis()->GetBinCenter(maxBinPi0);

    // 2) Identify approximate eta peak
    int binEtaMin = hMass->GetXaxis()->FindBin(etaMin);
    int binEtaMax = hMass->GetXaxis()->FindBin(etaMax);
    int maxBinEta = binEtaMin;
    double maxEtaContent = 0.0;
    for(int b=binEtaMin; b<=binEtaMax; ++b){
        double c = hMass->GetBinContent(b);
        if(c>maxEtaContent){
            maxEtaContent = c;
            maxBinEta = b;
        }
    }
    double meanEtaEstimate = hMass->GetXaxis()->GetBinCenter(maxBinEta);

    // Build a double‐Gaussian + pol4
    totalFit = new TF1("totalFit","gaus(0) + gaus(3) + pol4(6)", fitStart, fitEnd);
    totalFit->SetLineColor(kRed);

    // Pi0 initial guesses
    totalFit->SetParameter(0, maxPi0Content);
    totalFit->SetParameter(1, meanPi0Estimate);
    totalFit->SetParameter(2, 0.012); // sigma guess
    totalFit->SetParLimits(2, 0.008, 0.02);

    // Eta initial guesses
    totalFit->SetParameter(3, maxEtaContent);
    totalFit->SetParameter(4, meanEtaEstimate);
    totalFit->SetParameter(5, 0.03); // sigma guess
    totalFit->SetParLimits(5, 0.02, 0.05);

    // pol4 init
    for(int p=6; p<11; p++){
        totalFit->SetParameter(p,0.0);
    }

    std::cout << "[DEBUG] Starting TFitResultPtr on histogram: " << hMass->GetName() << std::endl;
    TFitResultPtr fitResult = hMass->Fit(totalFit,"SR");
    int fitStatus = fitResult;
    if(fitStatus!=0){
        std::cerr<<"[WARNING] Fit Status="<<fitStatus<<" => may not have converged properly.\n";
    }

    gaussPi0Fit = new TF1("gaussPi0Fit","gaus", fitStart, fitEnd);
    gaussPi0Fit->SetParameters(totalFit->GetParameter(0),
                               totalFit->GetParameter(1),
                               totalFit->GetParameter(2));
    gaussPi0Fit->SetLineColor(kBlue);
    gaussPi0Fit->SetLineStyle(2);

    gaussEtaFit = new TF1("gaussEtaFit","gaus", fitStart, fitEnd);
    gaussEtaFit->SetParameters(totalFit->GetParameter(3),
                               totalFit->GetParameter(4),
                               totalFit->GetParameter(5));
    gaussEtaFit->SetLineColor(kGreen+2);
    gaussEtaFit->SetLineStyle(2);

    polyFit = new TF1("polyFit","pol4", fitStart, fitEnd);
    for(int p=6; p<11; p++){
        polyFit->SetParameter(p-6, totalFit->GetParameter(p));
    }
    polyFit->SetLineColor(kOrange+7);
    polyFit->SetLineStyle(2);

    return fitResult;
}

//______________________________________________________________________________
// Draw text onto the canvas indicating final mass windows
void DrawInvMassCanvasText(double pTmin, double pTmax,
                           double meanPi0, double sigmaPi0,
                           double meanEta, double sigmaEta)
{
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.SetTextColor(kRed+3);
    latex.SetTextFont(62);

    latex.DrawLatex(0.20, 0.93, Form("pT bin: %.1f - %.1f GeV", pTmin, pTmax));
    latex.DrawLatex(0.20, 0.88, Form("#pi^{0} => mean=%.3f, #sigma=%.3f", meanPi0, sigmaPi0));
    latex.DrawLatex(0.20, 0.83, Form("#eta  => mean=%.3f, #sigma=%.3f", meanEta, sigmaEta));
}

//______________________________________________________________________________
// MAIN LOGIC
void analyzeSimulationLocally()
{
    bool doInvMass = true;
    
    std::cout << "[INFO] Starting analyzeSimulationLocally()...\n";

    // Input & output
    /*
     /sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run22/photon10/condorout_waveform/
     */
    std::string inputFile = "/Users/patsfan753/Desktop/caloana0130.root";
    /*
     /sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/
     */
    std::string outDir    = "/Users/patsfan753/Desktop/SimOut";
    ensureOutputDirectory(outDir);

    int nBins = (int)pTedges.size() - 1;
    std::cout << "[INFO] We have " << nBins << " pT bins.\n";

    // Histos for the di-cluster mass => [leading pT bin]
    std::vector<TH1F*> hInvMassLeading(nBins, nullptr);
    for(int i=0; i<nBins; i++){
        float low  = pTedges[i];
        float high = pTedges[i+1];
        std::string hName  = Form("hInvMassLeading_bin%d", i);
        std::string hTitle = Form("Dicluster InvMass; M_{#gamma#gamma} (GeV); Entries [%.1f,%.1f) GeV", low, high);
        hInvMassLeading[i] = new TH1F(hName.c_str(), hTitle.c_str(), 100, 0., 1.0);
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // (A) FIRST PASS: Fill di-cluster IM distributions using proper TLorentzVectors
    //                 *only* for pairs passing E>1 GeV and asym<0.5 cuts
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    {
        TFile* fIn = TFile::Open(inputFile.c_str(),"READ");
        if(!fIn || fIn->IsZombie()){
            std::cerr << "[ERROR] cannot open file: " << inputFile << std::endl;
            return;
        }
        TTree* slim = (TTree*) fIn->Get("slimtree");
        if(!slim){
            std::cerr<<"[ERROR] 'slimtree' not found. Closing.\n";
            fIn->Close();
            return;
        }

        slim->SetBranchAddress("ncluster_CLUSTERINFO_CEMC",       &ncluster);
        slim->SetBranchAddress("cluster_Et_CLUSTERINFO_CEMC",     cluster_Et);
        slim->SetBranchAddress("cluster_E_CLUSTERINFO_CEMC",      cluster_E);
        slim->SetBranchAddress("cluster_Eta_CLUSTERINFO_CEMC",    cluster_Eta);
        slim->SetBranchAddress("cluster_Phi_CLUSTERINFO_CEMC",    cluster_Phi);

        Long64_t nEntries = slim->GetEntries();
        std::cout << "[INFO] First pass => TTree has " << nEntries << " entries.\n";

        for(Long64_t iev=0; iev<nEntries; iev++){
            if( (iev % 200000 == 0) && iev>0 ){
                std::cout << "[DEBUG] First pass: processed " << iev
                          << " / " << nEntries << " events...\n";
            }
            slim->GetEntry(iev);
            if(ncluster < 2) continue;

            // Loop over all cluster pairs
            for(int i=0; i<ncluster; i++){
                for(int j=i+1; j<ncluster; j++){
                    float pt1  = cluster_Et[i];
                    float pt2  = cluster_Et[j];
                    float eta1 = cluster_Eta[i];
                    float eta2 = cluster_Eta[j];
                    float phi1 = cluster_Phi[i];
                    float phi2 = cluster_Phi[j];
                    float e1   = cluster_E[i];
                    float e2   = cluster_E[j];

                    // (1) Energy cut: require each cluster E>1 GeV
                    if(e1 < 1.0 || e2 < 1.0) continue;

                    // (2) Asymmetry cut: (|E1-E2|)/(E1+E2) < 0.5
                    float sumE = e1 + e2;
                    if(sumE < 1e-6) continue;  // avoid div-by-zero
                    float asym = std::fabs(e1 - e2)/sumE;
                    if(asym > 0.5) continue;

                    // leading cluster => find pT bin
                    float leadingPt = (pt1 >= pt2 ? pt1 : pt2);
                    int binIdx = findPtBin(leadingPt);
                    if(binIdx < 0) continue; // out of range

                    // Build Lorentz vectors
                    TLorentzVector v1, v2;
                    v1.SetPtEtaPhiE(pt1, eta1, phi1, e1);
                    v2.SetPtEtaPhiE(pt2, eta2, phi2, e2);

                    float invMass = (v1 + v2).M();
                    hInvMassLeading[binIdx]->Fill(invMass);
                }
            }
        }
        fIn->Close();
        std::cout << "[INFO] Done filling di-cluster IM distributions (E>1GeV, asym<0.5).\n";
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // (B) FIT each IM distribution => store results,
    //     then output each individually AND on a single multi‐pad canvas
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    std::vector<FitResultData> fitParams(nBins);
    for(int ib=0; ib<nBins; ib++){
        fitParams[ib].binIndex = ib;
        fitParams[ib].pTmin    = pTedges[ib];
        fitParams[ib].pTmax    = pTedges[ib+1];
        fitParams[ib].meanPi0  = 0.;
        fitParams[ib].sigmaPi0 = 0.;
        fitParams[ib].meanEta  = 0.;
        fitParams[ib].sigmaEta = 0.;
    }

    std::cout << "[INFO] Now performing fits for each pT bin and producing plots...\n";

    // Create big canvas with sub‐pads for all bins, e.g. 3 wide x N rows
    int nCols = 3;
    int nRows = (nBins + nCols - 1) / nCols;
    TCanvas cAll("cAll","Invariant Mass Distributions by pT Bin",1200,800);
    cAll.Divide(nCols, nRows);

    for(int i=0; i<nBins; i++){
        TH1F* hist = hInvMassLeading[i];
        if(!hist){
            std::cout << "[WARN] null hist pointer for bin " << i << std::endl;
            continue;
        }
        if(hist->GetEntries() < 10){
            std::cout << "[WARN] " << hist->GetName()
                      << " has <10 entries => skipping fit.\n";
            continue;
        }

        // Fit
        TF1 *totalF=nullptr, *gaussPi0=nullptr, *gaussEta=nullptr, *polyF=nullptr;
        double fStart=0., fEnd=1.;
        TFitResultPtr fitRes = PerformFitting(hist, totalF, gaussPi0, gaussEta, polyF, fStart, fEnd);

        // Extract parameters
        double mp0 = totalF->GetParameter(1);
        double sp0 = totalF->GetParameter(2);
        double met = totalF->GetParameter(4);
        double set = totalF->GetParameter(5);
        fitParams[i].meanPi0  = mp0;
        fitParams[i].sigmaPi0 = sp0;
        fitParams[i].meanEta  = met;
        fitParams[i].sigmaEta = set;

        // (i) Save an individual plot for this bin
        TCanvas cInd(Form("cInvMassBin%d", i),
                     Form("InvMass Leading pT = %.1f-%.1f GeV", pTedges[i], pTedges[i+1]),
                     800,600);
        hist->Draw("E");
        if(gaussPi0) gaussPi0->Draw("SAME");
        if(gaussEta) gaussEta->Draw("SAME");
        if(polyF)    polyF->Draw("SAME");
        if(totalF)   totalF->Draw("SAME");
        DrawInvMassCanvasText(pTedges[i], pTedges[i+1], mp0, sp0, met, set);

        std::string outName = Form("%s/InvMassLeading_%.1fto%.1f.png",
                                   outDir.c_str(), pTedges[i], pTedges[i+1]);
        cInd.SaveAs(outName.c_str());
        std::cout << "[INFO] Saved bin " << i
                  << " => " << outName << std::endl;

        // (ii) Also put this histogram on the multi‐pad canvas
        cAll.cd(i+1);
        hist->Draw("E");
        if(gaussPi0) gaussPi0->Draw("SAME");
        if(gaussEta) gaussEta->Draw("SAME");
        if(polyF)    polyF->Draw("SAME");
        if(totalF)   totalF->Draw("SAME");
        DrawInvMassCanvasText(pTedges[i], pTedges[i+1], mp0, sp0, met, set);
    }

    // Save the multi‐bin figure
    std::string outMulti = outDir + "/InvMassLeading_AllBins.png";
    cAll.SaveAs(outMulti.c_str());
    std::cout << "[INFO] Saved multi‐bin IM plot: " << outMulti << std::endl;

    // Optionally write final fit results to CSV
    {
        std::string csvName = outDir + "/InvMassFitResults.csv";
        std::ofstream ofs(csvName);
        if(!ofs.good()){
            std::cerr<<"[WARNING] Could not open CSV file: "<<csvName<<"\n";
        } else {
            ofs << "BinIndex,pTmin,pTmax,meanPi0,sigmaPi0,meanEta,sigmaEta\n";
            for(auto &fp : fitParams) {
                ofs << fp.binIndex << ","
                    << fp.pTmin << "," << fp.pTmax << ","
                    << fp.meanPi0 << "," << fp.sigmaPi0 << ","
                    << fp.meanEta << "," << fp.sigmaEta << "\n";
            }
            ofs.close();
            std::cout << "[INFO] Wrote fit results to " << csvName << std::endl;
        }
    }
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    
    if(doInvMass) {
        std::cout << "[INFO] doInvMass flag is true. Skipping second pass and further analysis.\n";
        return;
    }

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // (C) SECOND PASS: fraction-of-pi0/eta, mark clusters in +/- 2sigma,
    //     do original & enhanced purity, plus track how many clusters
    //     are flagged as meson by IM approach
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    nAllClustersInBin.resize(nBins, 0);
    nPi0EtaInBin.resize(nBins, 0);
    nPromptTagCorrect.resize(nBins, 0);
    nPromptTagFake.resize(nBins, 0);
    nPromptTagCorrectEnhanced.resize(nBins, 0);
    nPromptTagFakeEnhanced.resize(nBins, 0);
    nFlaggedAsMesonInBin.resize(nBins, 0);

    long long totalClustersProcessed=0;
    long long totalEventsProcessed=0;

    {
        TFile* fIn = TFile::Open(inputFile.c_str(),"READ");
        if(!fIn || fIn->IsZombie()){
            std::cerr<<"[ERROR] cannot open file: "<<inputFile<<"\n";
            return;
        }
        TTree* slim = (TTree*) fIn->Get("slimtree");
        if(!slim){
            std::cerr<<"[ERROR] 'slimtree' not found in second pass. Closing.\n";
            fIn->Close();
            return;
        }
        slim->SetBranchAddress("ncluster_CLUSTERINFO_CEMC",&ncluster);
        slim->SetBranchAddress("cluster_pid_CLUSTERINFO_CEMC",cluster_pid);
        slim->SetBranchAddress("cluster_Et_CLUSTERINFO_CEMC",  cluster_Et);
        slim->SetBranchAddress("cluster_E_CLUSTERINFO_CEMC",   cluster_E);
        slim->SetBranchAddress("cluster_Eta_CLUSTERINFO_CEMC", cluster_Eta);
        slim->SetBranchAddress("cluster_Phi_CLUSTERINFO_CEMC", cluster_Phi);
        slim->SetBranchAddress("cluster_iso_04_CLUSTERINFO_CEMC", cluster_iso_04);
        slim->SetBranchAddress("cluster_iso_04_emcal_CLUSTERINFO_CEMC", cluster_iso_04_emcal);

        slim->SetBranchAddress("nparticles",&nparticles);
        slim->SetBranchAddress("particle_pid", particle_pid);
        slim->SetBranchAddress("particle_photonclass", particle_photonclass);
        slim->SetBranchAddress("particle_Pt",  particle_Pt);
        slim->SetBranchAddress("particle_Eta", particle_Eta_);
        slim->SetBranchAddress("particle_Phi", particle_Phi_);

        Long64_t nEntries2 = slim->GetEntries();
        std::cout << "[INFO] Second pass => TTree has " << nEntries2 << " entries.\n";

        for(Long64_t iev=0; iev<nEntries2; iev++){
            if( (iev % 200000 == 0) && iev>0 ) {
                std::cout << "[DEBUG] Second pass: processed " << iev << " / "
                          << nEntries2 << " events so far...\n";
            }
            slim->GetEntry(iev);
            totalEventsProcessed++;

            if(ncluster<0 || ncluster>kMaxClusters){
                std::cerr << "[WARN] Invalid ncluster="<<ncluster<<" at event "<<iev<<". Skipping.\n";
                continue;
            }
            totalClustersProcessed += ncluster;

            // We'll track which clusters are "from meson" by the IM approach
            std::vector<bool> isFromMesonWindow(ncluster,false);

            // Pairwise check => mark clusters
            if(ncluster>=2){
                for(int i=0; i<ncluster; i++){
                    for(int j=i+1; j<ncluster; j++){
                        float pt1  = cluster_Et[i];
                        float pt2  = cluster_Et[j];
                        float eta1 = cluster_Eta[i];
                        float eta2 = cluster_Eta[j];
                        float phi1 = cluster_Phi[i];
                        float phi2 = cluster_Phi[j];
                        float e1   = cluster_E[i];
                        float e2   = cluster_E[j];

                        float leadingPt = (pt1>=pt2 ? pt1 : pt2);
                        int binL = findPtBin(leadingPt);
                        if(binL<0) continue;

                        double mPi0= fitParams[binL].meanPi0;
                        double sPi0= fitParams[binL].sigmaPi0;
                        double mEta= fitParams[binL].meanEta;
                        double sEta= fitParams[binL].sigmaEta;
                        if(sPi0<1e-6) sPi0=0.01;
                        if(sEta<1e-6) sEta=0.02;

                        // compute inv mass
                        TLorentzVector v1,v2;
                        v1.SetPtEtaPhiE(pt1,eta1,phi1,e1);
                        v2.SetPtEtaPhiE(pt2,eta2,phi2,e2);
                        float invM = (v1+v2).M();

                        bool inPi0Win = (fabs(invM - mPi0) < 2.*sPi0);
                        bool inEtaWin = (fabs(invM - mEta) < 2.*sEta);
                        if(inPi0Win || inEtaWin){
                            isFromMesonWindow[i]=true;
                            isFromMesonWindow[j]=true;
                        }
                    }
                }
            }

            // Single-cluster logic
            for(int ic=0; ic<ncluster; ic++){
                float cEt    = cluster_Et[ic];
                float cE     = cluster_E[ic];
                float cEta   = cluster_Eta[ic];
                float cPhi   = cluster_Phi[ic];
                float isoAll = cluster_iso_04[ic];
                float isoEM  = cluster_iso_04_emcal[ic];

                int binIdx = findPtBin(cEt);
                if(binIdx>=0){
                    // total clusters
                    nAllClustersInBin[binIdx]++;

                    // truly from pi0 or eta (based on cluster_pid==111 or 221)
                    if(cluster_pid[ic]==111 || cluster_pid[ic]==221){
                        nPi0EtaInBin[binIdx]++;
                    }
                    // if we flagged it by IM approach => increment
                    if(isFromMesonWindow[ic]){
                        nFlaggedAsMesonInBin[binIdx]++;
                    }
                }

                // Now the prompt-candidate logic
                bool isTaggedPrompt = passPromptCandidateCuts(cEt, isoAll, isoEM);
                if(!isTaggedPrompt) continue;

                // Is it truly prompt in MC?
                bool isRealPrompt=false;
                int bestMatch=-1;
                float bestDR=1e9;
                for(int ip=0; ip<nparticles; ip++){
                    if(particle_pid[ip]!=22) continue;
                    float dR= deltaR_func(cEta, particle_Eta_[ip], cPhi, particle_Phi_[ip]);
                    if(dR<0.05 && dR<bestDR){
                        bestDR=dR;
                        bestMatch=ip;
                    }
                }
                if(bestMatch>=0){
                    float truthPt= particle_Pt[bestMatch];
                    int photClass= particle_photonclass[bestMatch];
                    if(truthPt>0.f){
                        float ratio= cEt/truthPt;
                        if(ratio>0.5f && photClass==1){
                            isRealPrompt=true;
                        }
                    }
                }

                // Original purity
                if(binIdx>=0){
                    if(isRealPrompt) nPromptTagCorrect[binIdx]++;
                    else             nPromptTagFake[binIdx]++;
                }

                // Enhanced => skip clusters if flagged as meson
                if(isFromMesonWindow[ic]) {
                    continue;
                }
                if(binIdx>=0){
                    if(isRealPrompt) nPromptTagCorrectEnhanced[binIdx]++;
                    else             nPromptTagFakeEnhanced[binIdx]++;
                }
            }
        }
        fIn->Close();
        std::cout << "[INFO] Done with second pass.\n";
    }

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // (D) Build TGraphs for:
    //     (1) fraction of pi0/eta (truth) + fraction meson-by-IM
    //     (2) original vs enhanced purity
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    std::cout << "[INFO] Building final TGraphs...\n";

    // fraction of pi0/eta
    std::vector<double> vBinCenterFrac(nBins,0.),
                        vFracTruth(nBins,0.),      // fraction of *truth pi0/eta*
                        vErrFracTruth(nBins,0.);
    // fraction from IM approach
    std::vector<double> vFracIM(nBins,0.), vErrFracIM(nBins,0.);

    double maxFrac=0.;
    for(int i=0;i<nBins;i++){
        double cX= 0.5*(pTedges[i]+pTedges[i+1]);
        double tot= (double)nAllClustersInBin[i];
        double mes= (double)nPi0EtaInBin[i];
        double flagged= (double)nFlaggedAsMesonInBin[i];

        double fracTruth=0., errTruth=0.;
        double fracIM=0., errIM=0.;

        if(tot>0.){
            fracTruth = mes / tot;
            errTruth  = std::sqrt(fracTruth*(1.-fracTruth)/ tot);

            fracIM = flagged / tot;
            errIM  = std::sqrt(fracIM*(1.-fracIM)/ tot);
        }
        if(fracTruth>maxFrac) maxFrac=fracTruth;
        if(fracIM>maxFrac)    maxFrac=fracIM;

        vBinCenterFrac[i]= cX;
        vFracTruth[i]    = fracTruth;
        vErrFracTruth[i] = errTruth;
        vFracIM[i]       = fracIM;
        vErrFracIM[i]    = errIM;
    }

    // Purities
    std::vector<double> vBinCenterPur(nBins,0.), vPurOriginal(nBins,0.), vErrPurOriginal(nBins,0.);
    std::vector<double> vBinCenterPur2(nBins,0.), vPurEnhanced(nBins,0.), vErrPurEnhanced(nBins,0.);

    double maxPur=0.;
    for(int i=0;i<nBins;i++){
        double cX= 0.5*(pTedges[i]+pTedges[i+1]);
        long long c1= nPromptTagCorrect[i];
        long long f1= nPromptTagFake[i];
        double d1= (double)(c1+f1);
        double pur1=0., epur1=0.;
        if(d1>0.){
            pur1= (double)c1/d1;
            epur1= std::sqrt(pur1*(1.-pur1)/ d1);
        }
        if(pur1>maxPur) maxPur= pur1;

        long long c2= nPromptTagCorrectEnhanced[i];
        long long f2= nPromptTagFakeEnhanced[i];
        double d2= (double)(c2+f2);
        double pur2=0., epur2=0.;
        if(d2>0.){
            pur2= (double)c2/d2;
            epur2= std::sqrt(pur2*(1.-pur2)/ d2);
        }
        if(pur2>maxPur) maxPur= pur2;

        vBinCenterPur[i]     = cX;
        vPurOriginal[i]      = pur1;
        vErrPurOriginal[i]   = epur1;

        vBinCenterPur2[i]    = cX;
        vPurEnhanced[i]      = pur2;
        vErrPurEnhanced[i]   = epur2;
    }

    // (1) fraction of pi0/eta (truth) AND fraction from IM approach
    {
        // First, build the ratio arrays:
        std::vector<double> vRatio(nBins, 0.0);
        std::vector<double> vErrRatio(nBins, 0.0);

        for (int i = 0; i < nBins; i++) {
            double x = vFracIM[i];
            double ex = vErrFracIM[i];
            double y = vFracTruth[i];
            double ey = vErrFracTruth[i];
            if (y > 1e-12) {
                vRatio[i] = x / y;
                // standard propagation of errors for ratio = x/y
                // ratioErr = ratio * sqrt( (ex/x)^2 + (ey/y)^2 )
                double relErrSq = 0.0;
                if (x > 1e-12) relErrSq += std::pow(ex/x, 2);
                relErrSq += std::pow(ey/y, 2);
                vErrRatio[i] = vRatio[i] * std::sqrt(relErrSq);
            } else {
                vRatio[i]     = 0.0;
                vErrRatio[i]  = 0.0;
            }
        }

        // Create a taller canvas for main plot + ratio subplot
        TCanvas cFrac("cFrac","Fraction: truth vs. mass-window + ratio",800,800);

        // Top pad (70% height) for the overlayed fraction TGraphs
        TPad* padTop = new TPad("padTop","padTop", 0.0, 0.3, 1.0, 1.0);
        padTop->SetBottomMargin(0.02); // reduce gap
        padTop->Draw();
        padTop->cd();

        // Decide log scale or not
        bool useLogYFrac = (maxFrac < 0.02);
        if (useLogYFrac) padTop->SetLogy();

        TGraphErrors* grFracTruth = new TGraphErrors(nBins,
                                                     &vBinCenterFrac[0], &vFracTruth[0],
                                                     nullptr, &vErrFracTruth[0]);
        grFracTruth->SetTitle("");
        grFracTruth->SetMarkerStyle(20);
        grFracTruth->SetMarkerColor(kBlue+2);
        grFracTruth->SetLineColor(kBlue+2);

        TGraphErrors* grFracIM = new TGraphErrors(nBins,
                                                  &vBinCenterFrac[0], &vFracIM[0],
                                                  nullptr, &vErrFracIM[0]);
        grFracIM->SetMarkerStyle(21);
        grFracIM->SetMarkerColor(kOrange+1);
        grFracIM->SetLineColor(kOrange+1);

        // Draw the "Truth" first
        grFracTruth->Draw("AP");
        grFracTruth->GetHistogram()->GetXaxis()->SetRangeUser(pTedges.front()-0.5,
                                                              pTedges.back() + 0.5);
        grFracTruth->GetHistogram()->GetYaxis()->SetTitle("fraction");
        grFracTruth->GetHistogram()->SetTitle("Fraction of #pi^{0}/#eta (truth) vs. IM-flag");
        double minY = 1e-4 * maxFrac, maxY = 1.5 * maxFrac;
        if (minY < 1e-7) minY = 1e-7;
        if (!useLogYFrac) {
            grFracTruth->GetHistogram()->GetYaxis()->SetRangeUser(0., 1.2 * maxFrac);
        } else {
            grFracTruth->GetHistogram()->GetYaxis()->SetRangeUser(minY, maxY);
        }

        // Overplot the "IM flagged"
        grFracIM->Draw("P SAME");

        // Legend
        TLegend leg(0.55, 0.75, 0.88, 0.88);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.AddEntry(grFracTruth, "Truth #pi^{0}/#eta fraction", "pl");
        leg.AddEntry(grFracIM,    "Mass-window fraction (IM flagged)", "pl");
        leg.Draw();

        cFrac.cd(); // return to main canvas

        // Bottom pad (30% height) for ratio
        TPad* padBot = new TPad("padBot","padBot", 0.0, 0.0, 1.0, 0.3);
        padBot->SetTopMargin(0.02);
        padBot->SetBottomMargin(0.25);
        padBot->Draw();
        padBot->cd();

        TGraphErrors* grRatio = new TGraphErrors(nBins,
                                                 &vBinCenterFrac[0], &vRatio[0],
                                                 nullptr, &vErrRatio[0]);
        grRatio->SetMarkerStyle(20);
        grRatio->SetMarkerColor(kBlack);
        grRatio->SetLineColor(kBlack);
        grRatio->SetTitle(""); // We’ll set axis labels directly
        grRatio->Draw("AP");

        grRatio->GetHistogram()->GetXaxis()->SetRangeUser(pTedges.front()-0.5,
                                                          pTedges.back() + 0.5);
        grRatio->GetHistogram()->GetXaxis()->SetTitle("p_{T} (GeV)");
        grRatio->GetHistogram()->GetXaxis()->SetTitleSize(0.1);
        grRatio->GetHistogram()->GetXaxis()->SetLabelSize(0.08);
        grRatio->GetHistogram()->GetYaxis()->SetTitle("IM / Truth");
        grRatio->GetHistogram()->GetYaxis()->SetTitleSize(0.1);
        grRatio->GetHistogram()->GetYaxis()->SetLabelSize(0.08);

        // Set ratio y-range (optional)
        grRatio->GetHistogram()->GetYaxis()->SetRangeUser(0., 1.2);

        // optional horizontal line at ratio=1
        TLine* line1 = new TLine(pTedges.front()-0.5, 1.0, pTedges.back()+0.5, 1.0);
        line1->SetLineColor(kGray+2);
        line1->SetLineStyle(2);
        line1->Draw("SAME");

        // Finally save
        cFrac.SaveAs((outDir + "/FractionPi0EtaVsPt_overlayIM.png").c_str());
        std::cout << "[INFO] Saved overlap fraction + ratio plot: "
                  << outDir + "/FractionPi0EtaVsPt_overlayIM.png" << std::endl;
    }

    // (2) original vs enhanced purity
    double overallMaxPur = maxPur;
    TCanvas cPur("cPur","PromptPhotonPurity", 1200,600);
    cPur.Divide(2,1);

    // left pad => original purity
    cPur.cd(1);
    {
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);

        TGraphErrors* grPurOrg = new TGraphErrors(nBins,
                                                  &vBinCenterPur[0], &vPurOriginal[0],
                                                  nullptr, &vErrPurOriginal[0]);
        grPurOrg->SetTitle("Original Prompt Photon Purity; p_{T} (GeV); Purity");
        grPurOrg->SetMarkerStyle(20);
        grPurOrg->SetMarkerColor(kRed+1);
        grPurOrg->SetLineColor(kRed+1);

        bool useLogY = (overallMaxPur<0.02);
        if(useLogY){
            gPad->SetLogy();
            double minY=1e-7;
            if(overallMaxPur>1e-10) minY=1e-4*overallMaxPur;
            grPurOrg->GetYaxis()->SetRangeUser(minY,1.5*overallMaxPur);
        } else {
            grPurOrg->GetYaxis()->SetRangeUser(0., 1.2*overallMaxPur);
        }
        grPurOrg->GetXaxis()->SetRangeUser(pTedges.front()-0.5,
                                           pTedges.back()+0.5);
        grPurOrg->Draw("AP");
    }

    // right pad => enhanced purity
    cPur.cd(2);
    {
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);

        TGraphErrors* grPurEn = new TGraphErrors(nBins,
                                                 &vBinCenterPur2[0], &vPurEnhanced[0],
                                                 nullptr, &vErrPurEnhanced[0]);
        grPurEn->SetTitle("Enhanced Prompt Photon Purity (Excl. Meson); p_{T} (GeV); Purity");
        grPurEn->SetMarkerStyle(21);
        grPurEn->SetMarkerColor(kGreen+2);
        grPurEn->SetLineColor(kGreen+2);

        bool useLogY = (overallMaxPur<0.02);
        if(useLogY){
            gPad->SetLogy();
            double minY=1e-7;
            if(overallMaxPur>1e-10) minY=1e-4*overallMaxPur;
            grPurEn->GetYaxis()->SetRangeUser(minY, 1.5*overallMaxPur);
        } else {
            grPurEn->GetYaxis()->SetRangeUser(0., 1.2*overallMaxPur);
        }
        grPurEn->GetXaxis()->SetRangeUser(pTedges.front()-0.5,
                                          pTedges.back()+0.5);
        grPurEn->Draw("AP");
    }

    std::string outPur = outDir + "/PromptPhotonPurity_Compare.png";
    cPur.SaveAs(outPur.c_str());
    std::cout << "[INFO] Saved purity comparison plot: " << outPur << std::endl;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Summary prints
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    std::cout << "\n============================================\n";
    std::cout << "[INFO] Done with two-pass approach.\n";
    std::cout << "We have:\n";
    std::cout << "  * Fit di-cluster invariant mass in each pT bin => Wrote PNG + CSV\n";
    std::cout << "  * Overlaid fraction(#pi^0/#eta) vs fraction(mass-window) => saved a PNG\n";
    std::cout << "  * Computed original prompt-photon purity & enhanced purity => saved side-by-side PNG\n\n";

    std::cout << "Second Pass => totalEventsProcessed="<< totalEventsProcessed
              <<", totalClustersProcessed="<< totalClustersProcessed<<"\n";

    // Print fraction from truth
    std::cout << "\nBreakdown by pT bin: Truth fraction(#pi^0/#eta) vs IM-Flag fraction\n";
    for(int i=0; i<nBins; i++){
        double low=pTedges[i], hi=pTedges[i+1];
        long long tot= nAllClustersInBin[i];
        long long mes= nPi0EtaInBin[i];
        long long flagged = nFlaggedAsMesonInBin[i];
        double fracTruth= (tot>0? double(mes)/tot : 0.);
        double fracIM= (tot>0? double(flagged)/tot : 0.);
        std::cout<<"  pT["<<low<<","<<hi<<") => tot="<<tot
                 <<", truthPi0Eta="<<mes<<" => fraction="<<fracTruth
                 <<", IMflagged="<<flagged<<" => fraction="<<fracIM<<"\n";
    }

    // Print original purity
    std::cout<<"\nBreakdown by pT bin for Original PromptPhotonPurity:\n";
    for(int i=0;i<nBins;i++){
        double low=pTedges[i], hi=pTedges[i+1];
        long long c1=nPromptTagCorrect[i];
        long long f1=nPromptTagFake[i];
        long long denom=c1+f1;
        double pur=(denom>0? double(c1)/denom : 0.);
        std::cout<<"  pT["<<low<<","<<hi<<") => correct="<<c1
                 <<", fake="<<f1<<", purity="<<pur<<"\n";
    }

    // Print enhanced purity
    std::cout<<"\nBreakdown by pT bin for Enhanced PromptPhotonPurity:\n";
    for(int i=0;i<nBins;i++){
        double low=pTedges[i], hi=pTedges[i+1];
        long long c2=nPromptTagCorrectEnhanced[i];
        long long f2=nPromptTagFakeEnhanced[i];
        long long d2=c2+f2;
        double pur=(d2>0? double(c2)/d2 : 0.);
        std::cout<<"  pT["<<low<<","<<hi<<") => correct="<<c2
                 <<", fake="<<f2<<", purity="<<pur<<"\n";
    }

    std::cout<<"============================================\n\n";
    std::cout << "[INFO] All tasks completed successfully.\n";
}
