#ifndef ANALYZEJET10SLIMTREE_H
#define ANALYZEJET10SLIMTREE_H

#include <cmath>      // std::abs, std::sqrt
#include <vector>     // std::vector  <-- needed for pTedges
#include <TMath.h>    // TMath::Pi


// ---------------------------------------------------------------------
// Extern declarations for your truth-level arrays/variables.
// ---------------------------------------------------------------------
extern int   nparticles;
extern int   particle_pid[];               // PDG PID of each truth particle
extern float particle_E[];                 // truth energy
extern float particle_Eta_[];              // truth eta
extern float particle_Phi_[];              // truth phi
extern int   particle_photonclass[];       // e.g. 1=prompt,2=frag,3=decay
extern int   particle_photon_mother_pid[]; // mother PID if photon

// ---------------------------------------------------------------------
//  For custom pT bin edges. We only declare them extern here,
//  so we can define them in the .cpp or macro.
// ---------------------------------------------------------------------
extern std::vector<float> pTedges;

/**
 * Find which bin index 'pt' belongs to in pTedges.
 * If pt is outside the final edge, return -1.
 */
inline int findPtBin(float pt)
{
    // We loop over pTedges[i], pTedges[i+1].
    // If pT is in [pTedges[i], pTedges[i+1]), we return i.
    for(int i = 0; i < static_cast<int>(pTedges.size()) - 1; i++) {
        if(pt >= pTedges[i] && pt < pTedges[i+1]) {
            return i;
        }
    }
    return -1; // Not found => outside last edge
}

bool passIsoCombo(float isoAll, float isoEM, int iC)
{
    // For example:
    //  - iC=0 => isoAll < 2.0 && isoEM < 1.0
    //  - iC=1 => isoAll < 3.0 && isoEM < 2.0
    //  etc.
    // Implement as needed. Here's a trivial example:

    if(iC==0) {
        return (isoAll < 2.0 && isoEM < 1.0);
    }
    else if(iC==1) {
        return (isoAll < 4.0 && isoEM < 2.0);
    }
    else if(iC==2) {
        return (isoAll < 6.0 && isoEM < 3.0);
    }
    return false;
}

// ---------------------------------------------------------------------
// The matching logic in a namespace
// ---------------------------------------------------------------------
namespace SlimTreeMatch {

  // 1) A strongly-typed enum for the "top-level" classification
  enum class TruthMatchCategory
  {
    kPromptPhoton,     ///< truly prompt photon (photonclass=1, mother not pi0/eta)
    kFragPhoton,       ///< fragmentation photon (photonclass=2)
    kPi0orEtaDecay,    ///< pi^0 or eta (±111/±221) OR photon (photonclass=3) with mother=±111/±221
    kOther,            ///< hadron/lepton/unknown
    kNoMatch           ///< for "cluster→truth" usage if no match found; can also use for invalid
  };

  // 2) toString => for printing the above category
  inline const char* toString(TruthMatchCategory cat)
  {
    switch(cat) {
      case TruthMatchCategory::kPromptPhoton:  return "PromptPhoton";
      case TruthMatchCategory::kFragPhoton:    return "FragPhoton";
      case TruthMatchCategory::kPi0orEtaDecay: return "Pi0orEtaDecay";
      case TruthMatchCategory::kOther:         return "Other";
      case TruthMatchCategory::kNoMatch:       return "NoMatch";
    }
    return "Unknown";
  }

  // -------------------------------------------------------------------
  // 3) A sub-categorization function for "other" hadrons
  // -------------------------------------------------------------------
  inline int mapOtherPidToSubcat(int pidMatch)
  {
    if(pidMatch < 0 || pidMatch == 999999) {
      return 5; // "No-match" or invalid
    }
    if(pidMatch == 211) {
      return 6;  // "pion"
    }
    if(pidMatch == 321 || pidMatch == 311 ||
       pidMatch == 310 || pidMatch == 130) {
      return 7;  // "kaon"
    }
    if(pidMatch == 2212 || pidMatch == 2112) {
      return 8;  // "nucleon"
    }
    return 9; // catch-all
  }

  // -------------------------------------------------------------------
  // 4) deltaR function
  // -------------------------------------------------------------------
  inline float deltaR(float eta1, float phi1, float eta2, float phi2)
  {
    float dEta = eta1 - eta2;
    float dPhi = phi1 - phi2;
    while(dPhi >  TMath::Pi()) dPhi -= 2.f*TMath::Pi();
    while(dPhi < -TMath::Pi()) dPhi += 2.f*TMath::Pi();
    return std::sqrt(dEta*dEta + dPhi*dPhi);
  }

  // -------------------------------------------------------------------
  // 5) findBestTruthMatch (CLUSTER → TRUTH):
  //    - loop over 0..nparticles-1
  //    - require deltaR < 0.05, E_reco / E_truth > 0.5
  //    - choose the truth particle with the smallest deltaR
  //    - return its index, or -1 if none pass
  // -------------------------------------------------------------------
  inline int findBestTruthMatch(float cEta, float cPhi, float cE)
  {
    const float maxDR = 0.05;
    float bestDR      = 9999.f;
    int   bestIdx     = -1;

    for(int ip = 0; ip < nparticles; ip++){
      float etaTruth = particle_Eta_[ip];
      float phiTruth = particle_Phi_[ip];
      float eTruth   = particle_E[ip];

      // skip negligible energy
      if(eTruth < 1e-3) {
        continue;
      }

      // check deltaR
      float dR = deltaR(cEta, cPhi, etaTruth, phiTruth);
      if(dR > maxDR) {
        continue;
      }

      // check E ratio
      float ratio = cE / eTruth;
      if(ratio < 0.5f) {
        continue;
      }

      // update best
      if(dR < bestDR){
        bestDR  = dR;
        bestIdx = ip;
      }
    }
    return bestIdx;
  }


  // -------------------------------------------------------------------
  // 6) getTruthMatchCategoryForCluster (CLUSTER → TRUTH):
  //    - calls findBestTruthMatch
  //    - if none => kNoMatch
  //    - if found => classify by pid, photonclass, mother
  // -------------------------------------------------------------------
    inline TruthMatchCategory getTruthMatchCategoryForCluster(
        float cEta, float cPhi, float cE,
        /*OUT*/ int* pBestMatchIdx = nullptr)
    {
        int bestIdx = findBestTruthMatch(cEta, cPhi, cE);
        if(pBestMatchIdx) {
            *pBestMatchIdx = bestIdx;
        }
        // If no match => kNoMatch
        if(bestIdx < 0) {
            return TruthMatchCategory::kNoMatch;
        }

        int pid     = particle_pid[bestIdx];
        int abspid  = std::abs(pid);
        int momPid  = std::abs(particle_photon_mother_pid[bestIdx]);
        int phClass = particle_photonclass[bestIdx];

        // (A) If the best-matched truth particle is pi0 or eta => "Pi0orEtaDecay"
        if(abspid == 111 || abspid == 221) {
            return TruthMatchCategory::kPi0orEtaDecay;
        }

        // (B) If it's a photon => check mother + photonclass
        if(abspid == 22)
        {
            // If photonclass=3 ("decay"), we must see if mother=pi0 or eta
            if(phClass == 3)
            {
                if(momPid == 111 || momPid == 221) {
                    return TruthMatchCategory::kPi0orEtaDecay; // same as direct pi0/eta
                }
                return TruthMatchCategory::kOther; // decay from some other hadron
            }

            // If mother=±111 or ±221 => "kPi0orEtaDecay" too
            if(momPid == 111 || momPid == 221) {
                return TruthMatchCategory::kPi0orEtaDecay;
            }

            // Then check photonclass for 1=prompt, 2=frag, else => other
            if(phClass == 1) {
                return TruthMatchCategory::kPromptPhoton;
            }
            if(phClass == 2) {
                return TruthMatchCategory::kFragPhoton;
            }

            // Otherwise => "other"
            return TruthMatchCategory::kOther;
        }

        // (C) Else => hadron, lepton, etc. => "other"
        return TruthMatchCategory::kOther;
    }


   inline TruthMatchCategory getTruthMatchCategoryForParticle(int iPar)
    {
        if(iPar < 0 || iPar >= nparticles) {
            return TruthMatchCategory::kNoMatch; // invalid index
        }

        int pid     = particle_pid[iPar];
        int abspid  = std::abs(pid);
        int momPid  = std::abs(particle_photon_mother_pid[iPar]);
        int phClass = particle_photonclass[iPar];

        // If the particle itself is pi^0 or eta => "kPi0orEtaDecay"
        if(abspid == 111 || abspid == 221) {
            return TruthMatchCategory::kPi0orEtaDecay;
        }

        // If it's a photon => check mother + photonclass
        if(abspid == 22)
        {
            // Mother = pi^0 or eta => "kPi0orEtaDecay"
            if(momPid == 111 || momPid == 221) {
                return TruthMatchCategory::kPi0orEtaDecay;
            }

            // Then classify by photonclass
            if(phClass == 1) return TruthMatchCategory::kPromptPhoton; // prompt
            if(phClass == 2) return TruthMatchCategory::kFragPhoton;   // frag
            if(phClass == 3) {
                // photonclass=3 => "decay."
                // Already checked mother=111/221 above, so if we reach here,
                // mother is NOT pi0/eta => => "kOther"
                return TruthMatchCategory::kOther;
            }
            // If photonclass is something else => "kOther"
            return TruthMatchCategory::kOther;
        }

        // Else => hadron, lepton, etc. => "kOther"
        return TruthMatchCategory::kOther;
    }
} // end namespace SlimTreeMatch

#endif // ANALYZEJET10SLIMTREE_H
