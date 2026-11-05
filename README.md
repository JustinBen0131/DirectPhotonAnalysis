# Identifying neutral-meson backgrounds to isolated photons

A neutral pion or eta meson can decay into two photons. If one daughter passes an isolated-photon selection, it can be mistaken for a prompt photon. Justin Bennett independently designed and implemented the two-pass background-tagging study in this repository for sPHENIX proton–proton data. The idea is to look for a second calorimeter cluster and use the pair’s invariant mass to identify candidates consistent with a neutral-meson decay. Calibrating those mass windows from data lets them follow changes in the peak position and width with momentum and selection.

This repository documents the implemented tagging study and accompanying trigger work. The nominal PPG12 isolated-photon measurement used a boosted decision tree for identification and a data-driven ABCD purity estimate. The collaboration’s public Hard Probes presentation describes that measurement.

## The two passes

1. Fit the π⁰ and η mass peaks in binned two-cluster spectra with two Gaussian peaks and a polynomial background. Record the fitted means and widths. [Mass-fit implementation](https://github.com/JustinBen0131/DirectPhotonAnalysis/blob/48a918658eef750ac03a1a4ad7609d89da2bec3c/macro/AnalyzeTriggerGroupings.cpp#L1180-L1221)
2. Look up the appropriate windows using the trigger, cuts and momentum bin, then test whether the pair lies within three fitted standard deviations of either peak. The code specifies fallback windows when a usable fit is unavailable. [Pair classification](https://github.com/JustinBen0131/DirectPhotonAnalysis/blob/48a918658eef750ac03a1a4ad7609d89da2bec3c/src/caloTreeGen.cc#L2655-L2698)

## Related work

- [Bennett’s Run‑24 p+p selection script in the collaboration repository](https://github.com/sPHENIX-Collaboration/analysis/blob/15c6049534ca3e19016e366f753717905299d56a/runListGeneration/Run2024pp/GoldenRunList_ConductorFile.sh) opens with the run-quality criteria and options.
- [PPG01 neutral-pion flow and interactive cut studies](docs/neutral-pion-flow/README.md) describes a separate, earlier measurement and points to the original collaboration code.
- [Photon-trigger turn-ons](https://github.com/JustinBen0131/DirectPhotonAnalysis/blob/48a918658eef750ac03a1a4ad7609d89da2bec3c/macro/AnalyzeTriggerGroupings.cpp#L8234-L8299) · [Public Hard Probes presentation of the nominal measurement](https://indico.cern.ch/event/1428985/contributions/6998769/)
