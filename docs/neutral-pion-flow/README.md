# Neutral-pion elliptic flow and interactive cut studies

The sPHENIX collaboration measured how neutral pions were distributed relative to the collision geometry in Run‑23 gold–gold data. This distribution is summarized by the second harmonic, v₂. The pions are reconstructed from photon pairs, so choices about cluster energy, pair asymmetry, shower quality and separation can change both the mass peak and the extracted flow.

Justin Bennett and Apurva Narde co-led the analysis. Bennett built downstream C++/ROOT analysis tools and an interactive Python interface for comparing cut combinations, fit quality and uncertainties. The GUI lets an analyst overlay selections and inspect the mass distributions behind the points. The scientific result is the collaboration’s public 2024 measurement. This page collects links to the contemporary implementation in the collaboration repository.

## Follow one selection through the analysis

1. Enter energy, asymmetry, χ² and angular-separation cuts. The interface constructs their Cartesian product and overlays each combination. [GUI code](https://github.com/sPHENIX-Collaboration/analysis/blob/9fea6c50312fd3cb25ebfe697db24db015265a0c/pi0-vn/AnalyzePi0s_Final/src/InteractiveCutAnalysis_twoCentralities.py#L1381-L1420)
2. Inspect an individual combination and open its underlying invariant-mass histogram. [Histogram controls](https://github.com/sPHENIX-Collaboration/analysis/blob/9fea6c50312fd3cb25ebfe697db24db015265a0c/pi0-vn/AnalyzePi0s_Final/src/InteractiveCutAnalysis_twoCentralities.py#L922-L978)
3. Follow the final C++ calculation that combines flow estimates and propagates uncertainties. [Flow finalization](https://github.com/sPHENIX-Collaboration/analysis/blob/5e305a80a7eb07857a943f0db545fb80da9e1d94/pi0-vn/AnalyzePi0s_Final/src/src_final/FinalizeResults_IntegratedPT.cpp#L811-L825)

[Bennett’s December 2023 GUI commit](https://github.com/sPHENIX-Collaboration/analysis/commit/7aadcb266c2565f5684dc979747febb5bc0bbf3c) · [Public sPHENIX conference note](https://www.sphenix.bnl.gov/sPH-CONF-BULK-2024-01)
