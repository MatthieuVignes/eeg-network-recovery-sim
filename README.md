# eeg-network-recovery-sim

Simulation-based feasibility check for the directed-connectivity (network)
analysis in the "Auditory Verbal Hallucinations in Schizophrenia: Testing the 
Shared Resource and Source Monitoring Hypotheses" project: can a top-down vs. 
bottom-up shift in PAC → STG → IFG connectivity, locked to auditory verbal 
hallucination (AVH) events, be recovered after reducing a 64-channel EEG montage 
down to a small set of latent components?

## What this does

1. Simulates three latent Region Of Interest (ROI) sources: **P**rimary 
   **A**uditory **C**ortex, **S**uperior **T**emporal **G**yrus, **I**nferior 
   **F**rontal **G**yrus. They are coupled via a lag-delayed directional process, 
   under two competing hypotheses:
   - `bottom_up`: PAC → STG → IFG dominant (feedforward)
   - `top_down`: IFG → STG → PAC dominant (feedback)
2. Projects these sources through a simulated 64-channel setup (random
   leadfield-like mixing + sensor noise), reduces to 10 components via PCA,
   and reconstructs 3 ROI estimates via a linear MMSE ("Wiener") filter.
3. Estimates directed influence on the reconstructed signals via a lagged
   cross-correlation asymmetry index (efficient proxy for Granger-style
   directionality).
4. Computes:
   - **Power** to detect an AVH-locked top-down shift (conditions 4/6:
     inner-speech-cued, natural speech), across a grid of participant N
     (20/30/40/50) and AVH event rate (10/20/35% of 60 trials)
   - **False-positive rate** in control conditions (motor, visual, silence),
     where no shift is expected (calibration/sanity check)
   - **Cross-talk (leakage) matrix** quantifying how much the three ROI
     estimates contaminate one another after the 64→10-PC reduction.

## Requirements

R (≥ 4.0) with `ggplot2` and `scales`.

```r
install.packages(c("ggplot2", "scales"))
```

## Running

```bash
Rscript roi_connectivity_simulation.R
```

Runtime: ~2–3 minutes (1,000 simulations per grid cell). Outputs are written
to the working directory:

| File | Contents |
|---|---|
| `roi_connectivity_power_plot.png` | Power curves by N × AVH event rate |
| `roi_connectivity_power.csv` | Power/SE for the main AVH-locked contrast |
| `roi_connectivity_falsepositive.csv` | False-positive rate in control conditions |
| `roi_leakage_matrix.csv` | 3×3 ROI cross-talk/gain matrix |

## Result summary

Power to detect the AVH-locked directionality shift tops out just over 70% at N=50 
with the most favourable AVH event rate. This reflects a real but doable SNR cost 
of going through 64 noisy channels and a generic (non-anatomical) dimensionality
reduction. The false-positive rate in control conditions sits at a well-calibrated 5%.

*Practical note:* we will not plan the network analysis around raw 64-channel/PCA 
reduction. We will use state-of-the-art source reconstruction (beamforming or eLORETA 
with a template/individual head model) to a small number of anatomically defined ROIs, 
and treat Dynamic-Causal-Modelling (DCM)-for-Event-Related Potential (ERP) on that 
reduced ROI set as the primary confirmatory network analysis, with any full-scalp 
connectivity mapping as exploratory.

## Some warnings

- This is a **linear approximation**, not a DCM-for-ERP emulation. Our planned 
  analysis (nonlinear neural-mass models, Bayesian model inversion via SPM) is not 
  developed yet. We present a simplified proof-of-concept validation of the reduction 
  *pipeline*, not a stand-in for the final analysis.
- Coupling-strength and noise parameters are **assumed and tuned**, not derived 
  from pilot data or literature (no established benchmark exists for "directed 
  coupling strength between PAC/STG/IFG from scalp EEG", to the best of our knowledge). 
  They were chosen to produce an informative, non-saturated power curve. The 
  qualitative conclusion (channel-level connectivity is a hard but doable inferential 
  problem) not the specific percentages.
- Simulations to be re-run with updated parameters once pilot data give empirical 
  estimates of single-trial noise and realistic AVH event rates.

