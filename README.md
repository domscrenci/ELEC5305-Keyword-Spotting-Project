# ELEC5305 — Keyword Spotting (MFCC + LSTM)

This project improves MATLAB’s baseline keyword-spotting model by fixing false detections caused by sibilant “sss” sounds.  
The model uses MFCC + Δ + ΔΔ features with an LSTM network (`KWSBaseline1.mat`).  
The update focuses on adaptive thresholding, smoothing, and hysteresis for more accurate “yes” detection without retraining.

---

## Files

**KWS_For_Wav.m** — Main evaluation script.  
Loads a `.wav` file, extracts MFCCs, runs inference through the LSTM model, applies adaptive thresholds, and plots waveform, posterior, and detected events. this is the main file for testing single wavs

**KWS_Batch_Analysis.m** — Runs `KWS_Eval_final` across all test clips to compare multiple experiments or parameter sets. call the function with the name or location of the folder with all the test clips

**KWS_Batch_Analysis_CSV_Compiler.m** — Collects CSV outputs from batch runs and combines them into summary tables or heatmaps. run by calling with name of location of the folder ensure that the summary csv is also in the folder.

**KWS_Network_Checker.m** — Confirms that the LSTM model loads correctly and displays network structure and I/O sizes.

**KWSBaseline.mat / KWSBaseline1.mat** — Pre-trained binary LSTM keyword-spotting models. USE KWSBaseline1.mat when running the code

**KWS.mlx** — MATLAB Live Script version of the workflow; shows main pipeline interactively. This is a condesned version of the orignal matlab file and it will succsessfully run the orignal matlab without errors of the orignal.

**KeywordSpottingInNoiseUsingMFCCAndLSTMNetworksExample.mlx** — Original MathWorks example used as the baseline reference.

**Record_Yes.m** — Simple recorder to create new “yes” audio samples for testing.

**audio/** — Folder containing test recordings such as `yes__01.wav`, `yes__02.wav`, and `keywordTestSignal.wav`.

**README.md** — This file.

---

## How to Run

1. Open MATLAB and navigate to this project folder.  
2. Make sure the test clips are inside the `audio/` folder.  
3. Run the main script:    KWS_For_Wav for single wav evaluation or KWS_Batch_Analysis for running testing all wavs.
