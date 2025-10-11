# ELEC5305 — Live Keyword Spotting (MFCC + LSTM)

This repository documents a small but focused research project on live, single-word keyword spotting (KWS) using MFCC features (13 static coefficients plus first- and second-order deltas for a 39-dimensional stream) and an LSTM sequence model. 
The practical goal is to make the stock MathWorks live demo stable, predictable, and demonstrably more robust in real time, with a special emphasis on reducing sibilant-driven false positives (for example, the “ss” sound in “yes”). 
Rather than proposing a new architecture, the project concentrates on reproducible engineering improvements: cleaning up the live pipeline, tightening buffering and feature extraction, adding decision-time post-processing (smoothing, hysteresis, duration checks), and introducing simple signal-based vetoes that distinguish hiss from genuine speech.

The intended outcome is a live demo that behaves reliably on a microphone at 16 kHz, a compact set of evaluation routines that show before/after behaviour, and a brief written discussion that explains the design choices, the trade-offs, and the limitations.
All code is organised to be easy to read and easy to modify during a live demonstration or marking session.


├─ code

│  ├─ KeywordSpottingInNoiseUsingMFCCAndLSTMNetworksExample.mlx   ← ORIGNAL main live script (without fixes)

|  |─KeywordSpottingInNoiseUsingMFCCAndLSTMNetworksExampleUPDATED.mlx  ← UPDATED main live script (with fixes)

│  └─ generateKeywordFeatures.m                                   ← generated 16 kHz, 512/384 MFCC(13)+Δ+ΔΔ

├─ models

│  └─ KWSBaseline.mat                                            ← pre-trained binary model 

│  └─ KWSBaseline1.mat                                            ← pre-trained binary model 

└─ audio

   └─ keywordTestSignal.wav                                       ← small test clip 
   
