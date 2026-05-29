# newtonian-vs-gr-waveforms
Comparing Newtonian and GR peak merger frequencies for BH-NS binaries using SPH data and PyCBC SEOBNRv4 - FFT analysis
Newtonian vs General Relativity: Gravitational Waveform Comparison
A computational study comparing Newtonian and General Relativistic frameworks for modelling gravitational wave emission from black hole–neutron star (BH–NS) binary mergers.
Developed as a capstone project for the BRICS Astronomy & IDIA Data Analytics Training Programme.

Overview
Newtonian gravity works well for slowly moving, weakly gravitating systems — but breaks down near compact object mergers. This project quantifies that breakdown by comparing peak gravitational wave frequencies predicted by:

Newtonian framework — SPH simulation data from Lee & Kluźniak (1998–2000), loaded and analysed via FFT
General Relativistic framework — waveforms generated using PyCBC's SEOBNRv4 approximant

The key finding: GR predicts significantly higher peak merger frequencies than Newtonian models, confirming that relativistic effects are essential for accurate gravitational wave modelling near merger.

Scientific Background
The Newtonian simulation data comes from a three-part study by Lee & Kluźniak on BH–NS coalescence using Smoothed Particle Hydrodynamics (SPH):
PaperConfigurationEquation of StateLee & Kluźniak (1998)Tidally lockedStiff (Γ = 3)Lee & Kluźniak (1999)Tidally lockedSoft (Γ = 5/3)Lee (2000)IrrotationalStiff & Soft
Five waveform files are analysed, varying the adiabatic index Γ and spin configuration (tidally locked vs irrotational). The black hole is treated as a Newtonian point mass with an absorbing boundary at the Schwarzschild radius.
For the GR comparison, PyCBC's SEOBNRv4 approximant models the same BH–NS system (M_BH = 4.516 M☉, M_NS = 1.4 M☉) using full relativistic waveform templates.

Methodology

Load Newtonian waveforms — five .out files containing time, h+ and h× polarisations
Apply FFT — Fast Fourier Transform with Hanning window to extract frequency content
Identify peak frequencies — dominant GW frequency at merger for each configuration
Generate GR waveform — PyCBC SEOBNRv4 for same physical parameters
Compare peak frequencies — Newtonian vs GR across all five configurations


Key Result

GR peak frequencies are significantly higher than Newtonian predictions, reflecting the strong-field relativistic dynamics near merger that Newtonian gravity cannot capture. Newtonian models are adequate for early inspiral but fail at the merger itself.

The variation in Newtonian peak frequencies across different equations of state and spin configurations is small — confirming that Newtonian dynamics does not strongly encode neutron star internal structure or spin-orbit interactions the way GR does.

Tech Stack

Python
PyCBC (waveform generation — SEOBNRv4 approximant)
GWpy, LALSuite
NumPy, SciPy (FFT analysis)
Matplotlib, Seaborn


File Structure
├── waveform_comparison.py       # Main analysis script
├── data/
│   ├── waveform1BHNS.out        # Newtonian SPH data: Γ=3.0, Irrotational
│   ├── waveform2BHNS.out        # Newtonian SPH data: Γ=2.0, Irrotational
│   ├── waveform3BHNS.out        # Newtonian SPH data: Γ=5/3, Irrotational
│   ├── waveform4BHNS.out        # Newtonian SPH data: Γ=3.0, Tidally Locked
│   └── waveform5BHNS.out        # Newtonian SPH data: Γ=5/3, Tidally Locked
└── README.md

References

Lee, W. H. & Kluźniak, W. (1998). Newtonian Hydrodynamics of BH–NS Coalescence I. arXiv:astro-ph/9808185
Lee, W. H. & Kluźniak, W. (1999). Newtonian Hydrodynamics of BH–NS Coalescence II. arXiv:astro-ph/9904328
Lee, W. H. (2000). Newtonian Hydrodynamics of BH–NS Coalescence III. arXiv:astro-ph/0007206


Author
Shivika Lamba — MSc Astrophysics, Cardiff University
Research interests: gravitational-wave data analysis, waveform modelling, EMRI science, LISA
