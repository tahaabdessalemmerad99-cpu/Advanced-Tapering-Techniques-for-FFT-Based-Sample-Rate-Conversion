# Advanced-Tapering-Techniques-for-FFT-Based-Sample-Rate-Conversion
MATLAB implementation of FFT-based sample rate conversion and DDC taper analysis.

## MATLAB Scripts

### Code01 – Taper Generation and Impulse Response Analysis

Generates and evaluates the tapering windows used in the FFT-SRC framework (No Taper, Cosine, Hann, Blackman, Chebyshev, and DDC). The script compares their impulse responses and validates the taper implementations before using them in the SRC experiments.

---

### Code02 – Real Audio SRC Evaluation

Performs FFT-based sample rate conversion on real audio signals and compares the spectrograms obtained with different tapering methods.

---

### Code03 – Reference SRC Benchmark

Implements the reference FFT-SRC benchmark using a synthetic test signal and compares the performance of all tapering methods.

---

### Code04 – Low-Frequency and Between-Bin Analysis

Investigates the influence of DDC tapering on low-frequency spectral components and performs high-resolution FFT analysis to study spectral differences between FFT bins.

## MATLAB Scripts

|  Code  |              File                 |                   Description                  |
|--------|-----------------------------------|------------------------------------------------|
| Code01 | compare_tapers_impulse_response.m | Taper generation and impulse response analysis |
| Code02 | Python_like_Figure_AUDIO_ALL_Heatmap.m |    Real audio SRC evaluation              |
| Code03 | Python_like_Figure_ALL_WINDOWS_Heatmap.m |    Reference benchmark reproduction     |
| Code04 |            between_bins.m         |      Low-frequency and between-bin analysis    |
 


