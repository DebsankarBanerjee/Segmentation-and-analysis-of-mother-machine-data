# SAM: High-Throughput Automated Image Segmentation & Tracking

**Lead Developer: Deb Sankar Banerjee, PhD** *A modular MATLAB pipeline for high-precision metrology and feature extraction in time-lapse microscopy.*

## 🔬 Engineering Overview

SAM (Segmentation and Analysis of Mother-machine data) is a high-performance computational pipeline designed to extract precise morphological data from noisy, high-volume time-lapse imaging. Originally developed for microfluidic bacterial studies, the core architecture is built to handle:

* **Low Signal-to-Noise Ratios (SNR):** Robust to intensity fluctuations and optical aberrations.
* **High-Throughput Processing:** Optimized to process **~1 GB of image data per 15 minutes**.
* **Automated Metrology:** Adaptive segmentation algorithms to eliminate manual correction and bias.

## 🛠 Key Technical Features

* **Adaptive Segmentation:** Implements dynamic thresholding to identify micro-scale features (cells/objects) amidst fluctuating backgrounds and intensity gradients.
* **Noise & Artifact Mitigation:** Specialized logic to handle overlapping features, filamentation, and intrusion artifacts—critical for maintaining data integrity in stochastic environments.
* **Automated Tracking:** High-fidelity temporal tracking of object kinetics, measuring growth rates and spatial distribution over thousands of frames.
* **Scalable Architecture:** Designed for easy integration with Shell/Bash scripting for parallel processing on High-Performance Computing (HPC) clusters.

## 📊 Pipeline Workflow

1. **Pre-processing:** Background subtraction and intensity normalization to improve feature contrast.
2. **Segmentation:** Adaptive morphological operations to isolate objects of interest.
3. **Refinement:** Automated detection and treatment of "unwanted aberrations" (e.g., cell overlapping or intrusion).
4. **Analysis:** Quantitative extraction of size, growth rate, and temporal dynamics.
5. **Visualization:** Systematic generation of output file structures and graphical diagnostic reports.

## 🚀 Performance

* **Processing Speed:** Optimized for large-scale datasets, achieving a throughput of ~15 min/GB.
* **Accuracy:** Engineered to minimize faulty detections during high-frequency division/splitting events.

## 💻 Getting Started

### Prerequisites

* MATLAB (Tested on R2020a and later)
* Image Processing Toolbox

### Installation

```bash
git clone https://github.com/DebsankarBanerjee/Segmentation-and-analysis-of-mother-machine-data.git
cd Segmentation-and-analysis-of-mother-machine-data

```

### Usage

1. Place your raw `.tif` or image sequence in the `/data` directory.
2. Run the main processing script:

```matlab
run('SAM_Main.m') % [Ensure this matches your actual entry point script]

```
### Workflow and segmentation steps

![Descriptive Title](workflow_plot.png)

## 📖 Citation

If you use this tool in your research, please cite:

> *Banerjee, D. S., Stephenson, G., & Das, S. M. (2020). Segmentation and analysis of mother machine data: SAM. bioRxiv.* [[Link to Paper](https://doi.org/10.1101/2020.10.01.322685)]
