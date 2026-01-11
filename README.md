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

---

## ⚙️ Theory of Operation: The SAM Algorithm

The SAM pipeline transforms raw optical sensor data into precise physical metrology through a three-stage computational framework.

### 1. Pre-Processing & Data Decomposition

High-resolution images of the microfluidic device are decomposed into localized image stacks. Each stack represents a single sensing channel, isolating the region of interest (ROI) to reduce computational overhead and focus the segmentation logic on individual data streams.

---

### 2. The Core Processing Engine (SAM)

The engine utilizes a hybrid **1D Signal Projection + 2D Morphological Analysis** approach to detect and track features with high temporal resolution.

#### **A. Signal Enhancement & Dimensionality Reduction**

* **Intensity Projection:** Each frame undergoes a mid-line spatial projection. By averaging pixel intensities across the channel width, the 2D image is reduced to a 1D longitudinal intensity signal.
* **Feature Localization:** Local minima detection is performed on the 1D signal to identify potential boundaries (cell poles/division planes).

#### **B. Multi-Pass Segmentation**

* **Initial Masking:** A binary mask is generated using a user-defined global threshold (). This mask is then partitioned into **Connected Components (CC)** based on the previously detected intensity minima.
* **Adaptive Error Correction:** The system runs a series of consistency cross-checks to identify optical aberrations (e.g., overlapping features, signal drift, or "cell intrusion"). If inconsistencies are detected, specialized adaptive operations are triggered to refine the boundaries.
* **Final Metrology:** A secondary, refined segmentation (CC2) is performed. The system applies **elliptical fitting** to each feature to extract high-precision measurements of the major axis, providing a robust estimate of size and growth kinetics.

![**SAM Computational Metrology Pipeline and Feature Extraction Workflow**. The multi-pass architecture (Left) illustrates the transformation of raw optical data into high-fidelity segmented features. Key stages include: [A] Raw sensor input characterized by low Signal-to-Noise Ratio (SNR); [C] Initial segmentation mask generated via adaptive thresholding and 1D longitudinal intensity projection; and [D] Finalized high-precision segmentation after automated aberration correction and elliptical fitting for robust morphological characterization.](workflow.png)

#### **C. Temporal Tracking & Lineage Mapping**

* The algorithm records extrusion and division events in a **lineage-tractable data structure**, maintaining a continuous record of the system's evolution over time.

---

### 3. Post-Processing & Statistical Visualization

The final stage provides an automated suite for analyzing structured data. This includes generating division statistics and visualizing the kinetic profiles of thousands of features simultaneously, enabling rapid identification of system-level trends.

---




## 📖 Citation

If you use this tool in your research, please cite:

> *Banerjee, D. S., Stephenson, G., & Das, S. M. (2020). Segmentation and analysis of mother machine data: SAM. bioRxiv.* [[Link to Paper](https://doi.org/10.1101/2020.10.01.322685)]
