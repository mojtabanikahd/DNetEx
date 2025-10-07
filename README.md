# DNetEx Supplementary Code

This repository contains all code needed to reproduce the numerical experiments in:

> **Testing Sparse Differential Network in Gaussian Graphical Models with False Discovery Rate Control**  
> NeurIPS 2025 Submission

---

## 📂 Repository Structure
```text
├── SD1/                          # Synthetic data experiment - Scenario 1
│   ├── main.py                  # Runs the experiment and computes results
│   ├── plot-figures.ipynb       # Jupyter notebook for plotting results
│   ├── generate_dataset.py      # Contains dataset generation functions
│   ├── R_codes/                 # R scripts used in this experiment
│   ├── Data/                    # Directory for raw experimental data
│   └── Results/                 # Directory for saving generated plots
│
├── SD2/                          # Synthetic data experiment - Scenario 2
│   ├── main.py                  # Runs the experiment and computes results
│   ├── agg.ipynb                # Jupyter notebook for aggregating and plotting results
│   ├── generate_dataset.py      # Contains dataset generation functions
│   ├── R_codes/                 # R scripts used in this experiment
│   ├── Data/                    # Directory for raw experimental data
│   └── Results/                 # Directory for saving generated plots
│
├── fMRI_ABIDE/                   # Real data experiment using ABIDE fMRI dataset
│   ├── real_data_experiments.ipynb  # Entry-point notebook for running the experiment
│   ├── library.py               # Contains supporting functions for analysis
│   ├── R_codes/                 # R scripts used in this experiment
│   └── Results/                 # Directory for saving generated plots
│
├── requirements.txt             # List of required Python packages
└── README.md                    # Project overview and usage instructions
```

---

## 🛠️ Prerequisites

The experiments were conducted on **Windows 10** with the following software versions:

1. **Python 3.12.7**  
   _Run in a shell (e.g., PowerShell or Command Prompt)_  
   We recommend using Conda or pyenv for environment management. For example, with Conda:  
   ```bash
   conda create -n dnetex python=3.12.7
   conda activate dnetex
   ```

2. **R 4.3.1**  
   - Download and install via your web browser:
   ```
   https://cran.r-project.org/bin/windows/base/
   ```
   - On Windows, install Rtools43 to enable package compilation:
   ```
   https://cran.r-project.org/bin/windows/Rtools/
   ```
   - **Required R packages** (and all their dependencies) must be installed prior to running the analysis.  
   Run in the **R console** or via `Rscript`:
   ```r
   install.packages(c("GGMselect", "igraph", "DNetFinder"), dependencies = TRUE)
   ```
   

3. **System Dependencies**  
  - `git`, `curl` (for data download scripts)  
  - On Linux/Mac: standard build tools (`make`, `gcc`)

4. **Python Dependencies**  
   _Run in the previously activated Python environment (`dnetex`):_
   ```bash
   pip install -r requirements.txt
   ```

_Note_: Ensure that both your Python and R interpreters are on your system `PATH` so that the commands `python`, `pip`, `Rscript`, and `R` are recognized in your shell.

---

## 🚀 Quick Start  
   To reproduce our experiments and generate result plots, follow these steps:

### 🔹 Scenario 1: Synthetic Data (SCD1)

```bash
# Navigate to the Scenario 1 directory
cd SD1

# Run the main experiment script
python main.py

# Plot the results
# Open and run all cells of the following notebook using Jupyter:
#   plot-figures.ipynb
```
All plots will be saved under Results/.

 - **Outputs**:
    - `Data/*`: Logs of Empirical FDR and Statistical Power in Various Experiments.
    - `Results/fdr-control.pdf`: Comparison of Our Method and DNetFinder Using FDR Metric.
    - `Results/pr.pdf`: Comparison of Our Method and DNetFinder Using PR Metric.

### 🔹 Scenario 2: Synthetic Data (SCD2)

```bash
# Navigate to the Scenario 2 directory
cd SD2

# Run the main experiment script
python main.py

# Plot the results
# Open and run all cells of the following notebook using Jupyter:
#   agg.ipynb
```
All plots will be saved under Results/.

 - **Outputs**:
    - `Data/*`: Logs of Empirical FDR and Statistical Power in Various Experiments.
    - `Results/erdos.pdf`: The plot displays empirical FDR curves for p=100 and p=200.


### 🔹 Real data scenario: fMRI data (ABIDE)

```bash
# Navigate to the real data scenario directory
cd fMRI_ABIDE

# Plot the results
# Open and run all cells of the following notebook using Jupyter:
#   Real_data_experiment.ipynb
```
All plots will be saved under Results/.

 - **Outputs**:
    - Summaries of classification accuracy are provided in the cell output within the Validation section.
    - `Results/brain.pdf`: Differential brain connectivity between autism and control groups.

---
