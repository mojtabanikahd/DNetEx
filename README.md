# DNetEx

This repository contains all material needed to reproduce the numerical experiments in:

> **DNetEx: FDR-controlled differential network analysis for knowledge discovery from graphs**  
---

## 📂 Repository Structure
```text
├── SD1/                          # Synthetic data experiment - Scenario 1 (Section 4.1)
│   ├── main.py                  # Runs the experiment and computes results
│   ├── plot-figures.ipynb       # Jupyter notebook for plotting results
│   ├── r_wrappers.py            # rpy2 wrappers calling R helpers (datasets, DiffNetFDR (pmat/pcor), SPDtrace)
│   ├── R_codes/                 # R scripts used in this experiment
│   ├── Data/                    # Directory for raw experimental data
│   └── Results/                 # Directory for saving generated plots
│
├── SD2/                          # Synthetic data experiment - Scenario 2 (Section 4.1)
│   ├── main.py                  # Runs the experiment and computes results
│   ├── agg.ipynb                # Jupyter notebook for aggregating and plotting results
│   ├── r_wrappers.py            # rpy2 wrappers calling R helpers (datasets, DiffNetFDR (pmat/pcor), SPDtrace)
│   ├── R_codes/                 # R scripts used in this experiment
│   ├── Data/                    # Directory for raw experimental data
│   └── Results/                 # Directory for saving generated plots
│
├── SD_dimension_variation/       # Appendix A.1 - Sensitivity to dimensionality (p up to 1000)
│   ├── main.py                  # Runs the experiment and computes results
│   ├── plot-figures.ipynb       # Jupyter notebook for plotting results
│   ├── r_wrappers.py            # rpy2 wrappers calling R helpers
│   ├── R_codes/                 # R scripts used in this experiment
│   ├── Data/                    # Directory for raw experimental data
│   └── Results/                 # Directory for saving generated plots
│
├── SD_screening_sensitivity_analysis/          # Appendix A.2 - Effect of the screening parameter c
│   ├── main.py                  # Runs the experiment and computes results
│   ├── agg.ipynb                # Jupyter notebook for aggregating and plotting results
│   ├── r_wrappers.py            # rpy2 wrappers calling R helpers
│   ├── R_codes/                 # R scripts used in this experiment
│   ├── Data/                    # Directory for raw experimental data
│   └── Results/                 # Directory for saving generated plots
│
├── SD_sample_splitting_sensitivity_analysis/   # Appendix A.3 - Effect of the sample-splitting ratio
│   ├── main.py                  # Runs the experiment and computes results
│   ├── agg.ipynb                # Jupyter notebook for aggregating and plotting results
│   ├── r_wrappers.py            # rpy2 wrappers calling R helpers
│   ├── R_codes/                 # R scripts used in this experiment
│   ├── Data/                    # Directory for raw experimental data
│   └── Results/                 # Directory for saving generated plots
│
├── SD_sample_splitting_randomness_analysis/    # Appendix A.4 - Sensitivity to random sample splitting
│   ├── main.py                  # Runs the experiment and computes results
│   ├── plot-figures.ipynb       # Jupyter notebook for plotting results
│   ├── r_wrappers.py            # rpy2 wrappers calling R helpers
│   ├── R_codes/                 # R scripts used in this experiment
│   ├── Data/                    # Directory for raw experimental data
│   └── Results/                 # Directory for saving generated plots
│
├── fMRI_ABIDE/                   # Real data experiment using ABIDE fMRI dataset (Section 4.2)
│   ├── Real_data_experiment.ipynb  # Entry-point notebook for running the experiment
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
   install.packages(c("GGMselect", "igraph", "DNetFinder", "DiffNetFDR"), dependencies = TRUE)
   ```
   The `DiffNetFDR` package provides the competing baselines used throughout the experiments: `pmat` (Xia et al., 2015) and `pcor` (Liu et al., 2017).
   - **Required R package SPDtrace**: The ``SPDtrace`` package is mandatory and is used by SD1, SD2, and the other experiments (via the ``SPDtrace`` Python wrapper in ``r_wrappers.py`` through rpy2). You must install [SPDtrace](https://github.com/mojtabanikahd/SPDtrace.git) **in that same R** that Python binds to—not via ``pip``. For example:

   ```r
   install.packages(c("devtools", "remotes"), dependencies = TRUE)
   devtools::install_github("mojtabanikahd/SPDtrace")
   ```

   Follow the upstream instructions (e.g. ``remotes``) if ``devtools`` is not available.

   

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
    - `Results/sd1-fdr-control.pdf`: Comparison of Our Method with the pmat and pcor methods Using FDR Metric.
    - `Results/sd1-pr.pdf`: Comparison of Our Method with the pmat and pcor methods Using PR Metric.

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
    - `Results/SD2_erdos_combined.pdf`: The plot displays empirical FDR curves under varying graph densities for p=100 and p=200.


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

## 🧪 Appendix Experiments (Sensitivity & Scalability)

These experiments correspond to the additional analyses reported in Appendix A of the manuscript. They focus exclusively on DNetEx and characterize its scalability, robustness, and sensitivity to key parameters. Each follows the same workflow as the experiments above: run `main.py`, then open and run all cells of the listed notebook. All plots are saved under the respective `Results/` directory.

### 🔹 Appendix A.1: Sensitivity to Dimensionality

```bash
# Navigate to the experiment directory
cd SD_dimension_variation

# Run the main experiment script
python main.py

# Plot the results
# Open and run all cells of the following notebook using Jupyter:
#   plot-figures.ipynb
```

 - **Outputs**:
    - `Data/*`: Logs of Empirical FDR and Statistical Power across dimensions (p = 100, 200, 500, 1000).
    - `Results/dimension-fdr-control.pdf`: Empirical FDR curves of DNetEx under varying problem dimensions.
    - `Results/dimension-recall-control.pdf`: Empirical power curves of DNetEx under varying problem dimensions.

### 🔹 Appendix A.2: Effect of the Screening Parameter `c`

```bash
# Navigate to the experiment directory
cd SD_screening_sensitivity_analysis

# Run the main experiment script
python main.py

# Plot the results
# Open and run all cells of the following notebook using Jupyter:
#   agg.ipynb
```

 - **Outputs**:
    - `Data/*`: Logs of Empirical FDR, Statistical Power, and computational time for c = 1, 2, 4, 8, 16, 32.
    - `Results/screening_sensitivity_analysis.pdf`: Empirical FDR and power curves of DNetEx under different values of `c`.

### 🔹 Appendix A.3: Effect of the Sample-Splitting Ratio

```bash
# Navigate to the experiment directory
cd SD_sample_splitting_sensitivity_analysis

# Run the main experiment script
python main.py

# Plot the results
# Open and run all cells of the following notebook using Jupyter:
#   agg.ipynb
```

 - **Outputs**:
    - `Data/*`: Logs of Empirical FDR and Statistical Power for splitting ratios r = 0.1, 0.3, 0.5, 0.7, 0.9.
    - `Results/sample_splitting_ratio_sensitivity_analysis.pdf`: Empirical FDR and power curves of DNetEx under different sample-splitting ratios.

### 🔹 Appendix A.4: Sensitivity to Random Sample Splitting

```bash
# Navigate to the experiment directory
cd SD_sample_splitting_randomness_analysis

# Run the main experiment script
python main.py

# Plot the results
# Open and run all cells of the following notebook using Jupyter:
#   plot-figures.ipynb
```

 - **Outputs**:
    - `Data/*`: Logs of Empirical FDP and Statistical Power across repeated random splits of a fixed dataset.
    - `Results/rand-fdr-control.pdf`: Empirical FDR (mean ± standard deviation) over random sample splits.
    - `Results/rand-power.pdf`: Empirical power (mean ± standard deviation) over random sample splits.

---
