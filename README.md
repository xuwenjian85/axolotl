# AXOLOTL: an accurate method for detecting aberrant gene expression in rare diseases using co-expression constraints

AXOLOTL is a computational method for detecting **expression outliers** in bulk RNA-seq data. It complementary integrates two methods — **OUTRIDER**, **OUTSINGLE**, and a non-parametric outlier detection method **AXO** — to robustly identify genes with aberrant expression levels across a cohort of samples.

> **Audience**: This document is intended for bioinformaticians and clinical lab scientists who want to detect expression outliers from their RNA-seq count matrices. Detailed instructions on how to install and run the pipeline are included to lower the barrrier.

---

## Table of Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation](#installation)
  - [1. Install Miniforge (Conda)](#1-install-miniforge-conda)
  - [2. Create Conda Environments](#2-create-conda-environments)
  - [3. Configure Paths](#3-configure-paths)
- [Quick Start](#quick-start)
- [Input Data Format](#input-data-format)
- [Output Description](#output-description)
- [Advanced Usage: Benchmark All Methods](#advanced-usage-benchmark-all-methods)
- [Troubleshooting](#troubleshooting)
  - [Memory Issues](#memory-issues)
  - [R / Bioconductor Issues](#r--bioconductor-issues)
  - [Runtime Errors](#runtime-errors)
- [Project Structure](#project-structure)
- [Citation](#citation)

---

## Overview

Axolotl runs two outlier-detection methods(OUTRIDER and OUTSINGLE) and produces a new combined score:

| Method | Language | Description |
|--------|----------|-------------|
| **OUTRIDER** | R | Autoencoder-based normalization + negative binomial test for expression outliers |
| **OUTSINGLE** | Python | SVD-based z-score estimation for outlier detection |
| **AXO** (this work) | Python | Combines OUTSINGLE P-values, OUTRIDER P-values, and a novel "deviation" feature using Local Outlier Factor (LOF) |

**Optional for advanced users**
Additionally, **ABEILLE** is available as a fourth method for benchmarking comparisons in this repo.

---

## System Requirements

| Requirement | Minimum | Recommended |
|-------------|---------|-------------|
| **OS** | Linux (x86_64) | Linux (x86_64) |
| **RAM** | 32 GB (AXO, ≤1000 samples) | 64 GB (if benchmark with ABEILLE) |
| **Disk** | 10 GB free | 50 GB+ (for large cohorts) |
| **Conda** | Miniforge3 or Miniconda3 | Miniforge3 |
| **Python** | 3.10 / 3.12 | Managed by Conda |
| **R** | 4.0+ | Managed by Conda |

> **Note for clinical labs**: A standard Linux server with 32 GB RAM is sufficient for AXO to process up to 1,000 samples. If you plan to run the full benchmark including ABEILLE, at least 64 GB RAM is required.

---

## Installation

### 1. Install Miniforge (Conda)

If you do not have Conda installed, install Miniforge (recommended over Anaconda/Miniconda for faster package resolution):

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3
source $HOME/miniforge3/etc/profile.d/conda.sh
conda init bash   # or zsh, depending on your shell
```

Restart your terminal after `conda init`.

### 2. Create Conda Environments

Axolotl requires **three** separate Conda environments. Create them from the provided YAML files:

#### Environment 1: `py3.12` (Main Python environment for AXO, OUTSINGLE)

```bash
conda env create -f py3.12.yaml
```

This environment includes: Python 3.12, pandas, numpy, scikit-learn, scipy, gseapy, matplotlib, seaborn, statsmodels, and more.


#### Environment 2: `r4_outrider` (For OUTRIDER)

```bash
conda env create -f r4_outrider.yaml
#conda create -n r4_outrider -c conda-forge -c bioconda \
#    bioconductor-outrider \
#    r-base=4.3
```

> **Note**: The OUTRIDER R package requires Bioconductor. The `bioconductor-outrider` package from bioconda will pull in all necessary dependencies including `BiocParallel`, `DESeq2`, etc.

#### Environment 3: `py3.10_tf2.8` (For ABEILLE, Optional)

```bash
conda env create -f py3.10_tf2.8.yaml
```

This environment includes: Python 3.10, TensorFlow 2.8, pandas, pyyaml.


Verify all environments are created:

```bash
conda env list
# You should see: py3.12,  r4_outrider, py3.10_tf2.8
```

### 3. Configure Paths

Edit `config.yaml` to match your local setup:

```yaml
paths:
  conda: /path/to/your/miniforge3/etc/profile.d/conda.sh
  script1: /path/to/axolotl
  abeille: /path/to/axolotl

  abeille_env_py: py3.10_tf2.8
  abeille_env_r: r4.0_abeille

  outrider_env: r4_outrider
  outrider_script: 'outrider_pred_dim.R'

  outsingle_env: py3.12
  outsingle_script: "outsingle.py"
  outsingle_dir: "/path/to/axolotl/outsingle"

  mymethod_env: py3.12
  mymethod_script: "run_mymethod.py"
```

**Key fields to update in the `config.yaml` file**:

| Field | Description |
|-------|-------------|
| `conda` | Absolute path to your `conda.sh` (e.g., `$HOME/miniforge3/etc/profile.d/conda.sh`) |
| `script1` | Absolute path to this repository's root directory |
| `outsingle_dir` | Absolute path to the `outsingle/` subdirectory |

In case you are new to conda, you can find your conda path with:

```bash
echo "$(conda info --base)/etc/profile.d/conda.sh"
```

---

## Quick Start

The simplest way to run AXOLOTL on your own data:

```bash
bash axo.sh <input_counts.tsv> <output_directory>
```

**Example with the provided test data:**

```bash
bash axo.sh testdata/toy_small.tsv output
```

This will:
1. Run **OUTRIDER** on your count matrix
2. Run **OUTSINGLE** on your count matrix
3. Combine results using the **AXO** method
4. Save all outputs to the `output/` directory

---

## Input Data Format

The input is a **tab-separated values (TSV)** file containing a gene expression count matrix:

- **Row 1 (header)**: `Gene` followed by sample names (e.g., `Sample_1`, `Sample_2`, ...)
- **Column 1**: Gene identifiers (e.g., `GeneA`, `GeneB`, ...)
- **Values**: Raw read counts (integers)

Example (`toy_small.tsv`):

```
Gene    GTEx_1  GTEx_2  GTEx_3  GTEx_4
Transcript_5    7638    4057    7978    3913
Transcript_7    3395    3714    7069    3346
Transcript_9    2148    1572    3274    1623
```

**Requirements**:
- The first column MUST be named `Gene`
- Values should be raw (non-normalized) **integer** counts
- Genes with very low counts across all samples should be filtered before input to AXOLOTL
- At least **30 samples** for reliable outlier detection, Recommended: **100+**

---

## Output Description

After running `axo.sh <input> <output_dir>`, the output directory will contain:

| File | Description |
|------|-------------|
| `<base>_outrider.txt.gz` | OUTRIDER adjusted P-value matrix (genes × samples) |
| `<base>_outrider.txt.gz_raw.gz` | Raw OUTRIDER results (long format with zScore, l2fc, normcounts, pValue, padjust) |
| `<base>_outrider.txt.gz_zScore.gz` | OUTRIDER z-score matrix |
| `<base>_outrider.txt.gz_l2fc.gz` | OUTRIDER log2 fold change matrix |
| `<base>_outrider.txt.gz_normcounts.gz` | OUTRIDER normalized count matrix |
| `<base>_outrider.txt.gz_pValue.gz` | OUTRIDER raw P-value matrix |
| `<base>_outsingle.txt.gz` | OUTSINGLE P-value matrix (genes × samples) |
| `<base>_axo.txt.gz` | **AXO combined outlier score** (genes × samples, lower = more outlier-like) |

For example, if the input is testdata/toy_small.tsv, then **base = toy_small**, producing output like toy_small_outrider.txt.gz.

### Interpreting AXO Scores

The AXO score is the **negative outlier factor** from Local Outlier Factor (LOF):
- **More negative values** indicate stronger expression outliers
- A typical threshold: values < -1.5 suggest potential outliers
- For clinical use, combine AXO scores with OUTRIDER and OUTSINGLE P-values for consensus calling

---

## (Optional) Benchmark All Methods

For benchmarking or running all four methods (OUTRIDER, OUTSINGLE, ABEILLE, AXO):

### Step 1: Prepare a task configuration file

You may need to use multiple datasets and sample subsets, therefore it's recommended to prepare a task configuration file.   

Create a tab-separated config file (see `testdata/task_config/t00_TOY_s128_g1000.config` for example):

| Column | Description |
|--------|-------------|
| `task` | Task ID (integer, starting from 0) |
| `Dname` | Dataset name |
| `cts` | Absolute path to the count matrix |
| `samples` | Absolute path to the sample subset file |
| `MyMethod` | Output path for AXO results |
| `OUTRIDER` | Output path for OUTRIDER results |
| `ABEILLE` | Output path for ABEILLE results |
| `OUTSINGLE` | Output path for OUTSINGLE results |

### Step 2: Prepare a sample subset file

See `testdata/samples/t00_TOY_s128_g1000.txt` for the format — a matrix where each row is a task and columns are sample names to include.

### Step 3: Run individual methods

```bash
# Run OUTRIDER only
bash run_outrider.sh <task_config> <task_id>

# Run OUTSINGLE only
bash run_outsingle.sh <task_config> <task_id>

# Run ABEILLE only (requires GPU for faster performance)
bash run_abeille.sh <task_config> <task_id>

# Run AXO only (needs OUTRIDER and OUTSINGLE results to be there ready first)
bash run_axo.sh <task_config> <task_id>

# OR Run all methods sequentially in one command
bash run.sh <task_config>
```

---

## Troubleshooting

### Memory Issues 

#### OUTRIDER runs out of memory (16 GB Servers)

**Problem**: OUTRIDER's autoencoder fitting uses significant memory, especially with >500 samples.

**Solutions**:
1. **Reduce the encoding dimension** (`dim_q0`). The default is `N/3` where N is the number of samples. You can manually set a smaller value:
   - In `parse_axo_directly.py`, change `div_q0 = 3` to `div_q0 = 5` or higher
   
2. **Use `SerialParam()` instead of parallel processing** (this is already the default in the provided scripts, but verify your OUTRIDER version uses it):
   ```r
   ods <- OUTRIDER(ods, BPPARAM = BiocParallel::SerialParam())
   ```

3. **Monitor memory usage**:
   ```bash
   # In another terminal, monitor memory
   watch -n 5 free -h
   ```

4. **Reduce gene count**: Pre-filter lowly-expressed genes before running OUTRIDER.

#### General memory recommendations

| Use case | Recommended RAM |
|----------|-----------------|
| AXO, ≤1000 samples | 32 GB |
| Benchmark include ABEILLE | ≥64 GB |

If your dataset exceeds available RAM:
```bash
# Check available memory before running
free -h
```

---

### R / Bioconductor Issues

#### `Error: package 'OUTRIDER' is not available`

**Problem**: OUTRIDER is not properly installed in the R environment.

**Solution**:
```bash
conda activate r4_outrider
conda install -c bioconda bioconductor-outrider
```

#### `Error in OUTRIDER(ods, ...)`: Convergence issues

**Problem**: OUTRIDER autoencoder fails to converge.

**Solution**: Increase the number of iterations:
```r
ods <- OUTRIDER(ods, q = q0, BPPARAM = SerialParam(), iterations = 10)
```
This is set in `outrider_pred_dim.R` — uncomment the iterations parameter line.

#### `Bioconductor version mismatch`

**Problem**: R package version conflicts when create R enviroments.

**Solution**: Create a fresh environment with a specific R version:
```bash
conda create -n **r4_outrider_new** -c conda-forge -c bioconda \
    r-base=**4.3** bioconductor-outrider
```

---

### Runtime Errors

#### `FileNotFoundError` for config.yaml

**Problem**: The script cannot find `config.yaml`.

**Solution**: Always run scripts from the repository root directory, or ensure `config.yaml` is in the same directory as the Python scripts being executed.

#### `ModuleNotFoundError: No module named 'general'`

**Problem**: Python cannot find the `general.py` module.

**Solution**: Ensure you are running scripts from the repository directory:
```bash
# Correct:
bash axo.sh testdata/toy_small.tsv output

# Incorrect (running from a different directory may fail):
cd /tmp && bash /path/to/axolotl/axo.sh ...
```

#### OUTSINGLE intermediate files not found

**Problem**: `outsingle.py` fails because `outsingle/` directory scripts are missing.

**Solution**: Ensure `outsingle_dir` in `config.yaml` points to the correct `outsingle/` directory containing the OUTSINGLE helper scripts (`fast_zscore_estimation.py`, `optht_svd_zs.py`, etc.).

#### `Permission denied` on shell scripts

**Problem**: Shell scripts are not executable.

**Solution**:
```bash
chmod +x axo.sh run.sh run_axo.sh run_outrider.sh run_outsingle.sh run_abeille.sh
```

---

## Project Structure

```
axolotl/
├── README.md                  # This documentation
├── config.yaml                # Path configuration (edit for your system)
│
├── axo.sh                     # Main entry point for end users
├── run.sh                     # Run all methods (benchmark mode)
├── run_axo.sh                 # Run AXO method only (benchmark mode)
├── run_outrider.sh            # Run OUTRIDER only (benchmark mode)
├── run_outsingle.sh           # Run OUTSINGLE only (benchmark mode)
├── run_abeille.sh             # Run ABEILLE only (benchmark mode)
│
├── parse_axo_directly.py      # Main pipeline: runs OUTRIDER + OUTSINGLE + AXO
├── parse_task_axo.py          # AXO task parser (benchmark mode)
├── parse_task_outrider.py     # OUTRIDER task parser (benchmark mode)
├── parse_task_outsingle.py    # OUTSINGLE task parser (benchmark mode)
├── parse_task_abeille.py      # ABEILLE task parser (benchmark mode)
│
├── run_mymethod.py            # AXO core algorithm (LOF + devi feature)
├── outsingle.py               # OUTSINGLE helper script
├── abeille_in_one.py          # ABEILLE wrapper (VAE + R IdentifyAGE)
├── general.py                 # Shared utilities (postprocessing, AUPRC)
│
├── outrider_pred_dim.R        # OUTRIDER R script (with custom dim_q0)
├── outrider_pred.R            # OUTRIDER R script (default parameters)
│
├── outrider.sh                # OUTRIDER shell wrapper
├── outsingle.sh               # OUTSINGLE shell wrapper
├── mymethod.sh                # AXO shell wrapper
├── abeille.sh                 # ABEILLE shell wrapper
│
├── py3.12.yaml                # Conda environment for AXO/OUTSINGLE
├── py3.10_tf2.8.yaml          # Conda environment for ABEILLE
├── r4.0_abeille.yaml          # Conda environment for ABEILLE R part
├── r4_outrider.yaml           # Conda environment for OUTRIDER
│
├── outsingle/                 # OUTSIGNAL SVD helper scripts
├── fig_script/                # Figure generation scripts
│
└── testdata/                  # Example dataset
    ├── toy_small.tsv          # Test count matrix (128 samples × 1000 genes)
    ├── samples/               # Sample subset definitions
    │   └── t00_TOY_s128_g1000.txt
    └── task_config/           # Task configuration files
        └── t00_TOY_s128_g1000.config
```

---

## Citation

If you use Axolotl in your research, please cite:

> *(Citation to be added upon publication)*

### Referenced Methods

- **OUTRIDER**: Brechtmann F, et al. *OUTRIDER: a statistical method for detecting aberrantly expressed genes in RNA sequencing data*.
- **OUTSINGLE**: Salkovic E, et al. *OutSingle: a novel method of detecting and injecting outliers in RNA-Seq count data using the optimal hard  threshold for singular values*.
- **ABEILLE**: Labory J, et al. *ABEILLE: a novel method for ABerrant Expression Identification empLoying machine LEarning from  RNA-sequencing data*.

---

## License

This project is provided for academic research use. Please contact the authors for commercial licensing inquiries.
