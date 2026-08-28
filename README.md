# MONTE analysis

This repository contains the code and notebooks used to reproduce the analyses presented in the MONTE paper. If you are looking for tutorial or the main package, please check [MONTE](https://github.com/ylaboratory/MONTE).

---

## Contents

- [MONTE analysis](#monte-analysis)
  - [Contents](#contents)
  - [Folder structures](#folder-structures)
  - [Environment setup](#environment-setup)
    - [Python environment](#python-environment)
    - [Verify the installation](#verify-the-installation)
    - [R environment](#r-environment)
  - [Running the analyses](#running-the-analyses)
  - [Citation](#citation)
  - [Contact](#contact)

## Folder structures

> [!NOTE]
> All data used in this study, including the TCGA and PCAWG datasets and generated files, are available on [Zenodo](https://zenodo.org/records/22118662). After downloading the archive, unzip it, and move `data` into this folder. The scripts should then be able to access all required files.

This repository contains two main folders: `data` and `src`. The scripts expect the required datasets to be available in the `data` folder. The `src` folder is divided into `jupyter-notebooks`, which contains the primary analysis code, and `plots`, which contains the scripts used to generate the figures. The notebooks corresponding to the analyses in the paper are listed in the table below. Because the plotting scripts are labeled with their corresponding figure numbers, they are not described further here.

|Filename|Description|Corresponding section / figures in paper|
|---|---|---|
|`01_monte_pancancer_model_training`|Train the pan-cancer model and evaluate its performance on the test set|Pan-cancer model training|
|`02_monte_model_ablation_study`|Evaluate the contribution of each MONTE component|Figure 2A|
|`03_monte_model_probe_subset`|Assess how different probe subsets affect purity prediction|Figure 2C|
|`04_monte_correlation_between_purity_metrics`|Assess the correlation between predicted purity and other purity metrics|Figure 2B|
|`05_benchmark_cancer_cor`|Compare purity predictions from different methods|Figure 3A|
|`06_monte_pancancer_probe_correction`|Train the model and perform probe correction|Figure 5A|
|`07_benchmark_methylation_correction`|Compare probe-correction methods|Figure 5B|
|`08_monte_pcawg_fine_tuning`|Fine-tune the model on the PCAWG dataset and evaluate different training probe sizes|Figures 4A and 4B|
|`09_application_cancer_correction`|Fine-tune the pan-cancer model for a target cancer type|Cancer-specific models (released in the MONTE repository)|
|`10_application_BRCA_marker`|Identify BRCA markers after probe correction|Figure 6|

## Environment setup

### Python environment

The Python environment for this repository is managed using [`uv`](https://docs.astral.sh/uv/). Before proceeding, make sure that Git and `uv` are installed and available from your terminal.

MONTE currently requires Python 3.13 or later.

After cloning this repository, enter the project directory and synchronize the environment:

```bash
cd MONTE-analysis
uv sync --locked
```

> [!NOTE]
> MONTE is installed automatically from the GitHub repository when you initialize the environment. Because the repository includes pretrained models, installation may take some time.

This command creates a virtual environment in `.venv` and installs the exact dependency versions recorded in `uv.lock`, including the MONTE package.

If you use Jupyter through an editor such as VS Code, select the Python interpreter located in the `.venv` directory.

### Verify the installation

You can verify that MONTE was installed correctly by running:

```bash
uv run python -c "import monte; print('MONTE was imported successfully.')"
```

If this command completes without an error, the Python environment is ready.

For information about installing and using MONTE independently of this analysis repository, see the [MONTE repository](https://github.com/ylaboratory/MONTE).

### R environment

The R packages used by the R scripts are not managed by `uv`. Install the required packages in your R environment before running these scripts.


## Running the analyses

Download the TCGA and PCAWG methylation datasets from [Zenodo](https://zenodo.org/records/22118662), and move the `data` folder to the current directory.

The analyses can then be run using the notebooks and scripts provided in this repository.

## Citation

The MONTE preprint is available on [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.01.22.701164v2).

If you use MONTE or the code in this repository in your research, please cite:

```bibtex
@article {
    Kim2026.01.22.701164,
    author = {Kim, Mirae and Lee, Wei-Hao and Yao, Vicky},
    title = {MONTE: Methylation-based Observation Normalization and Tumor purity Estimation},
    elocation-id = {2026.01.22.701164},
    year = {2026},
    doi = {10.64898/2026.01.22.701164},
    publisher = {Cold Spring Harbor Laboratory},
    URL = {https://www.biorxiv.org/content/early/2026/01/23/2026.01.22.701164},
    eprint = {https://www.biorxiv.org/content/early/2026/01/23/2026.01.22.701164.full.pdf},
    journal = {bioRxiv}
}
```

## Contact

If you encounter a problem with this repository, please open a GitHub issue or contact Wei-Hao Lee at [wl61@rice.edu](mailto:wl61@rice.edu).