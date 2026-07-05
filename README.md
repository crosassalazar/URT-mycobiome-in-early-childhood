# URT MYCOBIOME IN EARLY CHILDHOOD

This repository contains the datasets and statistical analysis scripts used in the publication titled *"Taxonomy, Determinants, and Potential Interactions of the Upper Respiratory Tract Mycobiome in Early Childhood."*

## Repository Contents

### Data_dictionary
Contains descriptions of the variables included in the metadata.csv file.

### Datasets

Contains the following datasets used in the statistical analyses:

- `asvtable.csv`: ASV count table.
- `metadata.csv`: Sample metadata and variables used in the statistical analyses.
- `taxonomy.csv`: Taxonomic assignments corresponding to the ASVs in the ASV count table.

### Manuscript_logs
Contains original session information, processing logs, and workflow outputs.

### QIIME2_scripts
Contains QIIME2 scripts used for sequence processing and taxonomic assignment.

### R_scripts
Contains R scripts for data processing, statistical analyses, figure generation, and the master script (`00_run_all.R`).

## Folder Structure

To ensure that all scripts run correctly, place the repository folder on your Desktop and name it `ITS`.

Example:

```text
ITS/
├── Data_dictionary/
├── Datasets/
├── Manuscript_logs/
├── QIIME2_scripts/
├── R_scripts/
└── README.md
```

**Important:** The R scripts use file paths that assume the repository folder is located on the Desktop and named `ITS`. The folder structure should be preserved exactly as provided in this repository.

## Usage

1. Download or clone this repository.
2. Place the repository folder on your Desktop.
3. Ensure that the folder is named `ITS`.
4. Do not modify the folder structure.
5. For sequence processing, use the scripts contained in the `QIIME2_scripts` folder.
6. To reproduce the complete statistical analyses, run `00_run_all.R` from the `R_scripts` folder.

Individual R scripts may also be executed separately in R or RStudio.

## Software Requirements

- R
- RStudio (recommended)
- QIIME2

## Data Statement

The datasets in this repository are provided solely for the purpose of reproducing the statistical analyses and results presented in the associated publication.

The datasets must not be used:

- To identify or re-identify any study participant.
- To contact study participants or their families.
- For clinical decision-making.
- For purposes inconsistent with the informed consent document or applicable data-use regulations governing the original study.

Users are expected to maintain the confidentiality of study participants and comply with all applicable ethical, institutional, and legal requirements.

By using these datasets, users agree to use them solely for research, educational, and reproducibility purposes.

## License

The code contained in this repository is released under the MIT License.

Copyright (c) 2026 Christian Rosas-Salazar

The code is provided "as is," without warranty of any kind.

The datasets are released under the Creative Commons Attribution 4.0 International (CC BY 4.0) License.

Users are encouraged to cite the associated publication when using these materials.

## Citation

If you use the data or code from this repository, please cite the associated publication:

Rosas-Salazar C, et al. *Taxonomy, Determinants, and Potential Interactions of the Upper Respiratory Tract Mycobiome in Early Childhood*. Under review.

## Correspondence

Christian Rosas-Salazar, MD, MPH  
Vanderbilt University Medical Center  
Email: c.rosas.salazar@vumc.org
