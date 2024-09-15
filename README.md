# MicroRNA-Based Myocardial Infarction Diagnosis

This repository contains the code for the published scientific paper:

"An Exploration into the Diagnostic Capabilities of MicroRNAs for Myocardial Infarction Using Machine Learning"

[Link to the paper will be added]

## Project Overview

The main objective of this research is to develop a diagnostic tool for myocardial infarction using a signature of microRNAs. The study utilizes machine learning techniques to analyze two publicly available datasets.

## Data Sources

The analysis is based on two publicly available datasets:
1. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61741
2. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29532

Note: The datasets are not included in this repository. The R code is designed to fetch and process the data directly from their online sources.

## Repository Contents

- `R_code.R`: R script for data preprocessing, initial analysis, differentially expression analysis, and feature selection. 
- `Python_code.ipynb`: Jupyter notebook containing the machine learning model implementation

## Workflow

1. The R script (`R_code.R`) reads and processes the online datasets, generating several output files.
2. Two key output files are produced:
   - `common_mirs_exp.csv`
   - `GSE29532_mirs_expression.txt`
3. These files are used as inputs for the Jupyter notebook (`Python_code.ipynb`), which implements the machine learning models.

## Requirements

- R (version 4.3.2)
- Python (version 3.9.7)
- Jupyter Notebook

## Citation

If you use this code or find it helpful in your research, please cite our paper:

[Citation details will be added]


## Contact

For any queries or further information, please contact Mehrdad Samadishadlou at mehrdad.samadi90@gmail.com.
