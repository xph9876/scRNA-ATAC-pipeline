# scRNA-seq + scATAC-seq analysis pipeline
## Description
This pipeline is used to perform the cluster on scRNA-seq and scATAC-seq data using WNN

## Methods
This pipeline contains the following steps:
1. Normalize both scRNA-seq and scATAC-seq dataset
3. Perform dimension reduction and clustering
4. Automatic annotation for cell typs

## Dependencies
Most of the dependencies are listed in Conda yaml file __env.yaml__.
To install:
```bash
conda env create -n sc -f env.yaml
conda activate sc
```
## Usage
All the parameters are listed in **config.yaml**. After modification of the config file, run with snakemake
```bash
Snakemake -c 24
```

## File structures
1. **analyze.R**: Main analysis script in R
2. **utils.R**: Functions for analysis in R
3. **config.yaml**: Config file in YAML format
4. **env.yaml**: Conda environment in YAML format
6. **sample_results**: Example output
