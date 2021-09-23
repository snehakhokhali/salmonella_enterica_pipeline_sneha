# Genomics pipeline for reproducible data analysis of Salmonella enterica genomes - Bioinformatics Practical Course

[TOC]
## Table of Contents
### Task 1: About Salmonella

Salmonella enterica is a highly diverse Gram negative bacterial species. A few S.enterica serovars includongy Typhi, Paratyphi A, B or C are highly adapted to the human host as their reservoir. They are the casuative agents of enteric fever(also known as typhoid fever or paratyphoid fever). Enteric fever is an invasive, life-threatening, systemic disease resulting more deaths. It is endemic in developing world in regions that lack water nad adequate sanitation, faciliating the spread of these pathogens via the faecal-oral route. 

Source: https://www.frontiersin.org/articles/10.3389/fmicb.2014.00391/full
### Task 2: Pipelining

#### creating the conda environment
 * environment.yml witth content:
``` 
channels:
        - bioconda
        - conda-forge
dependencies:
        - snakemake-minimal
```
 * conda environment 
```
conda env create -n salmonella-pipeline -f environment.yaml
```
 * activate the environment 
```
conda activate salmonella-pipeline
```
 * check installation:
```
snakemake --help
```

