# NeuroArchitecture-FNC

This repository contains the analysis code and processed data accompanying the manuscript:

**From micro-architecture to networks: canonical neurobiological architecture of functional brain organization across the lifespan and in psychopathology**

*Manuscript submitted to Nature Biomedical Engineering.*

## Overview

This study investigates the correspondence between functional network connectivity (FNC) and 49 neurobiological markers (24 cortical cell types, 19 neurotransmitter markers, and 6 mitochondrial markers) across the human lifespan and in psychiatric disorders.

## Repository Structure

- `codes/` — MATLAB, R, and Python analysis scripts.
- `data/` — Processed neurobiological networks, normative FNC data, spatial covariates, and group-level disorder effect sizes.

## Analysis Workflow

1. **S1:** Reconstruct neurobiological similarity networks from NIfTI atlases (optional).
2. **S2:** Prepare lifespan and disorder analyses and calculate univariate spatial correspondence.
3. **S3:** Perform multivariable regression and relative importance analyses.
4. **S4–S5:** Postprocess statistical results and apply multiple-comparison correction.

Configure paths in `codes/audit_config.m` before running the analyses. The processed data allow users to start directly from S2 without reconstructing the original atlases.

## Requirements

- MATLAB and SPM12
- R (`relaimpo`, `R.matlab`, `tidyverse`, `foreach`, `doParallel`)
- Python and PCNtoolkit 0.30.post2 (for normative modeling)

## Data Availability

Processed group-level data are included. Individual-level neuroimaging and clinical data are not redistributed.

## Citation

*From micro-architecture to networks: canonical neurobiological architecture of functional brain organization across the lifespan and in psychopathology.* Nature Biomedical Engineering (manuscript under consideration).
