# NeuroArchitecture-FNC

This repository contains code and processed data associated with the manuscript:

**From micro-architecture to networks: canonical neurobiological architecture of functional brain organization across the lifespan and in psychopathology**

**Zenodo DOI:** https://doi.org/10.5281/zenodo.23573340


## Overview

This study examines the spatial correspondence between functional network connectivity (FNC) and 49 neurobiological markers, including 24 cortical cell types, 19 neurotransmitter receptors and transporters, and 6 mitochondrial markers, across the human lifespan and in psychopathology.

## Publicly Available Data

The `data/` directory contains five processed MATLAB files:

| File | Description |
|---|---|
| `molecular_nets_remove_dice_auto.mat` | Neurobiological similarity networks constructed from 49 cell-type and molecular atlas maps |
| `spa_dis.mat` | Dice similarity and spatial-distance matrices used as spatial covariates |
| `nm_age6_80_75.mat` | Normative FNC trajectories at 75 ages (6–80 years) |
| `dv_1_age7_80_74.mat` | One-year differences in normative FNC trajectories at 74 ages (7–80 years) |
| `disease_effect_sizes.mat` | Group-level FNC effect-size patterns for five psychiatric disorder groups and ten between-disorder comparisons |

These processed files support reproduction of the spatial correspondence and regression analyses without access to individual-level participant data.

## External Data Sources

The following source datasets and atlases are relevant for independently reconstructing the inputs:

| Data | Source |
|---|---|
| NeuroMark functional network atlas (53 ICNs) | [NeuroMark / GIFT](http://trendscenter.org/software/gift) |
| 24 cortical cell-type maps | [Feng et al. / Cellular and molecular associations with intrinsic brain organization](https://trendscenter.org/data/) |
| 19 neurotransmitter receptor/transporter maps | [Hansen et al. / Mapping neurotransmitter systems to the structural and functional organization of the human neocortex](https://github.com/netneurolab/hansen_receptors) |
| Six mitochondrial maps | [Mosharov et al. / A human brain map of mitochondrial respiratory capacity and diversity](https://neurovault.org/collections/16418/) |

The original, spatially prepared NIfTI maps used in S1 are not redistributed. Reconstructing the same inputs from the source atlases may require additional spatial processing, reference mapping, and compliance with the original data providers' access conditions.

Individual-level FNC data from healthy participants and patients are not redistributed in this repository due to the data-sharing agreements, access restrictions, and participant privacy requirements associated with the original datasets. Researchers wishing to access these data should apply through the respective data providers and comply with their data-use policies. To facilitate reproducibility, we provide the processed normative FNC trajectories, neurobiological similarity networks, and group-level disorder effect sizes required for the downstream analyses.

## Repository Structure

- `codes/` — MATLAB, R, and Python analysis scripts.
- `data/` — Publicly released processed analysis inputs.

## Analysis Workflow

1. **S1:** Reconstruct cell-type, neurotransmitter, and mitochondrial similarity networks from NIfTI maps (optional).
2. **S2:** Prepare lifespan and disorder analyses and calculate univariate spatial associations.
3. **S3:** Perform multivariable regression and relative importance analyses in R.
4. **S4–S5:** Postprocess statistical results and apply multiple-comparison correction.

**To reproduce the downstream analyses**, download this repository, configure `codes/audit_config.m`, and start from S2 using the supplied files in `data/`.

To reconstruct the neurobiological networks from source NIfTI maps, configure the external atlas directories and run S1 first. Normative-model reconstruction requires the original participant-level FNC and covariates, which must be obtained separately under the relevant data-use agreements.

## Software Requirements

The analyses were performed using the following software and package versions:

**MATLAB**
- MATLAB R2024b
- SPM12

**R**
- R 4.5.0 development version (2025-02-15, revision 87726)
- `relaimpo` 2.2.7
- `R.matlab`, `foreach`, `doParallel`, `parallel`

R analyses were executed within a Singularity container using the software versions listed above.

**Python**
- PCNtoolkit 0.30.post2
- NumPy, pandas, and SciPy


## Citation

If you use the code or processed data provided in this repository, please cite:

*From micro-architecture to networks: canonical neurobiological architecture of functional brain organization across the lifespan and in psychopathology.* Zenodo (2026). https://doi.org/10.5281/zenodo.23573340

**DOI:** [10.5281/zenodo.23573340](https://doi.org/10.5281/zenodo.23573340)


