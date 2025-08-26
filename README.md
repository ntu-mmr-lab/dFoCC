# dFoCC

[![DOI](https://zenodo.org/badge/105618413.svg)](https://doi.org/10.5281/zenodo.16946872)

## Requirements

- CCP4 (only tested in version 8.0)
- Phenix (only tested in version 1.20.1)
- Coot (only tested in version 0.9.8)
- Miniconda (or other distribution of Anaconda)
  - Python 3.11 (will be installed in the environment)
    - Biopython 1.85
    - CCTBX 2025.1
    - Gemmi 0.7
    - NumPy 1.26
    - pandas 2.2
    - SciPy 1.14
    - TQDM 4.67
    - Matplotlib 3.10

## Installation

- Clone this repository: `git clone https://github.com/ntu-mmr-lab/dFoCC.git`
- Create a conda environment: `conda create -n {env-name} -c conda-forge --file requirements.txt`
  - Remove the environment by `conda remove -n {env-name} --all`
- Make sure that your Python installation in conda is correct (`python -V` should print `3.11.11`)

## Preparation

Users must preprocess the original data with the following procedure to use this program.

### Step 1: Generate calculated dark-adapted amplitudes and phases (Fcdark and PHIcdark)

Fcdark and PHIcdark will have to be calculated and stored in the same .MTZ file as Fodark. It is recommended to use sfall within the CCP4 suite (executing the following command or using CCP4i GUI interface) to generate this file from the experimental dark-adapted dataset and its structural model.

```bash
sfall HKLIN $Fodark_file HKLOUT $FoFcdark_file <<EOF
title Calculate Dark Structure Factor
LABIN FP=Fodark SIGFP=SIGFodark FREE=FREE
labout FC=FCalcdark PHIC=PHICalcdark
MODE SFCALC XYZIN HKLIN
resolution $low_resolution $high_resolution
symmetry $symmetry
END
EOF
```

### Step 2: Scale Fotr dataset against Fodark

In order to prepare the observed difference electron density map (DEDo), against which dFoCC will refine coordinates, The Fotr dataset should be scaled against Fodark. This can be done by first combining the two datasets using cad or sftools from the CCP4 suite, followed by scaling with scaleit from the same package (executing the following command or using CCP4i GUI interface).

```bash
scaleit HKLIN $FoFcdark_file HKLOUT $FoFcdark_scaled_file <<EOF
title Scale Light Structure Factor
resolution 10.0 2.5
NOWT
converge NCYC 4
converge ABS 0.001
converge TOLR -7
REFINE SCALE WILSON
auto
LABIN FP=Fodark SIGFP=SIGFodark FPH1=Fotr SIGFPH1=SIGFotr
EOF
```

### Step 3: Construct observed isomorphous difference map (DEDo)

DEDo can be constructed with `phenix.fobs_minus_fobs_map` (isomorphous difference map in Phenix) from the two scaled observed datasets and phased by the dark-state structural model. The resulting .MTZ file contains the difference amplitude dataset (FoFo) and their associated phases (PHFc)

### Step 4: Prepare configuration file

`./configuration.phil` should be modified according to the comments

## Usage

- Demo Mode
  - Run the whole algorithm only with the initial structure
  - If `*_log_cycle.csv` exist, modify the initial structure by the parameters specified in it before running the algorithm
  - `python iterate_dFoCC.py ./configuration.phil mode=demo`
- Probe Mode
  - Create the first library of structures to see if the reaction coordinates are suitable
  - `python iterate_dFoCC.py ./configuration.phil mode=probe`
- Run Mode
  - iteratively run the algorithm for each residue
  - `python iterate_dFoCC.py ./configuration.phil`
