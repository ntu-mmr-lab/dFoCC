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
- Create a conda environment: `conda create -n $env_name -c conda-forge --file requirements.txt`
  - Remove the environment by `conda remove -n $env-name --all`
- Make sure that your Python installation in conda is correct (`python -V` should print `3.11.11`)

## Preparation

Users must preprocess the original data with the following procedure to use this program.

### Step 0: Convert mmCIF model and structure factor files in PDB to .PDB and .MTZ

The algorithm only accepts MTZ (structure factor dataset) and PDB (structural model) filetypes as input.

For structure factor and models on the Protein Data Bank (PDB), please download "Structure Factors (CIF)" and "Legacy PDB Format" respectively.

Raw structure factor needs to be converted into MTZ files before using the algorithm. We recommend using `gemmi cif2mtz` command provided by the Gemmi library, which is usually packaged with the CCP4 suite. Note that light-triggered state datasets usually contains raw and extrapolated structure factor blocks, and can be differentiated with the key "_diffrn.details". Users must explicitly provide the block index in the command.

```bash
gemmi cif2mtz 8Z1J-sf.cif $Fodark_file -B 1
gemmi cif2mtz 8Z3X-sf.cif $Fotr_file -B 2

```

In case the Protein Data Bank drops support for legacy PDB format downloads, users may convert the mmCIF model file into a PDB file using `gemmi convert` command.

```bash
# In case there is only a mmCIF model file provided on PDB
gemmi convert 8Z1J.cif $pr0max_dark_model
```

### Step 1: Generate calculated dark-adapted amplitudes and phases (Fcdark and PHIcdark)

Fcdark and PHIcdark will have to be calculated and stored in the same .MTZ file as Fodark. It is recommended using `sfall` within the CCP4 suite (executing the following command or using CCP4i GUI interface) to generate this file from the experimental dark-adapted dataset and its structural model.

```bash
sfall XYZIN $pr0max_dark_model HKLIN $Fodark_file HKLOUT $FoFcdark_file <<EOF
TITLE Calculate Dark Structure Factor
LABIN FP=FP SIGFP=SIGFP FREE=FreeR_flag
LABOUT FC=Fcdark PHIC=PHIcdark
MODE SFCALC XYZIN HKLIN
RESOLUTION $low_resolution $high_resolution
SYMMETRY $symmetry
END
EOF
```

### Step 2: Scale Fotr dataset against Fodark

In order to prepare the observed difference electron density map (DEDo), against which dFoCC will refine coordinates, The Fotr dataset should be scaled against Fodark. This can be done by first combining the two datasets using cad or sftools from the CCP4 suite, followed by scaling with scaleit from the same package (executing the following command or using CCP4i GUI interface).

```bash
# Combine dark and triggered structure factor and rename several columns
cad HKLIN1 $FoFcdark_file HKLIN2 $Fotr_file HKLOUT $FoFcdark_and_Fotr_file <<EOF
TITLE Combine Dark and Triggered Structure Factor
MONITOR BRIEF
LABIN file 1 E1=FreeR_flag E2=FP E3=SIGFP E4=Fcdark E5=PHIcdark
LABOUT file 1 E1=FreeR_flag E2=Fodark E3=SIGFodark E4=Fcdark E5=PHIcdark
CTYPIN file 1 E1=I E2=F E3=Q E4=F E5=P
RESOLUTION file 1 $low_resolution $high_resolution
LABIN file 2 E1=FP E2=SIGFP
LABOUT file 2 E1=Fotr E2=SIGFotr
CTYPIN file 2 E1=F E2=Q
RESOLUTION file 2 $low_resolution $high_resolution
EOF

# # Alternatively, combine them with `sftools`
# sftools <<EOF
# read $FoFcdark_file col FreeR_flag FP SIGFP Fcdark PHIcdark
# set label col FP
# Fodark
# set label col SIGFP
# SIGFodark
# read $Fotr_file col FP SIGFP
# set label col FP
# Fotr
# set label col SIGFP
# SIGFotr
# write $FoFcdark_and_Fotr_file
# end
# EOF

# Scale the Fotr dataset against Fodark
scaleit HKLIN $FoFcdark_and_Fotr_file HKLOUT $FoFcdark_and_Fotr_scaled_file <<EOF
TITLE Scale Triggered Structure Factor
RESOLUTION 10.0 2.5
NOWT
CONVERGE NCYC 4
CONVERGE ABS 0.001
CONVERGE TOLR -7
REFINE SCALE WILSON
AUTO
LABIN FP=Fodark SIGFP=SIGFodark FPH1=Fotr SIGFPH1=SIGFotr
EOF
```

### Step 3: Construct observed isomorphous difference structure factor (dFo)

The observed isomorphous difference structure factor (dFo) can be constructed with `phenix.fobs_minus_fobs_map` (isomorphous difference map in Phenix) from the two scaled observed datasets and phased by the dark-state structural model. The resulting .MTZ file contains the difference amplitude dataset (FoFo) and their associated phases (PHFc)

```bash
# Split the dark and the triggered dataset, so that Phenix could read them correctly
cad HKLIN1 $FoFcdark_and_Fotr_scaled_file HKLOUT $Fodark_split_file <<EOF
TITLE Split Fodark
MONITOR BRIEF
LABIN file 1 E1=FreeR_flag E2=Fodark E3=SIGFodark
LABOUT file 1 E1=FreeR_flag E2=Fodark E3=SIGFodark
CTYPIN file 1 E1=I E2=F E3=Q
RESOLUTION file 1 $low_resolution $high_resolution
EOF

cad HKLIN1 $FoFcdark_and_Fotr_scaled_file HKLOUT $Fotr_split_file <<EOF
TITLE Split Fotr
MONITOR BRIEF
LABIN file 1 E1=Fotr E2=SIGFotr
LABOUT file 1 E1=Fotr E2=SIGFotr
CTYPIN file 1 E1=F E2=Q
RESOLUTION file 1 $low_resolution $high_resolution
EOF

# # Alternatively, split them with `sftools`
# sftools <<EOF
# read $FoFcdark_and_Fotr_scaled_file
# write $Fodark_split_file col FreeR_flag Fodark Fcdark
# write $Fotr_split_file col Fotr SIGFotr
# END
# EOF

# Generate dFo structure factor
phenix.fobs_minus_fobs_map job_title="Generate dFo Structure Factor" \
  f_obs_1_file=$Fotr_split_file \
  f_obs_2_file=$Fodark_split_file \
  output_file=$dFo_file \
  f_obs_1_label=Fotr \
  f_obs_2_label=Fodark \
  phase_source=$pr0max_dark_model \
  high_res=$high_resolution low_res=$low_resolution \
  sigma_cutoff=3.0
```

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
