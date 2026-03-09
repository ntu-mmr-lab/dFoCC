# dFoCC Tutorial Folder

This is the tutorial folder of dFoCC to show the reproducibility of this method. It includes:

- preprocessed datasets (inside `./input/`) for the input;
- a preparation script (`./prepare_input.sh`);
- a sample configuration file (`./configuration.phil`) which substantially reduce the number of structure models generate in the structure library;
- and some of the output files (inside `./sample_results/`) to compare with.

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

- Install Miniconda according to the instruction in [https://www.anaconda.com/docs/getting-started/miniconda/install](https://www.anaconda.com/docs/getting-started/miniconda/install), or use other Conda distribution
- Clone this repository: `git clone https://github.com/ntu-mmr-lab/dFoCC.git --depth 1`
- Go to this directory: `cd dFoCC/tutorial/`
- Create a conda environment: `conda create -n $env_name -c conda-forge --file ../requirements.txt`
  - Replace `$env_name` with any valid name for a Conda environment
  - Remove the environment by `conda remove -n $env-name --all`
- Make sure that your Python installation in Conda is correct (`python -V` should print `3.11.11`)

## Preparation

The input files inside `./input/` have already been preprocessed, so the preparation steps stated in [../README.md](../README.md#preparation) are not necessary.

If users would like to test this algorithm with other datasets, they may follow the following step, with the help of the `./prepare_input.sh` script.

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

### Step 1-3: Generate all input files for dFoCC

The `./prepare_input.sh` script preforms [the 3 steps stated in `../README.md`](../README.md#step-1-generate-calculated-dark-adapted-amplitudes-and-phases-fcdark-and-phicdark) automatically.

Users should modify the input filenames and some parameters in the script before executing it.

### Step 4: Prepare configuration file

`./configuration.phil` should be modified according to the comments.

For tutorial, this has been configured to use the preprocessed input files. However, users might need to tweak the `n_jobs` variable to let the parallel jobs run smoothly.

Also, the `step_per_run` parameters are smaller than the default values, which shrinks the total number of poses generated, so that the executing time would be much shorter than a full run.

### Step 5: Prepare an extrapolated dataset for the final conventional refinement

An extrapolated structure factor file is required for the optional conventional refinement step. This step refines the whole model except for regions that are dealt with by dFoCC main algorithm.

Here, we provided a manually extrapolated (N = 11) structure factor file `./input/nestimation_11.0.mtz`.

To skip this step, comment out the `input.extrapolated_data` block in `configuration.phil`.

## Usage

Assume that users are still in the `dFoCC/tutorial/` directory.

- Run Mode
  - iteratively run the algorithm for each residue
  - `python ../iterate_dFoCC.py ./configuration.phil`
- A set of sample results can be found in `./sample_results/`
  - Users may verify their own test results by comparing the structures with Coot, or compare the logged CSV file with `diff`
- Alternative modes (not necessary for the tutorial)
  - Demo Mode
    - Run the whole algorithm only with the initial structure
    - If `*_log_cycle.csv` exist, it will modify the initial structure by the parameters specified in it before running the algorithm
    - `python ../iterate_dFoCC.py ./configuration.phil mode=demo`
  - Probe Mode
    - Create the first library of structures to see if the reaction coordinates are suitable
    - `python ../iterate_dFoCC.py ./configuration.phil mode=probe`
