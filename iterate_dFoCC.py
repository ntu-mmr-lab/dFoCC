#!/usr/bin/env python

# # Generation of sigma-weighted local difference maps and analysis (LOE, MMR)

# ## Imports

from os import makedirs, rename as rename
from os.path import isfile, isdir, dirname
from glob import glob

import subprocess
import re
from datetime import datetime

from numpy import nan, abs
import pandas as pd

from Bio import PDB
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
pdb_parser = PDBParser()
pdb_io = PDBIO()

# for beta analysis and future use
import gemmi
import scipy.stats

# Progress bar
# from tqdm.contrib.itertools import product
from tqdm.auto import tqdm

# for parallel processing
from multiprocessing import Pool

# for cli argument parsing
import argparse

# for wrapping functions as class
from dataclasses import dataclass

# for checking pdb files if they are the same
import hashlib

# not important: for type hinting
from typing import cast, Any
from Bio.PDB.Structure import Structure

# to prevent some warnings that are not important
import warnings
warnings.filterwarnings('ignore')

# local Imports
from dFoCC.utils.external import *
from dFoCC.reaction_coordinates import ReactionCoordinateSeries
from dFoCC.reaction_coordinates.aa_dihedrals import get_dihedrals, rotate_dihedrals
from dFoCC.reaction_coordinates.rotate_main_chain_general import rotate_main_chain_general, calculate_main_chain_general
from dFoCC.difference_maps.ccp4 import CCP4_DED_Utils
from dFoCC.difference_maps.unexplained_and_overmodeled import extract_unexplained_and_overmodeled_voxels
from plot_dFoCC import prepare_all_plots, prepare_linregress_plot

from configuration import *

@dataclass
class MapManipulation:
    rxn_coord_utils: ReactionCoordinateSeries
    ccp4_utils_with_config: CCP4_DED_Utils
    work_dir: str
    last_xyzin_name: str
    dark_xyzin_name: str
    dark_data_name: str
    light_data_name: str
    default_sigma_level: float
    calculated_sigma_level: float
    calculated_cutoff_density: float | None
    xyzregion_file: str
    sigma_mapin_name: str
    selection: str
    selection_radius: float
    resolution_factor: float
    keep_files: bool = False

    # ### Step 1: Build the PDB file with the modified moiety

    def modify_pdb(self, *paramsList: Parameters, prefix: str | None = None):
        """
        Create a PDB file with the moiety modified by several reaction coordinates,
        and record the properties of the altered moiety (e.g. dihedral angles, bond length)
        """
        struct: Structure = cast(Structure, pdb_parser.get_structure('test', self.last_xyzin_name))

        if prefix == None:
            prefix = generate_prefix(paramsList[0])

        # in case there are no parameters at all in the list (for demo mode),
        # the default structural logging and checking should also be done
        if len(paramsList) == 0:
            is_all_ok, obj = self.rxn_coord_utils.calculate_moiety_properties(struct, self.rxn_coord_utils.get_best_params())
        else:
            is_all_ok, obj = True, {}
        for params in paramsList:
            # modify_moiety(struct, params)
            self.rxn_coord_utils.modify_moiety(struct, params)
            is_this_ok, obj = self.rxn_coord_utils.calculate_moiety_properties(struct, params)
            if not is_this_ok:
                is_all_ok = False

        pdb_io.set_structure(struct)
        pdb_io.save(f"{self.work_dir}{prefix}tmp.pdb") # , select=RemoveHydrogen()

        # ok, obj = calc_moiety_properties(struct, BEST_STRUCTURE_PARAM if len(paramsList) == 0 else paramsList[-1])

        # Add the header part of the original pdb to the modified one
        subprocess.run(f"cat {self.work_dir}{XYZIN_NAME}.header {self.work_dir}{prefix}tmp.pdb > {self.work_dir}{prefix}.pdb",shell=True)
        subprocess.run(['rm', f'{self.work_dir}{prefix}tmp.pdb'])

        md5 = hashlib.md5()
        with open(f'{self.work_dir}{prefix}.pdb', 'rb') as f:
            while True:
                data = f.read(65536)
                if not data:
                    break
                md5.update(data)

        obj.update({ 'hash': int(md5.hexdigest(), 16) })

        return is_all_ok, prefix, obj

    # ### Step 2: Generate new maps from the modified structure
    def make_maps(self, prefix: str):
        log1, log2, log3 = '', '', ''

        try:
            log1 = self.ccp4_utils_with_config.make_calc_difference_sf(
                model_in=f'{self.work_dir}{prefix}.pdb',
                dark_data=self.dark_data_name,
                light_data=self.light_data_name,
                data_out=f'{self.work_dir}dFc.{prefix}.mtz',
                dark_model_in=self.dark_xyzin_name,
                keep_files=self.keep_files,
            )
            log1 += u'\u2500' * 30 + '\n'
            log1 += self.ccp4_utils_with_config.make_difference_map(
                data_in=f'{self.work_dir}dFc.{prefix}.mtz',
                map_out=f'{self.work_dir}dFc.{prefix}.map',
            )

            # remove all used files unless we need them (e.g. regenerate after each cycle)
            used_files = [
                # dark_data, # input file and still in use
                # light_data, # input file and still in use
                # f'{workdir}{prefix}.pdb', # generated model, still in use
                f'{self.work_dir}dFc.{prefix}.mtz', # difference structure factors between generated model and calculated dark
                # f'{workdir}dFc.{prefix}.map', # still in use
            ]
            if not self.keep_files:
                for used_file in used_files:
                    subprocess.run(['rm', used_file])

            if self.calculated_cutoff_density:
                rms = self.calculated_cutoff_density
            else:
                log2, rms = self.ccp4_utils_with_config.calculate_sigma_from_mean_density(map_in=f'{self.work_dir}dFc.{prefix}.map')
                log2 += u'\u2500' * 30 + '\n'
            log2 += self.ccp4_utils_with_config.generate_sigmadiff_maps(
                map_in=f'{self.work_dir}dFc.{prefix}.map',
                map_out=f'{self.work_dir}dFc.{prefix}.sigma{self.default_sigma_level}.map',
                rms=rms,
                sigma_level=self.calculated_sigma_level,
                keep_files=self.keep_files,
            )

            # generate_sigma_cutoff_map(
            #     map_in=f'{self.work_dir}dFc.{prefix}.map',
            #     map_out=f'{self.work_dir}dFc.{prefix}_sigma{self.default_sigma_level}.map',
            #     sigma_level=self.calculated_sigma_level,
            # )
            # log2 = ''

            # This part is moved from Step 3, to prevent storing large files
            self.ccp4_utils_with_config.extract_region(
                # xyzin_ref=f'{workdir}{XYZINREF_NAME}',
                xyzin_ref=self.dark_xyzin_name,
                xyzin=f'{self.work_dir}{prefix}.pdb',
                xyzout=f'{self.work_dir}{self.xyzregion_file}'.replace('*', prefix, 1)
            )

            # mask +/- sigma cleared calculated difference map by the moiety
            log3 = map_box(
                pdb_file=f'{self.work_dir}{self.xyzregion_file}'.replace('*', prefix, 1),
                ccp4_map_file=f'{self.work_dir}dFc.{prefix}.sigma{self.default_sigma_level}.map',
                output_file_name_prefix=f'{self.work_dir}dFc.{prefix}.masked_sigma{self.default_sigma_level}',
                selection=self.selection,
                mask_atoms=True,
                mask_atoms_atom_radius=self.selection_radius,
                resolution_factor=self.resolution_factor,
                wrapping=True,
            )

            # remove all used files unless we need them (e.g. regenerate after each cycle)
            used_files = [
                f'{self.work_dir}dFc.{prefix}.map',
                f'{self.work_dir}dFc.{prefix}.sigma{self.default_sigma_level}.map', # used in `map_box`
            ]
            if not self.keep_files:
                for used_file in used_files:
                    subprocess.run(['rm', used_file])
            # Checkpoint: check if the last output file exists
            if not isfile(f'{self.work_dir}dFc.{prefix}.masked_sigma{self.default_sigma_level}.ccp4'):
                with open(f"{self.work_dir}log_dFoCC.make_maps.map_box.txt", "w") as f:
                    f.write(log3)
                raise Exception(log3)

        except Exception as err:
            raise Exception(log1 + log2 + log3 + str(err))

        return log1, log2, log3

    # ### Step 3: Calculate the CC
    # The `calculate_cc` function below calculates the correlation coefficient between dFo and dFc maps inside the selected moiety
    def calculate_cc(self, prefix: str) -> tuple[str, str, float, float | None]:
        # mask the observed difference map by the moiety
        log4 = map_box(
            pdb_file=f'{self.work_dir}{self.xyzregion_file}'.replace('*', prefix, 1),
            ccp4_map_file=self.sigma_mapin_name,
            output_file_name_prefix=f'{self.work_dir}dFo.{prefix}.masked_sigma{self.default_sigma_level}',
            selection=self.selection,
            mask_atoms=True,
            mask_atoms_atom_radius=self.selection_radius,
            resolution_factor=self.resolution_factor,
            wrapping=True,
            bounds_match_this_file=f'{self.work_dir}dFc.{prefix}.masked_sigma{self.default_sigma_level}.ccp4',
        )

        # Checkpoint: check if the last output file exists
        if not isfile(f'{self.work_dir}dFc.{prefix}.masked_sigma{self.default_sigma_level}.ccp4') \
            or not isfile(f'{self.work_dir}dFo.{prefix}.masked_sigma{self.default_sigma_level}.ccp4'):
            with open(f"{self.work_dir}log_dFoCC.calculate_cc.map_box.txt", "w") as f:
                f.write(log4)
            raise Exception(log4)

        if IS_BETA_ANALYSIS:
            slope, cc = None, None
            # calculate the slope of linear regression (?) of the moiety for the generated difference map and the observed difference map
            map_dFc = gemmi.read_ccp4_map(f'{self.work_dir}dFc.{prefix}.masked_sigma{self.default_sigma_level}.ccp4')
            map_dFo = gemmi.read_ccp4_map(f'{self.work_dir}dFo.{prefix}.masked_sigma{self.default_sigma_level}.ccp4')
            # indices = pd.MultiIndex.from_product([range(s) for s in map_dFc.grid.shape], names=['a', 'b', 'c'])
            #
            # map_df_1 = pd.DataFrame(map_dFc.grid.array.flatten(), index=indices)
            # map_df_2 = pd.DataFrame(map_dFo.grid.array.flatten(), index=indices)
            # map_df_join = map_df_1.merge(map_df_2, on=indices).dropna()
            try:
                # linregress_result = scipy.stats.linregress(map_df_join['0_x'], map_df_join['0_y'])
                linregress_result = scipy.stats.linregress(map_dFo.grid.array.flatten(), map_dFc.grid.array.flatten()) # type: ignore
                slope, cc = linregress_result.slope, linregress_result.rvalue # type: ignore
            except ValueError:
                slope, cc = nan, nan

        else:
            # calculate the CC of the moiety for the generated difference map and the observed difference map
            log4 += overlapmap(
                mapin1=f'{self.work_dir}dFc.{prefix}.masked_sigma{self.default_sigma_level}.ccp4',
                mapin2=f'{self.work_dir}dFo.{prefix}.masked_sigma{self.default_sigma_level}.ccp4',
                script='''
                    CORRELATE SECTION
                    END
                '''
            )

            cc, slope = nan, None

            for line in log4.splitlines():
                match = re.search(r'Total correlation coefficient is : *(-?\d\.\d+)', line)
                if match:
                    cc = float(match.group(1))

        # remove all used files
        used_files = [
            # f'{self.work_dir}dFc.{prefix}_sigma{SIGMA_LEVEL}.map', # already removed in Step 2
            *glob(f'{self.work_dir}{self.xyzregion_file}'.replace('*', prefix, 1)),
            f'{self.work_dir}{prefix}.pdb', # the library of structures
            # f'{workdir}dFc.{prefix}_masked_sigma{SIGMA_LEVEL}.map', # if use CCP4 `MAPMASK` instead of `phenix.map_box`, remove this instead of the next one
            *glob(f'{self.work_dir}dFc.{prefix}.masked_sigma{self.default_sigma_level}.*'),
            # f'{workdir}dFo_sigma{SIGMA_LEVEL}{prefix}.map', # if use CCP4 `MAPMASK` instead of `phenix.map_box`, remove this instead of the next one
            *glob(f'{self.work_dir}dFo.{prefix}.masked_sigma{self.default_sigma_level}.*'),
        ]
        if not self.keep_files:
            for used_file in used_files:
                subprocess.run(['rm', used_file])

        return log4, prefix, cc, slope

    def generate_coot_script(self, best_prefix: str, rms: float | None = None):
        dFo_read_script = f'(make-and-draw-map "{self.light_data_name}" "FoFo" "PHFc" "" 0 1)' if self.light_data_name.endswith('.mtz') \
            else f'(handle-read-ccp4-map "{self.light_data_name}" 1)'
        with open(f'{self.work_dir}coot.scm', 'w') as coot_script_file:
            coot_script_file.write(f''';;;; Coot script file in scheme
(set-recentre-on-read-pdb 0)
(set-nomenclature-errors-on-read "ignore")

;; load the dark and the calculated light structure
(read-pdb "{self.dark_xyzin_name}")
(read-pdb "{self.work_dir}{best_prefix}.pdb")

;; read the dFo region-extracted sigma-cutoff map
;; and show almost all signal on the map
;; (handle-read-ccp4-map "{self.work_dir}dFo.{best_prefix}.masked_sigma{self.default_sigma_level}.ccp4" 1)
;; (set-last-map-contour-level-by-sigma 0.01)

;; read the entire dFo struct factor
;; and apply suitable sigma cutoff on it
{dFo_read_script}
(set-last-map-contour-level-by-sigma {self.default_sigma_level:.2f})
(set-last-map-colour 0.25 1.0 0.25)
; (export-map-fragment 2 (rotation-centre-position 0) (rotation-centre-position 1) (rotation-centre-position 2) 20.0 "{self.work_dir}dFo_for_pymol.map")

;; read the dFc region-extracted sigma-cutoff map
;; show almost all signal on the map, and change its colour into cyan-yellow pair
;; (handle-read-ccp4-map "{self.work_dir}dFc.{best_prefix}.masked_sigma{self.default_sigma_level}.ccp4" 1)
;; (set-last-map-contour-level-by-sigma 0.01)

;; read the entire dFo struct factor
;; apply suitable sigma cutoff on it
;; and change its colour into cyan-yellow pair
(make-and-draw-map "{self.work_dir}dFc.{best_prefix}.mtz" "{'FoFo' if self.ccp4_utils_with_config.should_use_phenix_for_dFc else 'dFc'}" "{'PHFc' if self.ccp4_utils_with_config.should_use_phenix_for_dFc else 'dPHIc'}" "" 0 1)
{f"(set-last-map-contour-level {rms*self.default_sigma_level:.4f})" if rms else f"(set-last-map-contour-level-by-sigma {self.calculated_sigma_level:.2f})"}
(set-last-map-colour 0.0 1.0 1.0)
; (export-map-fragment 3 (rotation-centre-position 0) (rotation-centre-position 1) (rotation-centre-position 2) 20.0 "{self.work_dir}dFc_for_pymol.map")
''')

# ;; # FIXME: for ligands that have no C_beta, is there any way to handle it? (hard coded for now)
# (set-go-to-atom-chain-residue-atom-name "{CHAIN_ID}" {RESIDUE_ID if isinstance(RESIDUE_ID, int) else RESIDUE_ID[1]} "{CENTER_ATOM}")
# (set-zoom 15.00)

def generate_phenix_refinement_script(xyzin_name: str, mtzin_name: str, output_dir: str, script_name: str, selection_for_conventional_refinement: str, col_names: dict[str, str | None], working_params: Any):
    makedirs(dirname(script_name), exist_ok=True)
    if col_names['Fext'] and col_names['Fext_sigma']:
        fext_labels = f"labels = {col_names['Fext']} {col_names['Fext_sigma']}"
    else:
        fext_labels = ''
    if col_names['Fext_free']:
        fext_free_labels = f"label = {col_names['Fext_free']}"
    else:
        fext_free_labels = ''
    with open(script_name, "w") as phenix_script_file:
        phenix_script_file.write(f"""refinement {{
  input {{
    pdb {{
      file_name = "{xyzin_name}"
    }}
    xray_data {{
      file_name = "{mtzin_name}"
      {fext_labels}
      high_resolution = {working_params.high_resolution}
      low_resolution = {working_params.low_resolution}
      r_free_flags {{
        file_name = "{mtzin_name}"
        {fext_free_labels}
        test_flag_value = 0
      }}
    }}
  }}
  output {{
    prefix = "{output_dir}refined"
    job_title = "conventional real-space refinement after dFoCC real-space refinement"
    write_def_file = False
  }}
  electron_density_maps {{
    map_coefficients {{
      map_type = "2mFo-DFc"
      mtz_label_amplitudes = "2FOFCWT"
      mtz_label_phases = "PH2FOFCWT"
      fill_missing_f_obs = True
    }}
    map_coefficients {{
      map_type = "2mFo-DFc"
      mtz_label_amplitudes = "2FOFCWT_no_fill"
      mtz_label_phases = "PH2FOFCWT_no_fill"
    }}
    map_coefficients {{
      map_type = "mFo-DFc"
      mtz_label_amplitudes = "FOFCWT"
      mtz_label_phases = "PHFOFCWT"
    }}
    map_coefficients {{
      map_type = "anomalous"
      mtz_label_amplitudes = "ANOM"
      mtz_label_phases = "PHANOM"
    }}
  }}
  refine {{
    strategy = *individual_sites individual_sites_real_space *rigid_body \\
               individual_adp group_adp tls occupancies group_anomalous
    sites {{
      rigid_body = {selection_for_conventional_refinement}
      individual = {selection_for_conventional_refinement}
      torsion_angles = {selection_for_conventional_refinement}
    }}
  }}
  main {{
    number_of_macro_cycles = 5
    nproc = {working_params.n_jobs}
  }}
  hydrogens {{
    refine = individual *riding Auto
  }}
  target_weights {{
    optimize_xyz_weight = True
    optimize_adp_weight = True
  }}
}}
""")

def main(
    rxn_coord_utils: ReactionCoordinateSeries,
    init_xyzin_name: str,
    demo_dir: str,
    probe_dir: str,
    n_value: float,
    selection: str,
    col_names: dict[str, str | None],
    log_cycle_name='log_cycle.csv',
    generate_work_dir_name=lambda _, cycle_num: 'cycle*/' if cycle_num=='*' else f'cycle{cycle_num:>03}',
    max_cycle=-1,
    is_beta_analysis=False, # FIX
    demo_parameters: list[Parameters] = [],
    is_demo_mode=False,
    is_probe_mode=False,
    is_verbose_mode=False,
):
    ccp4_utils = CCP4_DED_Utils(
        col_names=col_names,
        high_resolution=args.high_resolution,
        low_resolution=args.low_resolution,
        symmetry=args.space_group,
        grid_ccp4=GRID_CCP4,
        default_sigma_level=args.input.difference_data.sigma_level,
        atom_pattern=['ATOM  ', 'HETATM'],
        should_use_phenix_for_dFc=SHOULD_USE_PHENIX_FOR_DFC,
    )

    print(f'Starting: {datetime.now().replace(microsecond=0).isoformat(" ")}')
    log0 = ''
    light_input_data_name: str = args.input.difference_data.file_name
    light_map_name = f"{PREPARING_DIR}{light_input_data_name.rsplit('/', 1)[-1].replace('.mtz', '.map')}"
    sigma_mapin_name = SIGMA_MAPIN_NAME or light_map_name.replace('.map', '_sigma_cutoff.map')
    # the sigma cutoff dFo map will be generated first in the input directory
    if isfile(sigma_mapin_name) and not isfile(light_map_name):
        rms = None
        light_map_name = sigma_mapin_name
        if is_beta_analysis:
            raise NotImplementedError('Sigma cutoff map as input file is not supported in beta analysis mode\n(difficult to extract sigma from a sigma cutoff map)')
    else:
        print('Prepare the input dFo map')

        if not isdir(PREPARING_DIR):
            subprocess.run(['mkdir', PREPARING_DIR])

        if light_input_data_name.endswith('.mtz'):
            log0 = ccp4_utils.make_difference_map(
                data_in=light_input_data_name,
                map_out=light_map_name,
                is_dFo=True,
            )
            log0 += u'\u2500' * 30 + '\n'
        else:
            light_map_name = light_input_data_name
        log, rms = ccp4_utils.calculate_sigma_from_mean_density(map_in=light_map_name)
        log0 += log
        log0 += u'\u2500' * 30 + '\n'
        log0 += ccp4_utils.generate_sigmadiff_maps(
            map_in=light_map_name,
            map_out=sigma_mapin_name,
            rms=rms,
            sigma_level=args.input.difference_data.sigma_level,
            keep_files=is_demo_mode or is_verbose_mode,
        )
        log0 += u'\u2500' * 30 + '\n'
        for line in log0.splitlines():
            if re.search('error', line, flags=re.IGNORECASE):
                raise RuntimeError(log0)

    # sanity check: at least one atom is selected with the `SELECTION` query string
    log_atom_seletion = pdb_atom_selection(init_xyzin_name, inselection=selection)
    if (match := re.search(r'^  (\d+) atoms? selected', log_atom_seletion, re.MULTILINE)):
        if match.group(1) == '0':
            raise RuntimeError(f'No atom matching the `{selection = }` found in {init_xyzin_name}')
    else:
        raise RuntimeError(log_atom_seletion)
    log0 += log_atom_seletion

    # ------------------------------

    # ## Auto Recovery
    # 
    # If the csv exists, recover from the last cycle.
    # Otherwise, create it with a header.

    # make a list of dicts of cycle log for later use
    log_cycle_list: list[dict] = []

    if is_demo_mode or is_probe_mode:
        round_num = 1
        cycle_num = 1

        last_work_dir = ''
        best_prefix = init_xyzin_name.replace('.pdb', '')

    elif isfile(log_cycle_name):
        csv_data_frame = pd.read_csv(log_cycle_name)
        last_round_log = csv_data_frame.iloc[-1] # last log -> the latest finished cycle

        round_num, cycle_num = 1, 1

        # check if the structure is the same as the last cycle
        if last_round_log['prefix'] == generate_prefix(rxn_coord_utils.get_best_params()):
            # round_num, cycle_num = last_round_log['round'] + 1, 1
            round_num, cycle_num = last_round_log['round'] + 1, last_round_log['cycle'] + 1
        else:
            round_num, cycle_num = last_round_log['round'], last_round_log['cycle'] + 1

        # Check if the last cycle is the final one
        is_last_round = round_num > 1 and rxn_coord_utils.generate_params_list(round_num - 1).is_last_round
        # is_last_round = round_num > 1 and generate_loop_params(round_num - 1, INIT_PARAM_INTERVALS, STRUCT_PER_PARAM)[0]
        assert not is_last_round, "The dFoCC algorithm has already reached the end."

        last_work_dir = generate_work_dir_name(last_round_log['round'], last_round_log['cycle'])

        # convert the LOG_ROUND_NAME csv to a list of dicts for later use
        log_cycle_list = csv_data_frame.to_dict('records')
        best_prefix = last_round_log['prefix']

    else:
        # during a fresh rerun from Cycle 1
        # if there are executed directories (`round*_*/` or `cycle_*/`) left,
        # archive all of them into the `ARCHIVE_DIR`
        if archiving_dirs := glob(generate_work_dir_name('*', '*')):
            subprocess.run(['mkdir', ARCHIVE_EXE_DIR])
            for archiving_dir in archiving_dirs:
                subprocess.run(['mv', archiving_dir, f'{ARCHIVE_EXE_DIR}{archiving_dir}'])

        # # create a log of all structures and a brief log of each cycle
        # with open(LOG_ROUND_NAME, 'w', newline='') as f:
        #     # create the csv writer
        #     csv_writer = csv.writer(f, delimiter=('\t' if LOG_ROUND_NAME[-4:] == '.tsv' else ','))
        #     csv_writer.writerow(['round', 'cycle', 'prefix', 'CC'])

        round_num = 1
        cycle_num = 1

        last_work_dir = ''
        best_prefix = init_xyzin_name.replace('.pdb', '')

        # log the parameters of the dark structure before the main loop

        dark_struct: Structure = cast(Structure, pdb_parser.get_structure('', args.input.pdb.dark_file_name))
        is_ok, props_obj = rxn_coord_utils.calculate_moiety_properties(dark_struct, rxn_coord_utils.get_best_params())
        md5 = hashlib.md5()
        with open(args.input.pdb.dark_file_name, 'rb') as f:
            while True:
                data = f.read(65536)
                if not data:
                    break
                md5.update(data)

        log_cycle_list.append({
            'round': 0,
            'cycle': 0,
            'prefix': generate_prefix(rxn_coord_utils.get_best_params()),
            'ok': is_ok
        } | props_obj | {
            'hash': int(md5.hexdigest(), 16),
            'CC': float('nan'),
            'unexplained': float('nan'),
            'overmodeled': float('nan'),
            'dFo_voxels': float('nan'),
            'dFc_voxels': float('nan'),
            'ttl_voxels': float('nan'),
            'reason': 'dark',
        })
        # pd.DataFrame(log_round_list).to_csv(log_cycle_name, sep='\t' if log_cycle_name[-4:] == '.tsv' else ',', index=False)

    # **MAIN LOOP**
    while True:
        if is_demo_mode:
            print('Start making demo structure and maps for dFoCC')
            work_dir = demo_dir
        elif is_probe_mode:
            print('Start making a library of structures for dFoCC')
            work_dir = probe_dir
        else:
            # print(f"Start Round {round_num}_{cycle_num}")
            print(f"Start Cycle {cycle_num} (Round {round_num})")
            work_dir = generate_work_dir_name(round_num, cycle_num)

        if isdir(work_dir):
            subprocess.run(['mv', work_dir, f"{work_dir[:-1]}{ARCHIVE_ROUND_SUFFIX}"])
        subprocess.run(['mkdir', work_dir])

        # copy all necessary files to the working directory of the next cycle
        subprocess.run(['cp', f'{last_work_dir}{best_prefix}.pdb', f'{work_dir}{XYZIN_NAME}'])

        print(f'Work directory: {work_dir}')

        # ----------------------------------------------------------------

        # Test of the `generate_loop_params` function,
        # to check whether the function works in this cycle
        is_last_round, loop_through_params = rxn_coord_utils.generate_params_list(round_num)
        total_iterations = len(loop_through_params)

        print(f'Maximum structure made for this cycle: {total_iterations}')

        # ----------------------------------------------------------------

        # Grab the header part of the original pdb, by excluding everything not in the header
        # As biopython don't generate any header after modifying positions of atoms
        subprocess.run(
            f"""
                grep -v ATOM {work_dir}{XYZIN_NAME} |
                grep -v ANISO |
                grep -v HETATM |
                grep -v END |
                grep -v TER |
                grep -v MODEL |
                grep -v ENDMDL > {work_dir}{XYZIN_NAME}.header
            """,
            shell=True
        )

        # ----------------------------------------------------------------

        sigma_level = args.input.difference_data.sigma_level

        map_manipulation = MapManipulation(
            rxn_coord_utils,
            ccp4_utils,
            work_dir,
            keep_files=is_verbose_mode or is_demo_mode or is_probe_mode,
            last_xyzin_name=f'{work_dir}{XYZIN_NAME}',
            dark_xyzin_name=args.input.pdb.dark_file_name,
            dark_data_name=args.input.dark_data.file_name,
            light_data_name=args.input.difference_data.file_name,
            default_sigma_level=sigma_level,
            calculated_sigma_level=sigma_level if is_beta_analysis else sigma_level * n_value,
            calculated_cutoff_density=rms if is_beta_analysis else None,
            xyzregion_file=XYZREGION_FILE,
            sigma_mapin_name=sigma_mapin_name,
            selection=selection,
            selection_radius=args.selection_radius,
            resolution_factor=args.high_resolution / GRID_SLICE,
        )

        if is_demo_mode:
            # in case the code stops early, the log is still partially produced
            # so the log should be defined with empty string first
            log1, log2, log3, log4 = '', '', '', ''
            try:
                print("Step 1/3: Build the PDB file with the modified moiety (including appending the header)")
                # prefix = 'DEMO__' + '__'.join(map(generate_prefix, demo_parameters))

                # modify the structure parameter by parameter (?)
                prefix = f'DEMO_cycle000_initial_struct'
                is_ok,  _, log_struct = map_manipulation.modify_pdb(prefix=prefix)
                for idx, param in enumerate(demo_parameters):
                    prefix = f'DEMO_cycle{idx + 1:>03}_{generate_prefix(param)}'
                    is_ok,  _, log_struct = map_manipulation.modify_pdb(param, prefix=prefix)
                    map_manipulation.last_xyzin_name = f'{work_dir}{prefix}.pdb'

                if len(demo_parameters) == 0:
                    print(f'The initial structure {init_xyzin_name} is used to generate the calculated structure factor')
                else:
                    print(f'The initial structure {init_xyzin_name} is modified with the following parameters before generating the calculated structure factor:')
                    print(*demo_parameters, sep=', ')
                print(f"The structure has the following properties:")
                print(log_struct)
                print(f'and is{"" if is_ok else " not"} OK to be used')

                # probe mode (structure-only mode): early return
                if is_probe_mode:
                    break

                print("Step 2/3: Generate new maps from the modified structure")
                log1, log2, log3 = map_manipulation.make_maps(prefix=prefix)

                # if there is any word 'error' in the log, throw an error
                # however, phenix does output the statistical error value, thus it is ignored in this situation
                if not ccp4_utils.should_use_phenix_for_dFc:
                    for line in (log1+log2+log3).splitlines():
                        if re.search('error', line, flags=re.IGNORECASE):
                            raise RuntimeError(log1+log2+log3)

                print("Step 3/3: Calculate the CC")
                log4, prefix, cc, slope = map_manipulation.calculate_cc(prefix=prefix)

                for line in log4.splitlines():
                    if re.search('error', line, flags=re.IGNORECASE):
                        raise RuntimeError(log4)

                print(f"The structure and maps of '{prefix}' is generated,")
                print(f"with the correlation coefficient of {cc}{',' if is_beta_analysis else '.'}")
                if is_beta_analysis:
                    print(f"and the beta value of {slope}.")

                prepare_linregress_plot(
                    work_dir,
                    n_value,
                    args.input.difference_data.sigma_level,
                    light_map_name, # type: ignore # only unbound when dFo is given as map in beta analysis, which is forbidden
                    prefix,
                    'Linear Regression on dFc-dFo in Demo mode',
                )
                extract_unexplained_and_overmodeled_voxels(work_dir, prefix, sigma_level)
                map_manipulation.generate_coot_script(prefix, rms=rms if is_beta_analysis else None)

            finally:
                with open(f"{work_dir}log_dFoCC.demo.txt", "w") as f:
                    for log in [log0, log1, log2, log3, log4]: # log1
                        if log != '':
                            f.write(log)
                            f.write(u'\u2500' * 30)
                            f.write('\n')

            break

        # ----------------------------------------------------------------

        log_struct_list: list[dict] = []
        log_cc_list: list[dict] = []

        print(f"Step 1/{1 if is_probe_mode else 4}: Generate a library of structures with the modified moiety (including appending the header)")
        with Pool(processes=args.n_jobs) as pool:
            for ok, prefix, log_struct \
            in tqdm(pool.imap_unordered(map_manipulation.modify_pdb, loop_through_params), total=total_iterations):
                log_struct_list.append({ 'prefix': prefix, 'ok': ok } | log_struct)
                if not ok:
                    if is_verbose_mode or is_probe_mode:
                        subprocess.run(['mv', f'{work_dir}{prefix}.pdb', f'{work_dir}{prefix}.bad.pdb'])
                    else:
                        subprocess.run(['rm', f'{work_dir}{prefix}.pdb'])

        # log all structures first
        log_struct_data_frame = pd.DataFrame(log_struct_list).set_index('prefix')
        log_struct_data_frame.to_csv(f'{work_dir}{LOG_STRUCT_NAME}', sep=('\t' if LOG_STRUCT_NAME[-4:] == '.tsv' else ','))

        # probe mode: early return
        if is_probe_mode:
            break

        ok_struct_params = [log_struct['prefix'] for log_struct in log_struct_list if log_struct['ok']]
        ok_total_iterations = len(ok_struct_params)

        print("Step 2/4: Generate new maps from the modified structure")
        with Pool(processes=args.n_jobs) as pool:
            for _ in tqdm(pool.imap_unordered(map_manipulation.make_maps, ok_struct_params), total=ok_total_iterations):
                pass

        print("Step 3/4: Calculate the CC")
        with Pool(processes=args.n_jobs) as pool:
            for _, prefix, cc, slope \
            in tqdm(pool.imap_unordered(map_manipulation.calculate_cc, ok_struct_params), total=ok_total_iterations):
                if is_beta_analysis:
                    log_cc_list.append({ 'prefix': prefix, 'CC': cc, 'beta': slope })
                else:
                    log_cc_list.append({ 'prefix': prefix, 'CC': cc })

        log_cc_data_frame = pd.DataFrame(log_cc_list).set_index('prefix')
        log_struct_data_frame = log_struct_data_frame.join(log_cc_data_frame)
        log_struct_data_frame.to_csv(f'{work_dir}{LOG_STRUCT_NAME}', sep=('\t' if LOG_STRUCT_NAME[-4:] == '.tsv' else ','))

        print("Step 4/4: Find the map with the highest CC")
        # check if at least one CC value is calculated
        if cast(bool, log_struct_data_frame['CC'].isna().all()):
            raise RuntimeError('No legitimate CC value from any structure, possibly due to overestimation of the N value')

        last_prefix = generate_prefix(rxn_coord_utils.get_best_params())
        best_prefix = str(log_struct_data_frame['CC'].idxmax())
        best_log = log_struct_data_frame.loc[best_prefix]

        if is_beta_analysis and (log_struct_data_frame['beta'] >= 1.).any():
            best_CC = best_log['CC']
            distance_from_beta_1 = abs(log_struct_data_frame['beta'] - 1)
            idx_closest_to_beta_1, distance_closest_to_beta_1 = distance_from_beta_1.idxmin(), distance_from_beta_1.min()

            # FIXME: for now, the CC checking is turned off so that the occupancy is reported no matter the CC value
            best_prefix = str(idx_closest_to_beta_1)
            best_log = log_struct_data_frame.loc[best_prefix]
            # if distance_closest_to_beta_1 < 0.05 and best_CC - log_struct_data_frame['CC'][idx_closest_to_beta_1] < 0.05:
            #     best_prefix = str(idx_closest_to_beta_1)
            #     best_log = log_struct_data_frame.loc[best_prefix]

        # best case scenario: if the structure from the last cycle has the highest CC also in this cycle
        if last_prefix == best_prefix:
            is_next_round = True
            reason = 'strict maxima'
        # if both structure from this cycle and the last cycle have the same CC value,
        # use last structure as the best structure of this cycle instead
        elif log_struct_data_frame['CC'].max() == log_struct_data_frame['CC'][last_prefix]:
            best_prefix = last_prefix
            is_next_round = True
            reason = 'plateau'
        else:
            if len(log_cycle_list) == 0:
                is_next_round = False
                reason = ''
            else:
                csv_data_frame_old = pd.DataFrame(log_cycle_list)
                is_duplicated = (csv_data_frame_old['hash'] == best_log['hash']).any()
                is_highest = csv_data_frame_old['CC'].max() <= best_log['CC']
                if is_duplicated and is_highest:
                    is_next_round = True
                    reason = 'duplicated'
                else:
                    is_next_round = False
                    reason = ''

        # generate best log again with the correct prefix
        best_log = log_struct_data_frame.loc[best_prefix]
        # best_cc = float(log_struct_data_frame['CC'][log_struct_data_frame['ok']].max())

        print(f"The best result in this cycle is '{best_prefix}',")
        if is_beta_analysis:
            print(f"with the correlation coefficient of {best_log['CC']}, and the beta value of {best_log['beta']}")
            print(best_log)
        else:
            print(f"with the correlation coefficient of {best_log['CC']}.")
            print(log_struct_data_frame.sort_values(by=['CC'], ascending=False).head(5))

        # ----------------------------------------------------------------

        # record the work directory for this cycle
        last_work_dir = work_dir

        # regenerate the model and other files for the best set of coordinates.

        print("Regenerating best structure in this cycle")
        map_manipulation.keep_files = True
        progress_bar = tqdm(total=100)

        map_manipulation.modify_pdb(split_prefix(best_prefix))
        progress_bar.update(5) # 5%

        map_manipulation.make_maps(best_prefix)
        progress_bar.update(40) # 45%

        map_manipulation.calculate_cc(best_prefix)
        progress_bar.update(30) # 75%

        dFo_unexplained_count, dFc_overmodeled_count, dFo_count, dFc_count, total_count = extract_unexplained_and_overmodeled_voxels(work_dir, best_prefix, sigma_level)
        log_cycle_list.append({
            'round': round_num,
            'cycle': cycle_num,
            'prefix': best_prefix,
        } | best_log.to_dict() | {
            'unexplained': dFo_unexplained_count,
            'overmodeled': dFc_overmodeled_count,
            'dFo_voxels': dFo_count,
            'dFc_voxels': dFc_count,
            'ttl_voxels': total_count,
            'reason': reason,
        })
        pd.DataFrame(log_cycle_list).to_csv(log_cycle_name, sep='\t' if log_cycle_name[-4:] == '.tsv' else ',', index=False)
        progress_bar.update(10) # 85%

        prepare_all_plots(work_dir, rxn_coord_utils, n_value, pd.DataFrame(log_cycle_list), cycle_num)
        prepare_linregress_plot(
            work_dir,
            n_value,
            sigma_level,
            light_map_name, # type: ignore # only unbound when dFo is given as map in beta analysis, which is forbidden
            best_prefix,
            f'Linear Regression on dFc-dFo in Cycle {cycle_num}, Round {round_num}',
        )
        progress_bar.update(15) # 100%

        map_manipulation.generate_coot_script(best_prefix, rms=rms if is_beta_analysis else None)
        progress_bar.close()
        print(f'You may check the structure later with coot by `coot -s {work_dir}coot.scm`')

        # ----------------------------------------------------------------

        # make a file called `INTERRUPT` (`touch INTERRUPT`) to interrupt after finishing a cycle
        # this is the only way (for now) to stop the execution in the middle gracefully
        if isfile(args.input.interrupt_file_name):
            break

        # stop the algorithm when the cycle number reaches a certain value
        # no matter what number is set, the algorithm will still run at least once
        if max_cycle > 0 and cycle_num >= max_cycle:
            print(f'Cycle {cycle_num} has reached the max cycle number {max_cycle}')
            break

        # ----------------------------------------------------------------

        # If the requirement of starting the next round is matched,
        # we can use a smaller interval (half) for the next cycle
        if is_next_round:
            if is_last_round:
                break
            print("The best structure of the last cycle is also the best in this cycle")
            print("Change the interval of parameters for further fine tuning")
            print(f'Reason: {reason}')

            # refresh the intervals
            round_num += 1
            cycle_num += 1

        # Otherwise, use this result as the input and train for one more cycle
        # No need to change the interval, as the best structure is far away...
        else:
            print("Use the best result in this cycle as the input for the next cycle")

            # `cycle` no increment; no need to refresh the interval 
            cycle_num += 1

        print(f'Cycle done: {datetime.now().replace(microsecond=0).isoformat(" ")}')
        print()
        print(u'\u2500' * 30)
        print()

    # DEPRECATED: user should explicitly execute `./coot_open.py` or `coot -s cycle???/coot.scm` for visualized results
    # # for probe mode, there's no coot.scm, so don't show the result in coot
    # if isfile(f'{work_dir}coot.scm'):
    #     print('This is the end of dFoCC')
    #     print('Open coot to confirm...')
    #     subprocess.run(['coot', '-s', f'{work_dir}coot.scm'], stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    print(f'Finished: {datetime.now().replace(microsecond=0).isoformat(" ")}')
    print()
    print(u'\u2500' * 30)
    print()

# ------------------------------

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(prog='iterate_dFoCC')
    # arg_parser.add_argument('-d', '--demo', action='store_true', dest='is_demo_mode')
    # arg_parser.add_argument('mode', choices=['run', 'demo', 'probe'], default='run', const='run', nargs='?')
    # arg_parser.add_argument('-i', '--init', dest='init_xyzin_name', help=f"Use a specific file as the initial light structure instead") # , default=INIT_XYZIN_NAME
    # arg_parser.add_argument('-d', '--demo-dir', dest='demo_dir', default=DEMO_DIR, help="Use a specific name for demo directory")
    # arg_parser.add_argument('-p', '--probe-dir', dest='probe_dir', default=PROBE_DIR, help="Use a specific name for probe directory")
    # arg_parser.add_argument('-n', '--n-value', type=float, default=N_VALUE)
    # arg_parser.add_argument('-b', '--beta', action='store_true', default=IS_BETA_ANALYSIS, dest='is_beta_analysis', help="Calculate beta value of linear regression in addition to correlation coeficient")
    # arg_parser.add_argument('-s', '--structure-only', action='store_true', dest='is_structure_only_mode', help="Only generate the structure(s) and terminate the run")
    # arg_parser.add_argument('-v', '--verbose', action='store_true', default=IS_VERBOSE_MODE, dest='is_verbose_mode', help="Keep all the intermediate files (several GB per cycle)")
    # arg_parser.add_argument('-c', '--max-cycle', type=int, default=MAX_CYCLE, help="The algorithm will stop at specific cycle number; never end eariler if set to -1; no matter what the value is set to, the algorithm will be run at least once")
    arg_parser.add_argument('main', nargs='+')
    sys_args = arg_parser.parse_args()

    args_interpreter = default_phil.command_line_argument_interpreter()
    args_phil = args_interpreter.process_args(args=sys_args.main)

    working_phil = cast(phil.scope, default_phil.fetch(sources=args_phil))

    print('The parameters are set as the following:')
    print('```')
    working_phil.show()
    print('```')

    args: Any = working_phil.extract()

    init_xyzin_name = args.input.pdb.initial_file_name
    init_or_dark_xyzin_name = init_xyzin_name or args.input.pdb.dark_file_name
    if not init_or_dark_xyzin_name or not isfile(init_or_dark_xyzin_name):
        raise RuntimeError(f'{init_or_dark_xyzin_name} does not exist')

    previous_best_struct = init_or_dark_xyzin_name

    all_selection_in_dFoCC = ''

    col_names: dict[str, str | None] = {
        # the column names of the `DARK_DATA_NAME` mtz files
        'Fo_dark': args.input.dark_data.labels[0],
        'Fo_dark_sigma': args.input.dark_data.labels[1],
        'Fc_dark': args.input.dark_data.labels[2],
        'Fc_dark_phase': args.input.dark_data.labels[3],
        'free': args.input.dark_data.r_free_flags_label,
        # the column names of the option `EXTRAPOLATED_DATA_NAME` mtz files
        'Fext': args.input.extrapolated_data.labels[0] if args.input.extrapolated_data.labels else None,
        'Fext_sigma': args.input.extrapolated_data.labels[1] if args.input.extrapolated_data.labels else None,
        'Fext_free': args.input.extrapolated_data.r_free_flags_label,
    }

    if args.mode == 'demo':
        if isdir(DEMO_DIR):
            rename(DEMO_DIR, f"{DEMO_DIR[:-1]}{ARCHIVE_ROUND_SUFFIX}")
        makedirs(DEMO_DIR)
        if init_xyzin_name and isfile(init_xyzin_name):
            print('<demo>')
            print('Generate dFc map for the initial structure')
            main(
                ReactionCoordinateSeries(),
                init_xyzin_name=init_xyzin_name,
                demo_dir=f'{DEMO_DIR}no_mod/',
                probe_dir='',
                log_cycle_name='',
                generate_work_dir_name=lambda round, cycle: '',
                n_value=1.00, # to show all the blobs clearly
                selection='all',
                max_cycle=1,
                col_names=col_names,
                demo_parameters=[],
                is_beta_analysis=IS_BETA_ANALYSIS, # FIX
                is_demo_mode=True,
                is_probe_mode=False,
                is_verbose_mode=True,
            )
    elif args.mode == 'probe':
        makedirs(PROBE_DIR)

    for idx, refining_residue in enumerate(args.residue):
        if refining_residue.output_prefix:
            name = refining_residue.output_prefix
        else:
            name = f'resi{refining_residue.residue_id}' + (f'_alt{refining_residue.alt_conf_id}' if refining_residue.alt_conf_id else '')

        selection = prepare_extraction_selection(refining_residue)
        if all_selection_in_dFoCC:
            all_selection_in_dFoCC = all_selection_in_dFoCC + ' or ' + selection
        else:
            all_selection_in_dFoCC = selection

        rxn_coord_utils = prepare_reaction_coordinate_utils(refining_residue, init_or_dark_xyzin_name)

        def generate_work_dir_name(_, cycle_num: int):
            return f"{idx}_{name}_cycle*/" if cycle_num == '*' else f'{name}_cycle{cycle_num:>03}/'

        print(f'<{name}>')

        # demo mode:
        demo_parameters = []
        log_cycle_name = f'{name}_log_cycle.csv'
        if log_cycle_name and isfile(log_cycle_name):
            log_cycle_df = pd.read_csv(log_cycle_name).set_index('cycle').drop(index=0, errors='ignore')
            demo_parameters: list[Parameters] = log_cycle_df['prefix'].apply(split_prefix).to_list()
        elif args.mode == 'demo':
            print(f'Skipping {name} due to the nonexistence of the demo parameter file {log_cycle_name}')
            print()
            print(u'\u2500' * 30)
            print()
            continue

        try:
            main(
                rxn_coord_utils,
                init_xyzin_name=init_or_dark_xyzin_name,
                demo_dir=f'{DEMO_DIR}{idx}_{name}/',
                probe_dir=f'{PROBE_DIR}{idx}_{name}/',
                log_cycle_name=log_cycle_name,
                generate_work_dir_name=generate_work_dir_name,
                n_value=args.n_value,
                selection=selection,
                max_cycle=args.max_cycle,
                col_names=col_names,
                demo_parameters=demo_parameters,
                is_beta_analysis=IS_BETA_ANALYSIS, # FIX
                # is_demo_mode=args.mode == 'demo',
                # is_probe_mode=args.mode == 'probe' or args.is_structure_only_mode,
                # is_verbose_mode=args.is_verbose_mode or args.mode in ('demo', 'probe'),
                is_demo_mode=args.mode == 'demo',
                is_probe_mode=args.mode == 'probe',
                is_verbose_mode=args.verbose_mode or args.mode in ('demo', 'probe'),
            )
        except AssertionError as e:
            print(e)
            print()
            print(u'\u2500' * 30)
            print()

        if args.mode == 'run' and isfile(args.input.interrupt_file_name):
            print(f'Interrupted by the file {args.input.interrupt_file_name}')
            break

        # combining the structures
        if args.mode != 'probe':
            log_cycle_df = pd.read_csv(log_cycle_name).set_index('cycle').drop(index=0, errors='ignore')
            demo_parameters: list[Parameters] = log_cycle_df['prefix'].apply(split_prefix).to_list()
            last_cycle = log_cycle_df.iloc[-1]

            demo_dir = DEMO_DIR if args.mode == 'demo' else ''
            demo_dir += f'{idx}_{name}__combined/'

            if idx != 0:
                print('Combining the results from previous residues')

                main(
                    rxn_coord_utils,
                    init_xyzin_name=previous_best_struct,
                    demo_dir=demo_dir,
                    probe_dir='',
                    log_cycle_name=log_cycle_name,
                    generate_work_dir_name=generate_work_dir_name,
                    n_value=args.n_value,
                    selection=selection,
                    max_cycle=1,
                    col_names=col_names,
                    demo_parameters=demo_parameters,
                    is_beta_analysis=IS_BETA_ANALYSIS, # FIX
                    is_demo_mode=True,
                    is_probe_mode=True, # structure-only mode
                    is_verbose_mode=True,
                )

                previous_best_struct = f'{demo_dir}DEMO_cycle{last_cycle.name:>03}_{last_cycle["prefix"]}.pdb'

            elif args.mode == 'run':
                previous_best_struct = f'{generate_work_dir_name(last_cycle["round"], last_cycle.name)}{last_cycle["prefix"]}.pdb'
            else:
                previous_best_struct = f'{DEMO_DIR}{idx}_{name}/DEMO_cycle{last_cycle.name:>03}_{last_cycle["prefix"]}.pdb'

        if args.mode == 'run' and isfile(args.input.interrupt_file_name):
            print(f'Interrupted by the file {args.input.interrupt_file_name}')
            break

    if args.mode == 'run' and isfile(args.input.interrupt_file_name):
        print(f'Interrupted by the file {args.input.interrupt_file_name}')

    elif args.mode != 'probe':
        work_dir = DEMO_DIR if args.mode == 'demo' else ''
        work_dir += f'{len(args.residue)}_phenix_refine/'
        if isdir(work_dir):
            rename(work_dir, f"{work_dir[:-1]}{ARCHIVE_ROUND_SUFFIX}")

        extrapolated_data_name = args.input.extrapolated_data.file_name

        generate_phenix_refinement_script(
            xyzin_name=previous_best_struct,
            mtzin_name=extrapolated_data_name or '!!!MODIFY HERE!!!',
            output_dir=work_dir,
            script_name=f'{work_dir}refine_script.phil',
            selection_for_conventional_refinement=f'not ({all_selection_in_dFoCC})',
            col_names=col_names,
            working_params=args,
        )

        if args.mode == 'run' and extrapolated_data_name and isfile(extrapolated_data_name):
            print('Start conventional refinement with `phenix.refine`')
            print(f'Atom selection: not ({all_selection_in_dFoCC})')
            phenix_refine(f'{work_dir}refine_script.phil')
            print()
            print(u'\u2500' * 30)
            print()

            print('Generate dFc map for the refined structure')
            dFc_dir = f'{len(args.residue)}_phenix_refine_dFc/'
            main(
                ReactionCoordinateSeries(),
                init_xyzin_name=f'{work_dir}refined_001.pdb',
                demo_dir=dFc_dir,
                probe_dir='',
                log_cycle_name='',
                generate_work_dir_name=lambda round, cycle: '',
                n_value=1.00, # to show all the blobs clearly
                selection='all',
                max_cycle=1,
                col_names=col_names,
                demo_parameters=[],
                is_beta_analysis=IS_BETA_ANALYSIS, # FIX
                is_demo_mode=True,
                is_probe_mode=False,
                is_verbose_mode=True,
            )

            with open(f'{dFc_dir}coot.scm', 'a') as coot_script_file:
                coot_script_file.write(f'''
;; read the extrapolated struct factor
;; (auto-read-make-and-draw-maps "{extrapolated_data_name}")
(make-and-draw-map "{work_dir}refined__001.mtz" "2FOFCWT_no_fill" "PH2FOFCWT_no_fill" "" 0 0)
''')

        else:
            print('If you would like to start the conventional refinement')
            print(f'modify the config file `{work_dir}refine_script.phil` first')
            print(f'and then execute `phenix.refine {work_dir}refine_script.phil`')
            print(f'The final dFc of the refined structure can be generate by running `python ./iterate_dFoCC.py ./configuration.phil mode=demo initial_file_name={work_dir}refined_001.pdb`')

        print()
        print('Refinement terminates')
        print()
        print(u'\u2500' * 30)
        print()
