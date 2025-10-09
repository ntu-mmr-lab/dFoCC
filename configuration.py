#!/usr/bin/env python

from os import makedirs, rename
from os.path import isfile, isdir, dirname
import sys
import re
import argparse
from typing import Any, Literal, NotRequired, TypedDict, cast

# for archiving old directory
from time import strftime

from Bio.PDB.PDBParser import PDBParser
pdb_parser = PDBParser()

from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Bio.PDB.vectors import Vector

from libtbx import phil

import pandas as pd

from dFoCC.reaction_coordinates import ParamTuple as Parameters
from dFoCC.reaction_coordinates import ReactionCoordinateSeries

from dFoCC.reaction_coordinates.aa_dihedrals import get_dihedrals, rotate_dihedrals

# These reaction coordinates are hard coded for each moiety
# Use them with care
from dFoCC.reaction_coordinates.rotate_main_chain_general import rotate_main_chain_general, calculate_main_chain_general

default_phil = phil.parse("""
    residue
        .multiple = True
    {
        output_prefix = None
            .type = path
        chain_id = None
            .type = str
            .optional = False
        residue_id = None
            .type = int
            .optional = False
        alt_conf_id = None
            .type = str
            .optional = False
        main_chain_rotation
        {
            start.residue_id = None
                .type = int
            start.atom_id = N CA C
                .type = choice
            end.residue_id = None
                .type = int
            end.atom_id = N CA C
                .type = choice
        }
        initial_parameters_step = None
            .type = ints(size_min=1, size_max=5, value_min=1)
            .optional = False
        step_per_run = None
            .type = ints(size_min=1, size_max=5, value_min=1)
        selection = None
            .type = str
    }
    selection_radius = 1.56
        .type = float
    n_jobs = 1
        .type = int(value_min=1)
    max_cycle = -1
        .type = int
    verbose_mode = False
        .type = bool
    n_value = None
        .type = float(value_min=0.0)
        .optional = False
    high_resolution = None
        .type = float(value_min=0.0)
        .optional = False
    low_resolution = None
        .type = float(value_min=0.0)
        .optional = False
    space_group = None
        .type = str
        .optional = False
    input
        .optional = False
    {
        pdb {
            dark_file_name = None
                .type = path
                .optional = False
            initial_file_name = None
                .type = path
        }
        dark_data {
            file_name = None
                .type = path
                .optional = False
            labels = None
                .type = strings
                .optional = False
            r_free_flags_label = None
                .type = str
                .optional = False
        }
        difference_data {
            file_name = None
                .type = path
                .optional = False
            sigma_level = 3.0
                .type = float
        }
        extrapolated_data
        {
            file_name = None
                .type = path
            labels = None
                .type = strings
            r_free_flags_label = None
                .type = str
                .optional = False
        }
        interrupt_file_name = ./INTERRUPT
            .type = path
    }
    mode = *run demo probe
        .type = choice
        .optional = False
""")

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(prog='dry_run_dFoCC_config')
    arg_parser.add_argument('main', nargs='+')
    sys_args = arg_parser.parse_args()

    args_interpreter = default_phil.command_line_argument_interpreter()
    args_phil = args_interpreter.process_args(args=sys_args.main)

    args = cast(phil.scope, default_phil.fetch(sources=args_phil))

    print('The parameters are set as the following:')
    print('```')
    args.show()
    print('```')
else:
    args: Any = default_phil.extract()

# ----------------------------------------

# ## Configurations

# # how many processor can be used in the same time
# # notice that keeping at least 4 of the processor / threads is recommended
# # to keep the interaction with the system and other visualizing tools simultaneously
# NUMBER_OF_JOBS = working_params.n_jobs
#
# # `N_VALUE` is 200 divided by the assumed occupancy, and is used to scale the dFc map
# N_VALUE = working_params.n_value
#
# LOW_RESOLUTION = working_params.low_resolution
# HIGH_RESOLUTION = working_params.high_resolution
# SYMMETRY = working_params.space_group

# class RefiningResidue(TypedDict):
#     # name: NotRequired[str]
#     name: str | None
#     chain_id: str
#     residue_id: int
#     # alt_conf_id: NotRequired[str]
#     alt_conf_id: str | None
#     # main_chain_rotation: NotRequired[tuple[tuple[int, Literal['N', 'CA', 'C']], tuple[int, Literal['N', 'CA', 'C']]]]
#     main_chain_rotation: tuple[
#         tuple[int | None, Literal['N', 'CA', 'C'] | None],
#         tuple[int | None, Literal['N', 'CA', 'C'] | None],
#     ]
#     # selection: NotRequired[str]
#     selection: str | None
#     initial_parameters_step: Parameters
#
# # list of residue to be refined
# # `chain_id`, `residue_id`, `n_value` should be set
# # `name` is optional and used as the prefix of the directory name
# # `alt_conf_id` should be set if and only if the initial structure has alternative conformations
# REFINING_RESIDUE_LIST: list[RefiningResidue] = [
#     {
#         'name': res.output_prefix,
#         'chain_id': res.chain_id,
#         'residue_id': res.residue_id,
#         'alt_conf_id': res.alt_conf_id,
#         'initial_parameters_step': tuple(res.initial_parameters_step),
#         'main_chain_rotation': (
#             (res.main_chain_rotation.start.residue_id, res.main_chain_rotation.start.atom_id),
#             (res.main_chain_rotation.end.residue_id, res.main_chain_rotation.end.atom_id),
#         ),
#         'selection': res.selection
#     }
#     for res in working_params.residue
# ]
# # REFINING_RESIDUE_LIST: list[RefiningResidue] = [
# #     # { 'name': 'K258', 'chain_id': 'A', 'residue_id': 258, 'initial_parameters_step': (720, 720, 720, 720) },
# #     { 'name': 'H309', 'chain_id': 'A', 'residue_id': 309, 'initial_parameters_step': (720, 720) },
# #     {
# #         'name': 'D321',
# #         'chain_id': 'A',
# #         'residue_id': 321,
# #         'main_chain_rotation': ((320, 'C'), (322, 'CA')),
# #         'initial_parameters_step': (720, 720, 64)
# #     },
# #     {
# #         # 'name': 'E384_altC',
# #         'name': 'E384',
# #         'chain_id': 'A',
# #         'residue_id': 384,
# #         # 'alt_conf_id': 'C',
# #         'initial_parameters_step': (720, 720, 720)
# #     },
# #     { 'name': 'N395', 'chain_id': 'A', 'residue_id': 395, 'initial_parameters_step': (720, 720) },
# # ]

# CHAIN_ID = 'A'
# # it is ok to provide just the residue number for normal residue,
# # otherwise other ligands or moiety should be provided in this format: `'H_HEM', 201, ' '`
# RESIDUE_ID: int | tuple[str, int, str] = 'H_60F', 69, ' '
# # change it to `None` if there is only one conformation in the refining region
# ALT_CONF_ID: str | None = None
# # change it to 'CB' for normal residue, only used when opening coot
# CENTER_ATOM = 'CB'

# how many slice should be created in the length of the high resolution
GRID_SLICE = 3
# in most of the situation, it is not necessary to modify `GRID_CCP4`
# if you encounter an error during Step 3 / CCP4 overlapmap
# you may try to provide this in the format of `200 100 30` which is shown in `mapdump` of the dFo map
GRID_CCP4 = f'SAMPLE {GRID_SLICE}'
# # the sigma level cutoff that for the dFo and dFc map, to create a clear boundary of them
# SIGMA_LEVEL = working_params.difference_data.sigma_level
#
# # select the region to be compared
# # within `SELECTION_RADIUS` Å of the selected reference atom
# SELECTION_RADIUS = working_params.selection_radius
# # # select the reference atom
# # # the syntax can be found [here](https://phenix-online.org/documentation/reference/atom_selections.html)
# # SELECTION = 'resseq 69 and (name SG or name C1 or name O1 or name C2 or name C3)'

# ----------------------------------------

# ## Filename Definition

# ### the following variables are for the directory names created during the process (must append '/')
PREPARING_DIR = 'prepare/'
# PREPARING_DIR = f'preparing_sigma{SIGMA_LEVEL}_gridslice{GRID_SLICE}/'
DEMO_DIR = 'demo/'
PROBE_DIR = 'probe/'
ARCHIVE_EXE_DIR = f"data_bu_from_{strftime('%m%d_%H%M%S')}/"
ARCHIVE_ROUND_SUFFIX = f"_bu_at_{strftime('%m%d_%H%M%S')}/"

# # ### the following are for the input files (full relative / absolute path)
#
# # 1 REQUIRED input file
# # for region selection and the initial structure
# DARK_XYZIN_NAME = working_params.input.pdb.dark_file_name
#
# # 2 REQUIRED input files:
# # Fo(dark) and Fc(dark)
# DARK_DATA_NAME = working_params.input.dark_data.file_name
# # dFo structure factor file, run `phenix.fobs_minus_fobs_map` to create that file
# # dFo map (DED map) can also be used directly (absent dFo reflections are not handled)
# LIGHT_DATA_NAME = working_params.input.difference_data.file_name

# 3 optional input files:
# the dFo map file generated from the input dFo mtz
# this should only be set if the sigma cutoff dFo map file is provided as input
# otherwise set to `None` instead
SIGMA_MAPIN_NAME: str | None = None
# # the pdb file provided as the initial structure to be modified
# # by default (set as `None`), it will use the `DARK_XYZIN_NAME` file instead
# INIT_XYZIN_NAME: str | None = working_params.input.pdb.initial_file_name
# # the extrapolated structure factor file for the conventional refinement
# # of the residues that are not refined by dFoCC
# # if not specified (set as `None`), the conventional refinement step is skipped
# EXTRAPOLATED_DATA_NAME: str | None = working_params.input.extrapolated_data.file_name

# col_names: dict[str, str] = {
#     # the column names of the `DARK_DATA_NAME` mtz files
#     'Fo_dark': working_params.input.dark_data.labels[0],
#     'Fo_dark_sigma': working_params.input.dark_data.labels[1],
#     'Fc_dark': working_params.input.dark_data.labels[2],
#     'Fc_dark_phase': working_params.input.dark_data.labels[3],
#     # the column names of the option `EXTRAPOLATED_DATA_NAME` mtz files
#     'Fext': working_params.input.extrapolated_data.labels[0],
#     'Fext_sigma': working_params.input.extrapolated_data.labels[1],
#     'free': working_params.input.dark_data.r_free_flags_label,
# }

# ### the following are for the intermediate and output files (only filename)

# the name of the best structure of each round
XYZIN_NAME = 'previous_best_struct.pdb'

LOG_STRUCT_NAME = 'log_struct.csv' # .tsv

# only the map of the moiety will be generated
# to be more clear, name the related files with the name of the moiety
# the '*' will be substituted into `prefix` later on
XYZREGION_FILE = '*_combine_dark.pdb'

# # keep the intermediate files if necessary (several GB per cycle)
# # intermediate files in demo mode are kept regradless of this setting
# IS_VERBOSE_MODE = working_params.verbose_mode
#
# # you may set a max cycle number, so that the algorithm will be stopped earlier
# # set to -1 if you don't need this feature
# MAX_CYCLE: int = working_params.max_cycle
#
# # apart from setting a max cycle number,
# # you may make a file called `INTERRUPT` (`touch INTERRUPT`)
# # to terminate after finishing a cycle during the execution
# INTERRUPT_FILE_NAME = working_params.input.interrupt_file_name

# TODO: calculate beta value of linear regression instead of correlation coeficient
IS_BETA_ANALYSIS = False

# use phenix.fobs_minus_fobs_map instead of sftools for generating dFc
# it should not be used unless sftools is not applicable (e.g. symmetry issue)
# as the absent reflection in dFo will not be ignored in the dFc
SHOULD_USE_PHENIX_FOR_DFC = False

# plot anything related to CC with the following range of CC value
EXPECTED_CC_RANGE = (0.00, 1.00)

# ----------------------------------------

# ## Filename Utilities
def generate_prefix(param: Parameters):
    """
    The calculated file of each modified structure will be infixed with something like "_coord_16_-16_16_-16_16_-16_16"
    """
    return "rc_" + "_".join(str(i) for i in param)

def split_prefix(prefix: str) -> Parameters:
    return tuple(int(i) for i in re.findall(r'(-?\d+)', prefix))

# ----------------------------------------

# # the script that will remove all the hydrogens in the pdb
# class RemoveHydrogen(Select):
#     """
#     Remove hydrogen of the modified moiety
#     Just in case?
#     """
#     def accept_atom(self, atom: Atom):
#         return atom.get_full_id()[3] != RESIDUE_ID \
#             or atom.get_id()[0] != 'H'

# ----------------------------------------

# prepare the reaction coordinates and extraction selection
def prepare_extraction_selection(refining_residue: Any) -> str:
    alt_conf_id = refining_residue.alt_conf_id

    if refining_residue.selection:
        selection: str = refining_residue.selection
    else:
        selection = f'(resseq {refining_residue.residue_id} and not (name N or name CA or name C or name O))'
        # if refining_residue['main_chain_rotation']:
        #     res_from, res_to = refining_residue['main_chain_rotation']
        #     match res_from[1]:
        #         case 'N':
        #             selection += f' or (resseq {res_from[0]} and not name N)'
        #         case 'CA':
        #             selection += f' or (resseq {res_from[0]} and (name C or name O))'
        #         case 'C':
        #             pass
        #         case _:
        #             assert res_from[1] in ('N', 'CA', 'C'), 'Atom is not in the main chain or is invalid'
        #     match res_to[0] - res_from[0]:
        #         case 0:
        #             raise IndexError(f'Single residue main chain rotation (from {res_from} to {res_to}) is not supported')
        #         case 1:
        #             pass
        #         case 2:
        #             selection += f' or (resseq {res_from[0]+1})'
        #         case _:
        #             selection += f' or (resseq {res_from[0]+1}:{res_to[0]-1})'
        #     match res_to[1]:
        #         case 'N':
        #             pass
        #         case 'CA':
        #             selection += f' or (resseq {res_to[0]} and name N)'
        #         case 'C':
        #             selection += f' or (resseq {res_from[0]} and not (name C or name O))'
        #         case _:
        #             assert res_from[1] in ('N', 'CA', 'C'), 'Atom is not in the main chain or is invalid'
        selection = f'chain {refining_residue.chain_id} and ({selection})'
        if alt_conf_id:
            selection += f' and altloc {alt_conf_id}'

    return selection

class AminoAcidGeneralReactionCoordinates:
    def __init__(self, refining_residue, alt_conf_id: str, has_main_chain_rotation: bool):
        self.rxn_coord_utils = ReactionCoordinateSeries()
        self.refining_residue = refining_residue
        self.has_main_chain_rotation = has_main_chain_rotation
        self.alt_conf_id = alt_conf_id
        self.rxn_coord_utils.use_rc(
            rc_mod_fn=self.modify_amino_acid,
            init_step_size=tuple(refining_residue.initial_parameters_step),
            step_per_run=refining_residue.step_per_run,
        )
        self.rxn_coord_utils.append_properties(self.calculate_amino_acid)

    def modify_amino_acid(self, struct, params):
        chain: Chain = struct[0][self.refining_residue.chain_id]
        residue: Residue = cast(Residue, chain[self.refining_residue.residue_id])
        if self.has_main_chain_rotation:
            rotate_dihedrals(residue, *params[:-1], alt_conf_id=self.alt_conf_id)
            rotate_main_chain_general(
                chain,
                (self.refining_residue.main_chain_rotation.start.residue_id, self.refining_residue.main_chain_rotation.start.atom_id),
                (self.refining_residue.main_chain_rotation.end.residue_id, self.refining_residue.main_chain_rotation.end.atom_id),
                params[-1],
                alt_conf_id=self.alt_conf_id,
            )
        else:
            rotate_dihedrals(residue, *params, alt_conf_id=self.alt_conf_id)

    def calculate_amino_acid(self, struct, params):
        chain: Chain = struct[0][self.refining_residue.chain_id]
        residue: Residue = cast(Residue, chain[self.refining_residue.residue_id])
        prop_obj = {
            f'chi_{idx+1}({self.refining_residue.residue_id})': chi for idx, chi in enumerate(get_dihedrals(residue, self.alt_conf_id)[0])
        }
        if self.has_main_chain_rotation:
            dihedral_main_chain_name, dihedral_main_chain = calculate_main_chain_general(
                chain,
                (self.refining_residue.main_chain_rotation.start.residue_id, self.refining_residue.main_chain_rotation.start.atom_id),
                alt_conf_id=self.alt_conf_id,
            )
            prop_obj.update({ dihedral_main_chain_name: dihedral_main_chain })
        return True, prop_obj

def prepare_reaction_coordinate_utils(refining_residue: Any, init_or_dark_xyzin_name: str) -> ReactionCoordinateSeries:
    """Automatically generate the reaction coordinates from the given residue id"""

    alt_conf_id = refining_residue.alt_conf_id

    has_main_chain_rotation = refining_residue.main_chain_rotation.start.atom_id != None

    init_struct = cast(Structure, PDBParser().get_structure('', init_or_dark_xyzin_name))
    residue = cast(Residue, init_struct[0][refining_residue.chain_id][refining_residue.residue_id])
    num_of_chis = len(get_dihedrals(residue, alt_conf_id)[0])

    if refining_residue.initial_parameters_step == None:
        if refining_residue.step_per_run != None:
            raise Exception('`step_per_run` is given but `initial_parameters_step` is not given')
        match num_of_chis, has_main_chain_rotation:
            case 1, False:
                refining_residue.initial_parameters_step = (1, )
                refining_residue.step_per_run = (3600, )
            case 1, True:
                refining_residue.initial_parameters_step = (16, 10) # (16, 32)
                refining_residue.step_per_run = (225, 15) # (225, 5)
            case 2, False:
                refining_residue.initial_parameters_step = (60, 60)
                refining_residue.step_per_run = (61, 61) # (61, 61)
            case 2, True:
                refining_residue.initial_parameters_step = (144, 144, 32)
                refining_residue.step_per_run = (25, 25, 5)
            case 3, False:
                refining_residue.initial_parameters_step = (240, 240, 240)
                refining_residue.step_per_run = (15, 15, 15)
            case 3, True:
                refining_residue.initial_parameters_step = (400, 400, 400, 32)
                refining_residue.step_per_run = (9, 9, 9, 5)
            case 4, False:
                refining_residue.initial_parameters_step = (600, 600, 600, 600)
                refining_residue.step_per_run = (7, 7, 7, 7)
            case 4, True:
                refining_residue.initial_parameters_step = (720, 720, 720, 720, 32)
                refining_residue.step_per_run = (5, 5, 5, 5, 5)
            case 0, _:
                raise Exception(f'Unexpected residue: no regular chi angle found in the side chain')
            case _, _:
                raise Exception(f'Unexpected residue: {num_of_chis} regular chi angles found in the side chain')

    amino_acid_general_reaction_coordinates = AminoAcidGeneralReactionCoordinates(refining_residue, alt_conf_id, has_main_chain_rotation)

    return amino_acid_general_reaction_coordinates.rxn_coord_utils
