#!/bin/env python

import typing
import itertools

import numpy as np
import pandas as pd
from scipy.stats import linregress
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import gemmi

from configuration import *
from dFoCC.reaction_coordinates import ReactionCoordinateSeries

min_cc, max_cc = EXPECTED_CC_RANGE

def prepare_linregress_plot(work_dir: str, n_value: float | None, sigma_level: float, dFo_filename: str, best_prefix: str, title: str):
    # Figure 3: Correlation coeficient plot

    fig, ax = plt.subplots(num=5, clear=True)

    # ax.set_title(title)

    # dFo_filename = f'{PREPARING_DIR}{LIGHT_DATA_NAME.rsplit("/", 1)[-1].replace(".mtz", ".map")}' if LIGHT_DATA_NAME.endswith('.mtz') else LIGHT_DATA_NAME
    dFo_all = gemmi.read_ccp4_map(dFo_filename)
    dFo_sigma = dFo_all.grid.array.std() * sigma_level
    if IS_BETA_ANALYSIS:
        dFc_sigma = dFo_sigma
    else:
        if not n_value:
            raise TypeError()
        dFc_all = gemmi.read_ccp4_map(f'{work_dir}dFc.{best_prefix}.map')
        dFc_sigma = dFc_all.grid.array.std() * sigma_level * n_value

    dFo_map = gemmi.read_ccp4_map(f'{work_dir}dFo.{best_prefix}.masked_sigma{sigma_level:.1f}.ccp4')
    dFc_map = gemmi.read_ccp4_map(f'{work_dir}dFc.{best_prefix}.masked_sigma{sigma_level:.1f}.ccp4')

    dFo_axis = dFo_map.grid.array.flatten()
    dFc_axis = dFc_map.grid.array.flatten()

    ax.plot(dFo_axis, dFc_axis, 'o', color='steelblue', markersize=1)

    ax.fill_betweenx([-dFc_sigma, dFc_sigma], x1=dFo_sigma, x2=dFo_axis.max(), color='teal', alpha=0.3, edgecolor=None)
    ax.fill_between([-dFo_sigma, dFo_sigma], y1=dFc_sigma, y2=dFc_axis.max(), color='blue', alpha=0.3, edgecolor=None)
    ax.fill_betweenx([-dFc_sigma, dFc_sigma], x1=dFo_axis.min(), x2=-dFo_sigma, color='magenta', alpha=0.3, edgecolor=None)
    ax.fill_between([-dFo_sigma, dFo_sigma], y1=dFc_axis.min(), y2=-dFc_sigma, color='yellow', alpha=0.3, edgecolor=None)

    lr_res = linregress(dFo_axis, dFc_axis)
    label = f'dFc = {lr_res.slope:.4f} dFo {"+" if lr_res.intercept >= 0 else "-"} {np.abs(lr_res.intercept):.4f}, CC = {lr_res.rvalue:.4f}'
    ax.plot(dFo_axis, lr_res.slope * dFo_axis + lr_res.intercept, label=label, color='darkorange', linewidth=1)

    ax.legend()
    ax.set_xlabel(r'$\mathrm{DED_o (e\cdot Å^{-3})}$')
    ax.set_ylabel(r'$\mathrm{DED_c (e\cdot Å^{-3})}$')

    fig.savefig(f'{work_dir}linregress_cc.svg')

    fig, ax = plt.subplots(num=6, clear=True)

    dFo_axis_no_cross_pattern = dFo_axis[~((dFo_axis == 0) ^ (dFc_axis == 0))] # remove cross pattern but not origin
    dFc_axis_no_cross_pattern = dFc_axis[~((dFo_axis == 0) ^ (dFc_axis == 0))] # remove cross pattern but not origin

    ax.plot(dFo_axis_no_cross_pattern, dFc_axis_no_cross_pattern, 'o', color='steelblue', markersize=1)

    lr_res = linregress(dFo_axis_no_cross_pattern, dFc_axis_no_cross_pattern)
    label = f'dFc = {lr_res.slope:.4f} dFo {"+" if lr_res.intercept >= 0 else "-"} {np.abs(lr_res.intercept):.4f}, CC = {lr_res.rvalue:.4f}'
    ax.plot(dFo_axis_no_cross_pattern, lr_res.slope * dFo_axis_no_cross_pattern + lr_res.intercept, label=label, color='darkorange', linewidth=1)

    ax.legend()
    ax.set_xlabel(r'$\mathrm{DED_o (e\cdot Å^{-3})}$')
    ax.set_ylabel(r'$\mathrm{DED_c (e\cdot Å^{-3})}$')

    fig.savefig(f'{work_dir}linregress_no_cross_pattern_cc.svg')

def prepare_all_plots(work_dir: str, rxn_coord_utils: ReactionCoordinateSeries, n_value: float, log_cycle: pd.DataFrame, cycle_num: int | None = None):
    try:
        init_param_intervals = pd.Series(INIT_PARAM_INTERVALS) # type: ignore
    except:
        init_param_intervals = pd.Series(init_step_size for rc in rxn_coord_utils.rc_list for init_step_size in rc.init_step_size) # type: ignore

    log_cycle_slice = log_cycle[log_cycle['cycle'] <= cycle_num] if cycle_num else log_cycle
    content_cycle = log_cycle_slice.iloc[-1]
    # work_dir = generate_work_dir_name(content_cycle['round'], content_cycle['cycle'])
    log_struct_cycle = pd.read_csv(work_dir + LOG_STRUCT_NAME).set_index('prefix')

    # Figure 0: RC accumulation

    fig, ax = plt.subplots(num=0, clear=True)

    # ax.set_title('Accumulation of Reaction Coordinates')

    params_cycle_frame = pd.DataFrame(split_prefix(prefix) for prefix in log_cycle_slice['prefix'].tolist())
    params_acc_frame = params_cycle_frame.cumsum()
    params_acc_norm_frame = params_acc_frame / init_param_intervals

    params_acc_norm_frame.plot(ax=ax, marker='o')
    ax.set_xlabel('cycle #')
    ax.set_ylabel('relative reaction coordinate')
    ax.legend()

    fig.savefig(f'{work_dir}accumulation_rc.png')

    # Figure 1:

    fig, ax = plt.subplots(num=1, clear=True)

    # ax.set_title(f'Histogram of CC Values in Cycle {content_cycle["cycle"]}')
    # ax.set_xlim((log_struct['CC'].min() - .1 * log_cycle['CC'].std(), log_cycle['CC'].max() + .1 * log_cycle['CC'].std()))
    ax.set_xlim(min_cc, max_cc)
    ax.hist(log_struct_cycle['CC'], bins=32)
    ax.set_xlabel('CC')
    ax.set_ylabel('# of structures')

    fig.savefig(f'{work_dir}histogram_cc.png')

    # Figure 2: CC-RC based on the closest RC
    best_prefix = typing.cast(str, content_cycle['prefix']) # log_struct['CC'].idxmax()
    best_log = log_struct_cycle[log_struct_cycle == best_prefix]
    best_params = pd.Series(split_prefix(best_prefix))

    params_frame = pd.DataFrame(split_prefix(prefix) for prefix in log_struct_cycle.index.tolist())

    for rc_index_pair in itertools.combinations(range(len(init_param_intervals)), 2):

        fig, ax = plt.subplots(num=2, clear=True)

        # ax.set_title(f'Heatmap of CC-RC ({rc_index_pair[0]}, {rc_index_pair[1]}) in Cycle {content_cycle["cycle"]}')

        selected_rc_a = params_frame[rc_index_pair[0]].sort_values().unique()
        selected_rc_b = params_frame[rc_index_pair[1]].sort_values().unique()

        selected_heatmap_mat = np.full((len(selected_rc_b), len(selected_rc_a)), np.nan)

        for i in range(len(selected_rc_b)):
            for j in range(len(selected_rc_a)):
                params = best_params.copy()
                params[rc_index_pair[0]] = selected_rc_a[j]
                params[rc_index_pair[1]] = selected_rc_b[i]
                prefix = generate_prefix(tuple(params))
                selected_heatmap_mat[i][j] = log_struct_cycle.loc[prefix, 'CC']

        ax_mappable = ax.imshow(selected_heatmap_mat, vmin=min_cc, vmax=max_cc) # , cmap='cool'
        # ax_mappable = ax.imshow(selected_heatmap_mat, norm=colors.PowerNorm(gamma=5, vmin=min_cc, vmax=max_cc)) # , cmap='cool'

        # ax.set_xticks(range(len(selected_rc_a)), selected_rc_a)
        ax.set_xticks(range(0, len(selected_rc_a), 10), selected_rc_a[::10])
        # ax.set_yticks(range(len(selected_rc_b)), selected_rc_b)
        ax.set_yticks(range(0, len(selected_rc_b), 10), selected_rc_b[::10])

        ax.invert_yaxis()
        ax.set_xlabel(f'RC {rc_index_pair[0]}')
        ax.set_ylabel(f'RC {rc_index_pair[1]}')
        fig.colorbar(ax_mappable)

        fig.savefig(f'{work_dir}heatmap_cc_rc{rc_index_pair[0]}_rc{rc_index_pair[1]}.png')

    # # Figure 3
    #
    # prepare_linregress_plot(work_dir, n_value, best_prefix,  f'Linear Regression on dFc-dFo in Cycle {content_cycle["cycle"]}, Round {content_cycle["round"]}')

    # Figure 4: Unexplained dFo and Overmodeled dFc blobs

    fig, ax = plt.subplots(num=3, clear=True)

    # ax.set_title('# of unexplained or overmodeled voxels')

    log_cycle_slice[['unexplained', 'overmodeled']].plot(ax=ax, marker='o', color=['magenta', 'blue'])
    ax.set_xlabel('cycle #')
    ax.set_ylabel('# of voxels')
    ax.legend()

    fig.savefig(f'{work_dir}unexplained_or_overmodeled.png')

    # Figure 4: Unexplained dFo and Overmodeled dFc blobs

    fig, ax = plt.subplots(num=7, clear=True)

    # ax.set_title('# of unexplained or overmodeled voxels')

    unexplained_out_of_dFo = log_cycle_slice['unexplained'] / log_cycle_slice['dFo_voxels']
    overmodeled_out_of_dFc = log_cycle_slice['overmodeled'] / log_cycle_slice['dFc_voxels']
    ax.plot(log_cycle_slice.index, unexplained_out_of_dFo * 100, label='unexplained / dFo', marker='o', color='magenta')
    ax.plot(log_cycle_slice.index, overmodeled_out_of_dFc * 100, label='overmodeled / dFc', marker='o', color='blue')
    ax.set_xlabel('cycle #')
    ax.set_ylabel('outlier ratio (%)')
    ax.legend()

    fig.savefig(f'{work_dir}unexplained_or_overmodeled_ratio.png')

    # Figure 4: RC log

    fig, ax = plt.subplots(num=4, clear=True)

    # ax.set_title('Evolution of Reaction Coordinates')

    props_cycle_frame = log_cycle_slice.set_index('cycle').drop(
        columns=['round', 'prefix', 'ok', 'hash', 'CC', 'beta', 'reason', 'unexplained', 'overmodeled', 'dFo_voxels', 'dFc_voxels', 'ttl_voxels'],
        errors='ignore', # beta may not be existed
    )
    # dark_cycle = props_cycle_frame.iloc[0]
    # relative_props_cycle_frame = (props_cycle_frame - dark_cycle + 180) % 360 - 180

    props_cycle_frame.rename(columns=lambda col_name: f'$\{col_name}$').plot(ax=ax, marker='o')
    # relative_props_cycle_frame.plot(ax=ax, marker='o')
    ax.set_xlabel('cycle #')
    ax.set_ylabel('absolute reaction coordinate')
    ax.legend()

    fig.savefig(f'{work_dir}evolution_rc.png')

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(prog='plot_dFoCC')
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

    for residue in args.residue:
        output_prefix = residue.output_prefix

        log_cycle_name = output_prefix + '_log_cycle.csv'
        log_cycle = pd.read_csv(log_cycle_name)

        rxn_coord_utils = prepare_reaction_coordinate_utils(residue, init_or_dark_xyzin_name)

        for cycle_num in log_cycle['cycle']:
            if cycle_num == 0: continue
            working_dir = output_prefix + f'_cycle{cycle_num:>03}/'
            prepare_all_plots(
                work_dir=working_dir,
                rxn_coord_utils=rxn_coord_utils,
                n_value=args.n_value,
                log_cycle=log_cycle,
                cycle_num=cycle_num,
            )
            prepare_linregress_plot(
                work_dir=working_dir,
                n_value=args.n_value,
                sigma_level=args.input.difference_data.sigma_level,
                dFo_filename=f"{PREPARING_DIR}{args.input.difference_data.file_name.rsplit('/', 1)[-1].replace('.mtz', '.map')}", 
                best_prefix=log_cycle['prefix'][cycle_num],
                title='',
            )
            print(f'The plots of Cycle {cycle_num} / {log_cycle["cycle"].max()} is done\r', end='')
        print()
