import numpy as np
import gemmi

def extract_unexplained_and_overmodeled_voxels(work_dir: str, best_prefix: str, sigma_level: float) -> tuple[float, float, float, float, float]:
    """
    To extract and count the voxels of unexplained (dFo != 0 but dFc == 0) and overmodeled (dFo == 0 but dFc != 0) region, which would significantly lower the CC value

    While overmodeled dFc would not affect the beta value (slope),
    unexplained dFo would lower the beta value quite a bit,
    leading to the overestimation of occupancy in dFo-beta analysis

    @param work_dir: the working directory of the input and output map files of a selected region
    @param best_prefix: the infix of the filenames
    @param sigma_level: another harcoded infix of the filenames
    @return: (unexplained_voxels, overmodelled_voxels, dFo_nonzero_voxels, dFc_nonzero_voxels, total_voxels)
    """
    dFo_map = gemmi.read_ccp4_map(f'{work_dir}dFo.{best_prefix}.masked_sigma{sigma_level:.1f}.ccp4')
    dFo_map.setup(float('nan'))
    dFc_map = gemmi.read_ccp4_map(f'{work_dir}dFc.{best_prefix}.masked_sigma{sigma_level:.1f}.ccp4')
    dFc_map.setup(float('nan'))

    def generate_map(grid: gemmi.FloatGrid, filename: str):
        map = gemmi.Ccp4Map()
        map.grid = grid
        map.update_ccp4_header()
        map.write_ccp4_map(filename)

    dFo_unexplained_grid = dFo_map.grid.clone()
    unexplained_selection = (dFo_map.grid.array != 0) & (dFc_map.grid.array == 0)
    dFo_unexplained_grid.array[~unexplained_selection] = 0.
    generate_map(dFo_unexplained_grid, f'{work_dir}dFo.{best_prefix}.unexplained.map')

    dFc_overmodeled_grid = dFc_map.grid.clone()
    overmodeled_selection = (dFo_map.grid.array == 0) & (dFc_map.grid.array != 0)
    dFc_overmodeled_grid.array[~overmodeled_selection] = 0.
    generate_map(dFc_overmodeled_grid, f'{work_dir}dFc.{best_prefix}.overmodeled.map')

    dFo_all_selection = (dFo_map.grid.array != 0) & ~ np.isnan(dFo_map.grid.array)
    dFc_all_selection = (dFc_map.grid.array != 0) & ~ np.isnan(dFc_map.grid.array)

    return (
        unexplained_selection.sum(),
        overmodeled_selection.sum(),
        dFo_all_selection.sum(),
        dFc_all_selection.sum(),
        dFo_map.grid.point_count,
    )
    # print(f'There are {unexplained_selection.sum()} unexplained dFo voxel(s)')
    # print(f'There are {overmodeled_selection.sum()} overmodeled dFc voxel(s)')
