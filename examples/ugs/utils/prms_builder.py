import os
from pathlib import Path

import numpy as np
import pandas as pd

import gsflow
from gsflow import (
    PrmsParameters, ParameterRecord, PrmsData,
    PrmsModel, ControlFile, GsflowModel
)
from utils import gridutil

DEBUG = False


def fill_zeros_from_neighbors(array):
    """Fill zero-valued cells from non-zero 8-connected neighbors.

    Uses iterative propagation with priority order: cardinal
    directions first (up, down, left, right), then diagonals.
    Repeats until no zeros remain or no further fills are possible.

    Args:
        array: 2D numpy array with zeros to be filled.

    Returns:
        Copy of the array with zeros filled from neighbors.
    """
    filled = array.copy()
    nrows, ncols = filled.shape
    offsets = [
        (-1, 0), (1, 0), (0, -1), (0, 1),
        (-1, -1), (-1, 1), (1, -1), (1, 1),
    ]

    changed = True
    while changed:
        changed = False
        for di, dj in offsets:
            # Source slice bounds (neighbor region)
            si0 = max(0, -di)
            si1 = nrows - max(0, di)
            sj0 = max(0, -dj)
            sj1 = ncols - max(0, dj)

            # Destination slice bounds (target region)
            di0 = max(0, di)
            di1 = nrows - max(0, -di)
            dj0 = max(0, dj)
            dj1 = ncols - max(0, -dj)

            src = filled[si0:si1, sj0:sj1]
            dst = filled[di0:di1, dj0:dj1]

            mask = (dst == 0) & (src != 0)
            if mask.any():
                dst[mask] = src[mask]
                changed = True

        if not (filled == 0).any():
            break

    return filled


def _set_control_output_paths(control_obj, control_file):
    """Configure output file paths in the PRMS control file.

    Derives output filenames from the control file stem
    (e.g., 'model_control' -> 'model_control_output.csv').

    Args:
        control_obj: PRMS ControlFile object.
        control_file: Path to the control file (str or Path).
    """
    p = Path(control_file)
    file_map = {
        "csv_output_file": "_output.csv",
        "gsflow_output_file": "_gsflow_output.out",
        "param_file": "_par.params",
        "stat_var_file": "_stat_var.dat",
        "data_file": "_climate.dat",
    }
    for key, suffix in file_map.items():
        out = p.with_name(p.stem + suffix)
        control_obj.set_values(key, [str(out)])


def build_prms(mi):
    """Build and write the unstructured PRMS/GSFLOW model.

    Orchestrates the full PRMS build pipeline:
    1. Modify fine-grid PRMS parameters (stream segments,
       snow depletion curves)
    2. Map structured parameters to the unstructured grid
    3. Adjust cascade and HRU parameters for unstructured grid
    4. Assemble control file, parameters, and climate data
    5. Write all GSFLOW input files

    Args:
        mi: Model_info object with grid and model attributes.
    """
    modify_fine_prms(mi)
    param_list = build_prms_parameters(mi)
    param_list = modify_unstructured_prms(mi, param_list)
    param_obj = PrmsParameters(param_list)

    # Copy climate data from the fine model
    cdf = mi.fine_gsf.prms.data.data_df
    prmsdata = PrmsData(data_df=cdf)

    # Build control file from fine model records
    control_obj = ControlFile(
        records_list=[], name=mi.usg_control_file
    )
    for record in mi.fine_gsf.prms.control.records_list:
        control_obj.records_list.append(record)

    # Disable cascade routing for unstructured grid
    control_obj.set_values("cascade_flag", [0])
    control_obj.set_values("cascadegw_flag", [0])

    _set_control_output_paths(control_obj, mi.usg_control_file)

    # Assemble and write GSFLOW model
    prms = PrmsModel(
        control_obj, parameters=param_obj, data=prmsdata
    )
    gsf = GsflowModel(
        control=control_obj, prms=prms, mf=mi.model
    )
    gsf.write_input(
        basename=mi.usg_base_name,
        workspace=mi.usg_model_ws
    )
    mi.gsflow = gsf


def modify_fine_prms(mi):
    """Fix fine-grid PRMS parameters before remapping.

    Updates the structured model's hru_strmseg_down_id using
    the UZF irunbnd array (with zero-fill from neighbors) and
    replaces the snow depletion curve (snarea_curve) with
    corrected 22-point values.

    Args:
        mi: Model_info object with fine_gsf attribute.
    """
    # Derive stream segment IDs from UZF runoff boundary
    new_value = mi.fine_gsf.mf.uzf.irunbnd.array
    new_value = fill_zeros_from_neighbors(new_value)
    new_value = new_value.flatten()
    dims = [['ncascade', len(new_value)]]

    par1 = ParameterRecord(
        name='hru_strmseg_down_id',
        values=new_value,
        dimensions=dims,
        datatype=1
    )
    mi.fine_gsf.prms.parameters.remove_record(
        "hru_strmseg_down_id"
    )
    mi.fine_gsf.prms.parameters.parameters_list.append(par1)

    # Replace snow depletion curve (2 curves x 11 points)
    mi.fine_gsf.prms.parameters.set_values(
        "ndeplval", np.array([22])
    )
    snarea_values = np.array([
        0.05000000074506, 0.2399999946356,
        0.4000000059605, 0.5299999713898,
        0.6499999761581, 0.75,
        0.8199999928474, 0.8799999952316,
        0.9300000071526, 0.9900000095367,
        1.0,
        0.05000000074506, 0.25,
        0.4000000059605, 0.4799999892712,
        0.5400000214577, 0.5799999833107,
        0.6100000143051, 0.6399999856949,
        0.660000026226, 0.6800000071526,
        0.6999999880791,
    ], dtype=float)
    par1 = ParameterRecord(
        name='snarea_curve',
        values=snarea_values,
        dimensions=[['ndeplval', 22]],
        datatype=2
    )
    mi.fine_gsf.prms.parameters.remove_record("snarea_curve")
    mi.fine_gsf.prms.parameters.parameters_list.append(par1)


def modify_unstructured_prms(mi, param_list):
    """Adjust mapped parameters for the unstructured grid.

    Removes parameters not applicable to unstructured grids
    and rebuilds cascade-related HRU parameters (hru_up_id,
    hru_down_id, hru_pct_up) with correct dimensions.

    Args:
        mi: Model_info object with gridprops attribute.
        param_list: List of ParameterRecord from mapping step.

    Returns:
        Filtered and updated list of ParameterRecord objects.
    """
    nhrus_new = mi.gridprops['nodes']

    # Parameters incompatible with unstructured grids
    remove_params = [
        'gw_down_id', 'gw_up_id', 'gw_pct_up',
        'gw_strmseg_down_id',
        'ssr2gw_sq', 'ssstor_init',
        'soil_moist_init', 'soil_rechr_init',
        'soil_rechr_max',
    ]

    new_param_list = []
    for par in param_list:
        if par.name in remove_params:
            continue

        elif par.name == 'hru_up_id':
            # Sequential 1-based HRU IDs
            dims = [[par.dimensions_names[0], nhrus_new]]
            par1 = ParameterRecord(
                name='hru_up_id',
                values=1 + np.arange(nhrus_new),
                dimensions=dims,
                datatype=1
            )
            new_param_list.append(par1)

        elif par.name == 'hru_down_id':
            # No cascade routing: all zeros
            dims = [[par.dimensions_names[0], nhrus_new]]
            par1 = ParameterRecord(
                name='hru_down_id',
                values=np.zeros(nhrus_new),
                dimensions=dims,
                datatype=1
            )
            new_param_list.append(par1)

        elif par.name == 'hru_pct_up':
            # Full contribution from each HRU
            dims = [[par.dimensions_names[0], nhrus_new]]
            par1 = ParameterRecord(
                name='hru_pct_up',
                values=np.ones(nhrus_new),
                dimensions=dims,
                datatype=2
            )
            new_param_list.append(par1)

        else:
            new_param_list.append(par)

    return new_param_list


def map_data(mi, dfmap, par):
    """Map a structured-grid parameter to the unstructured grid.

    For float parameters (datatype=2): area-weighted average.
    For integer parameters (datatype=1): value from the source
    cell with the largest overlap area (majority rule).

    Args:
        mi: Model_info object with fine_gsf attribute.
        dfmap: DataFrame with source-to-target mapping and
            area-based weights.
        par: Parameter-like object with name, values, and
            datatype attributes.

    Returns:
        numpy array of mapped values for the unstructured grid.
    """
    old_nhrus = mi.fine_gsf.mf.nrow * mi.fine_gsf.mf.ncol
    structured_df = pd.DataFrame({
        'hru_id': np.arange(old_nhrus),
        'value': par.values
    })

    dfmap_ = dfmap.merge(
        structured_df,
        left_on='source', right_on='hru_id',
        how='left'
    )

    if par.datatype == 2:
        # Area-weighted sum for continuous (float) parameters
        dfmap_['value'] = dfmap_['value'] * dfmap_['weight']
        ddf = dfmap_[['target', 'value']].copy()
        ddf = (
            ddf.groupby('target').sum()
            .reset_index()
            .sort_values(by='target')
        )

    elif par.datatype == 1:
        # Majority rule for discrete (integer) parameters
        idx = dfmap_.groupby('target')['weight'].idxmax()
        ddf = (
            dfmap_.loc[idx]
            .reset_index(drop=True)
            .sort_values(by='target')
        )

    if DEBUG:
        gridutil.plot_grid(
            mi.oct_grid2d, ddf['value'], title=par.name
        )

    return ddf['value'].values


def build_prms_parameters(mi):
    """Map all PRMS parameters from structured to unstructured.

    Iterates through the fine model's parameter list and:
    - Updates dimension records (nhru, ngw, etc.) to match the
      unstructured node count
    - Maps 1D HRU-dimensioned parameters using area weights
    - Maps 2D parameters (e.g., nhru x nmonths) column-by-column
    - Passes non-HRU parameters through unchanged

    Args:
        mi: Model_info object with fine_gsf, mapping_df,
            gridprops, and model attributes.

    Returns:
        List of ParameterRecord objects for the unstructured
        model.
    """
    s_prms = mi.fine_gsf.prms
    parameters_list = s_prms.parameters.parameters_list

    nhrus = mi.fine_gsf.mf.nrow * mi.fine_gsf.mf.ncol
    new_nhru = mi.gridprops['nodes']
    new_nreach = mi.model.sfr.nstrm

    # Build area-based mapping weights
    dfmap = mi.mapping_df[mi.mapping_df['area'] > 0].copy()
    normalize_area = (
        dfmap[['target', 'area']]
        .groupby('target').sum()
        .reset_index()
        .rename(columns={'area': 'area_sum'})
    )
    dfmap = dfmap.merge(
        normalize_area,
        left_on='target', right_on='target', how='left'
    )
    dfmap['weight'] = dfmap['area'] / dfmap['area_sum']
    del dfmap['area_sum']

    parameter_records = []
    for par in parameters_list:

        # --- Dimension records ---
        if par.section == 'Dimensions':
            if par.name in [
                'nhru', 'nhrucell', 'ngw',
                'ngwcell', 'nssr'
            ]:
                new_value = [new_nhru]
            elif par.name in ['ncascade', 'ncascdgw']:
                new_value = [new_nhru]
            elif par.name == 'nreach':
                new_value = [new_nreach]
            else:
                new_value = par.values

            param_record = ParameterRecord(
                name=par.name,
                values=new_value,
                dimensions=par.dims,
                datatype=par.datatype
            )
            parameter_records.append(param_record)

        # --- Parameter records ---
        elif par.section == 'Parameters':
            dims = par.dims

            if len(dims) == 1:
                if dims[0] == nhrus:
                    # 1D HRU parameter: remap to unstructured
                    mapped = map_data(mi, dfmap, par)
                    param_record = ParameterRecord(
                        name=par.name,
                        values=mapped,
                        dimensions=[[
                            par.dimensions_names[0],
                            new_nhru
                        ]],
                        datatype=par.datatype
                    )
                else:
                    # Non-HRU 1D parameter: pass through
                    param_record = ParameterRecord(
                        name=par.name,
                        values=par.values,
                        dimensions=[[
                            par.dimensions_names[0],
                            par.dims[0]
                        ]],
                        datatype=par.datatype
                    )
                parameter_records.append(param_record)

            elif len(dims) == 2:
                if 'nhru' in par.dimensions_names:
                    # 2D parameter with HRU dimension:
                    # map each column independently
                    dim2 = par.dims[1]
                    old_vals = par.values.reshape(nhrus, dim2)
                    mapped_cols = []
                    for i in range(dim2):

                        class _Par:
                            pass

                        col_par = _Par()
                        col_par.name = par.name
                        col_par.values = old_vals[:, i]
                        col_par.datatype = par.datatype
                        mapped_cols.append(
                            map_data(mi, dfmap, col_par)
                        )

                    mapped = np.array(mapped_cols).flatten()
                    dim_2d = [
                        ['nhru', new_nhru],
                        [par.dimensions_names[1], dim2]
                    ]
                    param_record = ParameterRecord(
                        name=par.name,
                        values=mapped,
                        dimensions=dim_2d,
                        datatype=par.datatype
                    )
                    parameter_records.append(param_record)

    return parameter_records
