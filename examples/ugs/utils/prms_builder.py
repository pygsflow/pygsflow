import os
from pathlib import Path
import pandas as pd
import numpy as np
import flopy
import gsflow
import gridutil
from gsflow import (PrmsParameters, ParameterRecord, PrmsData, PrmsModel,
ControlFile, GsflowModel)
from gsflow.builder import (
    GenerateFishnet,
    ModflowBuilder,
    ControlFileBuilder,
    PrmsBuilder,
    FlowAccumulation
)
Debug = False

def fill_zeros_from_neighbors(array):
    """Fast fill zeros using 8-neighbor propagation with priority.

    Priority order: up, down, left, right, then diagonals. Propagates
    iteratively until no zeros can be filled.
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
            si0 = max(0, -di)
            si1 = nrows - max(0, di)
            sj0 = max(0, -dj)
            sj1 = ncols - max(0, dj)

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

def build_prms(mi):
    modify_fine_prms(mi)
    param_list = build_prms_parameters(mi)
    param_list = modify_unstrutured_prms(mi, param_list)
    param_obj = PrmsParameters(param_list)

    cdf = mi.fine_gsf.prms.data.data_df
    prmsdata = PrmsData(data_df=cdf)  
    control_obj = ControlFile(records_list=[], name=mi.usg_control_file) 
    control_records = mi.fine_gsf.prms.control.records_list
    for record in control_records:
        control_obj.records_list.append(record)

    control_obj.set_values("cascade_flag", [0])
    control_obj.set_values("cascadegw_flag", [0])
    
    p = Path(mi.usg_control_file)
    p = p.with_name(p.stem + "_output.csv")  # C:/data/file_v2.tar.gz
    control_obj.set_values("csv_output_file", [str(p)])
    
    p = Path(mi.usg_control_file)
    p = p.with_name(p.stem + "_gsflow_output.out") 
    control_obj.set_values("gsflow_output_file", [str(p)])

    p = Path(mi.usg_control_file)
    p = p.with_name(p.stem + "_par.params") 
    control_obj.set_values("param_file", [str(p)])

    p = Path(mi.usg_control_file)
    p = p.with_name(p.stem + "_stat_var.dat") 
    control_obj.set_values("stat_var_file", [str(p)])

    p = Path(mi.usg_control_file)
    p = p.with_name(p.stem + "_climate.dat") 
    control_obj.set_values("data_file", [str(p)])
    #prmsdata.file_name = [str(p)]



    # build the PrmsModel
    prms = PrmsModel(control_obj, parameters=param_obj, data=prmsdata)

    gsf = GsflowModel(control=control_obj, prms=prms, mf=mi.model)

    gsf.write_input(basename= mi.usg_base_name,
                    workspace=mi.usg_model_ws)

    mi.gsflow = gsf

def modify_fine_prms(mi):
    prms = mi.fine_gsf.prms
  
    new_value = mi.fine_gsf.mf.uzf.irunbnd.array
    #new_value = np.flipud(new_value)
    # fill zeros from 8-neighbors using fast vectorized propagation
    new_value = fill_zeros_from_neighbors(new_value)
                    
    new_value = new_value.flatten()
    dims = [['ncascade',len(new_value)]]

    # par1 = ParameterRecord(name= 'hru_down_id',
    #                     values= new_value,
    #                     dimensions=dims,
    #                     datatype= 1)
    # mi.fine_gsf.prms.parameters.remove_record("hru_down_id")
    # mi.fine_gsf.prms.parameters.parameters_list.append(par1)

    par1 = ParameterRecord(name= 'hru_strmseg_down_id',
                            values= new_value,
                            dimensions=dims,
                            datatype= 1)
    mi.fine_gsf.prms.parameters.remove_record("hru_strmseg_down_id")
    mi.fine_gsf.prms.parameters.parameters_list.append(par1)

    mi.fine_gsf.prms.parameters.set_values("ndeplval", np.array([22]))
    dims = [['ndeplval',22]]
    values = np.array([
        0.05000000074506,
        0.2399999946356,
        0.4000000059605,
        0.5299999713898,
        0.6499999761581,
        0.75,
        0.8199999928474,
        0.8799999952316,
        0.9300000071526,
        0.9900000095367,
        1.0,
        0.05000000074506,
        0.25,
        0.4000000059605,
        0.4799999892712,
        0.5400000214577,
        0.5799999833107,
        0.6100000143051,
        0.6399999856949,
        0.660000026226,
        0.6800000071526,
        0.6999999880791,
    ], dtype=float)
    par1 = ParameterRecord(name= 'snarea_curve',
                            values= values,
                            dimensions=dims,
                            datatype= 2)
    mi.fine_gsf.prms.parameters.remove_record("snarea_curve")
    mi.fine_gsf.prms.parameters.parameters_list.append(par1)

    

    # if 0:
    #     mi.fine_gsf.prms.parameters.remove_record("gw_down_id")
    #     mi.fine_gsf.prms.parameters.remove_record("gw_up_id")
    #     mi.fine_gsf.prms.parameters.remove_record("gw_pct_up")
    #     mi.fine_gsf.prms.parameters.remove_record("gw_strmseg_down_id")
    #     mi.fine_gsf.prms.control.set_values("cascade_flag", [1])
    #     mi.fine_gsf.prms.control.set_values("cascadegw_flag", [2])

    #     mi.fine_gsf.prms.parameters.set_values("ncascade", np.array([len(new_value)]))
    #     mi.fine_gsf.prms.parameters.set_values("ncascdgw", np.array([len(new_value)]))
    #     mi.fine_gsf.prms.parameters.set_values("cascade_flg", np.array([0]))
    #     mi.fine_gsf.prms.parameters.set_values("cascade_tol", np.array([0.0]))   
    # #mi.fine_gsf.prms.parameters.set_values("circle_switch", [0.0])

   
       

    #     par1 = ParameterRecord(name= 'hru_up_id',
    #                         values= np.zeros_like(new_value),
    #                         dimensions=dims,
    #                         datatype= 1)
    #     mi.fine_gsf.prms.parameters.remove_record("hru_up_id")
    #     mi.fine_gsf.prms.parameters.parameters_list.append(par1)

    #     par1 = ParameterRecord(name= 'hru_strmseg_down_id',
    #                         values= new_value,
    #                         dimensions=dims,
    #                         datatype= 1)
    #     mi.fine_gsf.prms.parameters.remove_record("hru_strmseg_down_id")
    #     mi.fine_gsf.prms.parameters.parameters_list.append(par1)
        
    #     par1 = ParameterRecord(name= 'hru_pct_up',
    #                         values= np.ones_like(new_value),
    #                         dimensions=dims,
    #                         datatype= 2)
    #     mi.fine_gsf.prms.parameters.remove_record("hru_pct_up")
    #     mi.fine_gsf.prms.parameters.parameters_list.append(par1)

def modify_unstrutured_prms(mi, param_list):

    nhrus_new = mi.gridprops['nodes']

    to_be_removed_params = [
        'gw_down_id',
        'gw_up_id',
        'gw_pct_up',
        'gw_strmseg_down_id',
        'ssr2gw_sq',
        'ssstor_init',
        'soil_moist_init',
        'soil_rechr_init',
        'soil_rechr_max',
    ]

    new_param_list = []
    for par in param_list:
        if par.name in ["hru_strmseg_down_id"]:
            pass
        if par.name in to_be_removed_params:
            continue

        elif par.name in ['cascade_flg']:
            par.values = np.array([0])
            new_param_list.append(par)

        elif par.name in ['cascade_tol']:
            par.values = np.array([0.0])
            new_param_list.append(par)

        elif par.name in ['hru_up_id']:
            dims = [[par.dimensions_names[0], nhrus_new]]
            par1 = ParameterRecord(name= 'hru_up_id',
                                values=  1+np.arange(nhrus_new),
                                dimensions=dims,
                                datatype= 1)
            new_param_list.append(par1)
        
        elif par.name in ['hru_down_id']:
            dims = [[par.dimensions_names[0], nhrus_new]]
            par1 = ParameterRecord(name= 'hru_down_id',
                                values= np.zeros(nhrus_new),
                                dimensions=dims,
                                datatype= 1)
            new_param_list.append(par1)
        

        elif par.name in ['hru_pct_up']:
            dims = [[par.dimensions_names[0], nhrus_new]]
            par1 = ParameterRecord(name= 'hru_pct_up',
                                values= np.ones(nhrus_new),
                                dimensions=dims,
                                datatype= 2)
            new_param_list.append(par1)
        
        elif par.name in ['gvr_cell_pct']:
            dims = [['nhrucell', nhrus_new]]
            par1 = ParameterRecord(name= 'gvr_cell_pct',
                                values= np.ones(nhrus_new),
                                dimensions=dims,
                                datatype= 2)
            new_param_list.append(par1)
        
        elif par.name in ['gvr_hru_id']:
            dims = [['nhrucell', nhrus_new]]
            par1 = ParameterRecord(name= 'gvr_hru_id',
                                values= 1+np.arange(nhrus_new),
                                dimensions=dims,
                                datatype= 1)
            new_param_list.append(par1)
        
        elif par.name in ['gvr_cell_id']:
            dims = [['nhrucell', nhrus_new]]
            par1 = ParameterRecord(name= 'gvr_cell_id',
                                values= 1+np.arange(nhrus_new),
                                dimensions=dims,
                                datatype= 1)
            new_param_list.append(par1)
        
        elif par.name in ['gvr_hru_pct']:
            dims = [['nhrucell', nhrus_new]]
            par1 = ParameterRecord(name= 'gvr_hru_pct',
                                values= nhrus_new*[1.0],
                                dimensions=dims,
                                datatype= 2)
            new_param_list.append(par1)

        elif par.name in ['smidx_coef']:
            par.values = np.array(nhrus_new*[0.005])
            new_param_list.append(par)

        elif par.name in ['smidx_exp']:
            par.values = np.array(nhrus_new*[0.3])
            new_param_list.append(par)
        
        elif par.name in ['snowinfil_max']:
            par.values = np.array(nhrus_new*[2.0])
            new_param_list.append(par)
        
        elif par.name in ['max_missing']:
            par.values = np.array([3])
            new_param_list.append(par)
        
        elif par.name in ['sat_threshold']:
            par.values = np.array(nhrus_new*[999.0])
            new_param_list.append(par)
        
        elif par.name in ['radadj_slope']:
            par.values = np.array(12*[0.0])
            new_param_list.append(par)
 
        else:
            new_param_list.append(par)
        

    return new_param_list

        
            
    #########################################################
    # if 0:
    #         # mi.fine_gsf.prms.parameters.remove_record("gw_down_id")
    #         # mi.fine_gsf.prms.parameters.remove_record("gw_up_id")
    #         # mi.fine_gsf.prms.parameters.remove_record("gw_pct_up")
    #         # mi.fine_gsf.prms.parameters.remove_record("gw_strmseg_down_id")
    #         mi.fine_gsf.prms.control.set_values("cascade_flag", [1])
    #         mi.fine_gsf.prms.control.set_values("cascadegw_flag", [2])

    #         mi.fine_gsf.prms.parameters.set_values("ncascade", np.array([len(new_value)]))
    #         mi.fine_gsf.prms.parameters.set_values("ncascdgw", np.array([len(new_value)]))
    #         mi.fine_gsf.prms.parameters.set_values("cascade_flg", np.array([0]))
    #         mi.fine_gsf.prms.parameters.set_values("cascade_tol", np.array([0.0]))   
    #     #mi.fine_gsf.prms.parameters.set_values("circle_switch", [0.0])

    
        

    #         # par1 = ParameterRecord(name= 'hru_up_id',
    #         #                     values= np.zeros_like(new_value),
    #         #                     dimensions=dims,
    #         #                     datatype= 1)
    #         # mi.fine_gsf.prms.parameters.remove_record("hru_up_id")
    #         # mi.fine_gsf.prms.parameters.parameters_list.append(par1)

    #         par1 = ParameterRecord(name= 'hru_strmseg_down_id',
    #                             values= new_value,
    #                             dimensions=dims,
    #                             datatype= 1)
    #         mi.fine_gsf.prms.parameters.remove_record("hru_strmseg_down_id")
    #         mi.fine_gsf.prms.parameters.parameters_list.append(par1)
            
    #         par1 = ParameterRecord(name= 'hru_pct_up',
    #                             values= np.ones_like(new_value),
    #                             dimensions=dims,
    #                             datatype= 2)
    #         mi.fine_gsf.prms.parameters.remove_record("hru_pct_up")
    #         mi.fine_gsf.prms.parameters.parameters_list.append(par1)

    #########################################################
   
def map_data(mi, dfmap, par):



    old_nhrus = mi.fine_gsf.mf.nrow * mi.fine_gsf.mf.ncol
    structured_df = pd.DataFrame(columns=['hru_id', 'value'])
    structured_df['hru_id'] = np.arange(old_nhrus)
    structured_df['value'] = par.values

    dfmap_ = dfmap.merge(
                        structured_df, left_on='source', right_on='hru_id',
                        how='left'
                    )
    if par.datatype == 2:  # float
        dfmap_['value'] = dfmap_['value'] * dfmap_['weight']
        ddf = dfmap_[['target', 'value']].copy()
        ddf.groupby('target').sum().reset_index(inplace=True)
        ddf.sort_values(by='target', inplace=True)      

    elif par.datatype == 1:  # integer
        ddf = dfmap_.loc[dfmap_.groupby('target')['weight'].idxmax()].reset_index(drop=True)
        ddf.sort_values(by='target', inplace=True)
        
        
    else:
        text = 1
    
    if Debug:
        gridutil.plot_grid(
            mi.oct_grid2d, ddf['value'], title=par.name
        )
    values = ddf['value'].values
    return values


def build_prms_parameters(mi):

    s_prms = mi.fine_gsf.prms
    
    parameters_list = s_prms.parameters.parameters_list

    uparams = []
    nhrus = mi.fine_gsf.mf.nrow * mi.fine_gsf.mf.ncol
    cell_size = mi.fine_gsf.mf.dis.delr[0] * mi.fine_gsf.mf.dis.delc[0]
    dfmap = (mi.mapping_df[mi.mapping_df['area'] > 0]).copy()    
    normalize_area = dfmap[['target','area']].groupby('target').sum()
    normalize_area.reset_index(inplace=True)
    normalize_area.rename(columns={'area': 'area_sum'}, inplace=True)
    dfmap = dfmap.merge(normalize_area, left_on='target', right_on='target', how='left')
    dfmap['weight'] = dfmap['area'] / dfmap['area_sum']
    del(dfmap['area_sum'])

    new_nhru = mi.gridprops['nodes']
    new_nreach = mi.model.sfr.nstrm
    # structured_df = pd.DataFrame(columns=['hru_id', 'value'])
    # structured_df['hru_id'] = np.arange(nhrus)

    # Debug = False
    parameter_records = []
    for par in parameters_list:      
        if par.section in ['Dimensions']:
            if par.name in ['nhru', 'nhrucell', 'ngw',
                           'ngwcell', 'nssr']:
                           new_value = [new_nhru]
            elif par.name in ['ncascade', 'ncascdgw']:
                new_value = [new_nhru]
            elif par.name in ['nsegment']:
                new_value = par.values
            elif par.name in ['nreach']:
                new_value = [new_nreach]
            else:
                new_value = par.values

            param_record = ParameterRecord(
                                            name=par.name,
                                            values= new_value,
                                            dimensions=par.dims,
                                            datatype=par.datatype                                            
                                        )
            parameter_records.append(param_record)
            

        elif par.section in ['Parameters']:
            dims = par.dims
            dim_names = par.dimensions_names
            if len(dims) == 1:
                if dims[0] == nhrus:
                    mapped_values = map_data(mi, dfmap, par)
                    param_record = ParameterRecord(
                            name=par.name,
                            values=mapped_values,
                            dimensions=[[par.dimensions_names[0], new_nhru]],
                            datatype=par.datatype                           
                        )
                    parameter_records.append(param_record)

                else:
                    param_record = ParameterRecord(
                                            name=par.name,
                                            values= par.values,
                                            dimensions= [[par.dimensions_names[0], par.dims[0]]],
                                            datatype=par.datatype                                            
                                        )
                    parameter_records.append(param_record)

            elif len(dims) == 2:
                if 'nhru' in par.dimensions_names:
                    dim2 = par.dims[1]
                    mapped_values = []
                    old_values = par.values.reshape(nhrus, dim2)
                    for i in range(dim2):
                        class p: pass
                        par_ = p()
                        par_.name = par.name
                        par_.values = old_values[:,i]                        
                        par_.datatype = par.datatype
                        mapped_values.append(map_data(mi, dfmap, par_))
                        
                    mapped_values = np.array(mapped_values).flatten()
                    dim_2d = [['nhru', new_nhru], [par.dimensions_names[1], dim2]]                     
                    param_record = ParameterRecord(
                            name=par.name,
                            values=mapped_values,
                            dimensions= dim_2d,
                            datatype=par.datatype                           
                        )
                    parameter_records.append(param_record)
                else:
                    xx = 1 # more than 2 dimensions

               
            else:
                xx = 1 # more than 2 dimensions


    
    return parameter_records
    