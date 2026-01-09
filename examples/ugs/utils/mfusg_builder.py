import os
import shutil   
import numpy as np
import gsflow
from gsflow.builder import Defaults, ModflowDefaults, PrmsDefaults
import flopy
from utils import mf_utils
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d


def build_mfusg(mi):
    """
    Build the MODFLOW-USG model.
    """   
    modelname = os.path.basename(mi.usg_model_fn)
    mi.model = flopy.mfusg.MfUsg(model_ws= mi.usg_model_ws, modelname=modelname, structured=False)

    two_layers = False

    if two_layers: # if two layers, convert unstructured layer to 2 layers
        nlay = 2
        # convert unstructured layer to 2 layers
        nodes = mi.gridprops['nnodes'] * nlay
        njag = mi.gridprops['njag'] * nlay +  mi.gridprops['nnodes'] 
        iac = np.tile(mi.gridprops['iac'], 2) + 1

        ja_top = []
        ja_bottom = []
        ja_one_layer = mi.gridprops['ja'].tolist()
        start = 0
        ivc_top = []
        ivc_bottom = []
        for num in mi.gridprops['iac']:
            endd = start + num         
            current_ja =ja_one_layer[start:endd]
            start = endd
            
            # first layer ja
            lower_cell_id = current_ja[0] +  mi.gridprops['nnodes']
            first_layer_ja = current_ja + [lower_cell_id]
            second_layer_ja = [jj + mi.gridprops['nnodes'] for jj in current_ja]
            second_layer_ja = second_layer_ja + [current_ja[0]] # upper cell
        

            ja_top.extend(first_layer_ja)
            ja_bottom.extend(second_layer_ja)


        ja = ja_top + ja_bottom        
            
        disu = flopy.mfusg.MfUsgDisU(
            mi.model, nlay = nlay, nodes =nodes, ivsd = 0,
            njag = njag,
            idsymrd = 0, # full matrix
            nodelay = [int(nodes/nlay),   int(nodes/nlay)],
            iac = iac,
            ja = ja,
            ic12 = 1,
            nper =2, perlen = [1, 5356], nstp = [1, 5356], tsmult = 1, steady = [True, False],
            **mi.gridprops, itmuni=4, #days
            lenuni=2, #meters
            )
    else:
        disu = flopy.mfusg.MfUsgDisU(
        mi.model, **mi.gridprops,
        itmuni=4,
        nper=2, perlen=[1, 5356], nstp=[1, 5356], 
        steady=[True, False])

    build_bas(mi)
    build_npf(mi)
    build_oc(mi)    
    build_sfr(mi)

    # Write input files
    #mi.model.write_input()

def build_disu(mi):
    """
    Build the DISU package.
    """    
    
    nodes = mi.mfusg.disu.nodes
    nlay = mi.mfusg.disu.nlay
    njag = mi.mfusg.disu.njag
    nodelay = [int(nodes/nlay),   int(nodes/nlay)]

   

    flopy.mfusg.MfUsgDisu( 
        mi.model,
        nodes=nodes,
        nlay=nlay,
        njag=njag,
        ivsd=0, # no vertical nesting
        nper=2,
        itmuni=4, #days
        lenuni=2, #meters
        idsymrd=0,
        laycbd=0,
        nodelay=None,
        top=1,
        bot=0,
        area=1.0,
        iac=None,
        ja=None,
        ivc=None,
        cl1=None,
        cl2=None,
        cl12=None,
        fahl=None,
        perlen=[1, 5356],
        nstp=[1, 5356],
        tsmult=1,
        steady=[True, False],
        extension="disu",
        unitnumber=None,
        filenames=None,
        start_datetime=None,
    )

def build_bas(mi):
    """
    Build the BAS package.
    """
    #disu.write_file()
    nodes = mi.model.disu.nodes
    nod0_5 = int(nodes/mi.model.disu.nlay)
    nlay = mi.model.disu.nlay
    
    ibound = np.ones((nlay, nod0_5))

    # --- Initial conditions ---
    strt = np.zeros((nlay, nod0_5))
    if mi.nlay > 1: #todo: not tested for multiple layers
        for k in range(mi.nlay):
            strt[k] = mi.model.disu.top.array[k]        
    else:
        strt[0] = mi.model.disu.top.array


    bas = flopy.mfusg.MfUsgBas(
        mi.model,
        ibound=ibound.tolist(),
        strt=strt.tolist(),
        ifrefm=True,
        ixsec=False,
        ichflg=False,
        stoper=None,
        iprintfv=False,
        iprinttime=False,
        structured=False,
        converge=False,
        richards=False,
        double_prec=False,
        double_out=False,
        double_io=False,
        ihm=0,
        sy_all=False,
        hnoflo=-999.99,
        extension="bas",
        unitnumber=None,
        filenames=None)

def build_npf(mi):
    
    lpf = flopy.mfusg.MfUsgLpf(
        mi.model,
        laytyp=1,
        layavg=4,
        chani=1.0,
        layvka=0,
        laywet=0,
        ipakcb=None,
        hdry=-1e30,
        iwdflg=0,
        wetfct=0.1,
        iwetit=1,
        ihdwet=0,
        ikcflag=0,
        anglex=0,
        hk=0.0175,
        hani=1.0,
        vka=0.0175,
        ss=1e-6,
        sy=0.2,
        vkcb=0.0,
        wetdry=-0.01,
        ksat=1.0,
        storagecoefficient=False,
        constantcv=False,
        thickstrt=False,
        nocvcorrection=False,
        novfc=False,
        bubblept=False,
        fullydry=False,
        alpha=0,
        beta=0,
        sr=0,
        brook=0,
        bp=0
       )
    
        
def get_usg_sfr():
    # from reach_data and segment_data, get all cells ceneters
    # make a line for each segment.
    # intersect the line with cells of the unstructured grid.
    # compute reach_data for each cell.
    pass

def plot_line(list_of_lines, **kw):
    fig, ax = plt.subplots()
    x,y = list_of_lines.xy
    ax.plot(x, y, **kw)
    plt.show()

def build_sfr(mi):
    """

    """
    mf_utils.intersect_line_with_grid(mi)
    
    old_seg_df = pd.DataFrame(mi.fine_gsf.mf.sfr.segment_data[0])
    old_reach_df = pd.DataFrame(mi.fine_gsf.mf.sfr.reach_data)
    unique_segs = old_seg_df["nseg"].unique()
    new_reach_df = []
    for iseg in unique_segs:
        current_reach_df = old_reach_df[old_reach_df["iseg"] == iseg]
        current_reach_df_usg = mi.usg_sfr_df[mi.usg_sfr_df["segment_id"] == iseg]
        ordered_rech = current_reach_df_usg.sort_values(by = 'distance_from_segment_upstream')
        total_seg_length = current_reach_df['rchlen'].sum()
        df_reach = pd.DataFrame(columns = current_reach_df.columns)

        df_reach['node'] = ordered_rech['node_id']
        df_reach['rchlen'] = total_seg_length * (ordered_rech['intersection_length']/ordered_rech['intersection_length'].sum())
        df_reach['ireach'] = range(len(df_reach))
        df_reach['ireach'] = df_reach['ireach'] + 1
        df_reach['iseg'] = iseg        
        for col in ['strthick', 'strhc1', 'thts', 'thti', 'eps', 'uhc']:
            df_reach[col] = current_reach_df[col].iloc[0]
        
        # To compute top of stream bed op and slope, we will use the first reach elevation
        fine_str_x = current_reach_df['rchlen'].cumsum().values
        fine_str_z = current_reach_df['strtop'].values
        usg_str_x = df_reach['rchlen'].cumsum().values

        uniform_slope = (fine_str_z[0] - fine_str_z[-1])/total_seg_length
        df_reach['strtop'] = fine_str_z[0] - uniform_slope * usg_str_x
        df_reach['slope'] = uniform_slope  
        new_reach_df.append(df_reach)
    new_reach_df = pd.concat(new_reach_df)
    sfr_old = mi.fine_gsf.mf.sfr
    new_seg_df = sfr_old.segment_data
    nstrm = len(new_reach_df)
    nss = new_reach_df['iseg'].max()
    del (new_reach_df['reachID'])
    del (new_reach_df['outreach'])
    del (new_reach_df['i'])
    del (new_reach_df['j'])
    del (new_reach_df['k'])
    new_reach_df['rchlen'] = new_reach_df['rchlen'].astype(np.float32)
    new_reach_df['strtop'] = new_reach_df['strtop'].astype(np.float32)
    new_reach = new_reach_df.to_records(index = False)

    sfr = flopy.modflow.ModflowSfr2(mi.model,
       nstrm=nstrm,
        nss=nss,
        nsfrpar=sfr_old.nsfrpar,
        nparseg=sfr_old.nparseg,
        const=sfr_old.const,
        dleak=sfr_old.dleak,
        ipakcb=sfr_old.ipakcb,
        istcb2=sfr_old.istcb2,
        isfropt=sfr_old.isfropt,
        nstrail=sfr_old.nstrail,
        isuzn= sfr_old.isuzn,
        nsfrsets=sfr_old.nsfrsets,
        irtflg= sfr_old.irtflg,
        numtim=sfr_old.numtim,
        weight=sfr_old.weight,
        flwtol=0.0001,
        reach_data= new_reach,
        segment_data = new_seg_df,
        channel_geometry_data=None,
        channel_flow_data=None,
        dataset_5=sfr_old.dataset_5,
        irdflag=0,
        iptflag=0,
        reachinput=sfr_old.reachinput,
        transroute=sfr_old.transroute,
        tabfiles=sfr_old.tabfiles,
        tabfiles_dict=sfr_old.tabfiles_dict,
        options=sfr_old.options
    )
    return sfr


def build_oc(mi):
    
    oc = flopy.modflow.ModflowOc(
    mi.model,
    stress_period_data={
        (0, 0): ["print budget", "print head", "save head", "save budget"]
    },
)

def add_uzf(mi):

    basename= mi.usg_base_name
    workspace=mi.usg_model_ws

    uzf_fn = os.path.join(workspace, basename + ".uzf")

    fidw = open(uzf_fn, 'w')
    fidw.write('# UZF package\n')
    fidw.write('3 1 1 0 0 0 10 20 0 1.0\n')

    # 
    fidw.write("Constant 1 #IUZFBND")

   
    fidw.write("INTERNAL  1.0  (FREE)  -1  #irunbnd\n")
    val = mi.gsflow.prms.parameters.get_values('hru_strmseg_down_id')
    for v in val:
        fidw.write(f"{v} ")
    fidw.write("\n")

    fidw.write("CONSTANT    1.0                          #vks  \n")
    fidw.write("CONSTANT    3.5                           #eps  \n")
    fidw.write("CONSTANT    0.25                           #extdp  \n")
    fidw.write("CONSTANT    1.000000E+00                           #extdp  \n")

    fidw.write("         1 #finf for stress period 1\n")

    fidw.write("INTERNAL  1.0  (FREE)  -1  #finf1\n") 
    for v in val:
        fidw.write(f"{1.0} ")
    fidw.write("\n")

    fidw.write("-1 #finf for stress period 2")

    fidw.close()

def add_rch(mi):

    basename= mi.usg_base_name
    workspace=mi.usg_model_ws

    rch_fn = os.path.join(workspace, basename + ".rch")

    fidw = open(rch_fn, 'w')
    fidw.write('# rch package\n')
    fidw.write('3 0 # 2a. NRCHOP IRCHCB\n')
    fidw.write('1 0 # 5. INRECH INIRCH\n')
      
    fidw.write("INTERNAL  1.0  (FREE)  -1  #  6. RECH\n")
    val = mi.gsflow.prms.parameters.get_values('hru_strmseg_down_id')
    for v in val:
        fidw.write(f"{0.0} ")
    fidw.write("\n")
    fidw.write("-1 #finf for stress period 2")

    fidw.close()

def add_evt(mi):

    basename= mi.usg_base_name
    workspace=mi.usg_model_ws

    evt_fn = os.path.join(workspace, basename + ".evt")

    fidw = open(evt_fn, 'w')
    fidw.write('# EVT package\n')
    fidw.write('3 0 #NEVTOP IEVTCB\n')
    fidw.write('1 1 1 0 #  Set 5 - INSURF INEVTR INEXDP INIEVT\n')
    fidw.write("INTERNAL  1.0  (FREE)  -1  #  surf\n")
    val = mi.gsflow.mf.disu.top.array
    for v in val:
        fidw.write(f"{v} ")
    fidw.write("\n")

    fidw.write("INTERNAL  1.0  (FREE)  -1  #  EVTR\n")
    for v in val:
        fidw.write(f"{0.0} ")
    fidw.write("\n")

    fidw.write("INTERNAL  1.0  (FREE)  -1  #  EXDP\n")
    for v in val:
        fidw.write(f"{1.0} ")
    fidw.write("\n")


    fidw.write("-1 -1 -1 0# stress period 2")



    fidw.close()







