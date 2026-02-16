import os

import numpy as np
import pandas as pd
import flopy
import flopy.mfusg

from utils import mf_utils


def build_mfusg(mi):
    """Build the MODFLOW-USG model and all packages.

    Creates the MfUsg model object, the DISU (unstructured
    discretization) package, and delegates BAS, LPF, OC, and
    SFR package construction to helper functions.

    Args:
        mi: Model_info object with grid properties, paths,
            and fine-model attributes.
    """
    modelname = os.path.basename(mi.usg_model_fn)
    mi.model = flopy.mfusg.MfUsg(
        model_ws=mi.usg_model_ws,
        modelname=modelname,
        structured=False
    )

    flopy.mfusg.MfUsgDisU(
        mi.model, **mi.gridprops,
        itmuni=4,
        nper=2,
        perlen=[1, 5356],
        nstp=[1, 5356],
        steady=[True, False]
    )

    build_bas(mi)
    build_lpf(mi)
    build_oc(mi)
    build_sfr(mi)


def build_bas(mi):
    """Build the BAS (Basic) package.

    Sets all cells as active (ibound=1) and initialises
    hydraulic heads to the top elevation of each layer.

    Args:
        mi: Model_info object with model and nlay attributes.
    """
    nodes = mi.model.disu.nodes
    nlay = mi.model.disu.nlay
    nod_per_lay = int(nodes / nlay)

    ibound = np.ones((nlay, nod_per_lay))

    # Initial heads set to layer top elevations
    strt = np.zeros((nlay, nod_per_lay))
    if mi.nlay > 1:
        for k in range(mi.nlay):
            strt[k] = mi.model.disu.top.array[k]
    else:
        strt[0] = mi.model.disu.top.array

    flopy.mfusg.MfUsgBas(
        mi.model,
        ibound=ibound.tolist(),
        strt=strt.tolist(),
        ifrefm=True,
        structured=False,
        hnoflo=-999.99
    )


def build_lpf(mi):
    """Build the LPF (Layer-Property Flow) package.

    Uses convertible layers (laytyp=1) with the upstream
    weighting scheme (layavg=4). Hydraulic conductivity and
    storage parameters are set to uniform default values.

    Args:
        mi: Model_info object with the MODFLOW-USG model.
    """
    flopy.mfusg.MfUsgLpf(
        mi.model,
        laytyp=1,
        layavg=4,
        chani=1.0,
        layvka=0,
        laywet=0,
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
    )


def build_oc(mi):
    """Build the OC (Output Control) package.

    Configures head and budget output for the first time
    step of the first stress period.

    Args:
        mi: Model_info object with the MODFLOW-USG model.
    """
    flopy.modflow.ModflowOc(
        mi.model,
        stress_period_data={
            (0, 0): [
                "print budget", "print head",
                "save head", "save budget"
            ]
        },
    )


def build_sfr(mi):
    """Build the SFR2 (Streamflow-Routing) package.

    Maps structured-grid SFR reaches onto the unstructured
    grid by intersecting segment lines with octree cells.
    Reach lengths are redistributed proportionally to the
    intersection lengths, and streambed tops are linearly
    interpolated along each segment.

    Args:
        mi: Model_info object with fine_gsf, sfr_df, and
            usg_sfr_df attributes.

    Returns:
        The constructed ModflowSfr2 package object.
    """
    mf_utils.intersect_line_with_grid(mi)

    old_seg_df = pd.DataFrame(
        mi.fine_gsf.mf.sfr.segment_data[0]
    )
    old_reach_df = pd.DataFrame(
        mi.fine_gsf.mf.sfr.reach_data
    )
    unique_segs = old_seg_df["nseg"].unique()

    new_reach_dfs = []
    for iseg in unique_segs:
        cur_reach = old_reach_df[
            old_reach_df["iseg"] == iseg
        ]
        cur_usg = mi.usg_sfr_df[
            mi.usg_sfr_df["segment_id"] == iseg
        ]

        # Order reaches by distance from segment upstream end
        ordered = cur_usg.sort_values(
            by='distance_from_segment_upstream'
        )

        total_len = cur_reach['rchlen'].sum()
        df = pd.DataFrame(columns=cur_reach.columns)

        df['node'] = ordered['node_id']

        # Distribute total segment length proportionally
        int_len = ordered['intersection_length']
        df['rchlen'] = total_len * (int_len / int_len.sum())

        df['ireach'] = range(1, len(df) + 1)
        df['iseg'] = iseg

        # Copy constant reach properties from first reach
        for col in [
            'strthick', 'strhc1', 'thts',
            'thti', 'eps', 'uhc'
        ]:
            df[col] = cur_reach[col].iloc[0]

        # Interpolate streambed top using uniform slope
        fine_z = cur_reach['strtop'].values
        slope = (fine_z[0] - fine_z[-1]) / total_len
        usg_x = df['rchlen'].cumsum().values
        df['strtop'] = fine_z[0] - slope * usg_x
        df['slope'] = slope

        new_reach_dfs.append(df)

    new_reach_df = pd.concat(new_reach_dfs)

    # Remove structured-grid-specific columns
    for col in ['reachID', 'outreach', 'i', 'j', 'k']:
        if col in new_reach_df.columns:
            del new_reach_df[col]

    new_reach_df['rchlen'] = (
        new_reach_df['rchlen'].astype(np.float32)
    )
    new_reach_df['strtop'] = (
        new_reach_df['strtop'].astype(np.float32)
    )
    new_reach = new_reach_df.to_records(index=False)

    sfr_old = mi.fine_gsf.mf.sfr
    nstrm = len(new_reach_df)
    nss = new_reach_df['iseg'].max()

    sfr = flopy.modflow.ModflowSfr2(
        mi.model,
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
        isuzn=sfr_old.isuzn,
        nsfrsets=sfr_old.nsfrsets,
        irtflg=sfr_old.irtflg,
        numtim=sfr_old.numtim,
        weight=sfr_old.weight,
        flwtol=0.0001,
        reach_data=new_reach,
        segment_data=sfr_old.segment_data,
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


def _write_array(fid, values, fmt="{} "):
    """Write an array of values as space-separated text.

    Args:
        fid: Open file handle.
        values: Iterable of values to write.
        fmt: Format string for each value.
    """
    for v in values:
        fid.write(fmt.format(v))
    fid.write("\n")


def add_uzf(mi):
    """Write the UZF (Unsaturated-Zone Flow) package file.

    Writes a simple UZF input file with uniform soil
    properties and runoff boundary conditions derived
    from hru_strmseg_down_id.

    Args:
        mi: Model_info object with gsflow and path attributes.
    """
    uzf_fn = os.path.join(
        mi.usg_model_ws, mi.usg_base_name + ".uzf"
    )
    irunbnd = mi.gsflow.prms.parameters.get_values(
        'hru_strmseg_down_id'
    )

    with open(uzf_fn, 'w') as f:
        f.write("# UZF package\n")
        f.write("3 1 1 0 0 0 10 20 0 1.0\n")
        f.write("Constant 1 #IUZFBND\n")

        # Runoff boundary: route to stream segments
        f.write("INTERNAL  1.0  (FREE)  -1  #irunbnd\n")
        _write_array(f, irunbnd)

        # Uniform soil properties
        f.write("CONSTANT    1.0       #vks\n")
        f.write("CONSTANT    3.5       #eps\n")
        f.write("CONSTANT    0.25      #extdp\n")
        f.write("CONSTANT    1.0E+00   #extwc\n")

        # Stress period 1: uniform infiltration
        f.write("         1 #finf for stress period 1\n")
        f.write("INTERNAL  1.0  (FREE)  -1  #finf1\n")
        _write_array(f, [1.0] * len(irunbnd))

        # Stress period 2: reuse previous
        f.write("-1 #finf for stress period 2")


def add_rch(mi):
    """Write the RCH (Recharge) package file.

    Writes a simple recharge file with zero recharge for
    all cells (recharge is handled by PRMS/UZF coupling).

    Args:
        mi: Model_info object with gsflow and path attributes.
    """
    rch_fn = os.path.join(
        mi.usg_model_ws, mi.usg_base_name + ".rch"
    )
    nhru = len(
        mi.gsflow.prms.parameters.get_values(
            'hru_strmseg_down_id'
        )
    )

    with open(rch_fn, 'w') as f:
        f.write("# RCH package\n")
        f.write("3 0 # 2a. NRCHOP IRCHCB\n")
        f.write("1 0 # 5. INRECH INIRCH\n")

        # Zero recharge (handled by PRMS coupling)
        f.write("INTERNAL  1.0  (FREE)  -1  #  6. RECH\n")
        _write_array(f, [0.0] * nhru)

        # Stress period 2: reuse previous
        f.write("-1 #finf for stress period 2")


def add_evt(mi):
    """Write the EVT (Evapotranspiration) package file.

    Sets the ET surface to the grid top elevations, zero
    maximum ET rate, and uniform extinction depth of 1.0.

    Args:
        mi: Model_info object with gsflow and path attributes.
    """
    evt_fn = os.path.join(
        mi.usg_model_ws, mi.usg_base_name + ".evt"
    )
    top = mi.gsflow.mf.disu.top.array

    with open(evt_fn, 'w') as f:
        f.write("# EVT package\n")
        f.write("3 0 #NEVTOP IEVTCB\n")
        f.write(
            "1 1 1 0 #  Set 5 - INSURF INEVTR INEXDP INIEVT\n"
        )

        # ET surface = grid top elevations
        f.write("INTERNAL  1.0  (FREE)  -1  #  surf\n")
        _write_array(f, top)

        # Zero ET rate (handled by PRMS coupling)
        f.write("INTERNAL  1.0  (FREE)  -1  #  EVTR\n")
        _write_array(f, [0.0] * len(top))

        # Uniform extinction depth
        f.write("INTERNAL  1.0  (FREE)  -1  #  EXDP\n")
        _write_array(f, [1.0] * len(top))

        # Stress period 2: reuse previous
        f.write("-1 -1 -1 0# stress period 2")
