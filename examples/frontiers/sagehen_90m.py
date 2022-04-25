import os
import utm
import platform
import flopy
import numpy as np
import shapefile
import matplotlib.pyplot as plt
from flopy.utils import Raster
from flopy.plot import styles
from gsflow import GsflowModel, PrmsModel, PrmsData
from gsflow.builder import (
    GenerateFishnet,
    ModflowBuilder,
    ControlFileBuilder,
    PrmsBuilder,
    FlowAccumulation
)
import gsflow.builder.builder_utils as bu
import pandas as pd

pd.options.mode.chained_assignment = None


def nash_sutcliffe_efficiency(qsim, qobs, flg):
    if flg:
      qsim = np.log(qsim)
      qobs = np.log(qobs)
    qsim[np.isinf(qsim)] = np.nan
    qobs[np.isinf(qobs)] = np.nan
    numerator = np.nansum((qsim - qobs) ** 2)
    denominator = np.nansum((qobs - np.nanmean(qobs)) ** 2)
    nse = 1 - (numerator / denominator)
    return nse


def build_lut(f, dtype=int):
    d = {}
    with open(f) as foo:
        for line in foo:
            temp = line.strip().split("#")[0]
            if not temp:
                continue
            else:
                l = temp.split(":")
                d[dtype(l[0])] = float(l[1])
    return d


if __name__ == "__main__":
    sample_grid = True
    # set file names here
    ws = os.path.abspath(os.path.dirname(__file__))
    iws = os.path.join(ws, "..", "data", "geospatial")
    ows = os.path.join(ws, "temp")
    if not os.path.exists(ows):
        os.mkdir(ows)

    dem_file = os.path.join(iws, 'dem.img')
    pour_point_file = os.path.join(iws, "model_points.shp")

    resampled_dem = os.path.join(ows, "sagehen_90m_med.txt")

    stream_threshold = 810000 # m3 of drainage area
    cellsize = 90

    # generate a "Fishnet"
    modelgrid = GenerateFishnet(dem_file, xcellsize=cellsize, ycellsize=cellsize)

    # resample DEM to the model grid using minimum elevation
    if sample_grid:
        raster = Raster.load(dem_file)
        dem = raster.resample_to_grid(
            modelgrid,
            band=raster.bands[0],
            method="median",
            multithread=True,
            thread_pool=12
        )
        np.savetxt(resampled_dem, dem, delimiter="  ")
    else:
        dem = np.genfromtxt(resampled_dem)

    hru_type = np.ones((modelgrid.nrow, modelgrid.ncol), dtype=int)
    fa = FlowAccumulation(dem,
                          modelgrid.xcellcenters,
                          modelgrid.ycellcenters,
                          hru_type=hru_type, verbose=True)

    flow_dir = fa.flow_directions(dijkstra=True)
    qx, qy = fa.get_vectors

    flow_acc = fa.flow_accumulation()

    # read in pour point from shapefile and set watershed boundary
    with shapefile.Reader(pour_point_file) as r:
        shape = r.shape(0)
        pour_point = shape.points
        pour_point[0][1] -= 20
        pour_point[0][0] -= 20

    watershed = fa.define_watershed(pour_point, modelgrid, fmt='xy')

    threshold = stream_threshold / (cellsize ** 2)
    strm_obj = fa.make_streams(flow_dir, flow_acc, threshold)
    cascades = fa.get_cascades(strm_obj)

    i = np.floor(3646 / modelgrid.ncol)
    j = 3646 % modelgrid.ncol

    mfbuild = ModflowBuilder(modelgrid, dem, "sagehen_90m")
    botm = dem - 100
    botm.shape = (1, modelgrid.nrow, modelgrid.ncol)
    ml = mfbuild.build_all(strm_obj.reach_data,
                           strm_obj.segment_data,
                           strm_obj.irunbnd,
                           finf=np.ones(dem.shape),
                           botm=botm,
                           ibound=watershed.astype(int),
                           iuzfbnd=watershed.astype(int)
                          )

    # update dis file to create a transient model
    flopy.modflow.ModflowDis(
        ml,
        nlay=ml.dis.nlay,
        nrow=ml.dis.nrow,
        ncol=ml.dis.ncol,
        nper=2,
        delr=ml.dis.delr,
        delc=ml.dis.delc,
        laycbd=ml.dis.laycbd,
        top=ml.dis.top,
        botm=ml.dis.botm,
        perlen=[1, 5356],
        nstp=[1, 5356],
        tsmult=[1, 1],
        steady=[True, False],
        itmuni=ml.dis.itmuni,
        lenuni=ml.dis.lenuni
    )

    # update a few SFR parameters for GSFLOW!
    ml.sfr.segment_data[0]["flow"] *= 0
    ml.sfr.segment_data[0]["roughch"] = 0.04
    ml.sfr.reach_data["strhc1"] = 0.1

    # tune some of the other MODFLOW parameters
    ml.upw.hk *= 1.75e-03  #  (default is 10)
    ml.upw.ss *= 1.0 #  (default is 1e-06)

    prms_outfile = os.path.join(ows, "sagehen.param")
    prmsbuild = PrmsBuilder(strm_obj, cascades, modelgrid, fa.get_dem_data().ravel(),
                            hru_type=watershed, hru_subbasin=watershed)
    param_obj = prmsbuild.build()
    lat, lon = utm.to_latlon(modelgrid.xcellcenters.ravel(), modelgrid.ycellcenters.ravel(), 10, "N")
    param_obj.set_values("hru_lat", lat)
    param_obj.set_values("hru_lon", lon)

    sample_rasters = True
    nhru = modelgrid.nrow * modelgrid.ncol
    # load in rasters and luts for parameterizing prms
    veg_type_raster = os.path.join(iws, "us_140evt_utm.img")
    veg_cov_raster = os.path.join(iws, "us_140evc_utm.img")
    awc_raster = os.path.join(iws, "awc.img")
    clay_raster = os.path.join(iws, "clay.img")
    ksat_raster = os.path.join(iws, "ksat.img")
    sand_raster = os.path.join(iws, "sand.img")
    impervious_raster = os.path.join(iws, "nlcd2011_imp_utm.img")
    prism = {"ppt_utm": [], "tmax_utm": [], "tmin_utm": []}
    for folder in prism.keys():
        for f in os.listdir(os.path.join(iws, "climate", folder)):
            if os.path.isfile(os.path.join(iws, "climate", folder, f)) and f.endswith(".img"):
                prism[folder].append(os.path.join(iws, "climate", folder, f))

    resampled_veg_type = os.path.join(ows, "veg_type_nearest_90.txt")
    resampled_veg_cov = os.path.join(ows, "veg_cov_nearest_90.txt")
    resampled_awc = os.path.join(ows, "awc_median_90.txt")
    resampled_clay = os.path.join(ows, "clay_median_90.txt")
    resampled_ksat = os.path.join(ows, "ksat_median_90.txt")
    resampled_sand = os.path.join(ows, "sand_median_90.txt")
    resampled_impervious = os.path.join(ows, "impervious_median_90.txt")

    resampled_ppt = os.path.join(ows, "ppt_bilinear_90.txt")
    resampled_tmax = os.path.join(ows, 'tmax_bilinear_90.txt')
    resampled_tmin = os.path.join(ows, 'tmin_bilinear_90.txt')

    covtype_remap = os.path.join(iws, "..", 'remaps', "landfire", "covtype.rmp")
    covden_sum_remap = os.path.join(iws, "..", 'remaps', "landfire", "covdensum.rmp")
    covden_win_remap = os.path.join(iws, "..", 'remaps', "landfire", "covdenwin.rmp")
    root_depth_remap = os.path.join(iws, "..", 'remaps', "landfire", 'rtdepth.rmp')
    snow_intcp_remap = os.path.join(iws, "..", "remaps", "landfire", "snow_intcp.rmp")
    srain_intcp_remap = os.path.join(iws, "..", "remaps", "landfire", "srain_intcp.rmp")

    climate_dataframe = os.path.join(iws, 'climate', "sagehen_climate.csv")
    climate_lapse_rates = os.path.join(iws, 'climate', 'sagehen_lapse_rates.csv')

    if sample_rasters:
        ibound = watershed.astype(int)
        raster = Raster.load(veg_type_raster)
        veg_type = raster.resample_to_grid(
            modelgrid,
            band=raster.bands[0],
            method="nearest",
            multithread=True,
            thread_pool=12
        )
        veg_type[ibound == 0] = 0
        veg_type = veg_type.astype(int)
        np.savetxt(resampled_veg_type, veg_type, fmt="%d")

        raster = Raster.load(veg_cov_raster)
        # todo: this might need to be a zonal statistic (mode)
        veg_cov = raster.resample_to_grid(
            modelgrid,
            band=raster.bands[0],
            method="nearest",
            multithread=True,
            thread_pool=12
        )
        veg_cov[ibound == 0] = 0
        veg_cov = veg_cov.astype(int)
        np.savetxt(resampled_veg_cov, veg_cov, fmt="%d")

        raster = Raster.load(awc_raster)
        awc = raster.resample_to_grid(
            modelgrid,
            band=raster.bands[0],
            method="median",
            multithread=True,
            thread_pool=12
        )
        awc[ibound == 0] = 0
        awc[awc == raster.nodatavals[0]] = np.nanmedian(awc)
        np.savetxt(resampled_awc, awc)

        raster = Raster.load(ksat_raster)
        ksat = raster.resample_to_grid(
            modelgrid,
            band=raster.bands[0],
            method="median",
            multithread=True,
            thread_pool=12
        )
        ksat[ibound == 0] = 0
        ksat[ksat == raster.nodatavals[0]] = np.nanmedian(ksat)
        np.savetxt(resampled_ksat, ksat)

        raster = Raster.load(sand_raster)
        sand = raster.resample_to_grid(
            modelgrid,
            band=raster.bands[0],
            method="median",
            multithread=True,
            thread_pool=12
        )
        sand[ibound == 0] = 0
        sand[sand == raster.nodatavals[0]] = np.nanmedian(sand)
        sand /= 100
        np.savetxt(resampled_sand, sand)

        raster = Raster.load(clay_raster)
        clay = raster.resample_to_grid(
            modelgrid,
            band=raster.bands[0],
            method="median",
            multithread=True,
            thread_pool=12
        )
        clay[ibound == 0] = 0
        clay[clay == raster.nodatavals[0]] = np.nanmedian(clay)
        clay /= 100
        np.savetxt(resampled_clay, clay)

        raster = Raster.load(impervious_raster)
        impervious = raster.resample_to_grid(
            modelgrid,
            band=raster.bands[0],
            method="median",
            multithread=True,
            thread_pool=12
        )
        impervious[ibound == 0] = 0
        impervious /= 100
        np.savetxt(resampled_impervious, impervious)

        ppt = []
        for rstr in prism["ppt_utm"]:
            raster = Raster.load(rstr)
            tppt = raster.resample_to_grid(
                modelgrid,
                band=raster.bands[0],
                method="linear",
                multithread=True,
                thread_pool=12
            )
            ppt.append(tppt.ravel())
        ppt = np.array(ppt)
        np.savetxt(resampled_ppt, ppt)

        tmin = []
        for rstr in prism["tmin_utm"]:
            raster = Raster.load(rstr)
            ttmin = raster.resample_to_grid(
                modelgrid,
                band=raster.bands[0],
                method="linear",
                multithread=True,
                thread_pool=12
            )
            tmin.append(ttmin.ravel())
        tmin = np.array(tmin)
        np.savetxt(resampled_tmin, tmin)

        tmax = []
        for rstr in prism["tmax_utm"]:
            raster = Raster.load(rstr)
            ttmax = raster.resample_to_grid(
                modelgrid,
                band=raster.bands[0],
                method="linear",
                multithread=True,
                thread_pool=12
            )
            tmax.append(ttmax.ravel())
        tmax = np.array(tmax)
        np.savetxt(resampled_tmax, tmax)

    else:
        veg_type = np.genfromtxt(resampled_veg_type, dtype=int)
        veg_cov = np.genfromtxt(resampled_veg_cov, dtype=int)
        awc = np.genfromtxt(resampled_awc)
        ksat = np.genfromtxt(resampled_ksat)
        sand = np.genfromtxt(resampled_sand)
        clay = np.genfromtxt(resampled_clay)
        impervious = np.genfromtxt(resampled_impervious)
        ppt = np.genfromtxt(resampled_ppt)
        tmax = np.genfromtxt(resampled_tmax)
        tmin = np.genfromtxt(resampled_tmin)
        ppt.shape = (12, nhru)
        tmax.shape = (12, nhru)
        tmin.shape = (12, nhru)

    covtype_lut = build_lut(covtype_remap)
    covden_sum_lut = build_lut(covden_sum_remap)
    covden_win_lut = build_lut(covden_win_remap)
    root_depth_lut = build_lut(root_depth_remap)
    snow_intcp_lut = build_lut(snow_intcp_remap)
    srain_intcp_lut = build_lut(srain_intcp_remap)

    # read in "climate dataframe"
    cdf = pd.read_csv(climate_dataframe)
    ldf = pd.read_csv(climate_lapse_rates)

    # build vegatative cover parameters
    covtype = bu.covtype(veg_type, covtype_lut)
    covden_sum = bu.covden_sum(veg_cov, covden_sum_lut)
    covden_win = bu.covden_win(covtype.values, covden_win_lut)
    rad_trncf = bu.rad_trncf(covden_win.values)
    snow_intcp = bu.snow_intcp(veg_type, snow_intcp_lut)
    srain_intcp = bu.srain_intcp(veg_type, srain_intcp_lut)
    wrain_intcp = bu.wrain_intcp(veg_type, snow_intcp_lut)

    # add veg to param_obj
    param_obj.add_record_object(covtype, True)
    param_obj.add_record_object(covden_sum, True)
    param_obj.add_record_object(covden_win, True)
    param_obj.add_record_object(rad_trncf, True)
    param_obj.add_record_object(snow_intcp, True)
    param_obj.add_record_object(srain_intcp, True)
    param_obj.add_record_object(wrain_intcp, True)

    # build soil parameters
    root_depth = bu.root_depth(veg_type, root_depth_lut)
    hru_aspect = bu.d8_to_hru_aspect(flow_dir)
    hru_slope = bu.d8_to_hru_slope(
        flow_dir,
        dem,
        modelgrid.xcellcenters,
        modelgrid.ycellcenters
    )

    soil_type = bu.soil_type(clay, sand)
    soil_moist_max = bu.soil_moist_max(awc, root_depth)
    soil_moist_init = bu.soil_moist_init(soil_moist_max.values)
    soil_rech_max = bu.soil_rech_max(awc, root_depth)
    ssr2gw_rate = bu.ssr2gw_rate(ksat, sand, soil_moist_max.values)
    ssr2gw_sq = bu.ssr2gw_exp(nhru)
    soil_rech_init = bu.soil_rech_init(soil_rech_max.values)
    slowcoef_lin = bu.slowcoef_lin(ksat, hru_aspect.values, cellsize, cellsize)

    slowcoef_sq = bu.slowcoef_sq(
        ksat,
        hru_aspect.values,
        sand,
        soil_moist_max.values,
        cellsize,
        cellsize
    )

    # add soil parameters to prms object
    param_obj.add_record_object(hru_slope, replace=True)
    param_obj.add_record_object(hru_aspect, replace=True)
    param_obj.add_record_object(soil_type, replace=True)
    param_obj.add_record_object(soil_moist_max, replace=True)
    param_obj.add_record_object(soil_moist_init, replace=True)
    param_obj.add_record_object(soil_rech_max, replace=True)
    param_obj.add_record_object(soil_rech_init, replace=True)
    param_obj.add_record_object(ssr2gw_rate, replace=True)
    param_obj.add_record_object(ssr2gw_sq, replace=True)
    param_obj.add_record_object(slowcoef_lin, replace=True)
    param_obj.add_record_object(slowcoef_sq, replace=True)

    # imperviousness parameters
    hru_percent_imperv = bu.hru_percent_imperv(impervious)
    carea_max = bu.carea_max(impervious)

    # add imperv to prms obj
    param_obj.add_record_object(hru_percent_imperv, replace=True)
    param_obj.add_record_object(carea_max, replace=True)

    # climate parameters
    param_obj.add_record(name="nobs", values=[1,])
    outlet_sta = modelgrid.intersect(pour_point[0][0], pour_point[0][1])
    outlet_sta = modelgrid.get_node([(0,) + outlet_sta])

    cdf = bu.add_prms_date_columns_to_df(cdf, "date")
    cdf.rename(
        columns={
            'precip': 'precip_0',
            'tmin': 'tmin_0',
            'tmax': 'tmax_0',
            'runoff': 'runoff_0',
            'date': 'Date'
        },
        inplace=True
    )
    # reorder dataframe to later build a prms Data object from it
    cdfcols = [
        "Year", "Month", "Day", "Hour", "Minute", "Second",
        "tmax_0", "tmin_0", "precip_0", "runoff_0", "Date"
    ]
    cdf = cdf[cdfcols]

    # start climate parameter calculations
    mean_ppt = bu.get_mean_monthly_from_df(cdf, 'precip_0')
    cdf["tmax_0"] = bu.fahrenheit_to_celsius(cdf["tmax_0"].values)
    cdf["tmin_0"] = bu.fahrenheit_to_celsius(cdf["tmin_0"].values)
    mean_tmax = bu.get_mean_monthly_from_df(cdf, "tmax_0", temperature=True)
    mean_tmin = bu.get_mean_monthly_from_df(cdf, "tmin_0", temperature=True)

    rain_adj = bu.rain_adj(ppt, mean_ppt)
    snow_adj = bu.snow_adj(ppt, mean_ppt)

    tmin_lapse = bu.tmin_lapse(ldf.tmin_lapse.values * (5 / 9))
    tmax_lapse = bu.tmax_lapse(ldf.tmax_lapse.values * (5 / 9))

    tmax_adj = bu.tmax_adj(nhru)
    tmin_adj = bu.tmin_adj(nhru)

    jh_coef = bu.calculate_jensen_haise(dem, mean_tmin, mean_tmax)

    # add climate parameters to param obj
    param_obj.add_record_object(rain_adj, replace=True)
    param_obj.add_record_object(snow_adj, replace=True)
    param_obj.add_record_object(tmin_lapse, replace=True)
    param_obj.add_record_object(tmax_lapse, replace=True)
    param_obj.add_record_object(tmax_adj, replace=True)
    param_obj.add_record_object(tmin_adj, replace=True)
    param_obj.add_record_object(jh_coef, replace=True)
    param_obj.add_record(
        "outlet_sta",
        values=[outlet_sta[0] + 1,],
        dimensions=[["one", 1]],
        datatype=1
    )
    param_obj.add_record(
        "id_obsrunoff",
        values=[outlet_sta[0] + 1, ],
        dimensions=[["one", 1]],
        datatype=1
    )

    param_obj.add_record(
        "tsta_elev",
        values=[1932.4,],
        dimensions=[["ntemp", 1]],
        datatype=2
    )

    # build the prms data file
    prmsdata = PrmsData(data_df=cdf)   
    control_obj = ControlFileBuilder().build("saghen_90m", param_obj, ml)

    # build the PrmsModel
    prms = PrmsModel(control_obj, parameters=param_obj, data=prmsdata)

    gsf = GsflowModel(control=control_obj, prms=prms, mf=ml)

    gsf.control.set_values("start_time", [1982, 10, 1, 0, 0, 0])
    gsf.control.add_record("end_time", values=[1996, 9, 31, 0, 0, 0])
    gsf.control.add_record("print_debug", values=[0, ])
    gsf.control.add_record("modflow_time_zero", values=[1982, 10, 1, 0, 0, 0])
    gsf.control.add_record("data_file", values=["sagehen_90m.data",])
    gsf.control.add_record("srunoff_module", values=["srunoff_smidx"])
    gsf.control.add_record("model_output_file", values=["gsflow_sagehen_90.out"])
    gsf.control.set_values("model_mode", values=["GSFLOW5"])
    gsf.control.set_values("subbasin_flag", values=[0,])
    gsf.control.set_values("parameter_check_flag", values=[0, ])
    gsf.control.add_record("statsON_OFF", values=[1])
    gsf.control.add_record("nstatVars", values=[6])
    gsf.control.add_record("statVar_element", values=["1", "1", "1", "1", "1", "1"])
    gsf.control.add_record("statVar_names",
                           values=["runoff",
                                   "basin_cfs",
                                   "basin_ssflow_cfs",
                                   "basin_gwflow_cfs",
                                   "basin_sroff_cfs",
                                   "basin_dunnian"])
    gsf.control.add_record("stat_var_file", values=["statvar_uncal.dat"])

    # Modify PRMS paramters for calibration
    # temp dist
    tmax_lapse = gsf.prms.parameters.get_values('tmax_lapse')
    tmin_lapse = gsf.prms.parameters.get_values('tmin_lapse')
    tmax_lapse = tmax_lapse + 1.2   #0.7
    tmin_lapse = tmin_lapse + 1.2   #0.7
    gsf.prms.parameters.set_values("tmax_lapse", values=tmax_lapse)
    gsf.prms.parameters.set_values("tmin_lapse", values=tmin_lapse)
    max_missing = gsf.prms.parameters.get_values('max_missing')
    max_missing = max_missing * 2
    gsf.prms.parameters.set_values("max_missing", values=max_missing)
    # snow
    tmax_allsnow = gsf.prms.parameters.get_values('tmax_allsnow')
    tmax_allsnow[:] = 0.7
    gsf.prms.parameters.set_values("tmax_allsnow", values=tmax_allsnow)
    value = [2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1]
    gsf.prms.parameters.add_record("tmax_allrain_offset", values=value, dimensions=[('nmonths', 12)])
    covden_win = gsf.prms.parameters.get_values('covden_win')
    rad_trncf = gsf.prms.parameters.get_values('rad_trncf')
    rad_trncf = 0.8 * covden_win  # correlated to covden_win
    gsf.prms.parameters.set_values("rad_trncf", values=rad_trncf)
    # ET
    soil_moist_max = gsf.prms.parameters.get_values('soil_moist_max')
    soil_moist_max = soil_moist_max * 3.0
    gsf.prms.parameters.set_values("soil_moist_max", values=soil_moist_max)
    value = [0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03]
    gsf.prms.parameters.add_record("jh_coef", values=value, dimensions=[('nmonths', 12)])
    # adding a 2D parameter like rain_adj then dimensions = [('nmonths', 12), ('nrhu', 5545454)]
    # runoff
    snowinfil_max = gsf.prms.parameters.get_values('snowinfil_max')
    snowinfil_max = snowinfil_max * 15.0
    gsf.prms.parameters.set_values("snowinfil_max", values=snowinfil_max)
    smidx_coef = gsf.prms.parameters.get_values('smidx_coef')
    smidx_coef = smidx_coef / 100.0
    smidx_exp = gsf.prms.parameters.get_values('smidx_exp')
    smidx_exp = smidx_exp / 100.0
    carea_max = gsf.prms.parameters.get_values('carea_max')
    carea_max = carea_max / 100.0
    gsf.prms.parameters.set_values("smidx_coef", values=smidx_coef)
    gsf.prms.parameters.set_values("smidx_exp", values=smidx_exp)
    gsf.prms.parameters.set_values("carea_max", values=carea_max)
    # interflow
    slowcoef_sq = gsf.prms.parameters.get_values('slowcoef_sq')
    slowcoef_sq = slowcoef_sq * 0.1
    gsf.prms.parameters.set_values("slowcoef_sq", values=slowcoef_sq)
    slowcoef_lin = gsf.prms.parameters.get_values('slowcoef_lin')
    slowcoef_lin = slowcoef_lin * 3.0
    gsf.prms.parameters.set_values("slowcoef_lin", values=slowcoef_lin)
    # Recharge
    ssr2gw_rate = gsf.prms.parameters.get_values('ssr2gw_rate')
    ssr2gw_rate = ssr2gw_rate * 500.0
    gsf.prms.parameters.set_values("ssr2gw_rate", values=ssr2gw_rate)
    sat_threshold = gsf.prms.parameters.get_values('sat_threshold')
    sat_threshold = sat_threshold / 3
    gsf.prms.parameters.set_values("sat_threshold", values=sat_threshold)

    # clean unused parameters
    par_to_remove = ["gw_up_id", "gw_down_id", "gw_strmseg_down_id", "gw_pct_up"]
    for par_ in par_to_remove:
        gsf.prms.parameters.remove_record(par_)

    gsf.write_input(basename="sagehen_90m", workspace=ows)

    gsf.prms.parameters.remove_record("adjmix_rain")

    exe_name = os.path.join("..", "..", "bin", "gsflow")
    if platform.system().lower() == "windows":
        exe_name += ".exe"

    # reload the model to assure that it has valid formatting
    gsf = GsflowModel.load_from_file(os.path.join(ows, "sagehen_90m_cont.control"))
    gsf.run_model(gsflow_exe=exe_name)

    # load PRMS output with simulated and measured streamflow and flow components
    stats = gsf.prms.get_StatVar()

    # get data from 10/1/1985 onward
    stats = stats[1096:]
    stats.reset_index(inplace=True, drop=True)

    # Calculate N-S and Log(N-S)
    nse_val = nash_sutcliffe_efficiency(stats.basin_cfs_1, stats.runoff_1, False)
    nse_val_log = nash_sutcliffe_efficiency(stats.basin_cfs_1, stats.runoff_1, True)
    nse_val_log_str = str(nse_val_log)
    print(nse_val,nse_val_log)

    gw_seepage = stats.basin_cfs_1.values.copy() - (
            stats.basin_ssflow_cfs_1.values.copy() +
            stats.basin_sroff_cfs_1.values.copy() +
            stats.basin_dunnian_1.values.copy()
    )

    with styles.USGSMap():
        fig, axis = plt.subplots(2,1,figsize=(10, 6))
        plt.rcParams.update({'font.size': 100})
        axis[0].plot(stats.Date, stats.basin_cfs_1, color='r', linewidth=2.2, label='simulated')
        axis[0].plot(stats.Date, stats.runoff_1, '--', color='b', linewidth=1.5, label='measured')
        handles, labels = axis[0].get_legend_handles_labels()
        axis[0].legend(handles, labels, bbox_to_anchor=(0.25, 0.65))
        axis[0].set_xlabel("Date")
        axis[0].set_ylabel("Streamflow, in cfs")
        axis[0].set_ylim(0, 400)

        plt.xlabel("Date")
        plt.ylabel("Streamflow, in cfs")
        plt.ylim(0, 400)

    with styles.USGSMap():

        axis[1].set_xlabel("Date")
        axis[1].set_ylabel("Flow Components, in cfs")
        axis[1].set_yscale("log")
        plt.xlabel("Date")
        plt.ylabel("Flow Components, in cfs")
        plt.yscale("log")
        plt.ylim(1.0e-3, 1.0e4)
        axis[1].plot(stats.Date, stats.basin_ssflow_cfs_1, color='r', linewidth=1.5, label='Interflow')
        axis[1].plot(stats.Date, gw_seepage, color='purple', linewidth=1.5, label='Groundwater seepage')
        axis[1].plot(stats.Date, stats.basin_sroff_cfs_1, color='y', linewidth=1.5, label='Hortonian runoff')
        axis[1].plot(stats.Date, stats.basin_dunnian_1, color='b', linewidth=1.5, label='Dunnian runoff')
        handles, labels = axis[1].get_legend_handles_labels()
        axis[1].legend(handles, labels, bbox_to_anchor=(0.25, 0.65))
        plt.tight_layout()
        plt.show()
