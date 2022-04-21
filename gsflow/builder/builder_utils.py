import os
import numpy as np
import pandas as pd
from ..prms import ParameterRecord


def build_lut(f, dtype=int):
    """
    Method to load a remap file into a lut

    Parameters
    ----------
    f : str
        file name
    dtype : type
        int or float type, defaults to int

    Returns
    -------
        dict
    """
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


def covtype(cov_resampled, lut):
    """
    Method to associate resampled vegitation type to PRMS covtype

    Parameters
    ----------
    cov_resampled : np.ndarray

    lut : dict

    Returns
    -------
        gsflow.prms.ParameterRecord object
    """
    values = np.array([lut[i] for i in cov_resampled.ravel()], dtype=int)
    record = ParameterRecord(
        "cov_type",
        values,
        dimensions=[
            ["nhru", values.size],
        ],
        datatype=1,
    )
    return record


def root_depth(cov_resampled, lut):
    """
    Method to get a rooting depth from cover type. This is an intermediate
    array used to calculate other parameters

    Parameters
    ----------
    cov_resampled : np.ndarray

    lut : d

    Returns
    -------
        np.ndarray of rooting depths
    """
    values = np.array([lut[i] for i in cov_resampled.ravel()], dtype=float)
    return values


def covden_sum(cov_resampled, lut):
    """
    Method to associate resampled vegitation coverage to PRMS covden_sum

    Parameters
    ----------
    cov_resampled : np.ndarray

    lut : dict

    Returns
    -------
        gsflow.prms.ParameterRecord object
    """
    covden = np.array([lut[i] for i in cov_resampled.ravel()], dtype=float)
    covden = covden / 100.0
    covden = ParameterRecord(
        "covden_sum",
        covden,
        dimensions=[
            ["nhru", covden.size],
        ],
        datatype=2,
    )
    return covden


def covden_win(cov_resampled, lut):
    """
    Method to associate resampled vegitation coverage to PRMS covden_win

    Parameters
    ----------
    cov_resampled : np.ndarray

    lut : dict

    Returns
    -------
        gsflow.prms.ParameterRecord object
    """
    covden = covden_sum(cov_resampled, lut)
    covden.name = "covden_win"
    return covden


def rad_trncf(covden_win):
    """
    Method to calculate rad_trncf from the covden_win parameter values

    Parameters
    ----------
    covden_win : np.ndarray

    Returns
    -------
        gsflow.prms.ParameterRecord object
    """

    values = 0.9917 * np.exp(-2.7557 * covden_win)
    values[values < 1e-5] = 0
    record = ParameterRecord(
        "rad_trncf",
        values,
        dimensions=[
            ["nhru", values.size],
        ],
        datatype=2,
    )
    return record


def snow_intcp(cov_resampled, lut):
    """
    Method to associate resampled vegitation type to PRMS snow_intcp

    Parameters
    ----------
    cov_resampled : np.ndarray

    lut : dict

    Returns
    -------
        gsflow.prms.ParameterRecord object
    """
    values = np.array([lut[i] for i in cov_resampled.ravel()], dtype=float)
    record = ParameterRecord(
        "snow_intcp",
        values,
        dimensions=[
            ["nhru", values.size],
        ],
        datatype=2,
    )
    return record


def srain_intcp(cov_resampled, lut):
    """
    Method to associate resampled vegitation type to PRMS srain_intcp

    Parameters
    ----------
    cov_resampled : np.ndarray

    lut : dict

    Returns
    -------
        gsflow.prms.ParameterRecord object
    """
    record = snow_intcp(cov_resampled, lut)
    record.name = "srain_intcp"
    return record


def wrain_intcp(cov_resampled, lut):
    """
    Method to associate resampled vegitation type to PRMS wrain_intcp

    Parameters
    ----------
    cov_resampled : np.ndarray

    lut : dict

    Returns
    -------
        gsflow.prms.ParameterRecord object
    """
    record = snow_intcp(cov_resampled, lut)
    record.name = "wrain_intcp"
    return record


def soil_type(clay, sand):
    """
    Method to determine soil type from resampled clay and sand percentages

    Parameters
    ----------
    clay : np.ndarray
    sand : np.ndarray

    Returns
    -------
        gsflow.prms.ParameterRecord
    """
    values = np.zeros(clay.size, dtype=int)
    values = np.where(clay.ravel() >= 0.40, 3, values)
    values = np.where(sand.ravel() >= 0.50, 1, values)
    values[values == 0] = 2

    record = ParameterRecord(
        "soil_type",
        values,
        dimensions=[
            ["nhru", values.size],
        ],
        datatype=1,
    )
    return record


def soil_moist_max(awc, root_depth):
    """
    Method to calculate the soil_moist_max parameter

    Parameters
    ----------
    awc : np.ndarray
        available water content
    root_depth : np.ndarray

    Returns
    -------
        gsflow.prms.ParameterRecord object

    """
    values = awc.ravel() * root_depth.ravel()
    values[values < 3] = 3
    values[values > 20] = 20
    record = ParameterRecord(
        "soil_moist_max",
        values,
        dimensions=[
            ["nhru", values.size],
        ],
        datatype=2,
    )
    return record


def soil_moist_init(soil_moist_max, factor=0.1):
    """
    Method to calculate the soil_rech_init parameter

    Parameters
    ----------
    soil_moist_max : np.ndarray
        maximum soil moisture water content
    factor : float
        scaling factor

    Returns
    -------
        gsflow.prms.ParameterRecord object

    """
    values = soil_moist_max * factor
    record = ParameterRecord(
        "soil_moist_init",
        values,
        dimensions=[
            ["nhru", values.size],
        ],
        datatype=2,
    )
    return record


def soil_rech_max(awc, root_depth, max_depth=18):
    """
    Method to calculate the soil_rech_max parameter

    Parameters
    ----------
    awc : np.ndarray
        available water content
    root_depth : np.ndarray
    max_depth : float
        maximum rooting depth to calculate recharge (default 18)

    Returns
    -------
        gsflow.prms.ParameterRecord object

    """
    values = np.where(
        root_depth.ravel() > max_depth,
        awc.ravel() * max_depth,
        awc.ravel() * root_depth.ravel(),
    )

    values[values < 1e-05] = 1e-05
    values[values > 20] = 20

    record = ParameterRecord(
        "soil_rechr_max",
        values,
        dimensions=[
            ["nhru", values.size],
        ],
        datatype=2,
    )
    return record


def soil_rech_init(soil_rech_max, factor=0.1):
    """
    Method to calculate the soil_rech_init parameter

    Parameters
    ----------
    soil_rech_max : np.ndarray
        maximum recharge water content
    factor : float
        scaling factor

    Returns
    -------
        gsflow.prms.ParameterRecord object

    """
    record = soil_moist_init(soil_rech_max, factor)
    record.name = "soil_rechr_init"
    return record


def ssr2gw_rate(ksat, sand, soil_moist_max):
    """
    Method to calculate the ssr2gw_rate from ksat (inches/day) and % sand

    Parameters
    ----------
    ksat : np.ndarray
    sand : np.ndarray
    soil_moist_max : np.ndarray

    Return
    ------
        gsflow.prms.ParameterRecord object
    """
    values = ksat.ravel() / (sand.ravel() * soil_moist_max.ravel())

    values[np.isnan(values)] = 1e-04
    values[values > 999] = 999
    values[values < 1e-04] = 1e-04

    record = ParameterRecord(
        "ssr2gw_rate",
        values,
        dimensions=[
            ["nssr", values.size],
        ],
        datatype=2,
    )
    return record


def ssr2gw_exp(nhru):
    """
    Method to calculate the ssr2gw_sq for PRMS

    Parameters
    ----------
    nhru : int

    Return
    ------
        gsflow.prms.ParameterRecord object
    """
    values = np.zeros((nhru,))
    record = ParameterRecord(
        "ssr2gw_sq",
        values,
        dimensions=[
            ["nssr", nhru],
        ],
        datatype=2,
    )
    return record


def slowcoef_lin(ksat, aspect, dx, dy):
    """
    Method to calculate the slowcoef_lin parameter in prms

    Parameters
    ----------
    ksat : np.ndarray
    aspect : np.ndarray
    dx : float
        grid dx dimension for a single cell
    dy : float
        grid dy dimension for a single cell

    Returns
    -------
        gsflow.prms.ParameterRecord object
    """
    hypot = np.sqrt(dx**2 + dy**2)
    hru_len = np.where((aspect == 90) | (aspect == 270), dx, dy)
    hru_len = np.where(
        (aspect != 0) & (aspect != 90) & (aspect != 180) & (aspect != 270),
        hypot,
        hru_len,
    )

    values = (
        0.1 * np.abs(ksat.ravel() * np.sin(aspect * np.pi / 180.0)) / hru_len
    )

    values[np.isnan(values)] = 0
    values[values < 1e-05] = 0

    record = ParameterRecord(
        "slowcoef_lin",
        values,
        dimensions=[
            ["nhru", values.size],
        ],
        datatype=2,
    )
    return record


def slowcoef_sq(ksat, aspect, sand, soil_moist_max, dx, dy):
    """
    Method to calculate the slowcoef_lin parameter in prms

    Parameters
    ----------
    ksat : np.ndarray
    aspect : np.ndarray
    sand : np.ndarray
        fraction of sand in soil
    soil_moist_max : np.ndarray
        soil moist max
    dx : float
        grid dx dimension for a single cell
    dy : float
        grid dy dimension for a single cell

    Returns
    -------
        gsflow.prms.ParameterRecord object
    """
    hypot = np.sqrt(dx**2 + dy**2)
    hru_len = np.where((aspect == 90) | (aspect == 270), dx, dy)
    hru_len = np.where(
        (aspect != 0) & (aspect != 90) & (aspect != 180) & (aspect != 270),
        hypot,
        hru_len,
    )

    values = 0.9 * np.abs(
        (ksat.ravel() * np.sin(aspect * np.pi / 180.0))
        / (hru_len * (sand.ravel() * soil_moist_max.ravel()))
    )

    values[np.isnan(values)] = 0
    values[values < 1e-05] = 0

    record = ParameterRecord(
        "slowcoef_sq",
        values,
        dimensions=[
            ["nhru", values.size],
        ],
        datatype=2,
    )
    return record


def hru_percent_imperv(resampled_imperv):
    """
    Method to set hru_percent_imperv parameter

    Parameters
    ----------
    resampled_imperv : np.ndarray

    Returns
    -------
        gsflow.prms.ParameterRecord object
    """
    record = ParameterRecord(
        "hru_percent_imperv",
        resampled_imperv.ravel(),
        dimensions=[
            ["nhru", resampled_imperv.size],
        ],
        datatype=2,
    )
    return record


def carea_max(resampled_imperv):
    """
    Method to set carea_max parameter

    Parameters
    ----------
    resampled_imperv : np.ndarray

    Returns
    -------
        gsflow.prms.ParameterRecord object
    """
    record = hru_percent_imperv(resampled_imperv)
    record.name = "carea_max"
    return record


def rain_adj(resampled_ppt, mean_monthly_sta):
    """
    Method to calculate rain_adj from a single station

    Parameters
    ----------
    resampled_ppt : np.ndarray

    mean_monthly_sta : np.ndarray

    Returns
    -------
        gsflow.prms.ParameterRecord object
    """
    mean_monthly_sta.shape = (mean_monthly_sta.size, 1)
    values = resampled_ppt / mean_monthly_sta
    values[values > 10] = 10
    values[values < 0.5] = 0.5

    record = ParameterRecord(
        "rain_adj",
        values.ravel(),
        dimensions=[
            ["nhru", values.shape[-1]],
            ["nmonths", mean_monthly_sta.size],
        ],
        datatype=2,
    )
    return record


def snow_adj(resampled_ppt, mean_monthly_sta):
    """
    Method to calculate rain_adj from a single station

    Parameters
    ----------
    resampled_ppt : np.ndarray

    mean_monthly_sta : np.ndarray

    Returns
    -------
        gsflow.prms.ParameterRecord object
    """
    record = rain_adj(resampled_ppt, mean_monthly_sta)
    record.name = "snow_adj"
    return record


def tmax_lapse(lapse_rates):
    """
    Method to calculate tmax_lapse parameter

    Parameters
    ----------
    lapse_rates : np.ndarray

    Returns
    -------
        gsflow.prms.ParameterRecord object
    """

    record = ParameterRecord(
        "tmax_lapse",
        lapse_rates.ravel(),
        dimensions=[["nmonths", len(lapse_rates)]],
        datatype=2,
    )
    return record


def tmin_lapse(lapse_rates):
    """
    Method to calculate tmax_lapse parameter

    Parameters
    ----------
    lapse_rates : np.ndarray

    Returns
    -------
        gsflow.prms.ParameterRecord object
    """
    record = tmax_lapse(lapse_rates)
    record.name = "tmin_lapse"
    return record


def tmax_adj(nhru, nmonths=12):
    """
    Method to calculate tmax_adj from a single station
    Parameters
    ----------
    nhru : int
    nmonths : np.ndarray
    Returns
    -------
        gsflow.prms.ParameterRecord object
    """
    values = np.zeros((nhru * nmonths))
    record = ParameterRecord(
        "tmax_adj",
        values,
        dimensions=[["nhru", nhru], ["nmonths", nmonths]],
        datatype=2,
    )
    return record


def tmin_adj(nhru, nmonths=12):
    """
    Method to calculate tmin_adj from a single station
    Parameters
    ----------
    nhru : int
    nmonths : np.ndarray
    Returns
    -------
        gsflow.prms.ParameterRecord object
    """
    values = np.zeros((nhru * nmonths))
    record = ParameterRecord(
        "tmin_adj",
        values,
        dimensions=[["nhru", nhru], ["nmonths", nmonths]],
        datatype=2,
    )
    return record


def calculate_jensen_haise(dem, tmin_mean, tmax_mean):
    """
    Method to calculate the Jensen Haise coefficient for PRMS

    Parameters
    ----------
    dem : np.ndarray
    tmin_mean : np.ndarray
        list of monthly mean values for tmin
    tmax_mean : np.ndarray
        array of monthly mean values for tmax

    Returns
    -------
        gsflow.prms.ParameterRecord object
    """
    nhru = dem.size
    dem = np.ravel(dem)
    idx = np.where(tmax_mean == np.max(tmax_mean))[0]
    tmax = tmax_mean[idx[0]]
    tmin = tmin_mean[idx[0]]
    jh_coef = 27.5 - 0.25 * (ea(tmax) - ea(tmin)) - (dem / 1000.0)
    jh_coef = ParameterRecord(
        "jh_coef_hru", jh_coef, dimensions=[["nhru", nhru]], datatype=2
    )
    return jh_coef


def ea(tempc):
    """
    Method to calculate ea for Jensen Haise equation

    Parameters
    ----------
    tempc : np.ndarray
        monthly mean temperatures

    Returns
    -------
        np.ndarray
    """
    tea = 6.1078 * np.exp((17.269 * tempc) / (tempc + 237.3))
    return tea


def d8_to_hru_aspect(flow_directions):
    """
    Method to get cell aspect from a D8 flow accumulation array

    Parameters
    ----------
    flow_directions : np.ndarray
        numpy array of flow directions

    Returns
    -------
        gsflow.prms.ParameterRecord object of aspect angles (0 - 360 deg)
    """
    d8_to_aspect = {
        -2: 0,
        -1: 0,
        32: 315,
        64: 0,
        128: 45,
        16: 270,
        1: 90,
        8: 225,
        4: 180,
        2: 135,
    }

    aspect = np.array(
        [d8_to_aspect[i] for i in flow_directions.ravel()], dtype=float
    )
    aspect = ParameterRecord(
        "hru_aspect",
        aspect,
        dimensions=[
            ["nhru", aspect.size],
        ],
        datatype=2,
    )
    return aspect


def d8_to_hru_slope(flow_directions, dem, xcenters, ycenters):
    """
    Method to calculate slopes from a D8 flow accumulation array

    Parameters
    ----------
    flow_directions : np.ndarray
        2d numpy array of flow directions
    dem : np.ndarray
        2d numpy array of dem elevations
    xcenters : np.ndarray
        modelgrid x-centers
    ycenters : np.ndarray
        modelgrid y-centers

    Returns
    -------
        gsflow.prms.ParameterRecord object
    """
    ioff = {
        32: -1,
        64: -1,
        128: -1,
        16: 0,
        1: 0,
        8: 1,
        4: 1,
        2: 1,
        -1: 0,
        -2: 0,
    }
    joff = {
        32: -1,
        64: 0,
        128: 1,
        16: -1,
        1: 1,
        8: -1,
        4: 0,
        2: 1,
        -1: 0,
        -2: 0,
    }

    slope = np.zeros(flow_directions.shape, dtype=float)
    for i, row in enumerate(flow_directions):
        for j, fdir in enumerate(row):
            i1 = i + ioff[fdir]
            j1 = j + joff[fdir]

            asq = (xcenters[i, j] - xcenters[i1, j1]) ** 2
            bsq = (ycenters[i, j] - ycenters[i1, j1]) ** 2
            dist = np.sqrt(asq + bsq)
            delta_elev = dem[i, j] - dem[i1, j1]
            if dist == 0:
                continue
            else:
                slope[i, j] = delta_elev / dist

    slope[slope < 1e-05] = 0
    slope[slope > 10] = 10

    slope = ParameterRecord(
        "hru_slope",
        slope.ravel(),
        dimensions=[
            ["nhru", slope.size],
        ],
        datatype=2,
    )
    return slope


def add_prms_date_columns_to_df(df, column):
    """
    Method to add date columns for PRMS to a pandas dataframe

    Parameters
    ----------
    df : pd.DataFrame
    column : str
        date column name


    Returns
    -------
        pd.Dataframe
    """
    df["Year"] = pd.DatetimeIndex(df[column]).year
    df["Month"] = pd.DatetimeIndex(df[column]).month
    df["Day"] = pd.DatetimeIndex(df[column]).day
    df["Hour"] = pd.DatetimeIndex(df[column]).hour
    df["Minute"] = pd.DatetimeIndex(df[column]).minute
    df["Second"] = pd.DatetimeIndex(df[column]).second

    return df


def get_mean_monthly_from_df(df, column, nodataval=-999, temperature=False):
    """
    Method to calculate mean monthly values from a dataframe

    Parameters
    ----------
    df : pd.Dataframe
    column : str
        data column name
    nodataval : float
        no data value (sets to nan when calculating
    temperature : bool
        boolean flag to indicate this is temperature (does not sum())

    Returns
    -------
        np.ndarray
    """
    means = []
    for month in range(1, 13):
        tdf = df[df.Month == month]
        tdf[tdf[column] == nodataval] = np.nan
        if temperature:
            tdf = tdf.groupby(by=["Year"], as_index=False)[column].mean()
        else:
            tdf = tdf.groupby(by=["Year"], as_index=False)[column].sum()
        means.append(tdf[column].mean())

    return np.array(means)


def fahrenheit_to_celsius(deg_f, nodataval=-999.0):
    """
    Method to convert fahrenheit to celsius

    Parameters
    ----------
    deg_f : float, np.ndarray
    nodataval : float
        value of nodata instances

    Returns
    -------
        float, np.ndarray
    """
    if isinstance(deg_f, (int, float)):
        deg_f = np.array([deg_f])
    elif isinstance(deg_f, list):
        deg_f = np.array(deg_f)

    deg_c = np.where(deg_f != nodataval, (deg_f - 32.0) * (5.0 / 9.0), -999.0)
    if deg_c.size == 1:
        return deg_c[0]
    else:
        return deg_c
