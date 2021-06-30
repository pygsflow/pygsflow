import numpy as np
import flopy as fp


def param2netcdf(f, model, parameter, **kwargs):
    """
    Method to add PRMS parameter to a netcdf
    file

    Parameters
    ----------
    f : str or fp.export.NetCdf
    model : fp.modflow.Modflow object
        fp.modflow.Modflow or gsflow.modflow.Modflow
    parameter : gsflow.prms.ParameterRecord object
        gsflow.prms.ParameterRecord object
    kwargs :
        keyword arguments

    """
    if parameter.section.lower() in ("dimensions",):
        return f

    shape = model.modelgrid.shape[1:]
    nhru = shape[0] * shape[1]

    if parameter.values.size == nhru:
        return param2d(f, model, parameter, **kwargs)

    elif parameter.values.size % nhru == 0:
        return param3d(f, model, parameter, **kwargs)

    else:
        return f


def param2d(f, model, parameter, **kwargs):
    """
    Method to add a 2D parameter to a netcdf file

    Parameters
    ----------
    f : str or fp.export.NetCdf
    model : fp.modflow.Modflow object
        fp.modflow.Modflow or gsflow.modflow.Modflow
    parameter : gsflow.prms.ParameterRecord object
        gsflow.prms.ParameterRecord object
    kwargs :
        keyword arguments

    """
    min_valid = kwargs.get("min_valid", -1.0e9)
    max_valid = kwargs.get("max_valid", 1.0e09)

    if isinstance(f, str):
        f = fp.export.NetCdf(f, model)

    var_name = parameter.name
    var_name = var_name.replace(" ", "_").lower()
    array = np.copy(parameter.values)
    nrow, ncol = model.modelgrid.shape[1:]
    array.shape = (nrow, ncol)

    with np.errstate(invalid="ignore"):
        if array.dtype not in (int, np.int32, np.int64):
            array[array <= min_valid] = np.NaN
            array[array >= max_valid] = np.NaN
            array[array == 0.0] = np.NaN
            mx = np.nanmax(array)
            mn = np.nanmax(array)

        else:
            mx = np.nanmax(array)
            mn = np.nanmin(array)

    array[np.isnan(array)] = f.fillvalue

    units = get_units(parameter.name)
    longname = get_longname(parameter.name)
    coordinates = "latitude longitude"

    attribs = {
        "long_name": longname,
        "units": units,
        "min": mn,
        "max": mx,
        "coordinates": coordinates,
    }

    var = f.create_variable(var_name, attribs, dimensions=("y", "x"))

    var[:] = array
    return f


def param3d(f, model, parameter, **kwargs):
    """
    Method to add a 3D parameter to a netcdf file

    Parameters
    ----------
    f : str or fp.export.NetCdf
    model : fp.modflow.Modflow object
        fp.modflow.Modflow or gsflow.modflow.Modflow
    parameter : gsflow.prms.ParameterRecord object
        gsflow.prms.ParameterRecord object
    kwargs :
        keyword arguments

    """
    from ..prms import Helper

    min_valid = kwargs.get("min_valid", -1.0e9)
    max_valid = kwargs.get("max_valid", 1.0e09)

    if isinstance(f, str):
        f = fp.export.NetCdf(f, model)

    if "month" not in f.nc.dimensions:
        f.nc.createDimension("month", 12)
        attribs = {
            "units": "",
            "standard_name": "month",
            "long_name": "PRMS month dimension",
            "positive": "down",
            "axis": "Z",
        }
        mo = f.create_variable("month", attribs, dimensions=("month",))
        mo[:] = np.arange(0, 12)

    var_name = parameter.name
    var_name = var_name.replace(" ", "_").lower()
    array = np.copy(parameter.values)
    nrow, ncol = model.modelgrid.shape[1:]
    array.shape = (-1, nrow, ncol)

    with np.errstate(invalid="ignore"):
        if array.dtype not in (int, np.int32, np.int64):
            array[array <= min_valid] = np.NaN
            array[array >= max_valid] = np.NaN
            array[array == 0.0] = np.NaN
            mx = np.nanmax(array)
            mn = np.nanmax(array)

        else:
            mx = np.nanmax(array)
            mn = np.nanmin(array)

    array[np.isnan(array)] = f.fillvalue

    units = get_units(parameter.name)
    longname = get_longname(parameter.name)
    coordinates = "month latitude longitude"

    attribs = {
        "long_name": longname,
        "units": units,
        "min": mn,
        "max": mx,
        "coordinates": coordinates,
    }

    var = f.create_variable(var_name, attribs, dimensions=("month", "y", "x"))

    var[:] = array
    return f


def get_longname(name):
    """
    Get parameter long name for netcdf

    Parameters
    ----------
    name : str
        parameter name

    Returns
    -------
        s : str

    """
    from ..prms import Helper

    help = Helper()

    if name in help.prms_parameter_names:
        return help.prms_parameter_names[name]["Descr"]
    elif name in help.prms_output_variables:
        return help.prms_output_variables[name]["Descr"]
    else:
        return name


def get_units(name):
    """
    Get parameter units for netcdf

    Parameters
    ----------
    name : str
        parameter name
    unit_base : str
        base unit (ex. ft, m, etc...)

    Returns
    -------
        s : str

    """
    from ..prms import Helper

    help = Helper()

    if name in help.prms_parameter_names:
        return help.prms_parameter_names[name]["Units"]

    elif name in help.prms_output_variables:
        return help.prms_output_variables[name]["Units"]

    else:
        return "unitless"
