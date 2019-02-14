import os
import warnings
warnings.simplefilter('always', UserWarning)


def get_file_abs(control_file=None, fn=None):
    """

    Parameters
    ----------
    control_file : str
        absoulute path to the control file
    fn : str
        path to file

    Returns
    -------
        abs_file : str
            absolute file path

    """
    fn = os.path.normpath(fn)
    if os.path.isabs(fn):
        return fn
    else:
        control_folder = os.path.dirname(control_file)
        abs_file = os.path.abspath(os.path.join(control_folder, fn))
        return abs_file


def _get_relative_path(control, fn):
    """
    If relative files are used, they should be relative to the control file

    Parameters
    ----------
    control
    fn

    Returns
    -------

    """
    control_file_abs = os.path.abspath(control)
    fn_abs = os.path.abspath(fn)
    # find common path
    rel_dir = os.path.relpath(os.path.dirname(fn), os.path.dirname(control_file_abs))
    rel_path = os.path.join(rel_dir + os.path.basename(fn))
    return rel_path


def find_parameter(name, parameters_list):
    """
    Utility method to loop through parameters to find values
    etc...

    Parameters
    ----------
    name : str
        parameter name
    parameters_list : list
        list of ParameterRecord objects or ControlRecord objects

    Returns
    -------
        rec : object
    """

    if len(parameters_list) > 0:
        for rec in parameters_list:
            if rec.name.lower() == name.lower():
                return rec
        return None

    else:
        err = "parameter_list is empty"
        warnings.warn(err, UserWarning)
        return None

