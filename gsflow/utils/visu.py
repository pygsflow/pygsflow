"""
Visulization of PRMS and MODFLOW input and output files using vtk file format

"""
import os
import warnings
warnings.simplefilter('always', UserWarning)

try:
    import pyevtk
except:
    pyevtk = None


def gsflow_to_vtk(control_file = None, only = []):
    if pyevtk is None:
        err = "pyevtk is not installed. Visualization files cannot be generated"
        warnings.warn(err, UserWarning)
        # or we can raise an import error here and direct the user to install pyevtk

    pass

def prms():
    pass

def modflow_to_vtk(name_file = None):
    pass

def shp_to_vtk(shapefile = None):
    pass


