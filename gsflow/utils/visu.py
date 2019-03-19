"""
Visulization of PRMS and MODFLOW input and output files using vtk file format

"""
import os
import warnings
warnings.simplefilter('always', UserWarning)

try:
    import pyevtk
except:
    err = "pyevtk is not installed. Visualization files cannot be generated"
    warnings.warn(err, UserWarning)


def gsflow_to_vtk(control_file = None, only = []):


    pass

def prms():
    pass

def modflow_to_vtk(name_file = None):
    pass

def shp_to_vtk(shapefile = None):
    pass


