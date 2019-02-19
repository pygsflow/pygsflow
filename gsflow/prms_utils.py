import os, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from dbfpy import dbf
"""
Author: Ayman H. Alzraiee
Date : 4/14/2017
Version: 0.00
"""

"""
@ comment JL

utils of this sort should not be grouped into a class.
As this project grows it is going to become unwieldly in a single class.

I suggest thinking about refactoring to IO utilities module (contains IO related functions)
data utilities modult (contains common functions), shapefile utilities module (...), 
plotting utilities module (plotting helpers that are seperate from plotting classes),
export utilities, import utilities, etc.

It makes more sense to group by their use cases.
"""

class utils(object):
    def __init__(self):
        pass

    def write_to_exiting_shapefile(self,f_in,f_out):
        """
        :param
            f_in: is file for a shapefile layer that contains the model grid
            f_out: is the augmented shaefile


        :return:
        """

        pass

    def read_param_from_shapefile(self, f_in):
        """

        :param f_in:
        :return: return a dictionary with parameters

        """
        pass
    def write_prms_project_excel(self, f_control):
        """

        :param f_control: is the name of the control file
        :return:
        """

        pass

    """
    @ comment JL
    
    modflow gridfile is a specific file created by 
    modflow 6. 
    
    Rename this to read_modflow_grid_from_shapefile?
    
    FloPy also has utilities to do this. Does it make
    sense to leverage those instead. Also pyshp seems 
    like a better package to use for this. dbfpy is not actively 
    maintained.
    """
    def read_modflow_gridfile(self, shpfile, id_name):
        """

        :param shpfile: a polygon shapefile for the grid that includes cells modflow id
        :return: id_name : field name for cell id

        """
        base_file = os.path.dirname(shpfile)
        file_parts = os.path.basename(shpfile).split('.')
        dbf_file = os.path.join(base_file, file_parts[0] + ".dbf")

        db = dbf.Dbf(dbf_file)
        mfid = [rec[id_name] for rec in db]

        db.addField('')

        db.close()
        pass

    def add_prms_par_to_shapefile(self, shpfile, id_name, data):
        pass
