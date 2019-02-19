import pandas as pd
import datetime
from pyprms import prms_py
from gispy import shp
import copy
import matplotlib.pyplot as plt
from Pytoolbox import matlab
from pyprms import prms_plots
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable


"""
@ comment JL

This code is deprecated by the output.prms_plot module
consider removing from project or scattering into built in methods
"""
class Plotter(object):
    def __init__(self, hru_shp = None, proj = None):
        self.hru_param = hru_shp
        self.proj = proj


    def plot2D(self, param):
        plt.figure()
        para = self.proj.prms_parameters['Parameters'][param]
        nrow = np.max(self.hru_param['HRU_ROW'].values)
        ncol = np.max(self.hru_param['HRU_COL'].values)
        self.nhrus = nrow * ncol
        if para[0] == 1:
            if int(para[2]) == self.nhrus:
                var2D = para[4]
                plt.imshow(np.reshape(var2D, (nrow, ncol)), interpolation='none')
                plt.title(param)
                plt.xlabel('X')
                plt.ylabel('Y')
                plt.colorbar()
                plt.show()

            else:
                pass  # generate an error
        elif para[0]  == 2: # two-Dimensional data
            no_dim1 =self.nhrus
            no_dim2 = 12
            if no_dim1 == self.nhrus:
                allmons = para[4].reshape(no_dim2, no_dim1 )
                allmons = allmons.transpose()
                for mon in np.arange(0, 12, 1):
                    plt.subplot(2, 6, (mon + 1))
                    var2D = allmons[:, mon]
                    im = plt.imshow(np.reshape(var2D, (nrow, ncol)), interpolation='none')
                    plt.title(param + '_' + str(mon + 1))
                    plt.xlabel('X')
                    plt.ylabel('Y')
                    ax = plt.gca()
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes("right", size="5%", pad=0.05)
                    plt.colorbar(im, cax=cax)
                    plt.tight_layout()
                    plt.show()

        pass

        pass
