import os, sys
import numpy as np
import matplotlib.pyplot as plt
import gsflopy as gsf
from mpl_toolkits.axes_grid1 import make_axes_locatable

"""
@ comment JL

this code is somewhat deprecated by the output.prms_plot()
module. Some of this code could be scattered to built in's though,
so I haven't removed it yet.

"""
class Prms_plots(object):
    def __init__(self, gsf):
        self.gsf = gsf
        dis = gsf.modflow.dis
        delc = dis.delc.array
        delr = dis.delr.array
        xv =  np.cumsum(delc)- delc # add origin here
        yv = np.cumsum(delr)- delr
        icol = np.arange(1,gsf.modflow.dis.nrow+1)
        irow = np.arange(1,gsf.modflow.dis.ncol+1)
        self.icol, self.irow = np.meshgrid(irow, icol)
        self.X, self.Y = np.meshgrid(xv,yv)
        self.mfid = (self.irow - 1) * gsf.modflow.dis.ncol + self.icol
        self.gvr_cell_id = self.gsf.prms.parameters['gvr_cell_id']
        self.nhrus = self.gsf.prms.parameters.dimensions['nhru']

    def plot2D(self, param):

        plt.figure()

        nrow = self.gsf.modflow.dis.nrow
        ncol = self.gsf.modflow.dis.ncol
        if len(self.gsf.prms.parameters[param].shape) == 1:
            if len(self.gsf.prms.parameters[param]) == self.nhrus:
                var2D = self.gsf.prms.parameters[param]
                hru_id = np.zeros((self.mfid.size), dtype=float) + np.nan
                hru_id[self.gvr_cell_id - 1] = var2D
                plt.imshow(np.reshape(hru_id, (nrow, ncol)), interpolation='none')
                plt.title(param)
                plt.xlabel('X')
                plt.ylabel('Y')
                plt.colorbar()
                plt.show()

            else:
                pass # generate an error

        elif len(self.gsf.prms.parameters[param].shape) == 2:

            if self.gsf.prms.parameters[param].shape[1] == self.nhrus:
                allmons = self.gsf.prms.parameters[param]

                for mon in np.arange(0,12,1):
                    plt.subplot(2, 6,(mon+1))
                    var2D = allmons[mon,:]
                    hru_id = np.zeros((self.mfid.size), dtype=float) + np.nan
                    hru_id[self.gvr_cell_id - 1] = var2D
                    im = plt.imshow(np.reshape(hru_id, (nrow, ncol)), interpolation='none')
                    plt.title(param + '_'+ str(mon+1))
                    plt.xlabel('X')
                    plt.ylabel('Y')
                    ax = plt.gca()
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes("right", size="5%", pad=0.05)
                    plt.colorbar(im, cax=cax)
                    plt.tight_layout()
                    plt.show()




        pass

    def plot_ts(self):
        pass


