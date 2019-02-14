import os, sys
import numpy as np
import matplotlib.pyplot as plt
import gsflopy as gsf
from mpl_toolkits.axes_grid1 import make_axes_locatable

"""
@ comment JL

lots (specifically, why are you creating a grid when flopy already does this?)
dis.sr.xgrid, dis.sr.ygrid is the same as your self.X and self.Y

mfid can easily be made using np.arange(dis.nrow * dis.ncol).reshape(dis.nrow, dis.ncol)

consistant input parameter names must be used across a project. Use gsflow or model....

what if the user wants to use a PRMS model with a custom discretization? I think we should
support that.

Maybe input for PrmsPlots should be __init__(self, prms, discretization=None)

As usual CamelCase all classes!!!!
Needs docstrings

I have timeseries plotting code that can be added to this....

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
        """
        @ comment JL
        
        this method needs some work. If a user wants to layer things
        it is not possible. 
        
        **kwargs need to be able to be passed.
        
        should return an axis object so the user can do more formatting
        
        user should be able to pass in an axis object to the plot2D function.
        if none ax = plt.gca()
        
        much to do.... 
        """

        nrow = self.gsf.modflow.dis.nrow
        ncol = self.gsf.modflow.dis.ncol
        if len(self.gsf.prms.parameters[param].shape) == 1:
            if len(self.gsf.prms.parameters[param]) == self.nhrus:
                var2D = self.gsf.prms.parameters[param]
                """
                @ comment JL
                
                using an uninstantiated array and then instatiating with np.nan
                is more efficient for large data sets
                
                a = np.empty((size))
                a[:] = np.nan
                
                it seems like you're using the mfid only for it's size
                
                instead of storing an array you could just store mfid.size in init
                and then that uses less memory, especially for large models (ex. russian river!)
                """
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
            """
            @ comment JL
            
            This whole plotting function should be reconsidered.
            
            This is something that belongs in plot util and attached to the 
            params class as a built in. 
            
            Plotting classes should be flexible and let the user decide
            which portions of the data they want to represent!
            
            As stated before, this is better suited as a built-in method
            attached to a specific class, ex. PrmsParameters() 
            """

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


