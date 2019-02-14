import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import numpy as np
from .discretization import PrmsDiscretization


class PrmsPlot(object):
    """
    Plotting class for PRMS data

    Parameters
    ----------
        prms_dis : PrmsDiscretization object, optional

        extent : tuple of floats, optional
            plot extent, if None PrmsPlot will get extent from PrmsDiscretization
    """
    def __init__(self, prms_dis=None, extent=None):

        if prms_dis is not None:
            if not isinstance(prms_dis, PrmsDiscretization):
                raise TypeError("prms_dis must be a PrmsDiscretization object")
            self.__maps = True

        else:
            self.__maps = False

        self.prms_dis = prms_dis

        self.extent = extent
        if self.extent is None:
            try:
                self.extent = prms_dis.extent
            except AttributeError:
                pass

    def plot_array(self, array, ax=None, masked_values=None, **kwargs):
        """
        General array plotter

        Parameters
        ----------
        array : array to plot.
            size must be equal to nhru

        ax : matplotlib.pyplot.axes

        kwargs : matplotlib keyword arguments

        Returns
        -------

        """
        if ax is None:
            plt.gca()

        if not isinstance(array, np.ndarray):
            array = np.array(array)

        if masked_values is not None:
            for mval in masked_values:
                array = np.ma.masked_equal(array, mval)

        return self.__plot_patches(array, ax=ax, **kwargs)

    def plot_parameter_map(self, parameter, ax=None, **kwargs):
        """

        Parameters
        ----------
        parameter
        ax
        kwargs

        Returns
        -------

        """
        if ax is None:
            ax = plt.gca()

        nhru = self.prms_dis.nhru


    def plot_parameter_timeseries(self, parameter, ax=None, **kwargs):
        """

        Parameters
        ----------
        parameter
        ax
        kwargs

        Returns
        -------

        """
        return

    def plot_data_timeseries(self, data, names, ax=None, **kwargs):
        """

        Parameters
        ----------
        data
        names
        ax
        kwargs

        Returns
        -------

        """
        return


    def plot_model_discretization(self, ax=None, **kwargs):
        """

        Parameters
        ----------
        ax
        kwargs

        Returns
        -------

        """
        if self.__maps:
            return self.prms_dis.plot_discretization(ax=ax, **kwargs)
        else:
            print("No discretization provided, cannot plot")

    def __plot_patches(self, array, ax, **kwargs):
        """
        General patch plotter, workhorse of the class

        Parameters
        ----------
        array : array to plot
        ax : matplotlib.pyplot.axes
        kwargs : matplotlib keyword arguments

        Returns
        -------
            matplotlib.pyplot.axes
        """
        if "color" not in kwargs:
            kwargs["facecolor"] = "None"

        if 'vmin' in kwargs:
            vmin = kwargs.pop('vmin')
        else:
            vmin = None

        if 'vmax' in kwargs:
            vmax = kwargs.pop('vmax')
        else:
            vmax = None

        patches = []
        for hru in self.prms_dis.xypts:
            polygon = Polygon(hru, False)
            patches.append(polygon)

        p = PatchCollection(patches)
        p.set_array(array)
        p.set_clim(vmin=vmin, vmax=vmax)
        p.set(**kwargs)

        ax.add_collection(p)

        extent = self.extent
        ax.set_xlim([extent[0], extent[1]])
        ax.set_ylim([extent[2], extent[3]])
        return ax