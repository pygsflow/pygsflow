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