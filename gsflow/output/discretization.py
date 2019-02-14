import os
import sys
import numpy as np
import shapefile


class PrmsDiscretization(object):
    """
    The PrmsDiscretization object is for plotting hru's and
    hru data from Prms on matplotlib plots.

    Could also be used for exporting shape files of data

    Parameters
    ----------


    """
    # todo: maybe make this xypts and then unravel
    def __init__(self, xypts):
        pass
        # self.xpts = xpts
        # self.ypts = ypts

        # if not isinstance(xpts, np.ndarray):
        #     self.xpts = np.array(xpts)

        # if not isinstance(ypts, np.ndarray):
        #     self.ypts = np.array(ypts)

    def get_hru_points(self, hru):
        """

        Parameters
        ----------
        hru

        Returns
        -------

        """
        pass

    @staticmethod
    def load_from_flopy(model):
        """

        Parameters
        ----------
        model

        Returns
        -------

        """
        pass

    @staticmethod
    def load_from_shapefile(shp):
        """
        Load method from a polygon shapefile. Shapefile
        must also have a hru field in the dbf file to sort
        the hydrologic reservior units properly.

        Parameters
        ----------
        shp : str or shapefile.Reader object

        Returns
        -------
            PrmsDiscretization object

        """
        if isinstance(shp, str):
            sf = shapefile.Reader(shp)
        else:
            sf = shp

        if sf.shapeType not in (shapefile.POLYGON,
                                shapefile.POLYGONZ):
            err = "A polygon shapefile must be supplied;" \
                  " current shape type is: {}".format(sf.shapeTypeName)
            raise TypeError("Shapefile must ")

        hru_field = False
        for ix, field in enumerate(sf.fields):
            if "hru" == field[0].lower():
                # index the hru field
                hru_field = ix
                break

        if not hru_field:
            err = "A hru field must be supplied; this is the hru order " \
                  "from 1 to n_hru"
        

        shapes = sf.shapes()
        print('break')

    def plot_discretization(self, ax=None, **kwargs):
        """
        Method to plot the PRMS discretization on a matplotlib.pyplot
        plot

        Parameters
        ----------
        ax : matplotlib.pyplot.axes
            if None, gets current working axes
        kwargs : matplotlib.pyplot keyword arguments

        Returns
        -------
            matplotlib.pyplot.axes object
        """
        pass
