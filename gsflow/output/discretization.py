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
    def __init__(self, xypts):
        self._xypts = xypts

    def get_hru_points(self, hru):
        """

        Parameters
        ----------
        hru

        Returns
        -------

        """
        return self._xypts(hru - 1)

    @staticmethod
    def load_from_flopy(model, xll=None, yll=None, rotation=None):
        """

        Parameters
        ----------
        model : flopy.modflow.Modflow or gsflow.modflow.Modflow object
        xll : float
            xoffset for modflow grid
        yll : float
            yoffset for modflow grid
        rotation : float
            rotation for modflow grid

        Returns
        -------
            PrmsDiscretization object
        """
        import flopy
        from gsflow.modflow import Modflow
        if not isinstance(model, flopy.modflow.Modflow) or \
                not isinstance(model, Modflow):
            raise ValueError("Model must be a flopy.modflow.Modflow or "
                             "gsflow.modflow.Modflow model")

        sr = model.sr
        if not isinstance(sr, flopy.utils.SpatialReference):
            raise AssertionError("Cannot find flopy discretization")

        if (xll, yll, rotation) != (None, None, None):
            if rotation is None:
                rotation = 0.
            sr.set_spatialreference(xll=xll, yll=yll, rotation=rotation)

        print('break')

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
            if "hru_id" == field[0].lower():
                # index the hru field
                hru_field = ix - 1
                break

        if not hru_field:
            err = "A hru_id field must be present in the shapefile; hru_id must be " \
                  "from 1 to n_hru"
            raise AssertionError(err)
        

        shapes = sf.shapes()

        xypts_dict = {}
        for ix, sh in enumerate(shapes):
            rec = sf.record(ix)
            hru_id = rec[hru_field]
            xypts_dict[hru_id] = sh

        xypts = []
        for key, values in sorted(xypts_dict.items()):
            xypts.append(values)

        return PrmsDiscretization(xypts)

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
