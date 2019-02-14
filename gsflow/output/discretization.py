import os
import sys
import numpy as np
import shapefile
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt


class PrmsDiscretization(object):
    """
    The PrmsDiscretization object is for plotting hru's and
    hru data from Prms on matplotlib plots.

    Could also be used for exporting shape files of data

    Parameters
    ----------
    xypts : list
        list of hru [(x1, y1)....(xn, yn)] points
        3 dimensional

    """
    def __init__(self, xypts):
        self._xypts = xypts

        xmin, xmax, ymin, ymax = (None, None, None, None)
        for hru in xypts:
            hru = np.array(hru)
            if xmin is None:
                xmin = np.min(hru.T[0])
                xmax = np.max(hru.T[0])
                ymin = np.min(hru.T[1])
                ymax = np.max(hru.T[1])

            else:
                if np.min(hru.T[0]) < xmin:
                    xmin = np.min(hru.T[0])
                if np.max(hru.T[0]) > xmax:
                    xmax = np.max(hru.T[0])
                if np.min(hru.T[1]) < ymin:
                    ymin = np.min(hru.T[1])
                if np.max(hru.T[1]) > ymax:
                    ymax = np.max(hru.T[1])

        self._extent = (xmin, xmax, ymin, ymax)

    @property
    def xypts(self):
        return self._xypts

    @property
    def extent(self):
        return self._extent

    def get_hru_points(self, hru):
        """

        Parameters
        ----------
        hru

        Returns
        -------

        """
        return self._xypts[hru - 1]

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

        xypts = []
        # create closed polygons for each hru from UL to LR
        for i in range(1, sr.xgrid.shape[0]):
            for j in range(1, sr.xgrid.shape[1]):
                t = [(sr.xgrid[i - 1, j - 1], sr.ygrid[i - 1, j - 1]),
                     (sr.xgrid[i - 1, j], sr.ygrid[i - 1, j]),
                     (sr.xgrid[i, j], sr.ygrid[i, j]),
                     (sr.xgrid[i, j - 1], sr.ygrid[i, j - 1]),
                     (sr.xgrid[i - 1, j - 1], sr.ygrid[i - 1, j - 1])]
                xypts.append(t)

        return PrmsDiscretization(xypts)

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
            xypts_dict[hru_id] = sh.points

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
            only for Polygon patches

        Returns
        -------
            matplotlib.pyplot.axes object

        """
        if ax is None:
            ax = plt.gca()

        if "color" not in kwargs:
            kwargs["facecolor"] = "None"

        patches = []
        for hru in self._xypts:
            polygon = Polygon(hru, False, **kwargs)
            patches.append(polygon)

        p = PatchCollection(patches, facecolors="None")
        ax.add_collection(p)

        extent = self.extent
        ax.set_xlim([extent[0], extent[1]])
        ax.set_ylim([extent[2], extent[3]])
        return ax