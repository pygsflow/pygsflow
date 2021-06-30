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

    Examples
    --------

    load from shapefile

    >>> dis = gsflow.output.PrmsDiscretization.load_from_shapefile("myshape.shp")

    load from flopy grid

    >>> gsf = gsflow.GsflowModel.load_from_file("mycontrol.control")
    >>> ml = gsflow.mf
    >>> dis = gsflow.output.PrmsDiscretization.load_from_flopy(model=ml)

    """

    def __init__(self, xypts):
        self._xypts = xypts
        self._nhru = len(xypts)
        self._xcenters = None
        self._ycenters = None

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
    def x_hru_centers(self):
        """
        Returns
        -------
            np.ndarray of x-centers for each hru

        """
        if self._xcenters is None:
            self._get_centers()

        return self._xcenters

    @property
    def y_hru_centers(self):
        """
        Returns
        -------
            np.ndarray of y-centers for each hru

        """
        if self._ycenters is None:
            self._get_centers()

        return self._ycenters

    @property
    def nhru(self):
        """
        Returns
        -------
            number of hrus

        """
        return self._nhru

    @property
    def xypts(self):
        """
        Returns
        -------
            np.ndarray of xy-points for each hru

        """
        return self._xypts

    @property
    def extent(self):
        """
        Returns
        -------
            tuple (xmin, xmax, ymin, ymax)
        """
        return self._extent

    def _get_centers(self):
        """
        Method to get the mean center of the grid cell.
        This method is limited and will eventually need to
        be updated for complex shapes.

        """
        xc = []
        yc = []
        for hru in self._xypts:
            hru = np.array(hru)
            xc.append(np.mean(hru.T[0]))
            yc.append(np.mean(hru.T[1]))

        self._xcenters = xc
        self._ycenters = yc

    def get_hru_points(self, hru):
        """
        Get the x, y coordinate points for a hru

        Parameters
        ----------
        hru : int
            hru number

        Returns
        -------
            list of x, y coordinate points

        """
        return self._xypts[hru - 1]

    @staticmethod
    def load_from_flopy(model, xll=None, yll=None, rotation=None):
        """
        Method to load discretization from a flopy model

        Parameters
        ----------
        model : flopy.modflow.Modflow or gsflow.modflow.Modflow object
        xll : float, optional
            xoffset for modflow grid
        yll : float, optional
            yoffset for modflow grid
        rotation : float, optional
            rotation for modflow grid

        Returns
        -------
            PrmsDiscretization object

        """
        import flopy
        from gsflow.modflow import Modflow

        if not isinstance(model, flopy.modflow.Modflow) or not isinstance(
            model, Modflow
        ):
            raise ValueError(
                "Model must be a flopy.modflow.Modflow or "
                "gsflow.modflow.Modflow model"
            )

        mg = model.modelgrid
        if not isinstance(mg, flopy.discretization.grid.Grid):
            raise AssertionError("Cannot find flopy discretization")

        if (xll, yll, rotation) != (None, None, None):
            if rotation is None:
                rotation = 0.0
            mg.set_coord_info(xll=xll, yll=yll, angrot=rotation)

        xypts = []
        # create closed polygons for each hru from UL to LR
        xv = mg.xvertices
        yv = mg.yvertices
        for i in range(1, mg.shape[1] + 1):
            for j in range(1, mg.shape[2] + 1):
                t = [
                    (xv[i - 1, j - 1], yv[i - 1, j - 1]),
                    (xv[i - 1, j], yv[i - 1, j]),
                    (xv[i, j], yv[i, j]),
                    (xv[i, j - 1], yv[i, j - 1]),
                    (xv[i - 1, j - 1], yv[i - 1, j - 1]),
                ]
                xypts.append(t)

        return PrmsDiscretization(xypts)

    @staticmethod
    def load_from_shapefile(shp, hru_id="hru_id"):
        """
        Load method from a polygon shapefile. Shapefile
        must also have a hru field in the dbf file to sort
        the hydrologic reservior units properly.

        Parameters
        ----------
        shp : str or shapefile.Reader object
        hru_id : str
            field of hru id in shapefile

        Returns
        -------
            PrmsDiscretization object

        """
        if isinstance(shp, str):
            sf = shapefile.Reader(shp)
        else:
            sf = shp

        if sf.shapeType not in (shapefile.POLYGON, shapefile.POLYGONZ):
            err = (
                "A polygon shapefile must be supplied;"
                " current shape type is: {}".format(sf.shapeTypeName)
            )
            raise TypeError("Shapefile must ")

        hru_field = False
        for ix, field in enumerate(sf.fields):
            if hru_id == field[0].lower():
                # index the hru field
                hru_field = ix - 1
                break

        if not hru_field:
            err = (
                "A hru_id field must be present in the shapefile; hru_id must be "
                "from 1 to n_hru"
            )
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

        patches = [Polygon(hru, False, **kwargs) for hru in self._xypts]

        p = PatchCollection(patches, facecolors="None")
        ax.add_collection(p)

        extent = self.extent
        ax.set_xlim([extent[0], extent[1]])
        ax.set_ylim([extent[2], extent[3]])
        return ax
