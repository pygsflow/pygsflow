"""
rasters.py includes code to work with raster datasets,
particularly to sample raster using a PrmsDiscretization
object and then use it as input data
"""
import numpy as np


class Rasters(object):
    """
    Raster object which allows the user
    to snap a raster to a grid.

    Parameters
    ----------
        raster : osgeo.gdal.Dataset object

    """

    def __init__(self, raster):
        self.raster = raster
        self._band_data = None
        self._xpoints = None
        self._ypoints = None
        self._xypoints = None

    @property
    def extent(self):
        """
        Returns
        -------
            (xmin, xmax, ymin, ymax)

        """
        return (
            self.raster.bounds.left,
            self.raster.bounds.right,
            self.raster.bounds.bottom,
            self.raster.bounds.top,
        )

    @property
    def xpoints(self):
        """
        Returns
        -------
            Cell centered x points for raster values

        """
        if self._band_data is not None:
            if self._xpoints is None:
                xmin, xmax = self.extent[0:2]
                ynum, xnum = self._band_data.shape
                t = np.linspace(xmin, xmax, xnum)
                self._xpoints = np.tile(t, (ynum, 1))

            return self._xpoints

    @property
    def ypoints(self):
        """
        Returns
        -------
            Cell centered y points for raster values

        """
        if self._band_data is not None:
            if self._ypoints is None:
                ymin, ymax = self.extent[2:4]
                ynum, xnum = self._band_data.shape
                t = np.linspace(ymin, ymax, ynum)
                self._ypoints = np.tile(t, (xnum, 1)).T

            return self._ypoints

    @property
    def xypoints(self):
        """
        Returns
        -------
            np.ndarray cell centered x, y points for raster values

        """
        if self._xypoints is None:
            if self.xpoints is None or self.ypoints is None:
                return

            self._xypoints = np.array([self.xpoints, self.ypoints])

        return self._xypoints

    @property
    def band_array(self):
        """
        Returns
        -------
            np.ndarray of the raster band

        """
        if self._band_data is not None:
            nodata = self.raster.nodata
            self._band_data[self._band_data == nodata] = np.nan
            return self._band_data

    def sample_discretization(self, prms_discretizaiton):
        """
        Method to sample the raster using the prms_discretization
        object cell centers

        Parameters
        ----------
        prms_discretizaiton : PrmsDiscretization object

        Returns
        -------
            np.ndarray

        """

        xy = list(
            zip(
                prms_discretizaiton.x_hru_centers,
                prms_discretizaiton.y_hru_centers,
            )
        )
        temp = np.array(list(self.raster.sample(xy))).flatten()
        temp[temp == self.raster.nodata] = np.nan
        return temp

    def set_raster_band(self, band):
        """
        Method to use GDAL to set the raster band
        for visualization

        Parameters
        ----------
        band : int
            raster band number

        """
        self._band_data = self.raster.read(band)

    @staticmethod
    def load(name):
        """

        Parameters
        ----------
        name : raster file name

        Returns
        -------
            Rasters object

        """
        import rasterio

        raster = rasterio.open(name)
        return Rasters(raster)
