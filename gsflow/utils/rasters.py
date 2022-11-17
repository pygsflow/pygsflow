"""
rasters.py includes code to work with raster datasets,
particularly to sample raster using a PrmsDiscretization
object and then use it as input data
"""
import queue
import threading
import numpy as np
import flopy.utils
import warnings

warnings.simplefilter("ignore", RuntimeWarning)

try:
    from numba import jit, prange

    ENABLE_JIT = True
except ImportError:
    from .numba import jit

    ENABLE_JIT = False


class Raster:
    """
    Deprecated Raster class, use flopy.utils.Raster. This class is used to resample
    rasters and the overload allows for multiprocessing using the "ray"
    python library

    Parameters
    ----------
    array : np.ndarray
        a three dimensional array of raster values with dimensions
        defined by (raster band, nrow, ncol)
    bands : tuple
        a tuple of raster bands
    crs : int, string, rasterio.crs.CRS object
        either a epsg code, a proj4 string, or a CRS object
    transform : affine.Affine object
        affine object, which is used to define geometry
    nodataval : float
        raster no data value
    rio_ds : DatasetReader object
        rasterIO dataset Reader object

    Notes
    -----


    Examples
    --------
    >>> from gsflow.utils import Raster
    >>>
    >>> rio = Raster.load("myraster.tif")

    """

    @staticmethod
    def load(raster):
        """
        Static method to load a raster file
        into the raster object

        Parameters
        ----------
        raster : str

        Returns
        -------
            Raster object

        """
        msg = "pyGSFLOW Rasters has been deprecated, calling flopy.utils.Raster()"
        warnings.warn(msg, DeprecationWarning)
        return flopy.utils.Raster.load(raster)
