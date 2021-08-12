from flopy.discretization import StructuredGrid
import numpy as np
from flopy.utils import Raster


class GenerateFishnet(StructuredGrid):
    """
    Class to build a fishnet model grid

    The GenerateFishnet class creates a flopy.discretization.StructuredGrid
    object from basic geospatial information

    Parameters:
    ----------
    bbox : shapefile, raster, [xmin, xmax, ymin, ymax]
        bounding box for modelgrid. The bounding box can be a the extent read
        from a shapfile, the extent read in from a raster, or a list of
        [xmin, xmax, ymin, ymax]
    xcellsize : float
        cell size in the x direction for fishnet
    ycellsize : float
        cell size in the y direction for fishnet
    buffer : int
        number of cells to buffer around the input geometry

    """

    def __init__(self, bbox, xcellsize, ycellsize, buffer=None):
        """
        bbox can be
        str - dem path, shapefile
        or a bounding box
        """

        # set cellsize
        self._bbox = bbox
        self.xcs = xcellsize
        self.ycs = ycellsize
        # set buffer
        self._buffer = buffer
        # get the geometry extent
        self._extent = self._get_extent()
        self._xmin, self._ymin, self._xmax, self._ymax = self._extent
        # calculate x and y lengths
        self._lx = self._xmax - self._xmin
        self._ly = self._ymax - self._ymin
        # get the number of rows and columns
        self._nrows, self._ncols = self._get_nrows_ncols()
        # build delr, delc
        delr = np.array([self.xcs] * self._ncols)
        delc = np.array([self.ycs] * self._nrows)
        # generate flopy grid
        super(GenerateFishnet, self).__init__(
            xoff=self._xmin, yoff=self._ymin, delr=delr, delc=delc
        )

    def _get_extent(self):
        """
        gets extent of input geometry

        Returns
        -------
        list : [xmin, xmax, ymin, ymax]
        """

        import shapefile

        if isinstance(self._bbox, str):
            if self._bbox.endswith(".shp"):
                # is shapefile
                shp = shapefile.Reader(self._bbox)
                extent = shp.bbox
            else:
                # assumes the path is for a raster
                raster = Raster.load(self._bbox)
                extent = raster.bounds
                extent = [extent[0], extent[2], extent[1], extent[3]]

        elif isinstance(self._bbox, (list, tuple, np.ndarray)):
            if len(self._bbox) != 4:
                raise AssertionError("Length of bbox must be 4")
            extent = self._bbox

        else:
            raise TypeError(
                "FishnetGenerator Error: Unrecognized bbox type"
                "can be shapefile path, raster path, or "
                "[xmin, xmax, ymin, ymax]"
            )

        if self._buffer is not None:
            xbuffer_len = self.xcs * self._buffer
            ybuffer_len = self.ycs * self._buffer
            xmin, ymin, xmax, ymax = extent
            xmin -= xbuffer_len
            ymin -= ybuffer_len
            xmax += xbuffer_len
            ymax += ybuffer_len
            extent = [xmin, ymin, xmax, ymax]

        return extent

    def _get_nrows_ncols(self):
        """
        determines the number of rows and columns for fishnet

        Returns
        -------
        tuple: (nrows, ncols)
        """

        ncols = int(np.ceil(abs(self._lx) / self.xcs))
        nrows = int(np.ceil(abs(self._ly) / self.ycs))
        return nrows, ncols
