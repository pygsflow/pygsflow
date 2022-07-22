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


class Raster(flopy.utils.Raster):
    """
    Overloaded flopy.utils.Raster class. This class is used to resample
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

    def __init__(
        self,
        array,
        bands,
        crs,
        transform,
        nodataval,
        driver="GTiff",
        rio_ds=None,
    ):
        super(Raster, self).__init__(
            array, bands, crs, transform, nodataval, driver, rio_ds
        )

    @property
    def bounds(self):
        """
        Returns a tuple of xmin, xmax, ymin, ymax boundaries
        """
        return super(Raster, self).bounds

    @property
    def bands(self):
        """
        Returns a tuple of raster bands
        """
        return super(Raster, self).bands

    @property
    def nodatavals(self):
        """
        Returns a Tuple of values used to define no data
        """
        return super(Raster, self).nodatavals

    @property
    def xcenters(self):
        """
        Returns a np.ndarray of raster x cell centers
        """
        return super(Raster, self).xcenters

    @property
    def ycenters(self):
        """
        Returns a np.ndarray of raster y cell centers
        """
        return super(Raster, self).ycenters

    def sample_point(self, *point, band=1):
        """
        Method to get nearest raster value at a user provided
        point

        Parameters
        ----------
        *point : point geometry representation
            accepted data types:
            x, y values : ex. sample_point(1, 3, band=1)
            tuple of x, y: ex sample_point((1, 3), band=1)
            shapely.geometry.Point
            geojson.Point
            flopy.geometry.Point

        band : int
            raster band to re-sample

        Returns
        -------
            value : float
        """
        return super(Raster, self).sample_point(*point, band)

    def sample_polygon(self, polygon, band, invert=False, **kwargs):
        """
        Method to get an unordered list of raster values that are located
        within a arbitrary polygon

        Parameters
        ----------
        polygon : list, geojson, shapely.geometry, shapefile.Shape
            sample_polygon method accepts any of these geometries:

            a list of (x, y) points, ex. [(x1, y1), ...]
            geojson Polygon object
            shapely Polygon object
            shapefile Polygon shape
            flopy Polygon shape

        band : int
            raster band to re-sample

        invert : bool
            Default value is False. If invert is True then the
            area inside the shapes will be masked out

        Returns
        -------
            np.ndarray of unordered raster values

        """
        return super(Raster, self).sample_polygon(
            polygon, band, invert, **kwargs
        )

    def resample_to_grid(
        self,
        modelgrid,
        band,
        method="nearest",
        multithread=False,
        thread_pool=2,
        extrapolate_edges=False,
        no_numba=False,
        use_oldstyle=False
    ):
        """
        Method to resample the raster data to a
        user supplied grid of x, y coordinates.

        x, y coordinate arrays should correspond
        to grid vertices

        Parameters
        ----------
        modelgrid : flopy.Grid object
            model grid to sample data from
        band : int
            raster band to re-sample
        method : str
            scipy interpolation methods

            ``linear`` for bi-linear interpolation

            ``nearest`` for nearest neighbor

            ``cubic`` for bi-cubic interpolation

            ``mean`` for mean sampling

            ``median`` for median sampling

            ``min`` for minimum sampling

            ``max`` for maximum sampling

        multithread : bool
            boolean flag indicating if multithreading should be used with
            the ``mean`` and ``median`` sampling methods
        thread_pool : int
            number of threads to use for mean and median sampling
        extrapolate_edges : bool
            boolean flag indicating if areas without data should be filled
            using the ``nearest`` interpolation method. This option
            has no effect when using the ``nearest`` interpolation method.
        no_numba : bool
            method to turn off numba based resampling, default is False
        use_oldstyle : bool
            method to force point in polygon intersection routine, default
            is fast resampling based on raster warp sampling methods

        Returns
        -------
            np.array
        """
        from scipy.interpolate import griddata
        from scipy.stats import mode

        method = method.lower()
        if method in ("linear", "nearest", "cubic"):
            xc = modelgrid.xcellcenters
            yc = modelgrid.ycellcenters

            data_shape = xc.shape
            xc = xc.flatten()
            yc = yc.flatten()
            # step 1: create grid from raster bounds
            rxc = self.xcenters
            ryc = self.ycenters

            # step 2: flatten grid
            rxc = rxc.flatten()
            ryc = ryc.flatten()

            # step 3: get array
            if method == "cubic":
                arr = self.get_array(band, masked=False)
            else:
                arr = self.get_array(band, masked=True)
            arr = arr.flatten()

            # step 3: use griddata interpolation to snap to grid
            data = griddata(
                (rxc, ryc),
                arr,
                (xc, yc),
                method=method,
            )

        elif method in ("median", "mean", "min", "max", "mode"):
            # these methods are slow and could use a speed up
            ncpl = modelgrid.ncpl
            data_shape = modelgrid.xcellcenters.shape
            if isinstance(ncpl, (list, np.ndarray)):
                ncpl = ncpl[0]

            data = np.zeros((ncpl,), dtype=float)
            success = False
            if modelgrid.grid_type == "structured":
                # 1st try rasterio resampling methods:
                if modelgrid.angrot != 0:
                    print(
                        "modelgrid rotation not equal to 0, cannot use "
                        "rapid resampling methods. Trying other methods"
                    )
                elif np.sum(modelgrid.delr - modelgrid.delr[0]) > 0:
                    print(
                        "modelgrid delr not constant, cannot use "
                        "rapid resampling methods. Trying other methods"
                    )
                elif np.sum(modelgrid.delc - modelgrid.delc[0]) > 0:
                    print(
                        "modelgrid delc not constant, cannot use "
                        "rapid resampling methods. Trying other methods"
                    )
                elif use_oldstyle:
                    pass
                else:
                    import rasterio
                    from rasterio.enums import Resampling

                    xmin, xmax, ymin, ymax = modelgrid.extent
                    rxmin, rxmax, rymin, rymax = self.bounds
                    x0off, x1off, y0off, y1off = 0, 0, 0, 0
                    if rxmin > xmin or rxmax < xmax or rymin > ymin or rymax < ymax:
                        print(
                            "modelgrid outside bounds of raster, "
                            "offsetting indicies"
                        )
                        xedges, yedges = modelgrid.xyedges
                        xedges += modelgrid.xoffset
                        yedges += modelgrid.yoffset
                        if xmax > rxmax:
                            for ix, val in enumerate(xedges[::-1]):
                                if rxmax > val:
                                    x1off = -1 * ix
                                    xmax = val
                                    break

                        if xmin < rxmin:
                            for ix, val in enumerate(xedges):
                                if rxmin < val:
                                    x0off = ix
                                    xmin = val
                                    break

                        if ymax > rymax:
                            for ix, val in enumerate(yedges):
                                if rymax > val:
                                    y0off = ix
                                    ymax = val
                                    break

                        if ymin < rymin:
                            for ix, val in enumerate(yedges[::-1]):
                                if rymin < val:
                                    y1off = -1 * ix
                                    ymin = val
                                    break

                    poly = [
                        [xmin, ymin], [xmin, ymax], [xmax, ymax], [xmax, ymin]
                    ]

                    self.crop(poly)

                    if method == "mean":
                        resampler = Resampling.average
                    elif method == "median":
                        resampler = Resampling.med
                    elif method == "min":
                        resampler = Resampling.min
                    elif method == "max":
                        resampler = Resampling.max
                    else:
                        resampler = Resampling.mode

                    ncol = modelgrid.ncol - x0off + x1off
                    nrow = modelgrid.nrow - y0off + y1off
                    tmp = self.structured_downscale(
                        band,
                        modelgrid.delr[0],
                        modelgrid.delc[0],
                        nrow,
                        ncol,
                        resampler
                    )
                    data = np.ones((modelgrid.nrow, modelgrid.ncol)) * np.nan
                    if x1off == 0:
                        x1off = None
                    if y1off == 0:
                        y1off = None
                    data[y0off:y1off, x0off:x1off] = tmp
                    success = True

            if ENABLE_JIT and not multithread and not no_numba and not success:
                arr = self.__arr_dict[band].astype(float)
                xcenters = self.xcenters
                ycenters = self.ycenters

                print("Using parallel numba for resampling")
                verts = modelgrid.verts
                iverts = modelgrid.iverts
                ray_count = np.zeros(xcenters.shape, dtype=int)
                mask = np.ones(xcenters.shape, dtype=bool)
                arr = np.ravel(arr)

                polygons = np.zeros((ncpl, 5, 2))
                tmp = verts[np.array(iverts).ravel()]
                tmp.shape = (ncpl, 4, 2)
                polygons[:, 0:4, :] = tmp
                polygons[:, -1, :] = polygons[:, 0, :]

                data = compiled_resampling(
                    polygons,
                    method,
                    arr,
                    xcenters,
                    ycenters,
                    self.nodatavals,
                    ray_count,
                    mask,
                    data,
                    ncpl,
                )

            elif multithread and not success:
                q = queue.Queue()
                container = threading.BoundedSemaphore(thread_pool)

                # determine the number of thread pairs required to
                # fill the grid
                nthreadpairs = int(ncpl / thread_pool)
                if ncpl % thread_pool != 0:
                    nthreadpairs += 1

                # iterate over the tread pairs
                for idx in range(nthreadpairs):
                    i0 = idx * thread_pool
                    nthreads = thread_pool
                    if i0 + thread_pool > ncpl:
                        nthreads = ncpl - i0
                    i1 = i0 + nthreads
                    threads = []
                    for node in range(i0, i1):
                        t = threading.Thread(
                            target=self.__threaded_resampling,
                            args=(modelgrid, node, band, method, container, q),
                        )
                        threads.append(t)

                    # start the threads
                    for thread in threads:
                        thread.daemon = True
                        thread.start()

                    # wait until all threads are terminated
                    for thread in threads:
                        thread.join()

                    for idx in range(nthreads):
                        node, val = q.get()
                        data[node] = val

            elif not success:
                for node in range(ncpl):
                    verts = modelgrid.get_cell_vertices(node)
                    try:
                        rstr_data = self.sample_polygon(
                            verts, band, convert=False
                        ).astype(float)
                    except TypeError:
                        rstr_data = self.sample_polygon(
                            verts, band
                        ).astype(float)
                    msk = np.in1d(rstr_data, self.nodatavals)
                    rstr_data[msk] = np.nan

                    if rstr_data.size == 0:
                        val = self.nodatavals[0]
                    else:
                        if method == "median":
                            val = np.nanmedian(rstr_data)
                        elif method == "mean":
                            val = np.nanmean(rstr_data)
                        elif method == "max":
                            val = np.nanmax(rstr_data)
                        elif method == "mode":
                            val = mode(
                                rstr_data, axis=None, nan_policy="omit"
                            ).mode
                            if len(val) == 0:
                                val = np.nan
                            else:
                                val = val[0]
                        else:
                            val = np.nanmin(rstr_data)

                    data[node] = val
        else:
            raise TypeError(f"{method} method not supported")

        if extrapolate_edges and method != "nearest":
            xc = modelgrid.xcellcenters
            yc = modelgrid.ycellcenters

            xc = xc.flatten()
            yc = yc.flatten()

            # step 1: create grid from raster bounds
            rxc = self.xcenters
            ryc = self.ycenters

            # step 2: flatten grid
            rxc = rxc.flatten()
            ryc = ryc.flatten()

            arr = self.get_array(band, masked=True).flatten()

            # filter out nan values from the original dataset
            if np.isnan(np.sum(arr)):
                idx = np.isfinite(arr)
                rxc = rxc[idx]
                ryc = ryc[idx]
                arr = arr[idx]

            extrapolate = griddata(
                (rxc, ryc),
                arr,
                (xc, yc),
                method="nearest",
            )
            data = np.where(np.isnan(data), extrapolate, data)

        # step 4: return grid to user in shape provided
        data.shape = data_shape

        # step 5: re-apply nodata values
        data[np.isnan(data)] = self.nodatavals[0]

        return data

    def structured_downscale(self, band, delr, delc, nrow, ncol, resampler):
        """
        Fast method for structured grid downscaling. Requires 0 rotation and
        that delc and delr are constant for accurate downscaling.

        Parameters
        ----------
        band : int
            raster band
        delr : int
            x resolution
        delc : int
            y resolution
        resampling : rasterio Resampling object

        Returns
        -------
            np.ndarray
        """
        import rasterio.warp as riow

        arr = self.get_array(band)
        nodata = self.nodatavals
        if nodata is None:
            nodata = -1e+30
        bounds = self.bounds
        bbox = (bounds[0], bounds[2], bounds[1], bounds[3])
        import affine
        newaff0, width, height = riow.calculate_default_transform(
             self._meta["crs"],
             self._meta["crs"],
             self._meta["width"],
             self._meta["height"],
             *bbox,
             resolution=(delr, delc)
        )

        newaff = affine.Affine(
            delr,
            self._meta["transform"].b,
            self._meta["transform"].c,
            self._meta["transform"].d,
            -1 * delc,
            self._meta["transform"].f
        )

        newarr = np.ones(shape=(nrow, ncol), dtype=float) * np.nan
        newarr, newaff = riow.reproject(
            arr,
            newarr,
            src_transform=self._meta["transform"],
            dst_transform=newaff,
            width=ncol,
            height=nrow,
            src_nodata=nodata[0],
            dst_nodata=nodata[0],
            src_crs=self._meta["crs"],
            dst_crs=self._meta["crs"],
            resampling=resampler
        )

        return newarr

    def crop(self, polygon, invert=False):
        """
        Method to crop a new raster object
        from the current raster object

        Parameters
        ----------
        polygon : list, geojson, shapely.geometry, shapefile.Shape
            crop method accepts any of these geometries:

            a list of (x, y) points, ex. [(x1, y1), ...]
            geojson Polygon object
            shapely Polygon object
            shapefile Polygon shape
            flopy Polygon shape

        invert : bool
            Default value is False. If invert is True then the
            area inside the shapes will be masked out

        """
        return super(Raster, self).crop(polygon, invert)

    def get_array(self, band, masked=True):
        """
        Method to get a numpy array corresponding to the
        provided raster band. Nodata vals are set to
        np.NaN

        Parameters
        ----------
        band : int
            band number from the raster
        masked : bool
            determines if nodatavals will be returned as np.nan to
            the user

        Returns
        -------
            np.ndarray

        """
        return super(Raster, self).get_array(band, masked)

    def write(self, name):
        """
        Method to write raster data to a .tif
        file

        Parameters
        ----------
        name : str
            output raster .tif file name

        """
        super(Raster, self).write(name)

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
        import rasterio

        dataset = rasterio.open(raster)
        array = dataset.read()
        bands = dataset.indexes
        meta = dataset.meta

        return Raster(
            array,
            bands,
            meta["crs"],
            meta["transform"],
            meta["nodata"],
            meta["driver"],
        )

    def plot(self, ax=None, contour=False, **kwargs):
        """
        Method to plot raster layers or contours.

        Parameters
        ----------
        ax : matplotlib.pyplot.axes
            optional matplotlib axes for plotting
        contour : bool
            flag to indicate creation of contour plot

        **kwargs :
            matplotlib keyword arguments
            see matplotlib documentation for valid
            arguments for plot and contour.

        Returns
        -------
            ax : matplotlib.pyplot.axes

        """
        return super(Raster, self).plot(ax, contour, **kwargs)

    def histogram(self, ax=None, **kwargs):
        """
        Method to plot a histogram of digital numbers

        Parameters
        ----------
        ax : matplotlib.pyplot.axes
            optional matplotlib axes for plotting

        **kwargs :
            matplotlib keyword arguments
            see matplotlib documentation for valid
            arguments for histogram

        Returns
        -------
            ax : matplotlib.pyplot.axes

        """
        return super(Raster, self).histogram(ax, **kwargs)


@jit(nopython=True, parallel=True)
def compiled_resampling(
    polygons, method, arr, xc, yc, nodatavals, ray_count, mask, data, ncpl
):
    """
    Multiprocessing method to resample to grid using mean, min, max, mode,
    or median resampling. Method uses numba compiled python to take advantage
    of parallel compiled methods for performing point in polygon calculations.
    Code is ~ 3 - 4 times faster than python multithreading and
    multiprocessing.

    Parameters
    ----------
    polygons : np.ndarray
        3 dimensional array of polygon vertices for the model
    method : str
        sampling methods

            ``mean`` for mean sampling

            ``median`` for median sampling

            ``min`` for minimum sampling

            ``max`` for maximum sampling

            ``mode`` for most frequent sampling
    arr : np.ndarray
        array of raster data
    xc : np.ndarray
        raster cell xcenter points
    yc : np.ndarray
        raster cell ycenter points
    nodatavals : np.ndarray
        1d array of no data values for raster
    ray_count : np.ndarray
        2d zero array of shape nrow, ncol
    mask : np.ndarray
        2d boolean array (True) of shape nrow, ncol
    data : np.ndarray
        2d zero array of shape nrow, ncol to set resampled data to
    ncpl : int
        number of cells per layer

    Returns
    -------
        data : np.ndarray of resampled raster data
    """
    # perform ray casting algorithm for point in polygon calculation

    for node in range(ncpl):
        polygon = polygons[node]
        ray_count[:, :] = 0
        num = len(polygon)
        j = num - 1
        for i in range(num):
            tmp = polygon[i][0] + (polygon[j][0] - polygon[i][0]) * (
                yc - polygon[i][1]
            ) / (polygon[j][1] - polygon[i][1])

            comp = np.where(
                ((polygon[i][1] > yc) ^ (polygon[j][1] > yc)) & (xc < tmp)
            )

            j = i
            if len(comp[0]) > 0:
                for ix, ii in enumerate(comp[0]):
                    ray_count[ii, comp[1][ix]] += 1

        mask[:, :] = True
        mask = np.where(ray_count % 2 == 0, False, True)

        fmask = np.ravel(mask)
        rstr_data = arr[fmask]

        for nval in nodatavals:
            rstr_data[rstr_data == nval] = np.nan

        if rstr_data.size == 0:
            val = nodatavals[0]

        else:
            if method == "median":
                val = np.nanmedian(rstr_data)
            elif method == "mean":
                val = np.nanmean(rstr_data)
            elif method == "max":
                val = np.nanmax(rstr_data)
            elif method == "mode":
                mval = numba_mode(rstr_data)[0]
                if len(mval) == 0:
                    val = np.nan
                else:
                    val = mval[0]
            else:
                val = np.nanmin(rstr_data)

        data[node] = val

    return data


@jit(nopython=True)
def numba_mode(a):
    """
    Numba method to calculate most frequent value in an array

    Parameters
    ----------
    a : np.ndarray
        1d array of values

    Returns
    -------
        mode, count
    """
    a = np.ravel(a)
    a = a[~np.isnan(a)]
    if a.size == 0:
        return np.array([np.nan]), np.array([0])

    scores = np.unique(np.ravel(a))
    oldmostfreq = np.array([0])
    oldcounts = np.array([0])

    for score in scores:
        template = a == score
        counts = np.sum(template)
        mostfrequent = np.where(counts > oldcounts, score, oldmostfreq)
        oldcounts = np.maximum(counts, oldcounts)
        if len(mostfrequent) > 0:
            oldmostfreq[0] = mostfrequent[0]

    return mostfrequent, oldcounts
