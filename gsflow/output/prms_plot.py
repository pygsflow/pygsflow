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

    Examples
    --------

    create a plot object and plot an array

    >>> dis = gsflow.output.PrmsDiscretization.load_from_shapefile("myshp.shp")
    >>> plot = gsflow.output.PrmsPlot(prms_dis=dis)
    >>> array = np.random((dis.nhru,))
    >>> ax = plot.plot_array(array)
    >>> plt.show()

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

        masked_values : list, optional
            list of values to mask from plotting

        kwargs : matplotlib keyword arguments

        Returns
        -------
            matplotlib.collections.PatchCollection

        """
        if not self.__maps:
            raise AssertionError(
                "PrmsPlot must be given a PrmsDiscretization"
                " object to plot parameter maps"
            )

        if len(array) != self.prms_dis.nhru:
            raise AssertionError("array size does not match nhru")

        if ax is None:
            ax = plt.gca()

        if not isinstance(array, np.ndarray):
            array = np.array(array)

        if masked_values is not None:
            for mval in masked_values:
                array = np.ma.masked_equal(array, mval)

        return self.__plot_patches(array, ax=ax, **kwargs)

    def contour_array(self, array, ax=None, masked_values=None, **kwargs):
        """
        Contour an array.

        Parameters
        ----------
        a : numpy.ndarray
            Array to plot.

        ax : matplotlib.pyplot.axes
            optional axes object

        masked_values : list, optional
            list of values to mask.

        **kwargs : dictionary
            keyword arguments passed to matplotlib.pyplot.pcolormesh

        Returns
        -------
        contour_set : matplotlib.tri.tricontour object

        """
        import matplotlib.tri as tri

        if not self.__maps:
            raise AssertionError(
                "PrmsPlot must be given a PrmsDiscretization"
                " object to plot parameter maps"
            )

        if len(array) != self.prms_dis.nhru:
            raise AssertionError("array size does not match nhru")

        if ax is None:
            ax = plt.gca()

        xcentergrid = np.array(self.prms_dis.x_hru_centers)
        ycentergrid = np.array(self.prms_dis.y_hru_centers)

        if not isinstance(array, np.ndarray):
            array = np.array(array)

        if masked_values is not None:
            for mval in masked_values:
                array = np.ma.masked_equal(array, mval)

        if "colors" in kwargs.keys():
            if "cmap" in kwargs.keys():
                kwargs.pop("cmap")

        triang = tri.Triangulation(xcentergrid, ycentergrid)

        mask = None
        try:
            amask = array.mask
            mask = [False for i in range(triang.triangles.shape[0])]
            for ipos, (n0, n1, n2) in enumerate(triang.triangles):
                if amask[n0] or amask[n1] or amask[n2]:
                    mask[ipos] = True
            triang.set_mask(mask)
        except:
            pass

        contour_set = ax.tricontour(triang, array, **kwargs)

        ax.set_xlim(self.extent[0], self.extent[1])
        ax.set_ylim(self.extent[2], self.extent[3])

        return contour_set

    def plot_parameter(self, parameter, ax=None, masked_values=None, **kwargs):
        """
        Method to plot a parameter from a ParameterRecord

        Parameters
        ----------
        parameter : PrmsParameter object

        ax : matplotlib.pyplot.axes

        masked_values : list, optional
            list of values to mask from plotting

        kwargs : matplotlib keyword arguments

        Returns
        -------
            matplotlib.collections.PatchCollection or
            matplotlib.figure.Figure

        """
        if not self.__maps:
            raise AssertionError(
                "PrmsPlot must be given a PrmsDiscretization"
                " object to plot parameter maps"
            )

        dims = parameter.dims
        if len(dims) == 1:

            if ax is None:
                ax = plt.gca()

            if dims[0] != self.prms_dis.nhru:
                raise AssertionError("Parameter dimensions do not match nhru")

            array = parameter.values
            return self.plot_array(
                array, ax=ax, masked_values=masked_values, **kwargs
            )

        elif len(dims) == 2:
            if dims[0] != self.prms_dis.nhru:
                raise AssertionError("Parameter dimensions do not match nhru")

            hru_dim = dims[0]
            nplot = dims[1]
            values = parameter.values
            values.shape = (nplot, hru_dim)
            kwargs["vmin"] = np.min(values)
            kwargs["vmax"] = np.max(values)

            if nplot == 12:
                fig = plt.figure(figsize=(22, 9))
                for i in range(1, 13):
                    ax = fig.add_subplot(3, 4, i)
                    ax.set_aspect("auto")
                    txt = "Month {}".format(i)
                    ax.set_title(txt)
                    ax.tick_params(axis="both", which="both", labelsize=8)
                    ax.tick_params(axis="y", which="both", labelrotation=0)
                    p = self.plot_array(
                        values[i - 1],
                        ax=ax,
                        masked_values=masked_values,
                        **kwargs
                    )

                fig.subplots_adjust(
                    left=0.10, right=0.85, wspace=0.3, hspace=0.3
                )
                cbar_ax = fig.add_axes([0.9, 0.1, 0.05, 0.7])
                fig.colorbar(p, cax=cbar_ax)
                return fig

        else:
            raise ValueError("Valid number of dimensions is 1 or 2")

    def contour_parameter(
        self, parameter, ax=None, masked_values=None, **kwargs
    ):
        """
        Contour an array.

        Parameters
        ----------
        parameter : PrmsParameter object
        ax : matplotlib.pyplot.axes
            optional axes object
        masked_values : iterable of floats, ints
            Values to mask.
        **kwargs : dictionary
            keyword arguments passed to matplotlib.pyplot.pcolormesh

        Returns
        -------
        contour_set : matplotlib.tri.tricontour object

        """
        if not self.__maps:
            raise AssertionError(
                "PrmsPlot must be given a PrmsDiscretization"
                " object to plot parameter maps"
            )

        dims = parameter.dims
        if len(dims) == 1:

            if ax is None:
                ax = plt.gca()

            if dims[0] != self.prms_dis.nhru:
                raise AssertionError("Parameter dimensions do not match nhru")

            array = parameter.values
            return self.contour_array(
                array, ax=ax, masked_values=masked_values, **kwargs
            )

        elif len(dims) == 2:
            if dims[0] != self.prms_dis.nhru:
                raise AssertionError("Parameter dimensions do not match nhru")

            hru_dim = dims[0]
            nplot = dims[1]
            values = parameter.values
            values.shape = (nplot, hru_dim)
            kwargs["vmin"] = np.min(values)
            kwargs["vmax"] = np.max(values)

            if nplot == 12:
                fig = plt.figure(figsize=(22, 9))
                for i in range(1, 13):
                    ax = fig.add_subplot(3, 4, i)
                    ax.set_aspect("auto")
                    txt = "Month {}".format(i)
                    ax.set_title(txt)
                    ax.tick_params(axis="both", which="both", labelsize=8)
                    ax.tick_params(axis="y", which="both", labelrotation=0)

                    p = self.contour_array(
                        values[i - 1],
                        ax=ax,
                        masked_values=masked_values,
                        **kwargs
                    )

                fig.subplots_adjust(
                    left=0.10, right=0.85, wspace=0.3, hspace=0.3
                )
                cbar_ax = fig.add_axes([0.9, 0.1, 0.05, 0.7])
                fig.colorbar(p, cax=cbar_ax)
                return fig

        else:
            raise ValueError("Valid number of dimensions is 1 or 2")

    def plot_data_timeseries(self, data, names, ax=None, **kwargs):
        """
        Not implemented yet

        Parameters
        ----------
        data
        names
        ax
        kwargs

        Returns
        -------

        """
        raise NotImplementedError(
            "plot_data_timeseries is not yet implemented"
        )

    def plot_model_discretization(self, ax=None, **kwargs):
        """
        Plots the model grid

        Parameters
        ----------
        ax : matplotlib.axes object, optional

        kwargs :
            matplotlib keyword arguments

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
        # if "color" not in kwargs:
        #    kwargs["facecolor"] = "None"

        if "vmin" in kwargs:
            vmin = kwargs.pop("vmin")
        else:
            vmin = np.min(array)

        if "vmax" in kwargs:
            vmax = kwargs.pop("vmax")
        else:
            vmax = np.max(array)

        patches = [Polygon(hru, True) for hru in self.prms_dis.xypts]

        p = PatchCollection(patches)
        p.set_array(array)
        p.set_clim(vmin=vmin, vmax=vmax)
        p.set(**kwargs)

        ax.add_collection(p)

        extent = self.extent
        ax.set_xlim([extent[0], extent[1]])
        ax.set_ylim([extent[2], extent[3]])
        return p
