import os
import numpy as np
from gsflow.utils import Raster
from gsflow.output import PrmsDiscretization


ws = os.path.abspath(os.path.dirname(__file__))


def test_raster_sampling_methods():
    import gsflow
    from gsflow.utils import Raster

    rws = os.path.join(ws, "..", "examples", "data", "geospatial")
    iws = os.path.join(ws, "..", "examples", "data", "sagehen", "modflow")
    raster_name = "dem.img"

    try:
        rio = Raster.load(os.path.join(rws, raster_name))
    except:
        return

    ml = gsflow.modflow.Modflow.load(
        "saghen.nam", version="mfnwt",
        model_ws=iws
    )
    xoff = 214110
    yoff = 4366620
    ml.modelgrid.set_coord_info(xoff, yoff)

    x0, x1, y0, y1 = rio.bounds

    x0 += 3000
    y0 += 3000
    x1 -= 3000
    y1 -= 3000
    shape = np.array([(x0, y0), (x0, y1), (x1, y1), (x1, y0), (x0, y0)])

    rio.crop(shape)

    methods = {
        "min": 2088.52343,
        "max": 2103.54882,
        "mean": 2097.05035,
        "median": 2097.36254,
        "mode": 2088.52343,
        "nearest": 2097.81079,
        "linear": 2097.81079,
        "cubic": 2097.81079,
    }

    for method, value in methods.items():
        data = rio.resample_to_grid(
            ml.modelgrid, band=rio.bands[0], method=method, no_numba=True
        )

        print(data[34, 37])
        if np.abs(data[34, 37] - value) > 1e-05:
            raise AssertionError(
                f"{method} resampling returning incorrect values"
            )


def test_raster_sampling_methods_numba():
    try:
        from numba import jit
    except ImportError:
        return

    import gsflow
    from gsflow.utils import Raster

    rws = os.path.join(ws, "..", "examples", "data", "geospatial")
    iws = os.path.join(ws, "..", "examples", "data", "sagehen", "modflow")
    raster_name = "dem.img"

    try:
        rio = Raster.load(os.path.join(rws, raster_name))
    except:
        return

    ml = gsflow.modflow.Modflow.load(
        "saghen.nam", version="mfnwt",
        model_ws=iws
    )
    xoff = 214110
    yoff = 4366620
    ml.modelgrid.set_coord_info(xoff, yoff)

    x0, x1, y0, y1 = rio.bounds

    x0 += 3000
    y0 += 3000
    x1 -= 3000
    y1 -= 3000
    shape = np.array([(x0, y0), (x0, y1), (x1, y1), (x1, y0), (x0, y0)])

    rio.crop(shape)

    methods = {
        "min": 2088.52343,
        "max": 2103.54882,
        "mean": 2097.05035,
        "median": 2097.36254,
        "mode": 2088.0,  # note some precision is lost in the C mode routine
        "nearest": 2097.81079,
        "linear": 2097.81079,
        "cubic": 2097.81079,
    }

    for method, value in methods.items():
        data = rio.resample_to_grid(
            ml.modelgrid, band=rio.bands[0], method=method
        )

        print(data[34, 37])
        if np.abs(data[34, 37] - value) > 1e-05:
            raise AssertionError(
                f"{method} resampling returning incorrect values"
            )


def test_raster_warp_resampling():
    import gsflow
    from gsflow.builder import GenerateFishnet
    from gsflow.utils import Raster

    rws = os.path.join(ws, "..", "examples", "data", "geospatial")
    raster_name = "dem.img"

    try:
        rio = Raster.load(os.path.join(rws, raster_name))
    except:
        return

    modelgrid = GenerateFishnet(
        os.path.join(rws, raster_name),
        xcellsize=90,
        ycellsize=90,
        buffer=-1
    )

    methods = {
        "min": 2044.927734,
        "max": 2057.975585,
        "mean": 2051.2239,
        "median": 2051.086669,
        "mode": 2057.9756,
        "nearest": 2050.995117,
        "linear": 2050.995117,
        "cubic": 2050.995117,
    }

    for method, value in methods.items():
        # rio = Raster.load(os.path.join(rws, raster_name))
        data = rio.resample_to_grid(
            modelgrid, band=rio.bands[0], method=method
        )

        print(data[34, 37], value)
        if np.abs(data[34, 37] - value) > 1e-04:
            raise AssertionError(
                f"{method} resampling returning incorrect values"
            )


if __name__ == "__main__":
    # test_raster_sampling_methods()
    # test_raster_sampling_methods_numba()
    test_raster_warp_resampling()
