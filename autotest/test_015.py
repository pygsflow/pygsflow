import os
import numpy as np
from gsflow.utils import Raster
from gsflow.output import PrmsDiscretization


ws = os.path.abspath(os.path.dirname(__file__))


def test_raster():
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
        "mode": 2044.92773,
        "nearest": 2050.995117,
        "linear": 2050.995117,
        "cubic": 2050.995117,
    }

    for method, value in methods.items():
        data = rio.resample_to_grid(
            modelgrid, band=rio.bands[0], method=method
        )

        print(data[34, 37], value)
        if np.abs(data[34, 37] - value) > 1e-04:
            raise AssertionError(
                f"{method} resampling returning incorrect values"
            )


if __name__ == "__main__":
    test_raster()
