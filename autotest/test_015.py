import os
import numpy as np
from gsflow.utils import Rasters
from gsflow.output import PrmsDiscretization


def test_raster_functions():

    # load the PrmsDiscrtization object
    ws = r"../examples/data/sagehen/shapefiles"
    shp = "hru_params.shp"
    dis = PrmsDiscretization.load_from_shapefile(os.path.join(ws, shp))

    # load a DEM raster of the area
    path = r'../examples/data/geospatial'
    name = "dem.img"
    try:
        x = Rasters.load(os.path.join(path, name))
    except ImportError:
        # trap for travis issues with
        return
    x.set_raster_band(1)
    dem = x.sample_discretization(dis)

    if not isinstance(dem, np.ndarray):
        raise AssertionError()

    if len(dem) != len(dis.x_hru_centers):
        raise AssertionError

    xp = x.xpoints
    yp = x.ypoints
    xyp = x.xypoints

    if xp is None:
        raise AssertionError()
    if yp is None:
        raise AssertionError()
    if xyp is None:
        raise AssertionError()


    array = x.band_array
    if not isinstance(array, np.ndarray):
        raise AssertionError()

if __name__ == "__main__":
    test_raster_functions()