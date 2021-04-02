import os
import numpy as np
from gsflow.utils import Rasters
from gsflow.output import PrmsDiscretization


ws = os.path.abspath(os.path.dirname(__file__))

def test_raster_functions():

    # load the PrmsDiscrtization object
    local_ws = os.path.join(ws, "..", "examples", "data", "sagehen", "shapefiles")
    shp = "hru_params.shp"
    dis = PrmsDiscretization.load_from_shapefile(os.path.join(local_ws, shp))

    # load a DEM raster of the area
    path = os.path.join(ws, '..', 'examples', 'data', 'geospatial')
    name = "dem.img"
    try:
        x = Rasters.load(os.path.join(path, name))
    except ImportError:
        # trap for travis issues with GDAL
        return

    x.set_raster_band(1)

    # sample the raster using the discretization object
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