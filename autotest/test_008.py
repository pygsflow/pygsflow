import os
import numpy as np
from gsflow.output import PrmsDiscretization


ws = os.path.abspath(os.path.dirname(__file__))


def test_empty_dis_object():
    dis = PrmsDiscretization([])
    assert isinstance(dis, PrmsDiscretization)
    assert dis.nhru == 0


def test_dis_object_from_shp():
    local_ws = os.path.join(ws, "..", "examples", "data", "sagehen", "shapefiles")
    shp = "hru_params.shp"

    dis = PrmsDiscretization.load_from_shapefile(os.path.join(local_ws, shp))
    assert isinstance(dis, PrmsDiscretization)
    assert dis.nhru == 6468
    assert dis.extent is not None


def test_dis_object_from_flopy():
    from gsflow.modflow import Modflow

    local_ws = os.path.join(ws, "..", "examples", "data", "sagehen", "modflow")
    nam = "saghen.nam"

    ml = Modflow.load(nam, model_ws=local_ws)
    dis = PrmsDiscretization.load_from_flopy(ml)
    assert isinstance(dis, PrmsDiscretization)
    assert dis.nhru == 6468
    assert dis.extent is not None


if __name__ == "__main__":
    test_empty_dis_object()
    test_dis_object_from_shp()
    test_dis_object_from_flopy()