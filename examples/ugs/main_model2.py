import os
#from tempfile import TemporaryDirectory
import shutil

# import matplotlib.pyplot as plt
# import numpy as np
# import pandas as pd

import gsflow, flopy
from flopy.utils.gridgen import Gridgen
#from gsflow import GsflowModel, PrmsModel, PrmsData

# import utils
# import prms_utils
from utils import oct_tree_grid
from utils.mfusg_builder import build_mfusg, add_uzf, add_rch, add_evt
from utils.prms_builder import build_prms
from pathlib import Path
# ================================
# (1) Global Variables
# ================================
class Model_info:
    pass
sample_grid = True  # make true if first time running
# hold all the model info

mi = Model_info()
mi.ws = os.path.abspath(os.path.dirname(__file__))
mi.iws = Path(mi.ws).parent/'data'/'geospatial'
mi.ows = os.path.join(mi.ws, "temp")
if not os.path.exists(mi.ows):
    os.mkdir(mi.ows)
mi.dem_file = os.path.join(mi.iws, 'dem.img')
mi.pour_point_file = os.path.join(mi.iws, "model_points.shp")
mi.resampled_dem = os.path.join(mi.ows, 
                                 "sagehen_50m_med.txt")

mi.usg_model_ws = Path(mi.ws)/'ugmodel'
mi.usg_base_name = "usg_sagehen"
mi.gridgen_exe = r"C:\Users\sregan\Workspace\bin\gridgen.exe"

mi.stream_threshold = 810000  # m3 of drainage area
mi.cellsize = 50
mi.fine_model_ws = (r"C:\Users\sregan\Workspace\git_repositories\pygsflow\examples\frontiers\temp")
mi.fine_model_fn = os.path.join(mi.fine_model_ws,
                                 "sagehen_50m.nam")
mi.fine_control_file = os.path.join(mi.fine_model_ws, 
                                     "sagehen_50m_cont.control")

mi.usg_model_fn = os.path.join(mi.usg_model_ws, f"{mi.usg_base_name}.nam")
mi.usg_control_file = os.path.join(mi.usg_model_ws, f"{mi.usg_base_name}.control")

if os.path.exists(mi.usg_model_ws):
    shutil.rmtree(mi.usg_model_ws)
os.mkdir(mi.usg_model_ws)

# ================================
# (2) Load fine model
# ================================
mi.fine_gsf = gsflow.GsflowModel.load_from_file(mi.fine_control_file)
mi.nlay = mi.fine_gsf.mf.dis.nlay
mi.botm = mi.fine_gsf.mf.dis.botm
# set the crs
mi.crs = flopy.utils.Raster.load(mi.dem_file).crs
mi.fine_gsf.mf.modelgrid.set_coord_info(crs=mi.crs)



## ================================
# (2) Unstructured Grid Generation
## ================================

oct_tree_grid.create_oct_tree_grid(mi, 
                                    stream_buffer_distance=100.0)

# 

build_mfusg(mi)
build_prms(mi)

add_uzf(mi)
add_rch(mi)
add_evt(mi)


#import gridutil
#gridutil.plot_grid(mi.oct_grid2d, mi.gsflow.prms.parameters.get_values('hru_strmseg_down_id'))

end = 1