import os
import shutil
from pathlib import Path

import flopy
import gsflow
from flopy.utils.flopy_io import which

from utils import oct_tree_grid
from utils.mfusg_builder import (
    build_mfusg, add_uzf, add_rch, add_evt
)
from utils.prms_builder import build_prms

# ================================
# Configuration
# ================================


class Model_info:
    """Container for model configuration and paths."""
    pass


mi = Model_info()

# Directory paths
mi.ws = os.path.abspath(os.path.dirname(__file__))
mi.iws = Path(mi.ws).parent / 'data' / 'geospatial'
mi.ows = os.path.join(mi.ws, "temp")
if not os.path.exists(mi.ows):
    os.mkdir(mi.ows)

# Input files
mi.dem_file = os.path.join(mi.iws, 'dem.img')
mi.pour_point_file = os.path.join(
    mi.iws, "model_points.shp"
)
mi.resampled_dem = os.path.join(
    mi.ows, "sagehen_50m_med.txt"
)

# Unstructured grid (USG) model paths
mi.usg_model_ws = Path(mi.ws) / 'ugmodel'
mi.usg_base_name = "usg_sagehen"
mi.usg_model_fn = os.path.join(
    mi.usg_model_ws, f"{mi.usg_base_name}.nam"
)
mi.usg_control_file = os.path.join(
    mi.usg_model_ws, f"{mi.usg_base_name}.control"
)

# Gridgen executable (must be on PATH)
mi.gridgen_exe = which("gridgen")

# Grid generation parameters
mi.stream_threshold = 810000  # drainage area (m^2)
mi.cellsize = 50

# Fine (structured) model paths
mi.fine_model_ws = str(
    Path(mi.ws).parent / 'frontiers' / 'temp'
)
mi.fine_model_fn = os.path.join(
    mi.fine_model_ws, "sagehen_50m.nam"
)
mi.fine_control_file = os.path.join(
    mi.fine_model_ws, "sagehen_50m_cont.control"
)

# Create a clean USG output directory
if os.path.exists(mi.usg_model_ws):
    shutil.rmtree(mi.usg_model_ws)
os.mkdir(mi.usg_model_ws)

# ================================
# Load fine (structured) model
# ================================
mi.fine_gsf = gsflow.GsflowModel.load_from_file(
    mi.fine_control_file
)
mi.nlay = mi.fine_gsf.mf.dis.nlay
mi.botm = mi.fine_gsf.mf.dis.botm

# Set coordinate reference system from DEM
mi.crs = flopy.utils.Raster.load(mi.dem_file).crs
mi.fine_gsf.mf.modelgrid.set_coord_info(crs=mi.crs)

# ================================
# Generate unstructured grid
# ================================
oct_tree_grid.create_oct_tree_grid(
    mi, stream_buffer_distance=100.0
)

# ================================
# Build MODFLOW-USG and PRMS models
# ================================
build_mfusg(mi)
build_prms(mi)

# Add stress packages
add_uzf(mi)
add_rch(mi)
add_evt(mi)
