import os
import numpy as np

ws = os.path.abspath(os.path.dirname(__file__))
cellsize = 90


def test_flow_accumulation_import():
    from gsflow.builder import FlowAccumulation


def test_flow_directions():
    from gsflow.builder import FlowAccumulation, GenerateFishnet
    try:
        import rasterio
    except ImportError:
        return

    data_ws = os.path.join(ws, "..", "examples", "data", "geospatial")
    resampled_ws = os.path.join(data_ws, "resampled_90m")
    dem_file = os.path.join(data_ws, "dem.img")
    resampled_dem = os.path.join(resampled_ws, "dem_90m_median.txt")

    modelgrid = GenerateFishnet(
        dem_file,
        xcellsize=cellsize,
        ycellsize=cellsize
    )
    dem_data = np.genfromtxt(resampled_dem)

    fa = FlowAccumulation(
        dem_data,
        modelgrid.xcellcenters,
        modelgrid.ycellcenters,
        hru_type=np.ones((modelgrid.nrow, modelgrid.ncol), dtype=int),
        verbose=False
    )
    flow_dir = fa.flow_directions(dijkstra=True)

    undef = (-1, -2)
    # check for undefined cells:
    for i in undef:
        chk = np.where(flow_dir == i)
        if len(chk[0]) > 0:
            raise AssertionError(
                "dijkstra flow direction failed and left undefined cells"
            )

    ddict = {
        1: 1235, 2: 775, 4: 811, 8: 374, 16: 399, 32: 549, 64: 846, 128: 1402
    }
    flow_dir = np.ravel(flow_dir)
    for direction, count in ddict.items():
        t = np.where(flow_dir == direction)
        cnt = len(t[0])
        if cnt != count:
            raise AssertionError(
                "dijkstra flow directions not properly assigned"
            )

    fa = FlowAccumulation(
        dem_data,
        modelgrid.xcellcenters,
        modelgrid.ycellcenters,
        hru_type=np.ones((modelgrid.nrow, modelgrid.ncol), dtype=int),
        verbose=False
    )
    flow_dir = fa.flow_directions(dijkstra=False)

    for i in undef:
        chk = np.where(flow_dir == i)
        if len(chk[0]) > 0:
            raise AssertionError(
                "topological flow direction failed and left undefined cells"
            )


def test_flow_accumulation():
    from gsflow.builder import FlowAccumulation, GenerateFishnet
    try:
        import rasterio
    except ImportError:
        return

    data_ws = os.path.join(ws, "..", "examples", "data", "geospatial")
    dem_file = os.path.join(data_ws, "dem.img")
    resampled_ws = os.path.join(data_ws, "resampled_90m")
    resampled_dem = os.path.join(resampled_ws, "dem_90m_median.txt")
    flow_directions = os.path.join(resampled_ws, "flow_directions_90m.txt")

    modelgrid = GenerateFishnet(
        dem_file,
        xcellsize=cellsize,
        ycellsize=cellsize
    )
    dem_data = np.genfromtxt(resampled_dem)
    flow_directions = np.genfromtxt(flow_directions, dtype=float)

    fa = FlowAccumulation(
        dem_data,
        modelgrid.xcellcenters,
        modelgrid.ycellcenters,
        hru_type=np.ones((modelgrid.nrow, modelgrid.ncol), dtype=int),
        flow_dir_array=flow_directions,
        verbose=False,
    )

    flow_acc = fa.flow_accumulation()
    if np.max(flow_acc) != 3567:
        raise AssertionError("flow accumulation did not accumulate properly")


def test_watershed_delineation():
    from gsflow.builder import FlowAccumulation, GenerateFishnet
    try:
        import rasterio
        import shapefile
    except ImportError:
        return

    data_ws = os.path.join(ws, "..", "examples", "data", "geospatial")
    dem_file = os.path.join(data_ws, "dem.img")
    resampled_ws = os.path.join(data_ws, "resampled_90m")
    resampled_dem = os.path.join(resampled_ws, "dem_90m_median.txt")
    flow_directions = os.path.join(resampled_ws, "flow_directions_90m.txt")
    pour_point_file = os.path.join(data_ws, "model_points.shp")

    modelgrid = GenerateFishnet(
        dem_file,
        xcellsize=cellsize,
        ycellsize=cellsize
    )
    dem_data = np.genfromtxt(resampled_dem)
    flow_directions = np.genfromtxt(flow_directions, dtype=float)

    fa = FlowAccumulation(
        dem_data,
        modelgrid.xcellcenters,
        modelgrid.ycellcenters,
        hru_type=np.ones((modelgrid.nrow, modelgrid.ncol), dtype=int),
        flow_dir_array=flow_directions,
        verbose=False,
    )

    # read in pour point from shapefile and set watershed boundary
    with shapefile.Reader(pour_point_file) as r:
        shape = r.shape(0)
        pour_point = shape.points
        pour_point[0][1] -= 20
        pour_point[0][0] -= 20

    watershed = fa.define_watershed(pour_point, modelgrid, fmt="xy")
    n = np.count_nonzero(watershed)

    if n != 3324:
        raise AssertionError("Watershed dilineation failed")


def test_stream_generation():
    from gsflow.builder import FlowAccumulation, GenerateFishnet
    try:
        import rasterio
        import shapefile
    except ImportError:
        return

    data_ws = os.path.join(ws, "..", "examples", "data", "geospatial")
    dem_file = os.path.join(data_ws, "dem.img")
    resampled_ws = os.path.join(data_ws, "resampled_90m")
    resampled_dem = os.path.join(resampled_ws, "dem_90m_median.txt")
    flow_directions = os.path.join(resampled_ws, "flow_directions_90m.txt")
    flow_accumulation = os.path.join(resampled_ws, "flow_accumulation_90m.txt")
    watershed = os.path.join(resampled_ws, "watershed_90m.txt")

    stream_threshold = 810000  # m3 of drainage area

    modelgrid = GenerateFishnet(
        dem_file,
        xcellsize=cellsize,
        ycellsize=cellsize
    )
    dem_data = np.genfromtxt(resampled_dem)
    flow_directions = np.genfromtxt(flow_directions, dtype=float)
    flow_accumulation = np.genfromtxt(flow_accumulation, dtype=float)
    watershed = np.genfromtxt(watershed, dtype=int)

    fa = FlowAccumulation(
        dem_data,
        modelgrid.xcellcenters,
        modelgrid.ycellcenters,
        hru_type=watershed,
        flow_dir_array=flow_directions,
        verbose=False,
    )

    threshold = stream_threshold / (cellsize ** 2)

    strm_object = fa.make_streams(
        flow_directions, flow_accumulation, threshold,
    )

    if np.max(strm_object.iseg) != 15:
        raise AssertionError("Stream generation failed")

    n = np.count_nonzero(strm_object.iseg)
    if n != 191:
        raise AssertionError("improper number of stream cells")


def test_cascade_generation():
    from gsflow.builder import FlowAccumulation, GenerateFishnet
    try:
        import rasterio
        import shapefile
    except ImportError:
        return

    data_ws = os.path.join(ws, "..", "examples", "data", "geospatial")
    dem_file = os.path.join(data_ws, "dem.img")
    resampled_ws = os.path.join(data_ws, "resampled_90m")
    resampled_dem = os.path.join(resampled_ws, "dem_90m_median.txt")
    flow_directions = os.path.join(resampled_ws, "flow_directions_90m.txt")
    pour_point_file = os.path.join(data_ws, "model_points.shp")
    watershed = os.path.join(resampled_ws, "watershed_90m.txt")
    streams = os.path.join(resampled_ws, "streams_90m.bin")

    modelgrid = GenerateFishnet(
        dem_file,
        xcellsize=cellsize,
        ycellsize=cellsize
    )
    dem_data = np.genfromtxt(resampled_dem)
    flow_directions = np.genfromtxt(flow_directions, dtype=float)
    watershed = np.genfromtxt(watershed, dtype=int)
    streams = FlowAccumulation.load_streams(streams)

    with shapefile.Reader(pour_point_file) as r:
        shape = r.shape(0)
        pour_point = shape.points
        pour_point[0][1] -= 20
        pour_point[0][0] -= 20

    fa = FlowAccumulation(
        dem_data,
        modelgrid.xcellcenters,
        modelgrid.ycellcenters,
        hru_type=watershed,
        flow_dir_array=flow_directions,
        verbose=False,
    )

    cascades = fa.get_cascades(streams, pour_point, modelgrid, fmt="xy")

    if cascades.ncascade != 3324:
        raise AssertionError(
            "cascade generation did not produced expected number of cascades"
        )

    n = np.count_nonzero(cascades.hru_strmseg_down_id)

    if n != 489:
        raise AssertionError("cascade generation failed")


if __name__ == "__main__":
    test_flow_accumulation_import()
    test_flow_directions()
    test_flow_accumulation()
    test_watershed_delineation()
    test_stream_generation()
    test_cascade_generation()