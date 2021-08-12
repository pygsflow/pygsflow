import os
import numpy as np

ws = os.path.abspath(os.path.dirname(__file__))


def test_fishnet_import():
    from gsflow.builder import GenerateFishnet


def test_fishnet_generation():
    from gsflow.builder import GenerateFishnet
    from flopy.discretization import StructuredGrid
    try:
        import rasterio
    except ImportError:
        return

    raster = os.path.join(
        ws, '..', 'examples', 'data', 'geospatial', 'dem.img'
    )
    modelgrid = GenerateFishnet(raster, xcellsize=50, ycellsize=100, buffer=0)

    if not isinstance(modelgrid, StructuredGrid):
        raise TypeError('GenerateFishnet should return SturcturedGrid object')

    if modelgrid.ncol != 149 or modelgrid.nrow != 69:
        raise AssertionError("DELC and DELR have not been properly calculated")


def test_import_modflow_builder():
    from gsflow.builder import ModflowBuilder


def test_modflow_builder():
    from gsflow.builder import ModflowBuilder, Defaults
    from gsflow.modflow import Modflow
    import flopy

    model_ws = os.path.join(ws, 'temp')
    model_name = "mfbuilder_model"
    defaults = Defaults().modflow.to_dict()

    LxLy = 100
    nrow = ncol = 10
    delc = delr = np.array([LxLy / nrow] * nrow)
    top = np.ones((nrow, ncol)) * 50
    ibound = np.ones((nrow, ncol), dtype=int)

    modelgrid = flopy.discretization.StructuredGrid(
        delc,
        delr,
        top=top,
        idomain=ibound,
        nlay=1,
        nrow=nrow,
        ncol=ncol
    )

    mf_builder = ModflowBuilder(modelgrid, top, model_name)
    mf_builder.build_dis()
    mf_builder.build_bas6(ibound)
    mf_builder.build_upw()
    mf_builder.build_oc()
    mf_builder.build_nwt()

    # create a fake stream network for sfr example,
    # assume each segment has 1 reach
    seg_connect = {i: i + 1 for i in range(1, 11)}
    seg_connect[10] = 0
    nsegments = nreaches = len(seg_connect)

    # segments connect in a diagonal across the model
    seg_ij = {i: (i - 1, i - 1) for i in seg_connect.keys()}
    rchlen = [np.sqrt(10 ** 2 + 10 ** 2)] * nreaches
    strtop = [top[v[0], v[1]] for k, v in seg_ij.items()]
    slope = [0.1] * nreaches

    # get an empty recarray for reach_data (nreaches = nseg in this example)
    reach_data = flopy.modflow.ModflowSfr2.get_empty_reach_data(nreaches)
    reach_data["iseg"] = list(sorted(seg_connect.keys()))
    reach_data['ireach'] = [1] * nreaches
    reach_data['i'] = [v[0] for k, v in sorted(seg_ij.items())]
    reach_data['j'] = [v[1] for k, v in sorted(seg_ij.items())]
    reach_data['rchlen'] = rchlen
    reach_data['strtop'] = strtop
    reach_data['slope'] = slope

    # apply sfr reach defaults
    reach_defaults = defaults['sfr']['reach']
    for key, val in reach_defaults.items():
        val = [val] * nreaches
        reach_data[key] = val

    # build the segment_data array
    segment_data = flopy.modflow.ModflowSfr2.get_empty_segment_data(nsegments)
    segment_data['nseg'] = list(sorted(seg_connect.keys()))
    segment_data['outseg'] = [v for k, v in sorted(seg_connect.items())]

    # apply segment defaults
    segment_defaults = defaults['sfr']['segment']
    for key, val in segment_defaults.items():
        val = [val] * nsegments
        segment_data[key] = val

    mf_builder.build_sfr(reach_data, segment_data)

    irunbnd = np.ones((nrow, ncol), dtype=int)
    for col in range(ncol):
        irunbnd[:, col] *= (col + 1)

    mf_builder.build_uzf(irunbnd)
    ml = mf_builder.model
    ml.change_model_ws(model_ws)
    ml.write_input()

    ml2 = Modflow.load(model_name + ".nam", model_ws=model_ws)
    if len(ml2.packagelist) < 7:
        raise AssertionError("Model was not built properly")


if __name__ == "__main__":
    test_fishnet_import()
    test_fishnet_generation()
    test_import_modflow_builder()
    test_modflow_builder()
