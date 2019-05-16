# test MODSIM class
import gsflow
import os
import pycrs
import shapefile


def test_modsim():
    ws = "../examples/data/sagehen/gsflow-modsim"
    control_file = "saghen_modsim_cont.control"
    shp = "./temp/test_modsim_modsim.shp"
    prj = "./temp/test_modsim_modsim.prj"
    proj4 = "+proj=utm +zone=10 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

    gsf = gsflow.GsflowModel.load_from_file(os.path.join(ws, control_file))
    modsim = gsflow.modsim.Modsim(gsf)

    if not modsim._ready:
        raise AssertionError

    sfr_topo = modsim.sfr_topology
    lak_topo = modsim.lake_topology

    modsim.write_modsim_shapefile(shp=shp, proj4=proj4)

    crs = pycrs.parse.from_proj4(proj4)

    sf = shapefile.Reader(shp)
    if not sf.numRecords == 16:
        raise AssertionError

    if not sf.shapeType == 3:
        raise AssertionError

    crs2 = pycrs.load.from_file(prj)

    if crs2.to_esri_wkt() != crs.to_esri_wkt():
        raise AssertionError


def test_gsflow_modsim_read_write():
    ws = "../examples/data/sagehen/gsflow-modsim"
    control_file = "saghen_modsim_cont.control"
    ws2 = "./temp"
    shp = "./temp/saghen_modsim_modsim.shp"
    prj = "./temp/saghen_modsim_modsim.prj"
    proj4 = "+proj=utm +zone=10 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

    gsf = gsflow.GsflowModel.load_from_file(os.path.join(ws, control_file))
    gsf.mf.modelgrid.set_coord_info(proj4=proj4)
    gsf.write_input(workspace=ws2)

    crs = pycrs.parse.from_proj4(proj4)

    sf = shapefile.Reader(shp)
    if not sf.numRecords == 16:
        raise AssertionError

    if not sf.shapeType == 3:
        raise AssertionError

    crs2 = pycrs.load.from_file(prj)

    if crs2.to_esri_wkt() != crs.to_esri_wkt():
        raise AssertionError


if __name__ == "__main__":
    test_modsim()
    test_gsflow_modsim_read_write()
