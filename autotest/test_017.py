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
    if not sf.numRecords == 17:
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
    if not sf.numRecords == 17:
        raise AssertionError

    if not sf.shapeType == 3:
        raise AssertionError

    crs2 = pycrs.load.from_file(prj)

    if crs2.to_esri_wkt() != crs.to_esri_wkt():
        raise AssertionError


def test_modsim_flag_spillway():
    ws = "../examples/data/sagehen_3lay_modsim/windows"
    ws2 = "./temp"
    control_file = "sagehen_modsim_3lay.control"
    proj4 = "+proj=utm +zone=10 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

    shp_iseg = "modsim_flg_iseg.shp"
    shp_elev = "modsim_flg_elev.shp"
    shp_flow = "modsim_flg_flow.shp"

    gsf = gsflow.GsflowModel.load_from_file(os.path.join(ws, control_file))
    gsf.mf.modelgrid.set_coord_info(proj4=proj4)
    gsf.modsim.write_modsim_shapefile(shp=os.path.join(ws2, shp_iseg),
                                      flag_spillway=[24,], nearest=False)

    gsf.modsim.write_modsim_shapefile(shp=os.path.join(ws2, shp_flow),
                                      flag_spillway='flow')

    gsf.modsim.write_modsim_shapefile(shp=os.path.join(ws2, shp_elev),
                                      flag_spillway='elev', nearest=False)

    sf = shapefile.Reader(os.path.join(ws2, shp_iseg))
    for record in sf.records():
        if record[0] == 24:
            if record[-1] != 1:
                raise AssertionError()

    sf = shapefile.Reader(os.path.join(ws2, shp_flow))
    for record in sf.records():
        if record[0] == 24:
            if record[-1] != 1:
                raise AssertionError()

    sf = shapefile.Reader(os.path.join(ws2, shp_elev))
    for record in sf.records():
        if record[0] == 24:
            if record[-1] != 1:
                raise AssertionError()


if __name__ == "__main__":
    # test_modsim()
    # test_gsflow_modsim_read_write()
    test_modsim_flag_spillway()