# test exporting capabilities
import gsflow
import os


ws = os.path.abspath(os.path.dirname(__file__))


def test_netcdf_export():
    local_ws = os.path.join(ws, "..", "examples", "data", "sagehen", "gsflow-modsim")
    control_file = "saghen_modsim_cont.control"

    gsf = gsflow.GsflowModel.load_from_file(os.path.join(local_ws, control_file))
    gsf.mf.modelgrid.set_coord_info(
        xoff=438979.0,
        yoff=3793007.75,
        proj4="+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
    )
    # gsf.prms.parameters.export_nc("test.nc", gsf.mf)
    gsf.export_nc(os.path.join(ws, "temp", "test.nc"))


if __name__ == "__main__":
    test_netcdf_export()