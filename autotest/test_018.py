# test exporting capabilities
import gsflow
import os


def test_netcdf_export():
    ws = "../examples/data/sagehen/gsflow-modsim"
    control_file = "saghen_modsim_cont.control"

    gsf = gsflow.GsflowModel.load_from_file(os.path.join(ws, control_file))
    gsf.mf.modelgrid.set_coord_info(xoff=438979.0, yoff=3793007.75,
                                    proj4="+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
    # gsf.prms.parameters.export_nc("test.nc", gsf.mf)
    gsf.export_nc(os.path.join("temp", "test.nc"))


if __name__ == "__main__":
    test_netcdf_export()