import gsflow
import os
import numpy as np
from gsflow import ControlFile, PrmsParameters, PrmsData


def test_load_write_model_prms_only():
    ws = "../examples/data/sagehen/prms/windows"
    control_file = "sagehen.control"
    gs = gsflow.GsflowModel.load_from_file(os.path.join(ws, control_file))
    assert isinstance(gs.control, ControlFile)
    assert isinstance(gs.prms.parameters, PrmsParameters)
    assert isinstance(gs.prms.data, PrmsData)

    ws2 = "./temp"
    gs.write_input(workspace=ws2)

    gs2 = gsflow.GsflowModel.load_from_file(os.path.join(ws2, control_file))

    assert len(gs2.control.record_names) == len(gs.control.record_names)
    assert len(gs2.prms.parameters.record_names) == len(gs.prms.parameters.record_names)


def test_load_write_gsflow_modflow():
    ws = "../examples/data/sagehen/modflow"
    nam = "saghen.nam"
    ml = gsflow.modflow.Modflow.load(nam, model_ws=ws)
    assert isinstance(ml, gsflow.modflow.Modflow)

    ws2 = "./temp"
    ml.change_model_ws(ws2)
    ml.write_input()

    ml2 = gsflow.modflow.Modflow.load(nam, model_ws=ws2)

    assert len(ml.packagelist) == len(ml2.packagelist)
    assert ml.nrow_ncol_nlay_nper == ml2.nrow_ncol_nlay_nper


def test_load_write_gsflow():
    # todo: this must be changed to write relative paths, not absolutes!
    ws = "../examples/data/sagehen/gsflow"
    control_file = "saghen_new_cont.control"

    gs = gsflow.GsflowModel.load_from_file(os.path.join(ws, control_file))
    assert isinstance(gs.control, ControlFile)
    assert isinstance(gs.prms.parameters, PrmsParameters)
    assert isinstance(gs.prms.data, PrmsData)
    assert isinstance(gs.mf, gsflow.modflow.Modflow)

    ws2 = "./temp"

    # change ws only ...
    gs.write_input(workspace=ws2)

    gs2 = gsflow.GsflowModel.load_from_file(os.path.join(ws2, control_file))
    assert len(gs2.control.record_names) == len(gs.control.record_names)
    assert len(gs2.prms.parameters.record_names) == len(gs.prms.parameters.record_names)
    assert len(gs2.mf.packagelist) == len(gs.mf.packagelist)
    assert gs2.mf.nrow_ncol_nlay_nper == gs.mf.nrow_ncol_nlay_nper

    ws = "../examples/data/sagehen/gsflow"
    control_file = "saghen_new_cont.control"

    gs = gsflow.GsflowModel.load_from_file(os.path.join(ws, control_file))
    assert not gs.modflow_only
    assert not gs.prms_only

    ws2 = "./temp"
    basename = "test2"
    gs.write_input(basename=basename, workspace=ws2)

    gs2 = gsflow.GsflowModel.load_from_file(os.path.join(ws2, "test2_cont.control"))
    assert len(gs2.control.record_names) == len(gs.control.record_names)
    assert len(gs2.prms.parameters.record_names) == len(gs.prms.parameters.record_names)
    assert len(gs2.mf.packagelist) == len(gs.mf.packagelist)
    assert gs2.mf.nrow_ncol_nlay_nper == gs.mf.nrow_ncol_nlay_nper





if __name__ == "__main__":
    #test_load_write_model_prms_only()
    #test_load_write_gsflow_modflow()
    test_load_write_gsflow()
