from gsflow import PrmsParameters
from gsflow.param_base import ParameterBase
from gsflow import ParameterRecord
from gsflow import ControlFile
from gsflow import ControlRecord
import numpy as np
import os


ws = os.path.abspath(os.path.dirname(__file__))


def test_empty_prms_parameter_object():
    pp = PrmsParameters([])

    assert isinstance(pp, PrmsParameters)
    assert isinstance(pp, ParameterBase)


def test_empty_control_file_object():
    cf = ControlFile([])

    assert isinstance(cf, ControlFile)
    assert isinstance(cf, ParameterBase)


def test_build_prms_parameter_object():
    file_name = "test.param"
    name = "ssr2gw_rate 0"
    values = np.random.randn(128)
    dimensions = [["nssr", 128]]
    datatype = 2
    pr = ParameterRecord(name, values=values, dimensions=dimensions,
                         datatype=datatype, file_name=file_name)

    pp = PrmsParameters([pr], header="this is a test")

    assert pp.parameter_files[0] == file_name
    assert pp.record_names[0] == name
    assert isinstance(pp.parameters_list[0], ParameterRecord)
    assert np.allclose(pp.parameters_list[0].values, pr.values)


def test_build_control_file_object():
    cf_name = "test.control"
    model_dir = os.path.join(ws, ".")
    name = "model_mode"
    datatype = 4
    values = ["GSFLOW"]
    cr = ControlRecord(name=name, values=values, datatype=datatype)

    control = ControlFile([cr], name=cf_name, model_dir=model_dir,
                          header="this is a test")

    assert control.control_file == os.path.join(model_dir, cf_name)
    assert control.record_names[0] == name
    assert control.records_list[0].values[0] == values[0]
    assert control.records_list[0].datatype == datatype
    assert isinstance(control.records_list[0], ControlRecord)


def test_load_prms_parameter_duplicates():
    iws = os.path.join(ws, "..", "examples", "data", "sagehen", "duplicates")

    param_files = [os.path.join(iws, "sagehen_duplicates.params"),]
    params = PrmsParameters.load_from_file(param_files)

    adjmix_rain = params.adjmix_rain.values
    if not np.all(adjmix_rain == adjmix_rain[0]):
        raise AssertionError("parameter replacement failed")


def test_load_control_file_duplicates():
    iws = os.path.join(ws, "..", "examples", "data", "sagehen", "duplicates")

    control_file = os.path.join(iws, "sagehen_duplicates.control")
    control = ControlFile.load_from_file(control_file)

    pc_flag = control.parameter_check_flag.values
    if not pc_flag[0] == -2:
        raise AssertionError("controlfile parameter replacement failed")


if __name__ == "__main__":
    test_empty_prms_parameter_object()
    test_empty_control_file_object()
    test_build_prms_parameter_object()
    test_build_control_file_object()
    test_load_prms_parameter_duplicates()
    test_load_control_file_duplicates()
