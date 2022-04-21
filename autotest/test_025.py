from gsflow import GsflowModel
from gsflow.prms import ParameterRecord
import os
import numpy as np


ws = os.path.abspath(os.path.dirname(__file__))
model_ws = os.path.join(ws, "..", "examples", "data", "sagehen", "gsflow")
control_file = os.path.join(model_ws, "saghen_new_cont.control")


def test_parameter_getattr():
    gsf = GsflowModel.load_from_file(control_file)
    params = gsf.prms.parameters

    ssr2gw_rate = params.ssr2gw_rate
    if not isinstance(ssr2gw_rate, ParameterRecord):
        raise AssertionError()

    names = params.record_names
    if not isinstance(names, list):
        raise AssertionError()


def test_parameter_setattr():
    gsf = GsflowModel.load_from_file(control_file)
    params = gsf.prms.parameters

    ssr2gw_rate = params.ssr2gw_rate
    t = np.ones(ssr2gw_rate.values.shape)

    params.ssr2gw_rate = t

    if not isinstance(params.ssr2gw_rate, ParameterRecord):
        raise AssertionError()

    if not np.allclose(params.ssr2gw_rate.values, t):
        raise AssertionError()


def test_parameter_math_float():
    gsf = GsflowModel.load_from_file(control_file)
    params = gsf.prms.parameters

    ssr2gw_rate = params.ssr2gw_rate
    valid = ssr2gw_rate.values.copy()

    ssr2gw_rate += 2
    if not np.allclose(ssr2gw_rate.values, valid + 2):
        raise AssertionError("Addition failed")

    ssr2gw_rate.values = valid.copy()
    ssr2gw_rate -= 2
    if not np.allclose(ssr2gw_rate.values, valid - 2):
        raise AssertionError("Subtraction failed")

    ssr2gw_rate.values = valid.copy()
    ssr2gw_rate *= 2
    if not np.allclose(ssr2gw_rate.values, valid * 2):
        raise AssertionError("Multiplication failed")

    ssr2gw_rate.values = valid.copy()
    ssr2gw_rate /= 2
    if not np.allclose(ssr2gw_rate.values, valid / 2):
        raise AssertionError("Division failed")

    ssr2gw_rate.values = valid.copy()
    ssr2gw_rate **= 2
    if not np.allclose(ssr2gw_rate.values, valid ** 2):
        raise AssertionError("Power failed")


def test_parameter_slice():
    gsf = GsflowModel.load_from_file(control_file)
    params = gsf.prms.parameters

    ssr2gw_rate = params.ssr2gw_rate
    valid = ssr2gw_rate.values.copy()
    x = np.sum(valid)
    ssr2gw_rate[0:150] += 1

    x1 = np.sum(ssr2gw_rate.values)

    if not np.isclose(x1, x + 150):
        raise AssertionError("Slicing failed")

    for i in range(150):
        if np.abs(ssr2gw_rate.values[i] - valid[i]) > 1.000001:
            raise AssertionError("Slicing failed")


def test_control_getattr():
    gsf = GsflowModel.load_from_file(control_file)
    control = gsf.control

    record = control.model_mode

    if record.values[0].lower() != "gsflow":
        raise AssertionError("ControlFile getattr failed")

    control.model_mode.values = ["MODFLOW",]

    if control.model_mode.values[0].lower() != "modflow":
        raise AssertionError("ControlFile set_values failed")


def test_control_setattr():
    gsf = GsflowModel.load_from_file(control_file)
    control = gsf.control

    record = control.model_mode

    if record.values[0].lower() != "gsflow":
        raise AssertionError("ControlFile getattr failed")

    control.model_mode = ["MODFLOW", ]

    if control.model_mode.values[0].lower() != "modflow":
        raise AssertionError("ControlFile set_values failed")


if __name__ == "__main__":
    test_parameter_getattr()
    test_parameter_setattr()
    test_parameter_math_float()
    test_parameter_slice()
    test_control_getattr()
    test_control_setattr()
