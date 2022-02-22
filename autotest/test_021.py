import os
import platform

ws = os.path.abspath(os.path.dirname(__file__))


def test_import_crt():
    from gsflow.crt import CRT


def test_read_write_crt():
    from gsflow.crt import CRT
    exe_name = os.path.join(ws, "..", "bin", "CRT_1.3.1")
    if platform.system().lower() == "windows":
        exe_name += ".exe"
    model_ws = os.path.join(ws, "..", "examples", "data", "crt")
    out_ws = os.path.join(ws, "temp")

    crt = CRT.load(model_ws)
    crt.write_input(out_ws)
    success, buff = crt.run_model(exe_name)

    if not success:
        raise AssertionError("CRT run failed")


def test_crt_parameter_output():
    from gsflow.crt import CRT
    from gsflow.prms import ParameterRecord

    model_ws = os.path.join(ws, "temp")
    crt = CRT.load(model_ws)
    params = crt.load_cascade_parameters_output()
    if len(params) != 6:
        raise AssertionError("cascade parameters not loaded properly")

    if not isinstance(params[0], ParameterRecord):
        raise TypeError("cascade parameter type is incorrect")

    if params[0].values[0] != 6266:
        raise ValueError("Incorrect value for ncascade dimension")


def test_crt_vis_output():
    from gsflow.crt import CRT
    import numpy as np

    vis_file = os.path.join(ws, 'temp', 'vis.txt')
    shp_name = os.path.join(ws, 'temp', 'vis.shp')
    recarray = CRT.get_vis_output(vis_file, shp_name)

    if not isinstance(recarray, np.recarray):
        raise TypeError("recarray not output from get_vis_output()")

    if not os.path.exists(shp_name):
        raise FileNotFoundError('shapefile not written')


if __name__ == "__main__":
    test_import_crt()
    test_read_write_crt()
    test_crt_parameter_output()
    test_crt_vis_output()