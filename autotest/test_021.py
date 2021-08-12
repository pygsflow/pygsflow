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


if __name__ == "__main__":
    test_import_crt()
    test_read_write_crt()