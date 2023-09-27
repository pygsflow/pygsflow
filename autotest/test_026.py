from gsflow import GsflowModel
import os


ws = os.path.abspath(os.path.dirname(__file__))
model_ws = os.path.join(ws, "..", "examples", "data", "sagehen", "gsflow")
control_file = os.path.join(model_ws, "saghen_new_cont.control")


def test_upslope_neighbors():
    gsf = GsflowModel.load_from_file(control_file)
    params = gsf.prms.parameters
    conn, pct = params.upslope_neighbors()
    if not isinstance(conn, dict):
        raise AssertionError()
    if not isinstance(pct, dict):
        raise AssertionError


def test_upslope_neighbors_with_hru():
    gsf = GsflowModel.load_from_file(control_file)
    params = gsf.prms.parameters
    conn, pct = params.upslope_neighbors(5267)
    assert conn == [5351, 5352, 5350, 5266, 5182]
    assert pct == [0.253, 0.017, 0.292, 0.192, 0.089]


if __name__ == '__main__':
    test_upslope_neighbors()
    test_upslope_neighbors_with_hru()

