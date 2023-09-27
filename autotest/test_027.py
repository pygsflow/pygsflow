import gsflow
import os
import pytest

ws = os.path.abspath(os.path.dirname(__file__))
model_ws = os.path.join("..", "examples", "data", "sagehen", "50m_tutorials")
control = os.path.join(model_ws, "sagehen_50m_initial.control")


def test_upslope_streams():
    gsf = gsflow.GsflowModel.load_from_file(control)
    params = gsf.prms.parameters
    strms = params.upslope_streams()
    if not isinstance(strms, dict):
        raise AssertionError("Upslope_streams returned incorrect type")


def test_upslope_streams_with_hru():
    gsf = gsflow.GsflowModel.load_from_file(control)
    params = gsf.prms.parameters
    strms = params.upslope_streams(stream_hru=10787)
    print(strms)
    if not isinstance(strms, dict):
        raise AssertionError("Upslope_streams returned incorrect type")
    if not len(strms) == 41:
        raise AssertionError("Upslope streams returned stream_conn_stack of incorrect length")
    assert strms[10787] == [10935]


def test_upslope_streams_no_cascades():
    gsf = gsflow.GsflowModel.load_from_file(control)
    params = gsf.prms.parameters
    params.hru_strmseg_down_id[:] = 0
    with pytest.raises(AssertionError):
        params.upslope_streams()


if __name__ == '__main__':
    test_upslope_streams()
    test_upslope_streams_with_hru()
    test_upslope_streams_no_cascades()
