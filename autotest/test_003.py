from gsflow import ControlRecord, ParameterRecord
from gsflow.record_base import RecordBase
import numpy as np


def test_empty_parameter_record():
    pr = ParameterRecord(name="test", values=[1])
    assert isinstance(pr, ParameterRecord)
    assert isinstance(pr, RecordBase)


def test_empty_control_record():
    cr = ControlRecord(name="test", values=[1])
    assert isinstance(cr, ControlRecord)
    assert isinstance(cr, RecordBase)


def test_build_parameter_record():
    name = "ssr2gw_rate 0"
    values = np.random.randn(128)
    dimensions = [["nssr", 128]]
    datatype = 2
    pr = ParameterRecord(name, values=values, dimensions=dimensions,
                         datatype=datatype)
    assert pr.ndim == 1
    assert pr.name == name
    assert np.allclose(values, pr.values)
    assert pr.datatype == datatype


def test_build_control_record():
    name = "model_mode"
    datatype = 4
    values = ["GSFLOW"]
    cr = ControlRecord(name=name, values=values, datatype=datatype)

    assert cr.name == name
    assert cr.nvalues == len(values)
    assert cr.datatype == datatype
    assert cr.values[0] == values[0]

if __name__ == "__main__":
    test_empty_parameter_record()
    test_empty_control_record()
    test_build_parameter_record()
    test_build_control_record()
