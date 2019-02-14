from gsflow import PrmsModel
from gsflow import PrmsData, PrmsParameters, ParameterRecord, ControlRecord, ControlFile
import os
import numpy as np
import pandas as pd


def test_build_prms_model():
    # buld control file object
    cf_name = "test.control"
    model_dir = "./"
    name = "model_mode"
    datatype = 4
    cr_values = ["GSFLOW"]
    cr = ControlRecord(name=name, values=cr_values, datatype=datatype)

    control = ControlFile([cr], name=cf_name, model_dir=model_dir,
                          header="this is a test")

    # build parameter object
    file_name = "test.param"
    name = "ssr2gw_rate 0"
    values = np.random.randn(128)
    dimensions = [["nssr", 128]]
    datatype = 2
    pr = ParameterRecord(name, values=values, dimensions=dimensions,
                         datatype=datatype, file_name=file_name)

    parameters = PrmsParameters([pr], header="this is a test")

    # build Data object
    data_df = pd.DataFrame()
    name = "test.data"
    model_dir = "./"
    header = "this is a test"
    data = PrmsData(data_df, name=name, model_dir=model_dir, header=header)

    prms = PrmsModel(control, parameters=parameters, data=data)

    assert isinstance(prms.control, ControlFile)
    assert isinstance(prms.parameters, PrmsParameters)
    assert isinstance(prms.data, PrmsData)
    assert prms.control_file == os.path.join(model_dir, cf_name)
    assert prms.parameters == parameters
    assert prms.data == data
    assert prms.control == control
    assert np.allclose(prms.parameters.parameters_list[0].values, values)
    assert prms.control.records_list[0].values[0] == cr_values[0]


if __name__ == "__main__":
    test_build_prms_model()
