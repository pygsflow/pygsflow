from gsflow import PrmsData, GsflowModel
from gsflow.prms import PrmsDay
import pandas as pd
import numpy as np
import os


ws = os.path.abspath(os.path.dirname(__file__))


def test_empty_prms_data():
    data = PrmsData({})

    assert isinstance(data, PrmsData)


def test_build_prms_data():
    data_df = pd.DataFrame()
    name = "test.data"
    model_dir = os.path.join(ws, ".")
    header = "this is a test"
    data = PrmsData(data_df, name=name, model_dir=model_dir, header=header)

    assert data.file_name == os.path.join(model_dir, name)
    assert data.header == header
    assert isinstance(data.data_df, pd.DataFrame)


def test_load_write_prms_day():
    model_dir = os.path.join(ws, "..", "examples", "data", "day")
    ows = os.path.join(ws, "temp")
    day_file = os.path.join(model_dir, "prcp.cbh")
    day = PrmsDay.load_from_file(day_file)
    df = day.dataframe
    if not isinstance(df, pd.DataFrame):
        raise AssertionError("day file not dataframe, error")

    day.change_file_ws(ows)
    day.write()

    day2 = PrmsDay.load_from_file(os.path.join(ows, "prcp.cbh"))
    df2 = day2.dataframe

    if not isinstance(df, pd.DataFrame):
        raise AssertionError("day file write or load error")

    for col in list(df2):
        if not np.allclose(df[col].values, df2[col].values):
            raise ValueError("Day file data is not consistent")


def test_load_write_prms_day_model():
    model_dir = os.path.join(ws, "..", "examples", "data", "day", "model")
    ows = os.path.join(ws, "temp")
    control_file = os.path.join(model_dir, "control_HW17_singlerun")
    gsf = GsflowModel.load_from_file(control_file, prms_only=True)

    gsf.write_input(workspace=ows)

    gsf2 = GsflowModel.load_from_file(
        os.path.join(ows, "control_HW17_singlerun"), prms_only=True
    )

    if not isinstance(gsf2.prms.day, dict):
        raise AssertionError("Day file not properly written to output")

    if "precip" not in gsf2.prms.day:
        raise AssertionError("PRMS variable type is incorrect for day file")

    df = gsf.prms.day['precip'].dataframe
    df2 = gsf2.prms.day["precip"].dataframe

    for col in list(df2):
        if not np.allclose(df[col].values, df2[col].values):
            raise ValueError("Day file data is not consistent")


if __name__ == "__main__":
    test_empty_prms_data()
    test_build_prms_data()
    test_load_write_prms_day()
    test_load_write_prms_day_model()