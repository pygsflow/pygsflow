from gsflow import PrmsData
import pandas as pd
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


if __name__ == "__main__":
    test_empty_prms_data()
    test_build_prms_data()
