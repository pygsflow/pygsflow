import pandas as pd
import os
from gsflow import GsflowModel
from gsflow.output import StatVar


def test_read_stat_var():
    ws = "./temp"
    name = "saghen_new_stat_var.dat"

    sv = StatVar(os.path.join(ws, name))
    assert isinstance(sv.stat_df, pd.DataFrame)


def test_read_stat_var_from_control():
    ws = "./temp"
    control_file = "saghen_new_cont.control"

    gs = GsflowModel.load_from_file(os.path.join(ws, control_file))
    df = gs.prms.get_StatVar()
    assert isinstance(df, pd.DataFrame)


if __name__ == "__main__":
    test_read_stat_var()
    test_read_stat_var_from_control()