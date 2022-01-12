import pandas as pd
import os
from gsflow import GsflowModel
from gsflow.output import StatVar


ws = os.path.abspath(os.path.dirname(__file__))


def test_read_stat_var():
    local_ws = os.path.join(ws, "temp")
    name = "saghen_new_stat_var.dat"

    sv = StatVar(os.path.join(local_ws, name))
    assert isinstance(sv.stat_df, pd.DataFrame)


def test_read_stat_var_from_control():
    local_ws = os.path.join(ws, "temp")
    control_file = "saghen_new_cont.control"

    gs = GsflowModel.load_from_file(os.path.join(local_ws, control_file))
    df = gs.prms.get_StatVar()
    assert isinstance(df, pd.DataFrame)


if __name__ == "__main__":
    test_read_stat_var()
    test_read_stat_var_from_control()