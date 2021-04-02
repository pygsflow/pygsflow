# test ModflowAg package
import os
from gsflow.modflow import Modflow, ModflowAwu, ModflowAg
import numpy as np


ws = os.path.abspath(os.path.dirname(__file__))


def test_ModflowAg_load_write():
    local_ws = os.path.join(ws, "..", "examples", "data", "sagehen", "modflow")
    agfile = "sagehen.awu"
    nper = 344

    ml = Modflow("agtest", model_ws=local_ws)

    ag = ModflowAwu.load(os.path.join(local_ws, agfile), ml,
                         nper=nper, ext_unit_dict={})

    if not isinstance(ag, ModflowAg):
        raise AssertionError("Override of Awu failed")

    ws2 = os.path.join(ws, "temp")
    ml.change_model_ws(ws2)
    ag.write_file()

    agfile2 = "agtest.ag"
    ml2 = Modflow("agtest2", model_ws=ws2)
    ag2 = ModflowAg.load(os.path.join(ws2, agfile2), ml2,
                          nper=nper, ext_unit_dict={})

    assert repr(ag.options) == repr(ag2.options)

    for ix, rec in enumerate(ag.time_series):
        assert rec == ag2.time_series[ix]

    for ix, rec in enumerate(ag.well_list):
        assert rec == ag2.well_list[ix]

    for per in range(nper):

        for ix, rec in enumerate(ag.irrdiversion[per]):
            rec2 = ag2.irrdiversion[per][ix]
            assert rec['segid'] == rec2['segid']
            assert rec['hru_id0'] == rec2['hru_id0']

        for ix, rec in enumerate(ag.irrwell[per]):
            rec2 = ag2.irrwell[per][ix]
            assert rec['wellid'] == rec2['wellid']
            assert rec["hru_id0"] == rec2["hru_id0"]

        for ix, rec in enumerate(ag.supwell[per]):
            rec2 = ag2.supwell[per][ix]
            assert rec['wellid'] == rec2['wellid']
            assert rec['segid0'] == rec2['segid0']

    if ag.plottable:
        raise AssertionError("ModflowAg should be non-plottable")

    if not ModflowAg._ftype() == "AG":
        raise AssertionError("ModflowAg ftype should be AG")


if __name__ == "__main__":
    test_ModflowAg_load_write()