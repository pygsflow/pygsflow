# test sfr renumbering schemes and other random utilities
import gsflow
import os
from gsflow.utils import SfrRenumber


def test_sfr_renumber():
    # simple test to ensure no crashes in the renumbering schemes
    # expand this later to test LAK, AG, and GAGE
    ws = "../examples/data/sagehen/gsflow"
    control_file = "saghen_new_cont.control"

    gsf = gsflow.GsflowModel.load_from_file(os.path.join(ws, control_file))
    ml = gsf.mf

    # renumber by topology
    sfrenum = SfrRenumber(model=ml)
    sfrenum.renumber_sfr()
    sfrenum.renumber_all()

    # renumber by dis
    sfrenum = SfrRenumber(model=ml, scheme="dis")
    sfrenum.renumber_sfr()
    sfrenum.renumber_all()

    # renumber by strtop
    sfrenum = SfrRenumber(model=ml, scheme="sfr")
    sfrenum.renumber_sfr()
    sfrenum.renumber_all()


if __name__ == "__main__":
    test_sfr_renumber()