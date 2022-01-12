# test exporting capabilities
import gsflow
import os
from gsflow.utils.vtk import Gsflowvtk, Mfvtk


ws = os.path.abspath(os.path.dirname(__file__))


def test_vtk():
    local_ws = os.path.join(ws, "..", "examples", "data", "sagehen", "gsflow")
    control_file = "saghen_new_cont.control"
    Gsflowvtk.gsflow_to_vtk(control_file=os.path.join(local_ws, control_file))


if __name__ == "__main__":
    test_vtk()