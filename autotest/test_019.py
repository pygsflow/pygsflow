# test exporting capabilities
import gsflow
import os
from gsflow.utils.vtk import Gsflowvtk, Mfvtk

def test_vtk():
    ws = "../examples/data/sagehen/gsflow"
    control_file = "saghen_new_cont.control"
    Gsflowvtk.gsflow_to_vtk(control_file=os.path.join(ws, control_file))




if __name__ == "__main__":
    test_vtk()