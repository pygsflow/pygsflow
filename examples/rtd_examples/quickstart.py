# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.5.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Quickstart for pyGSFLOW
#
# This tutorial shows how to load an existing model into pyGSFLOW and create
# a new (empty) GsflowModel object for building a new GSFLOW model

# ## Package import

import os
import gsflow
import flopy
import platform
import pandas as pd

# ## Load an existing GSFLOW model

model_ws = os.path.join("..", "..", "data", "sagehen", "gsflow")

control_file = os.path.join(model_ws, "saghen_new_cont.control")
gsf = gsflow.GsflowModel.load_from_file(control_file)

# ## Load and run a GSFLOW model from pyGSFLOW

control_file = os.path.join(model_ws, "saghen_new_cont.control")

exe_name = os.path.join("..", "..", "bin", "gsflow")
if platform.system().lower() == "windows":
    exe_name += ".exe"

gsf = gsflow.GsflowModel.load_from_file(control_file, gsflow_exe=exe_name)
gsf.run_model()

# ## Loading common output files

# ### stat var file

stat_var = gsf.prms.get_StatVar()

# ### gsflow.csv

csv_name = os.path.join(model_ws, "saghen_new_csv_output.csv")
df = pd.read_csv(csv_name)

# ### head file

head_name = os.path.join(model_ws, "saghen_new.hds")
head = flopy.utils.HeadFile(head_name)

# ### budget file

cbc_name = os.path.join(model_ws, "saghen_new.cbc")
cbc = flopy.utils.CellBudgetFile(cbc_name)

# ## Create an empty GsflowModel object for building a new model

control = gsflow.ControlFile(records_list=[])
gsf = gsflow.GsflowModel(control)
