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

# # Working with GSFLOW control file, prms parameters, and modflow
#
# This tutorial gives an overview on how to access, edit, remove, and add data
# using pyGSFLOW on a GSFLOW model.

# ## Working with GSFLOW control file parameters
#
# This section shows how to access, edit, remove, and add new control file
# parameters and packages to a GSFLOW model using pyGSFLOW

# Package import
import os
import gsflow
import numpy as np
import flopy

# ### Load a demonstration model

model_ws = os.path.join("..", "..", "data", "sagehen", "gsflow")

control_file = os.path.join(model_ws, "saghen_new_cont.control")
gsf = gsflow.GsflowModel.load_from_file(control_file)

# ### Accessing the ControlFile object

control = gsf.control

# ### Getting parameter values from the ControlFile object
# The `get_values()` method allows the user to get a list of parameter values

csv_out = control.get_values("csv_output_file")

# ### Adjusting parameter values
# The `set_values()` method allows the user to adjust control file parameter
# values

csv_out = "gsflow_example.csv"
control.set_values("csv_output_file", [csv_out,])

# ### Removing a parameter
# The `remove_record()` method will remove a parameter value from the ControlFile
# object

control.remove_record("csv_output_file")

# ### Adding a new parameter
# The `add_record()` method allows users to add new records to the ControlFile
# object.

csv_out = "gsflow_example.csv"
control.add_record("csv_output_file", [csv_out,])

# ## Working with PRMS parameters
#
# This section shows how to access, edit, remove, and add new PRMS
# parameters and packages to a GSFLOW model using pyGSFLOW

# ### Load a demonstration model

control_file = os.path.join(model_ws, "saghen_new_cont.control")
gsf = gsflow.GsflowModel.load_from_file(control_file)

# ### Accessing the PrmsParameters object
# The PrmsParameters object allows the user to access data from PRMS
# parameter files.

params = gsf.prms.parameters

# ### Getting a list all of the model parameters
# The `parameters_list` method returns all parameter names in the
# PrmsParameters object

param_names = params.parameters_list

# ### Getting parameter values from the PrmsParameters object
# The `get_values()` method returns a numpy array of parameter values
# for PRMS parameters

ssr2gw = params.get_values("ssr2gw_rate")

# ### Adjusting parameter values
# The `set_values()` method allows the user to adjust a PRMS parameter value

ssr2gw = params.get_values("ssr2gw_rate")
ssr2gw *= 0.8
params.set_values("ssr2gw_rate", ssr2gw)

# ### Removing a parameter
# The `remove_record()` method allows the user to remove a parameter from
# the PrmsParameters object.

params.remove_record("ssr2gw_rate")

# ### Adding a new parameter
# The `add_record()` method allows the user to add new parameters to the
# PrmsParameters object
#
# `add_record()` parameters include
#
#    - `name`: (str) parameter name
#
#    - `values`: (np.array, list) numpy array or list of parameter values
#
#    - `dimensions`: (list) 2 dimensional list that defines 1) the dimension
# name and 2) the dimension size
#
#    - `datatype`: (int) PRMS data type
#
#    - `filename`: (str) Optional parameter that allows the user to set which
# PRMS parameter file a given parameter is written to. Default is None which
# writes to the parameter file that contains the PRMS dimensions block.

nhru = gsf.mf.nrow * gsf.mf.ncol
ssr2gw = np.random.random(nhru)
params.add_record("ssr2gw_rate",
                  ssr2gw,
                  dimensions=[["nhru", nhru]],
                  datatype=3)

# ## Working with MODFLOW packages
#
# This tutorial shows how to access, edit, and add new MODFLOW parameters
# and packages to a GSFLOW model using pyGSFLOW
#
# Note: this is a minimal overview for working with the Modflow object.
# For a more information and examples showing FloPy's capabilities and data
# types please visit the flopy homepage
# [here](https://github.com/modflowpy/flopy)

# ### Load a demonstration model

control_file = os.path.join(model_ws, "saghen_new_cont.control")
gsf = gsflow.GsflowModel.load_from_file(control_file)

# ### Accessing the Modflow object
# The modflow object can be accessed using the mf parameter.

ml = gsf.mf

# ### Getting packages from the Modflow object

dis = gsf.mf.dis

# ### Getting parameter values from the Modflow object

top = ml.dis.top.array

# ### Adjusting parameter values

ml.dis.top *= 1.2

# ### Adding a new package to the Modflow object

spd = {i: [[0, 30, 30, -150.],] for i in range(ml.nper)}
wel = flopy.modflow.ModflowWel(ml, stress_period_data=spd)
