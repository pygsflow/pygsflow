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

# # Working with PRMS parameters
#
# This tutorial shows how to access, edit, remove, and add new PRMS
# parameters and packages to a GSFLOW model using pyGSFLOW

# Package import
import os
import gsflow
import numpy as np

# ## load a demonstration model

model_ws = os.path.join("..", "..", "data", "sagehen", "gsflow")

control_file = os.path.join(model_ws, "saghen_new_cont.control")
gsf = gsflow.GsflowModel.load_from_file(control_file)

# ## Accessing the PrmsParameters object
# The PrmsParameters object allows the user to access data from PRMS
# parameter files.

params = gsf.prms.parameters

# ## Getting a list all of the model parameters
# The `parameters_list` method returns all parameter names in the
# PrmsParameters object

param_names = params.parameters_list

# ## Getting parameter values from the PrmsParameters object
# The `get_values()` method returns a numpy array of parameter values
# for PRMS parameters

ssr2gw = params.get_values("ssr2gw_rate")

# ## Adjusting parameter values
# The `set_values()` method allows the user to adjust a PRMS parameter value

ssr2gw = params.get_values("ssr2gw_rate")
ssr2gw *= 0.8
params.set_values("ssr2gw_rate", ssr2gw)

# ## Removing a parameter
# The `remove_record()` method allows the user to remove a parameter from
# the PrmsParameters object.

params.remove_record("ssr2gw_rate")

# ## Adding a new parameter
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
