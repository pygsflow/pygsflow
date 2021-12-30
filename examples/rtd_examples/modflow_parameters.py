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

# # Working with MODFLOW packages
#
# This tutorial shows how to access, edit, and add new MODFLOW parameters
# and packages to a GSFLOW model using pyGSFLOW
#
# Note: this is a minimal overview for working with the Modflow object.
# For a more information and examples showing FloPy's capabilities and data
# types please visit the flopy homepage
# [here](https://github.com/modflowpy/flopy)

# Package import
import os
import gsflow
import flopy

# ## Load a demonstration model

model_ws = os.path.join("..", "..", "data", "sagehen", "gsflow")

control_file = os.path.join(model_ws, "saghen_new_cont.control")
gsf = gsflow.GsflowModel.load_from_file(control_file)

# ## Accessing the Modflow object
# The modflow object can be accessed using the mf parameter.

ml = gsf.mf

# ## Getting packages from the Modflow object

dis = gsf.mf.dis

# ## Getting parameter values from the Modflow object

top = ml.dis.top.array

# ## Adjusting parameter values

ml.dis.top *= 1.2

# ## Adding a new package to the Modflow object

spd = {i: [[0, 30, 30, -150.],] for i in range(ml.nper)}
wel = flopy.modflow.ModflowWel(ml, stress_period_data=spd)
