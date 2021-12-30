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

# # Working with GSFLOW control file parameters
#
# This tutorial shows how to access, edit, remove, and add new control file
# parameters and packages to a GSFLOW model using pyGSFLOW

# Package import
import os
import gsflow

# ## Load a demonstration model

model_ws = os.path.join("..", "..", "data", "sagehen", "gsflow")

control_file = os.path.join(model_ws, "saghen_new_cont.control")
gsf = gsflow.GsflowModel.load_from_file(control_file)

# ## Accessing the ControlFile object

control = gsf.control

# ## Getting parameter values from the ControlFile object
# The `get_values()` method allows the user to get a list of parameter values

csv_out = control.get_values("csv_output_file")

# ## Adjusting parameter values
# The `set_values()` method allows the user to adjust control file parameter
# values

csv_out = "gsflow_example.csv"
control.set_values("csv_output_file", [csv_out,])

# ## Removing a parameter
# The `remove_record()` method will remove a parameter value from the ControlFile
# object

control.remove_record("csv_output_file")

# ## Adding a new parameter
# The `add_record()` method allows users to add new records to the ControlFile
# object.

csv_out = "gsflow_example.csv"
control.add_record("csv_output_file", [csv_out,])
