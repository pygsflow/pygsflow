---
title: 'Integrated hydrologic model development and postprocessing for GSFLOW using pyGSFLOW'  
tags:
  - Python
  - hydrology
  - integrated hydrologic modeling
  - groundwater
  - surfacewater
  - GSFLOW
  - MODFLOW
  - PRMS  
authors:
  - name: Joshua D. Larsen^[corresponding author]
    ocrid: 0000-0002-1218-800X
    affiliation: 1
  - name: Ayman Alzraiee
    ocrid: 0000-0001-7576-3449
    affiliation: 1
  - name: Richard G. Niswonger
    ocrid: 0000-0001-6397-2403
    affiliation: 2  
affiliations:
  - name: U.S. Geological Survey, California Water Science Center, Sacramento, CA
    index: 1
  - name: U.S. Geological Survey, Integrated Modeling and Prediction Division, Water Mission Area, Menlo Park, CA 
    index: 2  
date: 18 February 2022
bibliography: paper.bib
---

# Overview
pyGSFLOW is a python package designed to create new GSFLOW integrated hydrologic 
models, read existing models, edit model input data, run GSFLOW models, process 
output, and visualize model data.

# Introduction
Hydrologic modeling has been steadily increasing in complexity over the years. 
The addition of new types of data sets and the need to represent the interaction 
between surface-water and groundwater systems has been driving this complexity. 
Ignoring landscape changes in space and time and applying simple boundary conditions 
within hydrologic models is no longer adequate to address watershed and basin scale 
issues [@Fatichi:2016]. Instead, integrated hydrologic models (IHMs), that couple 
governing equations for surface-water and groundwater flow, are used to represent 
feedback mechanisms between these systems.

Among these IHMs is GSFLOW simulation code [@Markstrom:2008] that simulates surface 
and subsurface hydrologic processes by integrating the Precipitation Runoff 
Modeling System (PRMS)  [@Markstrom:2015] and MODFLOW [@Harbaugh:2005; @Niswonger:2011] 
into a single code that simulates feedbacks between the two systems. Because 
modelers are moving toward simulating greater portions of the hydrologic cycle, 
larger datasets, from multiple sources are used to parameterize these models. 
Beyond the scope of example problems, most applied problems require custom workflows 
and code to process large datasets related to model inputs and outputs. Scripting 
languages like Python make it easier to process large data sets and provide standard 
methods that can be used for developing, editing, and properly formatting model 
input files and for analyzing model output data. These developments have led to 
major advancements in model reproducibility and improvements in model 
applicability [@Bakker:2016; @Gardner:2018; @Ng:2018].

# Statement of need
GSFLOW model development previously has been a piecemeal approach. Arcpy-GSFLOW 
scripts [@Gardner:2018] or GSFLOW-GRASS packages have been used to 
process surface-water input data into model files. PRMS-Python 
[@Volk:2019] could be used to edit most of the PRMS inputs to 
GSFLOW. Finally, FloPy [@Bakker:2016; @Bakker:2021] could be used to edit most 
of the MODFLOW inputs to GSFLOW. This approach is not tightly coupled and still 
requires manual edits and additional external scripts to edit, run models, and 
process output data. Because of the complexity of integrated hydrologic models 
and the need for model reproducibility, a single integrated scripting package 
will help standardize and streamline model development and calibration. 

# pyGSFLOW
pyGSFLOW [@Larsen:2021] is a Python package for creating new GSFLOW models, 
importing existing models, running GSFLOW models, processing model outputs, 
and visualizing model data. Instead of working directly with formatted model 
input files, the pyGSFLOW Application Programming Interface (API) allows the 
user to work with class-based methods to create GSFLOW [@Markstrom:2008], 
PRMS [@Markstrom:2015], MODSIM [@Labadie:2006] vectorized surface-water 
operations networks, and MODFLOW [@Harbaugh:2005] model packages and binds them 
into a single integrated model instance. pyGSFLOW leverages features from FloPy, 
an existing Python package for the MODFLOW suite of groundwater modeling software
[@Harbaugh:2005; @Niswonger:2011; @Panday:2013; @Langevin:2017] 
and extends the capabilities for integrated hydrologic models. pyGSFLOW relies 
on FloPy model and package objects and interfaces with these features to provide 
FloPy users with familiar code syntax and to ensure the long-term maintainability 
of the code base.

The pyGSFLOW package was developed for hydrologic modelers and researchers who 
are developing, calibrating, or applying models as part of their scientific 
investigation (e.g., prediction scenarios, stream capture) with GSFLOW. The code 
base currently is being used for the development of several river-basin scale 
hydrologic models in western basins, including the example shown below 
highlighting application to the Russian River basin and the Santa Rosa Plain 
Watershed (figure 1) [@Woolfenden:2014; @Gardner:2018]. 

The Santa Rosa Plain (SRP) model [@Woolfenden:2014] ) is an IHM that was developed 
as a tool to provide scientific information to water managers about future 
climate change scenarios. The SRP model applied four global-climate models and 
simulated relative change in groundwater storage and availability under each 
scenario. Prior to simulating future change, the model was calibrated to 
historical groundwater and surface-water conditions. Part of the calibration 
process involves identifying sensitive and insensitive parameters. 
Calibration and sensitivity analysis experiments on model parameters can provide 
insight into reducing model error when predicting results such as simulated 
streamflow (figure 1). Daily streamflow data from National Water Information 
System (NWIS) site 11466800 [@USGS:2022] provided observations 
for this experiment. Twenty-one iterations of the SRP model were run in series 
with pyGSFLOW to test the sensitivity of the snarea_curve (snow depletion curve), 
ssr2gw_rate (gravity reservoir to groundwater reservior routing coeficient), 
and gsflow_coef (linear groundwater discharge equation coeficient) with a 
simple “for loop.” The  ssr2gw_rate was identified as a much more sensitive 
parameter to model calibration than the snarea_curve and gsflow_coef 
parameters (figure 1). Insights like these allow researchers to focus their 
calibration efforts on the most sensitive parameters and fix insensitive parameters, 
thus reducing the time and complexity of the calibration process. Although actual 
calibration generally is not done directly with pyGSFLOW, it provides an easy to use 
interface to update parameters based on grid cell location or parameter zone that 
can be used for complex operations in conjuction with external calibration software.

![Root mean squared error in streamflow predictions at U.S Geological Survey 
streamgage 11466800 [@USGS:2022] for three PRMS parameters 
(gsflow_coef, snarea_curve, and ssr2gw_rate) during calibration experiments on 
the Santa Rosa Plain Integrated Hydrologic Model, Santa Rosa, 
California.](calibration_example_mwc.png)


The pyGSFLOW package also includes features to visualize input and output data 
spatially using Matplotlib [@Hunter:2007] plots and by exporting datasets to 
shapefile or the visualization toolkit (VTK) format [@Schroeder:2006]. By 
providing pyGSFLOW a shapefile or list of hydrologic response unit (HRU) geometries, 
the code is able to plot arrays and contour arrays of unique parameter values 
and is fully compatible with the FloPy plotting routines for MODFLOW. PRMS 
input parameter values can be layered over MODFLOW output and can potentially 
help identify trends and sensitive parameters controlling those trends in 
streamflow, recharge, and groundwater levels throughout the model. The ssr2gw_rate 
parameter, which scales the exchange between the PRMS gravity reservoir and the 
MODFLOW groundwater reservior in GSFLOW, can then be overlain on top of recharge 
arrays to inspect the input and output for correlated trends (figure 2). Figure 
2 shows that in the western part of the basin, both the ssr2gw_rate and the 
relative amount of areal recharge is slightly greater than the eastern part of 
the basin. The simulated recharge data also show the highest volume of recharge 
occurs along a few short losing stream sections. These insights can help the 
researcher adjust input parameters for both streamflow and groundwater level 
calibration.

![Sagehen Creek GSFLOW model, Truckee, California. Mean recharge for the entire 
simulation from MODFLOW is overlain with a spatial contour plot of PRMS 
parameter ssr2gw_rate which is a multiplier that scales the volume of recharge 
from PRMS to MODFLOW. Black fill indicates inactive model cells.](sagehen_plot.png)


The online documentation for pyGSFLOW [@Larsen:2021] contains API information 
for all major classes and methods and is updated with each new major release. 
In addition to the online documentation, sample Jupyter notebooks [@Kluyver:2016] 
are included in the repository to help users become familiar with the code.

# Package architecture
The pyGSFLOW package includes the gsflow module and 5 sub-packages (figure 3):

   - gsflow: the gsflow module contains the integrated modeling object 
   GsflowModel, which allows the user to build new GSFLOW models and import 
   existing models. This module calls classes and methods from the following 
   5 packages within pyGSFLOW.
   - prms: the prms sub-package contains classes and methods to build new PRMS 
   models, import existing PRMS models, edit model input data, and write PRMS 
   input data to file.
   - modflow: the modflow package contains modules and classes that interface 
   with FloPy and allow the user to create new MODFLOW packages, edit existing 
   packages, and write MODFLOW input data to file.
   - modsim: the modsim sub-package contains classes that translate MODFLOW 
   model stream and lake networks into vectorized shapefile representations 
   that can be used to define surface-water operation networks in MODSIM.
   - output: the output sub-package contains modules that allow the user to 
   define their surface-water model discretization and visualize output data 
   via matplotlib plots.
   - utils: includes general use utilities that are integrated into built in 
   functions in the gsflow module, and prms, modflow, and modsim sub-packages.   

![Hierarchical representation of the pyGSFLOW package. Each sub-package lists 
the model building classes within each package. The GsflowModel class interacts 
with each of these listed modules and the FloPy package.](Package_architecture.png)


# Conclusion
GSFLOW integrated hydrologic models simulate complex processes and interactions 
between surface-water and groundwater flow systems. Parameterizing these model 
processes requires large datasets from multiple sources to represent the hydrologic 
cycle. Previous approaches involved multiple disconnected scripts and packages, 
some of which rely on proprietary code, which makes reproducibility difficult. 
pyGSFLOW is a tightly coupled software package, that allows the user to import 
all parts of their model in one script, which helps to standardize and streamline 
model development, calibration, and output analysis.

# Acknowledgments
The authors welcome additions, suggestions, and assistance from the scientific 
community, and thank all past contributors for their work. Any use of trade, 
product, or firm names is for descriptive purposes only and does not imply 
endorsement by the U.S. Government.

# References



