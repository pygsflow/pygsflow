---
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
    affiliation: U.S. Geological Survey, California Water Science Center, Sacramento, CA
  - name: Ayman Alzraiee
    ocrid: 0000-0001-7576-3449
    affiliation: U.S. Geological Survey, California Water Science Center, Sacramento, CA
  - name: Richard G. Niswonger
    ocrid: 0000-0001-6397-2403
    affiliation: U.S. Geological Survey, Water Mission Area, Menlo Park, CA
date: 3/19/2021
bibliography: pygsflow_JOSS.bib
---

# Overview
Paste overview text here

# Introduction
Paste introduction text here

# Statement of need
Paste statement of need text here

# pyGSFLOW
Paste body text here

![Mean squared error in streamflow preditions for three PRMS parameters 
(gsflow_coef, snarea_curve, and ssr2gw_rate) during calibration experiments 
on the Santa Rosa Plain Integrated Hydrologic Model, Santa Rosa, 
California.\label{fig:calibration}](calibration_example.png)

\autoref{fig:calibration}

more body text

![Mean recharge for the entire simulation from MODFLOW is overlain with a 
spatial contour plot of PRMS ssr2gw_rate which is a multiplier that scales 
the volume of recharge from PRMS to MODFLOW. MODFLOW’s IBOUND array 
(black is inactive cells) is also plotted to distinguish active versus 
inactive model cells, Sagehen Creek GSFLOW model, 
Truckee, California.\label{fig:plotting}](sagehen_plot.png)

\autoref{fig:plotting}

more body text

# Package architecture
Paste package text here

![Hierarchical representation of the pyGSFLOW package. Each sub-package lists 
the model building classes within each package. The GSFlowModel class interacts 
with each of these listed modules and the FloPy 
package.\label{fig:arch}](Package_architecture.png)

\autoref{fig:arch}

# Conclusion
Paste conclusion text here

# References