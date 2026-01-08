import os
import numpy as np
import gsflow

def get_prms_watershed(gsf):
    #
    nrow = gsf.mf.dis.nrow
    ncol = gsf.mf.dis.ncol
    hru_type = gsf.prms.parameters.get_values("hru_type")
    hru_type = hru_type.reshape(nrow, ncol)
    mask = np.logical_not(hru_type == 0)
    hru_type[mask] = 1
    return hru_type