import numpy as np
from . import PrmsDefaults
from ..prms import ParameterRecord,PrmsParameters

class ParameterFileBuilder(object):
    """
    Class to load PRMS parameters from JSON files.

    Parameters
    ----------
    defaults : JSON file
        JSON file with all dimensions and parameters

    Examples
    --------

    Load parameters from file

    >>> from gsflow.builder import ParameterFileBuilder
    >>> pb = ParameterFileBuilder('path_to_json')
    >>> params = pb.build()
    
    Formatting json file examples
    --------
    {
        "parameter": {
            "dimensions": {
                "nhru": 10,
            },
            "parameters": {
                "elev_units": {
                    "dimension": "one",
                    "dtype": 1,
                    "record": 1
                },
                "hru_area": {
                    "dimension": "nhru",
                    "dtype": 2,
                    "record": [10]
                },
                "snarea_curve": {
                    "dimension": "ndeplval",
                    "dtype": 2,
                    "record": [
                        0.05,0.24,0.4,0.53,0.65,0.74,0.82,0.88,0.93,0.97,1.0
                    ]
                },
                "adjust_rain": {
                    "dimension": [["nrain","nmonths"]],
                    "dtype": 2,
                    "record": [
                        [-0.4]
                    ]
                }
            }
        }
    }

    """
    def __init__(self, defaults=None):
        if defaults is None:
            raise TypeError(F"Missing input default file")
        elif isinstance(defaults, PrmsDefaults):
            self._defaults = defaults.to_dict()
            if len([*self._defaults])>1:
                raise TypeError("Must only have one main key {'key':{'parameters': ...,'dimensions':...}}")
            else:
                self._key = [*self._defaults][0]
        else:
            raise TypeError(
                "Must be a PrmsDefaults object"
            )
    def build(self):
        """
        Method to build a parameters object

        Returns
        -------
            gsflow.prms.prms_parameter.PrmsParameters
        """
        f_keys = self._defaults[self._key]        
        if 'dimensions' not in f_keys:
            raise TypeError("Missing 'dimensions' key in json file")
        elif 'parameters' not in f_keys:
            raise TypeError("Missing 'parameters' key in json file")
        else:
            ## starting with dimensions
            dim = self._defaults[self._key]['dimensions']
            init_dim = [ParameterRecord(name=i, values=[dim[i]]) for i in [*dim]]
            params = PrmsParameters(init_dim)
            
            ## continuing with parameters
            par = self._defaults[self._key]['parameters']
            for p in [*par]:
                print(p)
                datatype = par[p]['dtype']
                dimension = par[p]['dimension']
                values = par[p]['record']
                if dimension == 'one':
                    dimension1 = [[par[p]['dimension'], dim[par[p]['dimension']]]]
                    values1=[values]
                    rec = ParameterRecord(name=p,
                                          dimensions=dimension1,
                                          datatype=datatype,
                                          values=values1,width='')
                elif (isinstance(dimension, str)) or ((len(dimension)==1) and (isinstance(dimension[0][1], int))):
                    l1 = np.array(values).shape[0]
                    if (isinstance(dimension,str)): #For defaults
                        dimension1 = [[par[p]['dimension'], dim[par[p]['dimension']]]]
                        if l1==1: # Defaults not modified or with single value modification
                            values1 = np.array(values*dimension1[0][1]).flatten()
                        elif l1 == dimension1[0][1]: # If it has the shape ready
                            values1 = np.array(values).flatten()
                        else:
                            print(p,"MODIFY!!! 1")
                    elif (len(dimension)==1) and (isinstance(dimension[0][1], int)): # For building in notebook
                        dimension1 = dimension
                        if l1 ==1:
                            values1 = np.array(values*dimension1[0][1]).flatten()
                        elif l1 == dimension1[0][1]:
                            values1 = np.array(values).flatten()
                        else:
                            raise TypeError(F"Wrong values length")                  
                    rec = ParameterRecord(name=p,
                                          dimensions=dimension1,
                                          datatype=datatype,
                                          values=values1,width='')
                    values1 = None
                elif len(dimension[0]) == 2:
                    try:
                        dimension2 = [[i,dim[i]] for i in dimension[0]]
                    except:
                        dimension2 = [[i[0],dim[i[0]]] for i in dimension]
                    
                    #if given the list of data
                    if (len(values)  == dimension2[0][1]*dimension2[1][1]) and (dimension2[0][1] != 1):
                        values2 = np.array(values).flatten()
                    # If defaults used
                    else:
                        l1,l2 = np.array(values).shape # lenghts of parameters
                        if (dimension2[0][1] == l1) and (dimension2[1][1] == l2):# if given the shape
                            values2=np.array(values).flatten()
                        elif (dimension2[0][1] == l1) and (1==l2):
                            values2 = np.array([i*dimension2[1][1] for i in values]).flatten() # To fill
                        elif (1==l1) and (1==l2):
                            values2 = np.array(values[0]*dimension2[0][1]*dimension2[1][1])
                        else:
                            raise TypeError(F"Error formatting")

                    rec = ParameterRecord(name=p,
                                          dimensions=dimension2,
                                          datatype=datatype,
                                          values=values2,width='')
                params.add_record_object(rec)
        return params