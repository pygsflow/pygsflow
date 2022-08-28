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

    """
    def __init__(self, defaults=None):
        if defaults is None:
            self._defaults = PrmsDefaults('prms_defaults.json').to_dict()
            self._key = [*self._defaults][0]
        #elif isinstance(defaults, Defaults):
        #    self._defaults = defaults.control.to_dict()
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
            # starting with dimensions
            dim = self._defaults[self._key]['dimensions']
            init_dim = [ParameterRecord(name=i, values=[dim[i]]) for i in [*dim]]
            params = PrmsParameters(init_dim)
            # continuing with parameters
            
            par = self._defaults[self._key]['parameters']
            for p in [*par]:
                datatype = par[p]['dtype']
                dimension = par[p]['dimension']
                values = par[p]['record']
                if len(dimension[0]) == 1:
                    dimension1 = [[par[p]['dimension'], dim[par[p]['dimension']]]]
                    values1 = np.ones(dim[par[p]['dimension']])*values
                    rec = ParameterRecord(name=p,
                                                 dimensions=dimension1,
                                                 datatype=datatype,
                                                 values=values1,width='')
                elif len(dimension[0]) == 2:
                    dimension2 = [[i,dim[i]] for i in dimension[0]]
                    if len(values) != dimension2[1][1]:
                        raise TypeError(F"Too many 'dimensions' in the parameter '{p}', maximum is 2")
                    else:
                        values2 = np.array([i*dimension2[0][1] for i in values]).flatten()
                        rec = ParameterRecord(name=p,
                                                     dimensions=dimension2,
                                                     datatype=datatype,
                                                     values=values2,width='')
                params.add_record_object(rec)
        return params