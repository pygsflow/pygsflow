import numpy as np
from gsflow.builder import Defaults, PrmsDefaults
from gsflow import PrmsParameters, ParameterRecord


class PrmsBuilder(object):
    """
    Class for building PrmsParameter object using built in Defaults

    PrmsBuilder builds a PrmsParameter object using default values
    from gsflow.builder.Defaults or a user supplied Defaults object.

    The user can then edit the PrmsParameter object to customize their
    model runs and add climate, land cover, and soils information.

    Parameters
    ----------
    stream_data_obj : gsflow.builder._StreamsObj
        stream data from the flow accumulation stream creation method
    crt_param_obj : gsflow.prms.PrmsParameters object
        PrmsParameters object loaded from Crt() output
    modelgrid : flopy.discretization.StructuredGrid object
    hru_type : np.ndarray
        hru_type array
    defaults : gsflow.builder.Defaults object
        optional gsflow.builder.Defaults object to build the PrmsParameters
        object with.

    """

    def __init__(
        self,
        stream_data_obj,
        crt_param_obj,
        modelgrid,
        hru_type=None,
        defaults=None,
    ):
        if defaults is None:
            self._defaults = Defaults().prms.to_dict()
        elif isinstance(defaults, Defaults):
            self._defaults = defaults.prms.to_dict()
        elif isinstance(defaults, PrmsDefaults):
            self._defaults = defaults.to_dict()
        else:
            raise TypeError("Defaults must be Default or PrmsDefault object")

        self.stream_data_obj = stream_data_obj
        self.crt_param_obj = crt_param_obj
        self.modelgrid = modelgrid
        self.hru_type = hru_type

    def build(self, name=None):
        """
        Method to build a PrmsParameters object that can be used to build
        a GSFLOW or PRMS model

        Parameters
        ----------
        name : str
            name of the parameter file

        Returns
        -------
            gsflow.prms.parameters.PrmsParameters object

        """
        param_list = []

        # parameter file defaults
        paramfile_defaults = self._defaults["parameter"]

        # load dimension defaluts
        dimension_defaults = paramfile_defaults["dimensions"]

        # get crt parameters??
        crt_params = PrmsParameters(self.crt_param_obj)  # update this

        # get the nhru from the model grid
        nhru = self.modelgrid.nrow * self.modelgrid.ncol

        # set segment
        dimension_defaults[
            "nsegment"
        ] = (
            self.stream_data_obj.iseg.max()
        )  # I think this needs to be number of reaches!!!
        dimension_defaults["ngw"] = nhru
        dimension_defaults["nhru"] = nhru  # may change this
        dimension_defaults["nssr"] = nhru
        dimension_defaults[
            "ncascade"
        ] = nhru  # is this right??? I'm not so sure that this is. We should be able to get this from CRT...
        dimension_defaults["ncascadegw"] = crt_params.get_values("ncascdgw")[0]
        dimension_defaults["ndeplval"] = dimension_defaults["ndepl"] * 11
        dimension_defaults["nsub"] = np.count_nonzero(np.unique(self.hru_type))
        dimension_defaults["nhrucell"] = nhru
        dimension_defaults["ndeplval"] = paramfile_defaults["dimensions"][
            "ndeplval"
        ]

        for key, value in dimension_defaults.items():

            param_record = ParameterRecord(
                name=key, values=[value], datatype=1, file_name=name
            )
            param_list.append(param_record)

        param_defaults = paramfile_defaults["parameters"]

        for key, vals in param_defaults.items():
            dtype = vals["dtype"]
            dimension = vals["dimension"]
            record = vals["record"]

            if isinstance(dimension, list):
                if isinstance(
                    dimension[0], list
                ):  # check this, may need to be changed in .json
                    dimension = dimension[0]
            elif isinstance(dimension, str):
                dimension = [dimension]
            else:
                raise Exception("data type not supported")

            dim = [[nm, dimension_defaults[nm]] for nm in dimension]

            if dim[0][0] == "one":
                record = [record]
            else:
                record = record * dim[0][1]
            record = np.array(record).ravel()

            try:
                param_record = ParameterRecord(
                    name=key,
                    values=record,
                    dimensions=dim,
                    datatype=dtype,
                    file_name=name,
                )
            except ValueError:
                print(print(key, dim, record, type(dim), type(record)))
                raise ValueError()
            # #
            param_list.append(param_record)

        param_dict = {}
        cell_area = self.modelgrid.xcs * self.modelgrid.ycs
        hru_area = np.full(nhru, cell_area)
        param_dict["hru_area"] = {"record": hru_area, "dtype": 2}

        hru_lon = self.modelgrid.xcellcenters.ravel()
        hru_lat = self.modelgrid.ycellcenters.ravel()
        param_dict["hru_lat"] = {"record": hru_lat, "dtype": 2}
        param_dict["hru_lon"] = {"record": hru_lon, "dtype": 2}

        hru_slope = self.stream_data_obj.slope.ravel()
        param_dict["hru_slope"] = {"record": hru_slope, "dtype": 2}

        # hru type - get it from the FA??? modified ... - nhru
        # @ comment JL: I think we should have an option to return the
        #  FA hru type and then we can
        hru_type = (
            self.hru_type.ravel()
        )  # this may need to come from somewhere else or modify
        param_dict["hru_type"] = {"record": hru_type, "dtype": 1}

        # # hru_elev = dem elev (should it be the sink filled..???). -
        # passing for now
        # todo: look at returning a sink filled dem elevation in FlowAcc. for this
        aspect = self.stream_data_obj.aspect.ravel()
        param_dict["aspect"] = {"record": aspect, "dtype": 2}

        # hru_aspect - nhru
        dim = [["nhru", nhru]]
        for key, records in param_dict.items():
            record = records["record"]
            dtype = records["dtype"]

            param_record = ParameterRecord(
                name=key,
                values=record,
                dimensions=dim,
                datatype=dtype,
                file_name=name,
            )

            param_list.append(param_record)

        return PrmsParameters(param_list)
