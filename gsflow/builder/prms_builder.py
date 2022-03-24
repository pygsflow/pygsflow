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
        stream data from the flow accumulation make_streams() method
    cascade_obj : gsflow.builder._Cascades object
        cascade data object from the flow accumulation get_cascades() method
    modelgrid : flopy.discretization.StructuredGrid object
    hru_type : np.ndarray
        hru_type array
    hru_subbasin : np.ndarray
        hru_subbasin array
    defaults : gsflow.builder.Defaults object
        optional gsflow.builder.Defaults object to build the PrmsParameters
        object with.

    """

    def __init__(
        self,
        stream_data_obj,
        cascades_obj,
        modelgrid,
        dem,
        hru_type=None,
        hru_subbasin=None,
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
        self.cascades_obj = cascades_obj
        self.modelgrid = modelgrid
        self.dem = dem
        self.hru_type = hru_type
        self.hru_subbasin = hru_subbasin

    def build(self, name=None, area_conv=2.47105e-4):
        """
        Method to build a PrmsParameters object that can be used to build
        a GSFLOW or PRMS model

        Parameters
        ----------
        name : str
            name of the parameter file
        area_conv : float
            conversion factor from model units sq. to arces

        Returns
        -------
            gsflow.prms.parameters.PrmsParameters object

        """
        param_list = []

        # parameter file defaults
        paramfile_defaults = self._defaults["parameter"]

        # load dimension defaluts
        dimension_defaults = paramfile_defaults["dimensions"]

        # get the nhru from the model grid
        nhru = self.modelgrid.nrow * self.modelgrid.ncol

        # set segment and reach
        dimension_defaults["nsegment"] = self.stream_data_obj.iseg.max()
        dimension_defaults["nreach"] = self.stream_data_obj.reach_data.size

        dimension_defaults["ngw"] = nhru
        dimension_defaults["ngwcell"] = nhru
        dimension_defaults["nhru"] = nhru
        dimension_defaults["nhrucell"] = nhru
        dimension_defaults["nssr"] = nhru
        dimension_defaults["ncascade"] = self.cascades_obj.ncascade
        dimension_defaults["ncascadgw"] = self.cascades_obj.ncascade
        dimension_defaults["ndeplval"] = dimension_defaults["ndepl"] * 11
        dimension_defaults["nsub"] = np.count_nonzero(np.unique(self.hru_type))
        dimension_defaults["nhrucell"] = nhru
        dimension_defaults["ndeplval"] = paramfile_defaults["dimensions"][
            "ndeplval"
        ]
        try:
            dimension_defaults["ngwcell"] = self.modelgrid.nnodes
        except TypeError:
            dimension_defaults["ngwcell"] = self.modelgrid.ncpl

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
        cell_area = (self.modelgrid.xcs * self.modelgrid.ycs) * area_conv
        hru_area = np.full(nhru, cell_area)
        param_dict["hru_area"] = {"record": hru_area, "dtype": 2}

        hru_lon = self.modelgrid.xcellcenters.ravel()
        hru_lat = self.modelgrid.ycellcenters.ravel()
        param_dict["hru_lat"] = {"record": hru_lat, "dtype": 2}
        param_dict["hru_lon"] = {"record": hru_lon, "dtype": 2}

        hru_slope = self.stream_data_obj.slope.ravel()
        hru_slope[hru_slope < 1e-04] = 0
        param_dict["hru_slope"] = {"record": hru_slope, "dtype": 2}

        if self.hru_type is None:
            hru_type = np.ones((nhru,), dtype=int)
        else:
            hru_type = self.hru_type.ravel()
        param_dict["hru_type"] = {"record": hru_type, "dtype": 1}

        # # hru_elev = dem elev (should it be the sink filled..???). -
        param_dict["hru_elev"] = {"record": self.dem.ravel(), "dtype": 2}

        # todo: look at returning a sink filled dem elevation in FlowAcc. for this
        aspect = self.stream_data_obj.aspect.ravel()
        param_dict["hru_aspect"] = {"record": aspect, "dtype": 2}

        if self.hru_subbasin is None:
            hru_subbasin = np.ones((nhru,), dtype=int)
        else:
            hru_subbasin = self.hru_subbasin.ravel()
        param_dict["hru_subbasin"] = {"record": hru_subbasin, "dtype": 1}

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

        params = PrmsParameters(param_list)

        # add cascade parameters
        for rname in ("hru_up_id", "hru_down_id", "hru_strmseg_down_id"):
            values = self.cascades_obj.__dict__[rname]
            params.add_record(
                name=rname,
                values=values,
                dimensions=[["ncascade", self.cascades_obj.ncascade]],
                datatype=1,
                file_name=name,
                replace=True,
            )

            params.add_record(
                name=rname.replace("hru", "gw"),
                values=values,
                dimensions=[["ncascadgw", self.cascades_obj.ncascade]],
                datatype=1,
                file_name=name,
                replace=True,
            )

        params.add_record(
            name="hru_pct_up",
            values=self.cascades_obj.hru_pct_up,
            dimensions=[["ncascade", self.cascades_obj.ncascade]],
            datatype=2,
            file_name=name,
            replace=True,
        )

        params.add_record(
            name="gw_pct_up",
            values=self.cascades_obj.hru_pct_up,
            dimensions=[["ncascadgw", self.cascades_obj.ncascade]],
            datatype=2,
            file_name=name,
            replace=True,
        )

        if len(np.unique(self.hru_subbasin)) == 2:
            params.add_record(
                "subbasin_down",
                values=[
                    0,
                ],
                dimensions=[["nsub", 1]],
                datatype=1,
                file_name=name,
            )

        # add gvr linkage
        for t in ("cell", "hru"):
            params.add_record(
                f"gvr_{t}_id",
                values=list(range(1, nhru + 1)),
                dimensions=[["nhrucell", nhru]],
                datatype=1,
                file_name=name,
            )

        return params
