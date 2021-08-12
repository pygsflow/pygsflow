from . import Defaults, ModflowDefaults
import gsflow
import numpy as np
import flopy


class ModflowBuilder(object):
    """
    Class for building Modflow model objects using built in Defaults.

    ModflowBuilder builds a steady state model using default values
    from gsflow.builder.Defaults or a user supplied Defaults object.

    The user can then edit the Modflow package objects to customize
    their model runs. (ex. create transient model), etc...

    supported packages include : DIS, BAS, UPW, UZF, SFR, NWT, OC

    Parameters
    ----------
    modelgrid : gsflow.builder.FishnetGenerator
        Structured grid object from FishnetGenerator or from
        flopy.discretization.StructuredGrid
    dem_data : np.ndarray
        numpy array of dimension (nrow, ncol) of DEM elevations
    model_name : str
        model name ex. "my_test_model"
    defaults : gsflow.builder.Defaults
        optional parameter, user can supply a gsflow.builder.Defaults
        instance to ModflowBuilder to use a custom set of default values

    """

    def __init__(self, modelgrid, dem_data, model_name, defaults=None):
        exe_name = "mfnwt.exe"
        self._ml = gsflow.modflow.Modflow(model_name, exe_name=exe_name)

        assert (modelgrid.nrow, modelgrid.ncol) == dem_data.shape
        self._modelgrid = modelgrid
        self._dem_data = dem_data

        if defaults is None:
            self._defaults = Defaults().modflow.to_dict()
        elif isinstance(defaults, Defaults):
            self._defaults = defaults.modflow.to_dict()
        elif isinstance(defaults, ModflowDefaults):
            self._defaults = defaults.to_dict()
        else:
            raise TypeError(
                "Defaults must be Default or ModflowDefault object"
            )

    @property
    def model(self):
        """
        Returns the gsflow.modflow.Modflow model object

        """
        return self._ml

    def build_all(
        self,
        reach_data,
        segment_data,
        irunbnd,
        finf=None,
        botm=None,
        ibound=None,
        iuzfbnd=None,
    ):
        """
        Method to build all supported modflow packages

        Parameters
        ----------
        reach_data : np.recarray
            flopy's ModflowSfr2 reach data parameter
        segment_data : np.recarray
            flopy's ModflowSfr2 segment data parameter
        irunbnd : np.ndarray
            flopy's ModflowUZF1 irunbnd parameter
            (runoff connection to streams)
        finf : np.ndarray
            UZF1's finf array which describes precipitation for recharge
        botm : np.ndarray
            bottom elevation for single layer model
        ibound : np.ndarray
            ibound array of active model cells
        iuzfbnd : np.ndarray
            uzf ibound array of active model cells

        Returns
        -------
            gsflow.modflow.Modflow object

        """
        self.build_dis(botm=botm)
        self.build_bas6(ibound=ibound)
        self.build_upw()
        self.build_nwt()
        self.build_oc()
        self.build_uzf(irunbnd, finf=finf, iuzfbnd=iuzfbnd)
        self.build_sfr(reach_data, segment_data)
        return self._ml

    def build_dis(self, botm=None):
        """
        Method to build the dis package using defaults

        Parameters
        ----------
        botm : float, int, np.ndarray
            Model botm elevations for discretization file. If botm is None
            then botm elevation is set 50 length units below DEM elevation

        Returns
        -------
            flopy.modflow.ModflowDis object

        """
        if botm is None:
            botm = self._dem_data - 50.0
        else:
            botm = botm

        dis_defaults = self._defaults["dis"]

        dis = flopy.modflow.ModflowDis(
            self._ml,
            nrow=self._modelgrid.nrow,
            ncol=self._modelgrid.ncol,
            delc=self._modelgrid.delc,
            delr=self._modelgrid.delr,
            top=self._dem_data,
            botm=botm,
            **dis_defaults
        )
        return dis

    def build_bas6(self, ibound=None):
        """
        Method to build the BAS6 package

        Parameters
        ----------
        ibound : int, np.ndarray
            array of active modflow cells within the model, >0 for active, 0
            for inactive

        Returns
        -------
            flopy.modflow.ModflowBas object

        """
        if ibound is None:
            ibound = np.ones(
                (self._modelgrid.nrow, self._modelgrid.ncol), dtype=int
            )
        bas_defaults = self._defaults["bas"]
        bas = flopy.modflow.ModflowBas(
            self._ml, ibound=ibound, strt=self._dem_data, **bas_defaults
        )
        return bas

    def build_upw(self):
        """
        Method to build a default version of the UPW package

        Returns
        -------
            flopy.modflow.ModflowUpw

        """
        upw_defaults = self._defaults["upw"]
        upw = flopy.modflow.ModflowUpw(self._ml, **upw_defaults)
        return upw

    def build_sfr(self, reach_data, segment_data):
        """
        Method to build a default version of the SFR package

        Parameters
        ----------
        reach_data : np.recarray
            reach data recarray for ModflowSfr2
        segment_data : np.recarray
            segment data recarray for ModflowSfr2

        Returns
        -------
            flopy.modflow.ModflowSfr2

        """
        # get package defaults and build the sfr package
        sfr_defaults = self._defaults["sfr"]["pkg"]
        nreaches = len(reach_data)
        nsegments = len(segment_data)

        sfr = flopy.modflow.ModflowSfr2(
            self._ml,
            nstrm=nreaches,
            nss=nsegments,
            reach_data=reach_data,
            segment_data=segment_data,
            **sfr_defaults
        )
        return sfr

    def build_uzf(self, irunbnd, finf=None, iuzfbnd=None):
        """
        Method to build a UZF package object using built in defaults

        Parameters
        ----------
        irunbnd : np.ndarray
            flopy's ModflowUZF1 irunbnd parameter
            (runoff connection to streams)
        finf : np.ndarray
            optional finf array of precipitation
        iuzfbnd : np.ndarray
            uzf ibound array of active unsaturated zone model cells

        Returns
        -------
            flopy.modflow.ModflowUzf

        """
        uzf_defaults = self._defaults["uzf"]
        if finf is None:
            finf = 1e-08
        if iuzfbnd is None:
            iuzfbnd = 1

        uzf = flopy.modflow.ModflowUzf1(
            self._ml,
            irunbnd=irunbnd,
            iuzfbnd=iuzfbnd,
            finf=finf,
            **uzf_defaults
        )
        return uzf

    def build_nwt(self):
        """
        Method to build a SIMPLE nwt solver instance

        Returns
        -------
            flopy.modflow.ModflowNwt

        """
        nwt_defaults = self._defaults["nwt"]
        nwt = flopy.modflow.ModflowNwt(self._ml, **nwt_defaults)
        return nwt

    def build_oc(self):
        """
        Method to build a simple one stress period OC object

        Returns
        -------
            flopy.modflow.ModflowOc
        """
        # build output control
        oc_defaults = self._defaults["oc"]
        oc = flopy.modflow.ModflowOc(self._ml, **oc_defaults)
        return oc
