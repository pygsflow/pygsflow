from .. import ControlFile
from . import Defaults, ControlFileDefaults


class ControlFileBuilder(object):
    """
    Class for building Modflow model objects using built in Defaults.

    ControlFile Builder builds a GSFLOW control file using default values
    from gsflow.builder.Defaults or a user supplied Defaults object.

    The user can then edit the ControlFile object to customize their model
    runs.

    Parameters
    ----------
    defaults : gsflow.builder.Defaults
        optional parameter, user can supply a gsflow.builder.Defaults instance
        to the ControlFileBuilder to use a custom set of parameters

    """

    def __init__(self, defaults=None):
        if defaults is None:
            self._defaults = Defaults().control.to_dict()
        elif isinstance(defaults, Defaults):
            self._defaults = defaults.control.to_dict()
        elif isinstance(defaults, ControlFileDefaults):
            self._defaults = defaults.to_dict()
        else:
            raise TypeError(
                "Defaults must be Defaults or ControlFileDefaults object"
            )

    def build(
        self, name="gsflow_builder", parameter_obj=None, modflow_obj=None
    ):
        """
        Method to build a pyGSFLOW ControlFile object

        Parameters
        ----------
        name : str, optional
            name of the control file/model
        parameter_obj : gsflow.prms.parameters.PrmsParameters
            optional parameter object, if supplied this adds parameter
            file names to the gsflow control file
        modflow_obj : gsflow.modflow.Modflow object
            optional modflow object, if supplied this adds the modflow
            nam file to the control file.

        Returns
        -------
            gsflow.control.ControlFile
        """
        from ..prms import PrmsParameters
        from ..modflow import Modflow

        control = ControlFile(records_list=[], name=name)

        defaults = self._defaults["control"]
        for key, value in defaults.items():
            record = value["record"]
            if isinstance(record, (float, str, int)):
                record = [record]
            control.add_record(key, record)

        if isinstance(parameter_obj, PrmsParameters):
            parameter_files = parameter_obj.parameter_files
            record = []
            for f in parameter_files:
                if f is None:
                    record.append(f"{name}.param")
                else:
                    record.append(f)
            control.add_record("param_file", record)

        if isinstance(modflow_obj, Modflow):
            record = [
                modflow_obj.namefile,
            ]
            control.add_record("modflow_name", record)

        return control
