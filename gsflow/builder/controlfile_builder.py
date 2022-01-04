from .. import ControlFile, ControlRecord
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

    def build(self, name="gsflow_builder"):
        """
        Method to build a pyGSFLOW ControlFile object

        Parameters
        ----------
        name : str, optional
            name of the control file/model

        Returns
        -------
            gsflow.control.ControlFile
        """
        control = ControlFile(
            records_list=[],
            name=name
        )

        defaults = self._defaults['control']
        for key, value in defaults.items():
            record = value['record']
            if isinstance(record, (float, str, int)):
                record = [record]
            control.add_record(key, record)

        return control