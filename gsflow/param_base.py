import numpy as np
from .utils import gsflow_io
import copy
import warnings

warnings.simplefilter("always", PendingDeprecationWarning)
warnings.simplefilter("always", UserWarning)


class ParameterBase(object):
    """
    Base class for parameter based files within GSFLOW

    Parameters
    ----------
    records_list : list
        list of ParameterRecord objects
    name : str, optional
        parameter file name
    model_dir : str, optional
        parameter file directory
    header : str, optional
        header for the parameter file

    """

    def __init__(self, records_list, name=None, model_dir=None, header=None):

        if name is not None:
            self.name = name

        if model_dir is not None:
            self.model_dir = model_dir

        if header is not None:
            self.header = header

        self._records_list = []
        if records_list is not None:
            self._records_list = copy.deepcopy(records_list)

    @property
    def record_names(self):
        rn = []
        if self._records_list is not None:
            if len(self._records_list) > 0:
                for rec in self._records_list:
                    rn.append(rec.name.lower())
        return rn

    def get_record(self, name, rectype):
        """
        Method to get records

        Parameters
        ----------
        name : str
            name of record
        rectype : class object
            ParameterRecord or ControlRecord

        Returns
        -------
            RecordBase object

        """
        record = gsflow_io.find_parameter(name, self._records_list)

        if isinstance(record, rectype):
            return record
        else:
            err = "The record does not exist..."
            warnings.warn(err, UserWarning)
            return None

    def get_values(self, name):
        """
        Method to get record values

        Parameters
        ----------
        name : str
            name of record

        Returns
        -------
            np.ndarray

        """
        record = gsflow_io.find_parameter(name, self._records_list)

        if record is None:
            err = "The record does not exist..."
            warnings.warn(err, UserWarning)
            return None

        elif isinstance(record.values, np.ndarray):
            return record.values

        else:
            err = "The values does not exist..."
            warnings.warn(err, UserWarning)
            return None

    def set_values(self, name, values):
        """
        Method to set new values to a record

        Parameters
        ----------
        name : str
            record name
        values : list
            list of values

        """

        record = gsflow_io.find_parameter(name, self._records_list)

        if record is None:
            err = "The record does not exist {}".format(name)
            warnings.warn(err, UserWarning)
            return None

        elif isinstance(record.values, np.ndarray):
            record.values = values
            return record.values

        else:
            err = "The record does not exist..."
            warnings.warn(err, UserWarning)
            return None

    def _check_before_add(self, name, values):
        """
        Internal method to check a record before adding to _records_list

        Parameters
        ----------
        name : str
            record name
        values : list
            list of values

        """
        if isinstance(name, str):
            name = name.lower()
        else:
            raise ValueError("Record name must be string")
        if isinstance(values, list) or isinstance(values, np.ndarray):
            pass
        else:
            raise ValueError("Value must be a list or np.ndarray")

        if name in self.record_names:
            err = (
                "The record already exists, skipping add_record: {}...".format(
                    name
                )
            )
            warnings.warn(err, UserWarning)
            return False

        return True

    def add_record(self, recobj, where=None, after=None):
        """
        Generalized method to add a record to a record list

        Parameters
        ----------
        recobj : RecordBase object
            ParameterRecord or ControlRecord
        where : int
            index location to insert record
        after : int
            index location - 1 to insert record

        """
        if after:
            index = -1
            for index, rec in enumerate(self.record_names):
                if rec == after:
                    break
            self._records_list.insert(index + 1, recobj)
            return

        elif where:
            self._records_list.insert(where, recobj)
        else:
            self._records_list.append(recobj)

    def remove_record(self, name):
        """
        Method to remove a record

        Parameters
        ----------
        name : str
            parameter name

        """
        if isinstance(name, str):
            name = name.lower()
        else:
            raise ValueError("Record name must be a string")

        if name not in self.record_names:
            warnings.warn(
                "The record does not exist: {}".format(name), UserWarning
            )
            return

        for index, nm in enumerate(self.record_names):
            if nm == name:
                self.record_names.pop(index)
                self._records_list.pop(index)
