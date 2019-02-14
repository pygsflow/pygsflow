import os
import sys
import numpy as np
from .utils import io
from .utils import GsConstant
import copy
import warnings
warnings.simplefilter('always', PendingDeprecationWarning)
warnings.simplefilter('always', UserWarning)


class ParameterBase(object):
    """
    Base class for parameter based files within GSFLOW

    Parameters
    ----------
    records_list
    name
    model_dir
    header
    """
    def __init__(self, records_list, name=None, model_dir=None, header=None):
        self._record_names = []

        if name is not None:
            self.name = name

        if model_dir is not None:
            self.model_dir = model_dir

        if header is not None:
            self.header = header

        self._record_names = []
        self._records_list = []
        if records_list is not None:
            self._records_list = copy.deepcopy(records_list)

    @property
    def record_names(self):
        rn = []
        if self._records_list is not None:
            if len(self._records_list) > 0:
                for rec in self._records_list:
                    rn.append(rec.name)
        return rn

    def get_record(self, name, rectype):
        """

        Parameters
        ----------
        name : str
        rectype : class object
            ParameterRecord or ControlRecord

        Returns
        -------

        """
        record = io.find_parameter(name, self._records_list)

        if isinstance(record, rectype):
            return record
        else:
            err = "The record does not exist..."
            warnings.warn(err, UserWarning)
            return None

    def get_values(self, name):
        """

        Parameters
        ----------
        name

        Returns
        -------

        """
        record = io.find_parameter(name, self._records_list)

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

        Parameters
        ----------
        name
        values

        Returns
        -------

        """

        record = io.find_parameter(name, self._records_list)

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

        Parameters
        ----------
        name
        values

        Returns
        -------

        """
        if isinstance(name, str):
            pass
        else:
            raise ValueError("Record name must be string")
        if isinstance(values, list) or isinstance(values, np.ndarray):
            pass
        else:
            raise ValueError("Value must be a list or np.ndarray")

        if name in self._record_names:
            err = "The record already exists, skipping add_record: {}...".format(name)
            warnings.warn(err, UserWarning)
            return False

        return True

    def add_record(self, recobj, where=None, after=None):
        """
        Generalized method to add a record to a record list

        Parameters
        ----------
        recobj :
        where :
        after :

        """
        if after:
            index = -1
            for index, rec in enumerate(self._record_names):
                if rec == after:
                    break
            self._record_names.insert(index + 1, recobj.name)
            self._records_list.insert(index + 1, recobj)
            return

        if where:
            self._record_names.insert(where, recobj.name)
            self._records_list.insert(where, recobj)
        else:
            self._record_names.append(recobj.name)
            self._records_list.append(recobj)

    def remove_record(self, name):
        """

        Parameters
        ----------
        name : str

        Returns
        -------

        """
        if isinstance(name, str):
            pass
        else:
            raise ValueError("Record name must be a string")

        if name not in self._record_names:
            warnings.warn("The record does not exist: {}".format(name),
                          UserWarning)
            return

        for index, nm in enumerate(self._record_names):
            if nm == name:
                self._record_names.pop(index)
                self._records_list.pop(index)

    def get_record_names(self):
        """

        Returns
        -------
            list of record names

        """
        err = "get_record_names() is Deprecated," \
              " calling record_names"
        warnings.warn(err, PendingDeprecationWarning)
        return self.record_names
