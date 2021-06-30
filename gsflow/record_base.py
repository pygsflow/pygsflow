import numpy as np
from .utils import GsConstant


class RecordBase(object):
    """
    Base class for record objects. Class stores common base code
    that is used across other Record objects

    Parameters
    ----------
    name : str
        record name
    values : list, np.ndarray
        record values
    datatype : int
        record datatype

    """

    def __init__(self, name, values=None, datatype=None):

        if isinstance(name, str):
            self.name = name
        else:
            raise ValueError(
                "Error: Name of the record ( {} ) is not a string".format(name)
            )

        self._values = None
        self._check_values(values)

        self.datatype = datatype
        if datatype is None:
            self._infer_dtype()

        self._force_dtype()

    def _infer_dtype(self):
        """
        Method to infer the dtype from data available

        """
        print(
            "Warning: {} data type will be infered from data supplied".format(
                self.name
            )
        )
        if "float" in self._values.dtype.name:
            self.datatype = 2
        elif "int" in self._values.dtype.name:
            self.datatype = 1
        elif "str" in self._values.dtype.name:
            self.datatype = 4
        else:
            raise ValueError(
                "Value type is not recognized...{}", self.values.dtype
            )

    def _check_values(self, new_values):
        """
        Method to check if values can be numpy array

        Parameters
        ----------
        new_values : list

        """
        if isinstance(new_values, np.ndarray):
            self._values = new_values
        elif isinstance(new_values, list):
            try:
                new_values = np.array(new_values)
                self._values = new_values
            except:
                raise ValueError("Cannot convert the list to 1D numpy array ")
        else:
            raise ValueError("Values must be list or 1D numpy array ")

    @property
    def values(self):
        return self._values

    @property
    def nvalues(self):
        if isinstance(self._values, np.ndarray):
            return self._values.size
        elif isinstance(self._values, list):
            return len(self._values)
        else:
            return 0

    def _check_dtype(self):
        """
        Method that checks if a valid dtype has been suplied

        """
        if "float" in self._values.dtype.name:
            self.datatype = 2
        elif "int" in self._values.dtype.name:
            self.datatype = 1
        elif "str" in self._values.dtype.name:
            self.datatype = 4
        elif "unicode" in self._values.dtype.name:
            self.datatype = 4
        else:
            raise ValueError(
                "Value type is not recognized...{}".format(self._values.dtype)
            )

    def _force_dtype(self):
        """
        Method to force the dtype on an array

        """
        dtype = GsConstant.PRMS_DATA_TYPES[self.datatype]
        if dtype == "integer":
            self._values = self._values.astype(int)
        elif dtype in ("real", "double"):
            self._values = self._values.astype(float)
        elif dtype == "string":
            self._values = self._values.astype(str)
        else:
            raise ValueError("Error : Cannot recognize data type")

    def export_nc(self, f, modflow, **kwargs):
        pass

    def from_dict(self):
        pass

    def to_dict(self):
        pass
