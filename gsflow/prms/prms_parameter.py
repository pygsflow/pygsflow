from __future__ import absolute_import, division, print_function
import os
import numpy as np
from ..record_base import RecordBase
from ..param_base import ParameterBase


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


class PrmsParameters(ParameterBase):
    """
    Class to hold PRMS parameter information.
    This class enables reading, writing, and editing PRMS parameter files.

    Parameters
    ----------
    parameters_list : list
        list of ParameterRecord objects
    headers : str, optional
        file header

    Examples
    --------

    Load parameters from file

    >>> import gsflow
    >>> params = gsflow.prms.PrmsParameters.load_from_file(["myparams1.txt", "myparams2.txt"])

    create a new PrmsParameters object

    >>> params = gsflow.prms.PrmsParameters([parameter_record1, parameter_record2,])

    """

    def __init__(self, parameters_list, header=None):
        super(PrmsParameters, self).__init__(parameters_list, header=header)
        # default header

        if header is None:
            self.header = [
                "Parameter files generated by gsflow python package"
            ]

        self.__parameter_files = []

    def __getattr__(self, item):

        if item in self.record_names:
            return self.get_record(item)
        else:
            try:
                return super(PrmsParameters).__getattribute__(item)
            except AttributeError:
                raise AttributeError(
                    f"PrmsParameters does not have record {item}"
                )

    def __setattr__(self, key, value):
        if key in (
            "name",
            "model_dir",
            "header",
            "_records_list",
            "_PrmsParameters__parameter_files",
        ):
            super().__setattr__(key, value)
        elif key in self.record_names:
            self.set_values(key, value)
        else:
            raise AttributeError(f"PrmsParameters does not have record {key}")

    @property
    def parameters_list(self):
        """
        Returns a list of parameter records

        """
        return self._records_list

    @property
    def parameter_files(self):
        """
        Returns a list of parameter file names

        """
        all_files = []

        for rec in self.parameters_list:
            if rec.section == "Dimensions":
                if rec.file_name in all_files:
                    continue
                else:
                    all_files.insert(0, rec.file_name)
            else:
                if not (rec.file_name in all_files):
                    all_files.append(rec.file_name)

        self.__parameter_files = all_files
        return self.__parameter_files

    def reset_filenames(self, f):
        """
        Method to reset all filenames in the PrmsParameter object
        to a single user supplied filename

        Parameters
        ----------
        f : str

        Returns
        -------

        """

    @staticmethod
    def load_from_file(param_files):
        """
        Create a PrmsParameters object from parameter file(s)

        Parameters
        ----------
        param_files : list
            list of parameter files

        Returns
        -------
            PrmsParameters

        """
        if isinstance(param_files, str):
            param_files = [param_files]

        elif isinstance(param_files, list):
            pass

        # loop over parameter files and load them successively
        all_dims = {}
        headers = []
        parameters_list = []
        for ifile, file in enumerate(param_files):
            print("------------------------------------")
            print(
                "Reading parameter file : {}".format(os.path.split(file)[-1])
            )
            print("------------------------------------")
            if not (os.path.isfile(file)):
                raise FileNotFoundError("Invalid file name {}".format(file))

            with open(file, "r") as fid:

                EndOfFile = False
                _read_comments = True
                in_dim_section = False

                while True:
                    record = fid.readline()
                    if record == "\n":
                        continue
                    elif record == "":
                        break

                    else:
                        # read comments
                        record = record.strip()
                        if _read_comments:
                            if (
                                "####" in record
                                or ("** Dimensions **" in record)
                                or ("** Parameters **" in record)
                            ):
                                _read_comments = False
                            else:
                                headers.append([record, file])
                                continue

                        # read records information Not comments

                        if "** Parameters **" in record:
                            in_dim_section = False

                            continue

                        elif "** Dimensions **" in record:
                            in_dim_section = True
                            continue

                        if "####" in record:
                            continue

                        if in_dim_section:
                            # Reading Dimensions Section
                            field_name = record.strip()
                            if is_number(field_name):
                                raise ValueError(
                                    "The parameter name is a number. Check the dimensions of {} "
                                    "".format(parameters_list[-1].name)
                                )
                            value = int(fid.readline().strip())
                            curr_record = ParameterRecord(
                                name=field_name, values=[value], file_name=file
                            )
                            parameters_list.append(curr_record)
                            all_dims[field_name] = value
                        else:
                            # read Parameters section
                            field_name = record.strip().split()[0]
                            if is_number(field_name):
                                raise ValueError(
                                    "The parameter name is a number. Check file {}. The dimensions "
                                    "of {} might be wrong"
                                    "".format(file, parameters_list[-1].name)
                                )
                            ndim = int(fid.readline().strip())
                            dim_nms = []
                            for dim_ in range(ndim):
                                dim_nms.append(fid.readline().strip())

                            nvalues = int(fid.readline().strip())
                            datatype = int(fid.readline().strip())
                            value = []

                            # read values sequentially
                            val_count = 0
                            while val_count < nvalues:
                                val_ = fid.readline().strip()
                                if "*" in val_:
                                    comp_val = val_.split("*")
                                    value = value + [comp_val[1]] * int(
                                        comp_val[0]
                                    )
                                    val_count = val_count + int(comp_val[0])

                                else:
                                    for va_ in val_.split():
                                        if datatype != 4:
                                            if not (is_number(va_)):
                                                break
                                        value.append(va_)
                                        val_count = val_count + 1

                            par_dim = []
                            for dn in dim_nms:
                                par_dim.append([dn, all_dims[dn]])
                            curr_record = ParameterRecord(
                                name=field_name,
                                values=value,
                                dimensions=par_dim,
                                datatype=datatype,
                                file_name=file,
                            )

                            parameters_list.append(curr_record)

        return PrmsParameters(parameters_list=parameters_list, header=headers)

    def export_nc(self, f, modflow, **kwargs):
        """
        Method to export a PrmsParameters object
        to netcdf (.nc) file

        Parameters
        ----------
        f : str or fp.export.NetCdf
            filename to write the parameter to (*.nc)
        modflow : object
            fp.modflow.Modflow or gsflow.modflow.Modflow object

        Notes
        -----
        NetCdf export relies on flopy, so at the moment will
        only work for GSFLOW models where PRMS has the same
        discretization as the modflow grid

        """
        for parameter in self.parameters_list:
            f = parameter.export_nc(f, modflow, **kwargs)

        return f

    def get_record(self, name):
        """
        Method to get a parameter record by name

        Parameters
        ----------
        name : str
            parameter name

        Returns
        -------
            ParameterRecord object

        """
        return super(PrmsParameters, self).get_record(name, ParameterRecord)

    def get_values(self, name):
        """
        Method to get values from a parameter

        Parameters
        ----------
        name : str
            parameter name

        Returns
        -------
            np.array

        """
        return super(PrmsParameters, self).get_values(name)

    def set_values(self, name, values):
        """
        Method to set values to an existing parameter

        Parameters
        ----------
        name : str
            parameter name
        values : list
            list of values

        """
        return super(PrmsParameters, self).set_values(name, values)

    def add_record(
        self,
        name,
        values=None,
        dimensions=None,
        datatype=None,
        width=10,
        file_name=None,
        where=None,
        after=None,
        replace=False,
    ):
        """
        Method to add a new parameter record to the PrmsParameters object

        Parameters
        ----------
        name : str
            parameter name
        values : list
            parameter values
        dimensions : list
            parameter dimensions
        datatype : int
            datatype flag
        width : int, optional
        file_name : str
            filename parameter will be written to
        where : int, optional
            index location to insert parameter
        after : int, optional
            index location - 1 to insert parameter
        replace : bool
            replace an existing parameter when True, default is False

        """

        add = self._check_before_add(name, values, replace)

        if not add and replace:
            self.remove_record(name)
            add = True

        if add:
            if file_name is None:
                file_name = self.parameter_files[-1]
                where = None
                after = None

            new_record = ParameterRecord(
                name=name,
                values=values,
                dimensions=dimensions,
                datatype=datatype,
                width=width,
                file_name=file_name,
            )

            super(PrmsParameters, self).add_record(
                new_record, where=where, after=after
            )

    def remove_record(self, name):
        """
        Method to remove a parameter record

        Parameters
        ----------
        name: str
            parameter name

        """
        super(PrmsParameters, self).remove_record(name)

    def add_record_object(self, record_obj, replace=True):
        """
        Method to add a ParameterRecord object

        record_obj : ParameterRecord object
            ParameterRecord object
        replace : bool
            boolean flag that allows record replacement, default is True
        """
        add = self._check_before_add(
            record_obj.name, record_obj.values, replace
        )

        if not add and replace:
            self.remove_record(record_obj.name)
            add = True

        if add:
            if record_obj.file_name is None:
                record_obj.file_name = self.parameter_files[-1]

            super(PrmsParameters, self).add_record(record_obj)

    def write(self, name=None):
        """
        Method to write the PrmsParameters object to PRMS parameter files

        Parameters
        ----------
        name : str, optional
            file name

        """
        if name is not None:
            if len(self.parameter_files) > 1:
                raise NotImplementedError(
                    "Cannot change the name of " "more than one parameter file"
                )
            else:
                file_list = [name]

        else:
            file_list = self.parameter_files

        for ifile, filename in enumerate(file_list):
            with open(filename, "w") as fid:
                if ifile == 0:
                    txt = "Generated by Gsflow python Package....\n"
                    txt = txt + "Version : --.--"
                    fid.write(txt)
                    fid.write("\n** Dimensions **")
                # write dimension
                for record in self.parameters_list:
                    # write dimension first
                    if record.section != "Dimensions":
                        continue

                    if name is not None:
                        record.write(fid)

                    else:
                        if ifile == 0 and record.section == "Dimensions":
                            if os.path.normpath(
                                record.file_name
                            ) == os.path.normpath(filename):
                                record.write(fid)
                ##
                # write param
                if ifile == 0:
                    fid.write("\n")
                    fid.write("** Parameters **")
                for record in self.parameters_list:
                    # write dimension first
                    if record.section == "Dimensions":
                        continue

                    if name is not None:
                        record.write(fid)

                    else:
                        if os.path.normpath(
                            record.file_name
                        ) == os.path.normpath(filename):
                            record.write(fid)

                fid.write("\n")


class ParameterRecord(RecordBase):
    """
    ParameterRecord is a class for storing parameters

    Parameters
    ----------
    name : str
        parameter name
    values : list
        parameter values
    dimensions : list
        dimensions of the record
    datatype : int
        datatype flag of the record
    width : int, optional
    file_name : str
        parameter file to write the parameter to

    """

    def __init__(
        self,
        name=None,
        values=None,
        dimensions=None,
        datatype=None,
        width=10,
        file_name=None,
    ):

        super(ParameterRecord, self).__init__(name, values, datatype)

        self.file_name = file_name

        self._dimensions = dimensions
        if self._dimensions:
            self.dimensions_names = [nm[0] for nm in dimensions]
            self.dims = [nm[1] for nm in dimensions]

        self.width = width

        # the parameter belong to dimension sections or parameters section?
        if self._dimensions is None:
            if len(self._values) > 1:
                raise ValueError(
                    "Error : Values in the Dimension section must be scalar"
                )
            self.dimensions_names = None
            self.dims = None
            self.datatype = 1  # always integer
            self.ndim = 1
            self.section = "Dimensions"

        else:
            # the record belongs to parameters section
            self.ndim = len(self.dimensions_names)
            if self.ndim < 1:
                raise ValueError("No dimension names are specified")

            if self._values.ndim > 1:
                pass
            else:
                pass

            self.section = "Parameters"

            # number of values
            _nvalues = np.prod(np.array(self.dims))
            if _nvalues != self.nvalues:
                err = (
                    "Summation of values in all dimensions "
                    "is not equal to number of values"
                )
                raise ValueError(err)

    def __getitem__(self, item):
        return self.values[item]

    def __setitem__(self, key, value):
        self.values[key] = value

    def __add__(self, other):
        self._check_if_compatible()
        self.values += other
        return self

    def __mul__(self, other):
        self._check_if_compatible()
        self.values *= other
        return self

    def __sub__(self, other):
        self._check_if_compatible()
        self.values -= other

        return self

    def __truediv__(self, other):
        self._check_if_compatible()
        self.values /= other
        return self

    def __pow__(self, power):
        self._check_if_compatible()
        self.values **= power
        return self

    def _check_if_compatible(self):
        if self.datatype not in (1, 2, 3):
            raise AssertionError(
                "Mathematical operation cannot be performed on strings"
            )
        elif self.datatype == 1:
            raise AssertionError(
                "Mathematical operation not supported for integer dtypes"
            )

    @property
    def values(self):
        """
        np.ndarray of record values
        """
        return self._values

    @values.setter
    def values(self, new_values):
        if isinstance(new_values, ParameterRecord):
            new_values = new_values.values

        self._check_values(new_values)
        if (
            np.prod(np.array(self.dims)) != self.nvalues
        ) and self.section != "Dimensions":
            err = (
                " The number of values is not "
                "compatible with the dimensions"
            )
            raise ValueError(err)

        # change data type
        self._check_dtype()

    def export_nc(self, f, modflow, **kwargs):
        """
        Method to export netcdf

        Parameters
        ----------
        f : str or fp.export.NetCdf
            filename to write the parameter to (*.nc)
        modflow : object
            fp.modflow.Modflow or gsflow.modflow.Modflow object

        kwargs : dict
            keyword arguments

        Notes
        -----
        NetCdf export relies on flopy, so at the moment will
        only work for GSFLOW models where PRMS has the same
        discretization as the modflow grid

        """
        from ..utils.netcdf import param2netcdf

        f = param2netcdf(f, modflow, self, **kwargs)
        return f

    def _write_dimension(self, fid):
        """
        Write method for dimensions ParameterRecord

        Parameters
        ----------
        fid : File object

        """
        fid.write("\n")
        fid.write("####\n")
        fid.write(self.name)
        # write values
        for val in self.values:
            fid.write("\n")
            fid.write(str(val))

    def _write_parameter(self, fid):
        """
        Write method for parameters ParameterRecord

        Parameters
        ----------
        fid : File object

        """
        fid.write("\n")
        fid.write("####\n")
        fid.write(self.name)
        fid.write(" ")
        fid.write("{}\n".format(self.width))
        # write number of dimension
        fid.write(str(self.ndim))
        # write dimension names

        for nm in self.dimensions_names:
            fid.write("\n")
            fid.write(nm)

        # write nvalues
        fid.write("\n")
        fid.write(str(self.nvalues))

        # write datatype
        fid.write("\n")
        fid.write(str(self.datatype))

        # write values
        for val in self.values:
            fid.write("\n")
            fid.write(str(val))

    def write(self, fid):
        """
        Method to write to an open file

        Parameters
        ----------
        fid : File object

        """
        if self.section == "Dimensions":
            self._write_dimension(fid)
        else:
            self._write_parameter(fid)

    def __repr__(self):
        try:
            return self.name
        except:
            return "Parameter Record"

    def __str__(self):
        if self._dimensions is None:  # print record for dimension
            s = [
                "\n",
                "####",
                "\n",
                self.name,
                "\n",
                str(self.values[0]),
                "\n",
                "####",
            ]
            s = "".join(s)
            return s
        else:  # print regular param.
            s = [
                "\n",
                "####",
                "\n",
                self.name,
                " ",
                str(self.width),
                "\n",
                str(self.ndim),
                "\n",
            ]
            for dim_nam in self.dimensions_names:
                s += [dim_nam, "\n"]

            s += [str(self.nvalues), "\n", str(self.datatype)]

            # write values
            for i, val in enumerate(self.values):
                if i > 3:
                    s.append(".\n.\n.")
                    break
                s += ["\n", str(val)]

            s += ["\n", "####"]
            s = "".join(s)
            return s

    def from_dict(self):
        pass

    def to_dict(self):
        pass
