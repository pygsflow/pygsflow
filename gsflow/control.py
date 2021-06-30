from __future__ import absolute_import, division, print_function
import os
from .record_base import RecordBase
from .param_base import ParameterBase
from .utils import GsConstant
from .utils import gsflow_io
import warnings

warnings.simplefilter("always", PendingDeprecationWarning)


class ControlFile(ParameterBase):
    """
    Class to hold information about control file, also it reads and writes data.

    Parameters
    ----------
    records_list : list
        list of ControlRecord objects
    name : str
        file name
    model_dir : str
        model working directory
    header : list
        file header
    abs_path : bool
        optional flag to store control file path variables to abs paths,
        default is True

    Examples
    --------

    load from file

    >>> control = gsflow.ControlFile.load_from_file("gsflow.control")


    build a new empty object

    >>> control = gsflow.ControlFile(records_list=[])

    """

    def __init__(
        self,
        records_list,
        name="Control",
        model_dir="",
        header=None,
        abs_path=True,
    ):
        super(ControlFile, self).__init__(
            records_list, name=name, model_dir=model_dir, header=header
        )

        if header is None:
            self.header = ["Control File"]

        self.control_file = os.path.join(self.model_dir, self.name)

        if abs_path:
            self._make_pths_abs()

    @property
    def records_list(self):
        """
        Returns
        -------
            a list of ControlRecord objects
        """
        return self._records_list

    @staticmethod
    def load_from_file(control_file, abs_path=True):
        """
        Method to load and create a ControlFile object
        from a pre-built gsflow control file.

        Parameters
        ----------
        control_file : str
            control file path & name

        set_abs : bool
            optional flag to set control file variables to abs paths, default
            is True

        Returns
        -------
            ControlFile object

        """
        if not (os.path.isfile(control_file)):
            raise FileNotFoundError("Invalid file name ....")

        with open(control_file, "r") as fid:
            headers = []
            records_list = []
            EndOfFile = False
            _read_comments = True
            while True:
                if EndOfFile:
                    break
                record = fid.readline().strip()
                # read comments
                if _read_comments:
                    if "####" in record:
                        _read_comments = False
                        continue
                    headers.append(record)
                    continue

                # read records information
                field_name = record
                nvalues = int(fid.readline().strip())
                data_type = int(fid.readline().strip())
                values = []

                # loop over values
                while True:

                    record = fid.readline()
                    if record == "\n":
                        continue  # empty line
                    elif not record:  # end of the file
                        EndOfFile = True
                        curr_record = ControlRecord(
                            name=field_name, values=values, datatype=data_type
                        )
                        records_list.append(curr_record)
                        break

                    else:
                        record = record.strip()
                        if "####" in record:
                            curr_record = ControlRecord(
                                name=field_name,
                                values=values,
                                datatype=data_type,
                            )
                            records_list.append(curr_record)
                            break
                        else:
                            values.append(record)

        model_dir, name = os.path.split(control_file)

        return ControlFile(
            records_list,
            name=name,
            model_dir=model_dir,
            header=headers,
            abs_path=abs_path,
        )

    def _make_pths_abs(self):
        """
        Makes all file paths in control absoulte paths

        """
        for file in GsConstant.GSFLOW_FILES:

            if file in self.record_names:
                gs_fn = self.get_values(file)
                flist = []
                for ff in gs_fn:
                    abs_file = gsflow_io.get_file_abs(
                        control_file=self.control_file, fn=ff
                    )
                    flist.append(abs_file)
                self.set_values(file, flist)

    def _generate_attributes(self):
        for rec in self.records_list:
            setattr(self, rec.name, rec)

    def get_record(self, name):
        """
        Get a complete record object

        Parameters
        ----------
        name : str
            ControlRecord name

        Returns
        -------
            ControlRecord object

        """
        return super(ControlFile, self).get_record(name, ControlRecord)

    def get_values(self, name):
        """
        Get a record's values

        Parameters
        ----------
        name : str
            record name
        Returns
        -------
            np.ndarray or list

        """
        return super(ControlFile, self).get_values(name)

    def set_values(self, name, values):
        """
        Method to set new values to a control record

        Parameters
        ----------
        name : str
            control record name
        values : list
            list of values

        """

        super(ControlFile, self).set_values(name, values)

    def add_record(self, name=None, values=None, where=None, after=None):
        """
        Convience method to add a record to the control file

        Parameters
        ----------
        name : str
            record name
        values : list
            list of values
        where : int
            index location to insert record
        after : int
            index location - 1 to insert record

        """

        add = self._check_before_add(name=name, values=values)

        if add:
            new_record = ControlRecord(name=name, values=values)
            super(ControlFile, self).add_record(
                new_record, where=where, after=after
            )

    def remove_record(self, name):
        """
        Convenience method to remove a record from a control file

        Parameters
        ----------
        name : str
            control record name

        """
        super(ControlFile, self).remove_record(name)

    def write(self, name=None):
        """
        Method to write the control file

        Parameters
        ----------
        name : str, optional
            control file name

        """
        if name is None:
            filename = self.control_file
        else:
            filename = os.path.join(self.model_dir, name)

        with open(filename, "w") as fid:
            for iline, header in enumerate(self.header):
                if iline == 0:
                    txt = header.strip()
                else:
                    txt = "\n" + header.strip()
                fid.write(txt)

            for record in self.records_list:
                record.write(fid)
            fid.write("\n")


class ControlRecord(RecordBase):
    """
    ControlRecord is the object used for creating and editing
    control file record objects

    Parameters
    ----------
    name : str
        record name
    values : list
        list of values
    datatype : int
        integer datatype flag
    nvalues : int
        number of values in record

    Examples

    create a modflow_nam ControlRecord

    >>> rec = ControlRecord("modflow_nam", "gsflow_test.nam")

    """

    def __init__(self, name=None, values=None, datatype=None):

        super(ControlRecord, self).__init__(name, values, datatype)

    @property
    def values(self):
        return self._values

    @values.setter
    def values(self, values):
        self._values = values
        self._check_values(values)

        if len(values) != self.nvalues:
            print("Warning: the number of values is modefied")
        self._check_dtype()
        self._force_dtype()

    """
    @ comment JL
    
    these write methods can be cleaned up by using
    a list and then .join()
    
    Furthermore; we could template these and fill using
    a format string... Think about best method.
    
    Applies to PrmsParameters class too!
    """

    def write(self, fid):
        """
        Write method for a control record

        Parameters
        ----------
        fid : File object

        """
        fid.write("\n")
        fid.write("####")
        fid.write("\n")
        fid.write(self.name)

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

    def __repr__(self):
        try:
            return self.name
        except:
            return "control Record"

    def __str__(self):
        print_str = ""
        print_str = print_str + "\n"
        print_str = print_str + "####"
        print_str = print_str + "\n"
        print_str = print_str + self.name
        print_str = print_str + "\n"
        print_str = print_str + str(self.nvalues)
        print_str = print_str + "\n"
        print_str = print_str + str(self.datatype)

        # write values
        for i, val in enumerate(self.values):
            if i > 3:
                print_str = print_str + ".\n.\n."
                break
            print_str = print_str + "\n"
            print_str = print_str + str(val)

        print_str = print_str + "\n"
        print_str = print_str + "####"
        return print_str
