import os
import pandas as pd
from ..utils import GsConstant


class PrmsDay(object):
    """
    PRMS python class to dynamically read, write, and edit prms day files.
    Plotting via pandas built in methods.

    Parameters:
    -----------
    f : str
        prms day file
    variable_name : str
        optional variable name, required if not loading a file
    dataframe : pd.DataFrame
        optional dataframe of cell by hru data if not loading from file

    """

    date_header = ["year", "month", "day", "hh", "mm", "sec"]

    nhru_declaration = None

    def __init__(self, f, variable_name=None, dataframe=None):

        self.__day_file = f
        self.__header = ""
        self.__name, self.__ws = self.__get_fname_and_ws(f)
        self.__init_ws = self.__get_fname_and_ws(f)[-1]
        self.__init_name = self.__get_fname_and_ws(f)[0]
        self.__day_variable_declaration = variable_name
        self.__orad_flag = False
        self.__data_startline = None
        self.__pd_hru_header = []
        self.__load_metadata()
        self.__dataframe = dataframe

    @property
    def name(self):
        return self.__name

    def change_file_ws(self, output_ws):
        if isinstance(output_ws, str):
            self.__ws = output_ws
        else:
            raise TypeError("path must be a string")

    @property
    def variable_name(self):
        """
        Method to get the PRMS climate variable name

        """
        return self.__day_variable_declaration

    def change_file_name(self, output_file_name):
        if isinstance(output_file_name, str):
            self.__name = output_file_name
        else:
            raise TypeError("Filename must be a string")

    def __load_metadata(self):
        """
        Method loads & retains __header and sets the day_variable_declaration,
        nhru_declaration, and orad_flag if necessary

        Returns
        -------
            df: (pandas dataframe) day file data
        """

        with open(self.__day_file) as f:
            for idx, l in enumerate(f):
                if idx == 0:
                    prms_header = "& modified by PRMS python"
                    self.__header = "{} {}\n".format(l.rstrip(), prms_header)

                elif l.startswith(GsConstant.DAY_FILE_VARIABLES):
                    h = l.rstrip().split(" ")
                    self.__day_variable_declaration = h[0]
                    self._nhru_declaration = int(h[-1])

                elif (
                    self.__day_variable_declaration is not None
                    and l.startswith("orad")
                ):
                    s = l.rstrip().split(" ")
                    orad = int(s[-1])
                    if orad == 1:
                        self.__orad_flag = True

                elif l.startswith("####"):
                    self.__data_startline = idx + 1
                    break

    def __load_day_file(self):
        """
        Method to load a prms day file into a python pandas dataframe.

        Returns:
            df: (pandas dataframe) day file data
        """

        if self.__day_variable_declaration is None:
            raise NotImplementedError(
                f"PrmsDay variable type not implemented for {self.__name}"
            )

        self.__pd_hru_header = [
            str(i + 1) for i in range(self._nhru_declaration)
        ]

        if self.__orad_flag:
            self.__pd_hru_header.append("orad")

        missing_value = -999.0
        df = pd.read_csv(
            self.__day_file,
            header=None,
            skiprows=self.__data_startline,
            delim_whitespace=True,
            na_values=[missing_value],
        )

        df.columns = PrmsDay.date_header + self.__pd_hru_header

        date = pd.Series(
            pd.to_datetime(
                df.year * 10000 + df.month * 100 + df.day, format="%Y%m%d"
            ),
            index=df.index,
        )

        df.index = pd.to_datetime(date)
        df.drop(PrmsDay.date_header, axis=1, inplace=True)
        df.columns.name = "input variables"
        df.index.name = "date"

        return df

    def __get_fname_and_ws(self, dpath):
        """
        Get the specific workspace (directory) for a day file

        Parmaeters
        ----------
        dpath : (string) data file name and path

        Returns
        -------
        fname, fpath : (tuple) data file name, data file directory
        """
        fpath, fname = os.path.split(dpath)
        return fname, fpath

    @property
    def dataframe(self):
        """
        Method to get a pandas dataframe of the input data

        """
        if self.__dataframe is None:
            self.__dataframe = self.__load_day_file()
        return self.__dataframe

    @dataframe.setter
    def dataframe(self, val):
        self.__dataframe = val

    def write(self):
        """
        User write method for PRMS day files

        Parameters:
        -----------
            ws: (str) directory to write day file
        """
        out_path = os.path.join(self.__ws, self.__name)
        in_path = os.path.join(self.__init_ws, self.__init_name)

        if not os.path.exists(self.__ws):
            os.mkdir(self.__ws)

        if out_path == in_path:
            df = self.dataframe

        if self.__dataframe is None:
            # No modifications were made, therefore we can directly write
            # from input file to the output file.
            with open(out_path, "w") as f_out:
                with open(in_path) as f_in:
                    for line in f_in:
                        f_out.write(line)

        else:
            # reconstruct original datafile format
            df = self.dataframe
            df["year"] = df.index.year
            df["month"] = df.index.month
            df["day"] = df.index.day
            df["hh"] = df["mm"] = df["sec"] = 0

            df = df[PrmsDay.date_header + self.__pd_hru_header]

            with open(out_path, "w") as f:
                f.write(self.__header)
                f.write(
                    "{} {}\n".format(
                        self.__day_variable_declaration, self._nhru_declaration
                    )
                )
                if self.__orad_flag:
                    f.write("orad 1\n")

                f.write("#" * 40 + "\n")

                df.to_csv(f, sep=" ", header=None, index=False, na_rep=-999)

            del df
            self.__dataframe = None

    @staticmethod
    def load_from_file(f):
        """
        Method to load a Day (climate by hru) file from file

        Parameters
        ----------
        f : str
            filename

        Returns
        -------
        PrmsDay object
        """
        return PrmsDay(f)
