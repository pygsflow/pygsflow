import os
import pandas as pd
import datetime
from ..utils import GsConstant
import warnings

warnings.simplefilter("always", PendingDeprecationWarning)


class PrmsData(object):
    """
    PrmsData is a class to load/edit/write PrmsData files

    Parameters
    ----------
        data_df : pd.Dataframe
        name : str
            DataFile file name to write
        model_dir : str
            model directory to write file to
        header : str
            prms data header

    Examples
    --------

    >>> data = gsflow.prms.PrmsData.load_from_file("mydatafile")

    """

    data_names = [
        "tmax",
        "tmin",
        "precip",
        "runoff",
        "pan_evap",
        "solrad",
        "from_data",
        "rain_day",
    ]

    def __init__(self, data_df, name="Data", model_dir="", header="\n"):
        self.name = name
        self.model_dir = model_dir
        self.data_df = data_df
        self.header = header

    @property
    def file_name(self):
        """
        Returns
        -------
            the data file path
        """
        return os.path.join(self.model_dir, self.name)

    @staticmethod
    def load_from_file(data_file):
        """
        Method to load a data file into a PrmsData object

        Parameters
        ----------
        data_file : str
            data file path and name

        Returns
        -------
            PrmsData object

        """

        model_dir, name = os.path.split(data_file)

        fid = open(data_file, "r")
        headers = fid.readline().strip()
        columns = []
        while True:
            line = fid.readline()
            if line.strip() == "" or line.strip()[0:2] == "//":
                continue

            if "####" in line:
                break

            if any(item in line for item in PrmsData.data_names):
                val_nm = line.strip().split()
                for val in range(int(val_nm[1])):
                    columns.append(val_nm[0] + "_" + str(val))

        columns = GsConstant.COLUMN_HEADER + columns
        data_pd = pd.read_csv(fid, delim_whitespace=True, names=columns)
        Dates = []
        for index, irow in data_pd.iterrows():
            dt = datetime.datetime(
                year=int(irow["Year"]),
                month=int(irow["Month"]),
                day=int(irow["Day"]),
                hour=int(irow["Hour"]),
                minute=int(irow["Minute"]),
                second=int(irow["Second"]),
            )
            Dates.append(dt)

        data_pd["Date"] = Dates
        fid.close()
        return PrmsData(
            data_df=data_pd, model_dir=model_dir, name=name, header=headers
        )

    def write(self, name=None):
        """
        Method to write PrmsData input to a PRMS Data file

        Parameters
        ----------
        name : str
            Data file file name, if none it uses the class stored name

        """
        if name is None:
            data_file = os.path.join(self.model_dir, self.name)

        else:
            data_file = os.path.join(self.model_dir, name)

        with open(data_file, "w") as fid:
            fid.write(self.header)
            fid.write("\n")
            columns = self.data_df.columns
            climate_data = []
            climate_count = {}
            climate_unique = []
            for col in columns:
                nm = col.split("_")[0]
                if nm in PrmsData.data_names:
                    climate_data.append(nm)
                    if not (nm in climate_unique):
                        climate_unique.append(nm)
                    if nm in climate_count.keys():
                        climate_count[nm] = climate_count[nm] + 1
                    else:
                        climate_count[nm] = 1

            # write headers
            for clim_name in climate_unique:
                line = clim_name + " " + str(climate_count[clim_name]) + "\n"
                fid.write(line)
            fid.write(
                "#########################################################################\n"
            )
            pd_to_write = self.data_df.copy()
            pd_to_write = pd_to_write.drop(["Date"], axis=1)
            pd_to_write.to_csv(
                fid, index=False, sep=" ", line_terminator="\n", header=False
            )
