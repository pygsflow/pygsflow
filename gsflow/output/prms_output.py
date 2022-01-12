from ..utils import GsConstant
import pandas as pd
import datetime
import os
import warnings

warnings.simplefilter("always", PendingDeprecationWarning)


class StatVar(object):
    """
    Class to read in the statvar file into a pandas Dataframe

    Parameters
    ----------
    statvar_file : str
        name of the stat var file
    statvar_names : list, optional
        statvar_names will be determined from file if None
    statvar_elements : list, optional
        statvar_elements will be determined from file if None

    Examples
    --------

    load from control file

    >>> control = gsflow.ControlFile.load_from_file("gsflow.gsflow")
    >>> stats = StatVar.load_from_control_object(control)

    load from statistics file

    >>> stats = StatVar("mystatvar.txt")

    """

    def __init__(
        self, statvar_file, statvar_names=None, statvar_elements=None
    ):

        self.statvar_file = statvar_file
        self.statvar_names = statvar_names
        self.statvar_elements = statvar_elements
        self.stat_df = None

        self._load_statvar_file()

    @staticmethod
    def load_from_control_object(control):
        """
        Load the stats var from a ControlFile object

        Parameters
        ----------
        control : ControlFile object

        Returns
        -------
            Statistics object

        """
        statvar_file = control.get_values("stat_var_file")[0]
        ws, fil = os.path.split(statvar_file)
        if not ws:
            statvar_file = os.path.join(control.model_dir, fil)
        else:
            statvar_file = os.path.join(control.model_dir, statvar_file)

        statvar_flg = control.get_values("statsON_OFF")[0]
        if statvar_flg == 0:
            print("There is no statvar output since statsON_OFF = 0 ")
            return None

        statvar_names = control.get_values("statVar_names")
        statvar_elements = control.get_values("statVar_element")

        return StatVar(statvar_file, statvar_names, statvar_elements)

    def _load_statvar_file(self):
        """
        Loads the statvar file into memory as a pandas array

        """
        # check to see if statvar elements and names were supplied
        sve = False
        svn = False
        if self.statvar_elements is None:
            sve = True
            self.statvar_elements = []

        if self.statvar_names is None:
            svn = True
            self.statvar_names = []

        print("Loading the statvar output file .....")
        with open(self.statvar_file, "r") as fid:
            nvals = int(fid.readline().strip())
            var_names = []
            var_element = []
            for header in range(nvals):
                nm, elem = fid.readline().strip().split()

                if svn:
                    self.statvar_names.append(nm)
                if sve:
                    self.statvar_elements.append(elem)

                nm = nm + "_" + elem
                var_names.append(nm)
                var_element.append(int(elem))

            columns = ["ID"] + GsConstant.COLUMN_HEADER + var_names
            stat_df = pd.read_csv(fid, delim_whitespace=True, names=columns)
            Dates = []
            for index, irow in stat_df.iterrows():
                dt = datetime.datetime(
                    year=int(irow["Year"]),
                    month=int(irow["Month"]),
                    day=int(irow["Day"]),
                    hour=int(irow["Hour"]),
                    minute=int(irow["Minute"]),
                    second=int(irow["Second"]),
                )
                Dates.append(dt)

            stat_df["Date"] = Dates
            self.stat_df = stat_df

        print("Finished Load the statvar output file .....")

    """
    @ comment JL
    I think that plot should have the option to select an element/element(s)

    although the built in (pandas) method works, it makes it difficult to 
    layer two or three items on it. A custom plotting routine may be better
    suited to accomplish this... 
    """

    def plot(self):
        """
        Built in plotting method for the statvar object. Uses
        pandas built ins.

        """
        for i, name in enumerate(self.statvar_names):
            nm = name + "_" + self.statvar_elements[i]
            self.stat_df.plot(x="Date", y=nm)
