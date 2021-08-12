import numpy as np
import platform
import os
from ..utils.gsflow_io import multi_line_strip


class CRT(object):
    """
    Class method to build, edit, write, and run a model for the
    Cascade Routing Tool

    Parameters
    ----------
    hru_type : np.ndarray
        hru_type array for PRMS in structured (row, col) format
    elev : np.ndarray
        elevation array for each hru
    outflow_hrus: list, tuple, np.ndarray
        iteratable of row, column locations of outflow hrus
    stream_cells: list, tuple, np.ndarray
        nx5 list of [[row, col, seg, reach, on_off],] that
        describe the stream cells
    hruflg : bool
        boolean flag to indicate whether the hrus are numbered differently
        than the grid-cells (ex. upper left corner is 1). Default is False
    strmflg : bool
        boolean flag that indicates how CRT interacts with streams.
        Default is True, which means selected stream reaches will be
        included in the cascade computation
    flowflg : bool
        boolean flag that determines how outflow fractions are calculated.
        Default is True, outflow is partitioned on the basis of relative slope
    visflg : bool
        boolean flag that determines if data will be written for later
        visualization. Default is False
    iprn : bool
        boolean flag that determines how output is provided. If True output
        is written in outputstat.txt, If False, file is not written. Defaults
        to True.
    ifill : bool
        boolean flag that determines if the CRT fill procedure is used.
        Default is True
    dpit : float
        real value less than or equal to 0.1 to which HRU elevation will
        be adjusted during an iteration of the CRT fill procedure.
    outitmax : int
        Integer value equal to number of CRT fill procedure iterations
    model_ws : str
        directory path to write model output to
    hru_ids : dict
        dictionary of hru_id : cell_id relationships for all active hrus.
        Required parameter when hruflg=True
    x_loc : np.ndarray
        numpy array of x-cell centers for all cells in the model.
        Required parameter when visflg=True
    y_loc : np.ndarray
        numpy array of y-cell centers for all cells in the model.
        Required parameter when visflg=True

    """

    def __init__(
        self,
        hru_type,
        elev,
        outflow_hrus,
        stream_cells,
        hruflg=False,
        strmflg=True,
        flowflg=True,
        visflg=False,
        iprn=True,
        ifill=True,
        dpit=0.1,
        outitmax=10000,
        model_ws=".",
        hru_ids=None,
        x_loc=None,
        y_loc=None,
    ):

        self.__sws = os.path.abspath(os.path.dirname(__file__))
        self.__exe_name = os.path.join(
            self.__sws, "..", "..", "bin", "CRT_1.3.1"
        )
        self.hru_type = np.array(hru_type, dtype=int)
        self._numactive = np.count_nonzero(self.hru_type)
        self.hruflg = hruflg
        self.strmflg = strmflg
        self.flowflg = flowflg
        self.visflg = visflg
        self.iprn = iprn
        self.ifill = ifill
        self.dpit = dpit
        self.outitmax = outitmax
        self.elev = self._check_shape(elev, "elev")
        self.outflow_hrus = self._check_and_expand(outflow_hrus)
        self.stream_cells = self._check_and_expand(stream_cells)
        self.hru_ids = hru_ids
        if hru_ids is not None:
            if not isinstance(hru_ids, dict):
                raise TypeError("hru_ids must be a dictionary of pairs")
            elif len(hru_ids) != self._numactive:
                raise AssertionError(
                    "hru_ids not same size as number of active cells"
                )
            else:
                pass

        self.x_loc = x_loc
        if x_loc is not None:
            self.x_loc = self._check_shape(x_loc, "x_loc")
        self.y_loc = y_loc
        if y_loc is not None:
            self.y_loc = self._check_shape(y_loc, "y_loc")

        self.model_ws = model_ws

    @property
    def shape(self):
        """
        Returns the hru_type shape

        """
        return self.hru_type.shape

    @property
    def nrow(self):
        """
        Returns number of rows in model

        """
        return self.shape[0]

    @property
    def ncol(self):
        """
        Returns the number of columns in the model

        """
        return self.shape[1]

    @property
    def numoutflowhrus(self):
        """
        Returns the number of outflow_hrus

        """
        if self.outflow_hrus is not None:
            return len(self.outflow_hrus)

    @property
    def nreach(self):
        """
        Returns the number of stream locations (stream reaches)

        """
        if self.stream_cells is not None:
            return len(self.stream_cells)

    def _check_shape(self, array, name):
        """
        Method to check an arrays shape

        Parameters
        ----------
        array : np.array
        name : str
            variable name

        Returns
        -------
            np.array
        """
        array = np.array(array, dtype=float)
        if array.shape != self.shape:
            raise AssertionError(
                "{} shape does not match hru_type shape".format(name)
            )

        return array

    def _check_and_expand(self, array):
        """
        Method to check and expand the size of an array

        Parameters
        ----------
        array : np.array

        Returns
        -------
            np.array
        """
        if not isinstance(array, np.ndarray):
            array = np.array(array)

        if len(array.shape) == 1:
            array = np.expand_dims(array, axis=0)

        return np.asarray(array, dtype=int)

    def _write_hru_casc(self, model_ws):
        """
        Method to write HRU_CASC file

        Parameters
        ----------
        model_ws : str
            optional path to write model files to

        """
        with open(os.path.join(model_ws, "HRU_CASC.DAT"), "w") as foo:
            foo.write(
                "{:d} {:d} {:d} {:d} {:d} {:d} {} {:d}\n".format(
                    int(self.hruflg),
                    int(self.strmflg),
                    int(self.flowflg),
                    int(self.visflg),
                    int(self.iprn),
                    int(self.ifill),
                    self.dpit,
                    int(self.outitmax),
                )
            )
            np.savetxt(foo, self.hru_type, fmt="%d", delimiter=" ")

    def _write_land_elev(self, model_ws):
        """
        Method to write LAND_ELEV file

        Parameters
        ----------
        model_ws : str
            optional path to write model files to

        """
        if self.elev is None:
            print(
                "WARNING: skipping LAND_ELEV.DAT, CRT will not run without it"
            )
            return

        with open(os.path.join(model_ws, "LAND_ELEV.DAT"), "w") as foo:
            foo.write("{} {}\n".format(self.nrow, self.ncol))
            np.savetxt(foo, self.elev, fmt="%d", delimiter=" ")

    def _write_outflow_hrus(self, model_ws):
        """
        Method to write OUTFLOW_HRU file

        Parameters
        ----------
        model_ws : str
            optional path to write model files to

        """
        if self.outflow_hrus is None:
            print(
                "WARNING: skipping OUTFLOW_HRU.DAT, CRT will not run without it"
            )
            return

        with open(os.path.join(model_ws, "OUTFLOW_HRU.DAT"), "w") as foo:
            foo.write("{:d}\n".format(self.numoutflowhrus))
            for ix, loc in enumerate(self.outflow_hrus):
                foo.write(
                    "{:d} {:d} {:d}\n".format(ix + 1, int(loc[0]), int(loc[1]))
                )

    def _write_stream_cells(self, model_ws):
        """
        Method to write STREAM_LOCATIONS file

        Parameters
        ----------
        model_ws : str
            optional path to write model files to

        """
        if self.outflow_hrus is None:
            print(
                "WARNING: skipping STREAM_CELLS.DAT, CRT will not run without it"
            )
            return

        with open(os.path.join(model_ws, "STREAM_CELLS.DAT"), "w") as foo:
            foo.write("{:d}\n".format(self.nreach))
            for row in self.stream_cells:
                foo.write("{:d} {:d} {:d} {:d} {:d}\n".format(*row))

    def _write_hru_identifiers(self, model_ws):
        """
        Method to write HRU_ID file

        Parameters
        ----------
        model_ws : str
            optional path to write model files to

        """
        if not self.hruflg:
            return
        else:
            if self.hru_ids is None:
                raise AssertionError("hru_ids must be supplied when hruflg=1")

            with open(os.path.join(model_ws, "HRU_ID.DAT"), "w") as foo:
                foo.write("{:d}\n".format(self._numactive))
                for key, value in self.hru_ids.items():
                    foo.write("{:d} {:d}\n".format(int(key), int(value)))

    def _write_xy(self, model_ws):
        """
        Method to write XY file

        Parameters
        ----------
        model_ws : str
            optional path to write model files to

        """
        if not self.visflg:
            return
        else:
            if self.x_loc is None or self.y_loc is None:
                raise AssertionError(
                    "x_loc and y_loc must be supplied when visflg=1"
                )

            with open(os.path.join(model_ws, "XY.DAT"), "w") as foo:
                hru_id = 1
                for i, row in enumerate(self.x_loc):
                    for j, val in enumerate(row):
                        foo.write(
                            "{:d} {} {}\n".format(
                                hru_id, self.x_loc[i, j], self.y_loc[i, j]
                            )
                        )

    def write_input(self, model_ws=None):
        """
        Method to write CRT data to input files

        Parameters
        ----------
        model_ws : str
            optional path to write model files to

        """
        if model_ws is None:
            model_ws = self.model_ws
        else:
            print("Changing model_ws to {}".format(model_ws))
            self.model_ws = model_ws

        self._write_hru_casc(model_ws)
        self._write_land_elev(model_ws)
        self._write_outflow_hrus(model_ws)
        self._write_stream_cells(model_ws)
        self._write_hru_identifiers(model_ws)
        self._write_xy(model_ws)

    def run_model(self, exe_name=None, model_ws=None):
        """
        Method to run CRT from python

        Parameters
        ----------
        exe_name : str
            optional CRT binary executable name
        model_ws : str
            optional, location of CRT model files

        Returns
        -------
        success, buff
        """
        from flopy.mbase import run_model

        if exe_name is None:
            exe_name = self.__exe_name

        if platform.system().lower() == "windows":
            if not exe_name.endswith(".exe"):
                exe_name += ".exe"

        if model_ws is None:
            model_ws = self.model_ws

        return run_model(
            exe_name,
            None,
            model_ws,
            normal_msg="cascades successfully generated",
        )

    @staticmethod
    def load(model_ws):
        """
        Method to load an existing CRT model into the CRT object for
        editing, writing, or running.

        Parameters
        ----------
        model_ws : str
            directory where CRT input files are located

        Returns
        -------
            CRT object

        """

        optional_files = ("HRU_ID.DAT", "XY.DAT")

        if not os.path.exists(os.path.join(model_ws, "HRU_CASC.DAT")):
            raise FileNotFoundError(os.path.join(model_ws, "HRU_CASC.DAT"))

        with open(os.path.join(model_ws, "HRU_CASC.DAT")) as foo:
            line = multi_line_strip(foo)
            t = line.strip().split()
            hruflg, strmflg, flowflg, visflg, iprn, ifill = [
                bool(int(i)) for i in t[:6]
            ]
            dpit = float(t[6])
            outitmax = int(t[7])
            hru_type = np.genfromtxt(foo, dtype=int)

        if not os.path.exists(os.path.join(model_ws, "LAND_ELEV.DAT")):
            raise FileNotFoundError(os.path.join(model_ws, "LAND_ELEV.DAT"))

        with open(os.path.join(model_ws, "LAND_ELEV.DAT")) as foo:
            line = multi_line_strip(foo)
            elev = np.genfromtxt(foo, dtype=float)

        if not os.path.exists(os.path.join(model_ws, "OUTFLOW_HRU.DAT")):
            raise FileNotFoundError(os.path.join(model_ws, "OUTFLOW_HRU.DAT"))

        with open(os.path.join(model_ws, "OUTFLOW_HRU.DAT")) as foo:
            line = multi_line_strip(foo)
            nouthru = int(line.split()[0])
            outflow_hrus = [[]] * nouthru

            i0 = 0
            while i0 < nouthru:
                line = multi_line_strip(foo)
                t = line.split()
                oid = int(t[0]) - 1
                rec = [int(t[1]), int(t[2])]
                outflow_hrus[oid] = rec
                i0 += 1

        if not os.path.exists(os.path.join(model_ws, "STREAM_CELLS.DAT")):
            raise FileNotFoundError(os.path.join(model_ws, "STREAM_CELLS.DAT"))

        stream_cells = []
        with open(os.path.join(model_ws, "STREAM_CELLS.DAT")) as foo:
            line = multi_line_strip(foo)
            nreach = int(line.split()[0])

            i0 = 0
            while i0 < nreach:
                line = multi_line_strip(foo)
                t = line.split()
                rec = [int(i) for i in t[:5]]
                stream_cells.append(rec)
                i0 += 1

        if not os.path.exists(os.path.join(model_ws, "HRU_ID.DAT")):
            hru_ids = None
        else:
            hru_ids = {}
            with open(os.path.join(model_ws, "HRU_ID.DAT")) as foo:
                line = multi_line_strip(foo)
                numactive = int(line.split()[0])

                i0 = 0
                while i0 < numactive:
                    line = multi_line_strip(foo)
                    t = line.split()
                    hru_id, cell_id = int(t[0]), int(t[1])
                    hru_ids[hru_id] = cell_id
                    i0 += 1

        if not os.path.exists(os.path.join(model_ws, "XY.DAT")):
            x_loc = None
            y_loc = None
        else:
            x_loc = np.zeros(hru_type.size, dtype=float)
            y_loc = np.zeros(hru_type.size, dtype=float)
            with open(os.path.join(model_ws, "XY.DAT")) as foo:
                i0 = 0
                while i0 < hru_type.size:
                    line = multi_line_strip(foo)
                    t = line.split()
                    hru_id = int(t[0]) - 1
                    x_loc[hru_id] = float(t[1])
                    y_loc[hru_id] = float(t[2])
                    i0 += 1

            x_loc.shape = hru_type.shape
            y_loc.shape = hru_type.shape

        return CRT(
            hru_type,
            elev,
            outflow_hrus,
            stream_cells,
            hruflg=hruflg,
            strmflg=strmflg,
            flowflg=flowflg,
            visflg=visflg,
            iprn=iprn,
            ifill=ifill,
            dpit=dpit,
            outitmax=outitmax,
            model_ws=model_ws,
            hru_ids=hru_ids,
            x_loc=x_loc,
            y_loc=y_loc,
        )
