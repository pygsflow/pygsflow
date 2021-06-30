import os
import numpy as np
import flopy
from flopy.modflow.mfag import _read_block_21_25_or_29
from flopy.utils.recarray_utils import create_empty_recarray
from gsflow.utils.gsflow_io import multi_line_strip
from flopy.utils.optionblock import OptionBlock
from collections import OrderedDict


class ModflowAg(flopy.modflow.ModflowAg):
    """
    The ModflowAg class is used to build read, write, and edit data
    from the MODFLOW-NWT AG package.

    Parameters
    ----------
    model : gsflow.modflow.Modflow object
        model object
    options : flopy.utils.OptionBlock object
        option block object
    time_series : np.recarray
        numpy recarray for the time series block
    well_list : np.recarray
        recarray of the well_list block
    irrdiversion : dict {per: np.recarray}
        dictionary of the irrdiversion block
    irrwell : dict {per: np.recarray}
        dictionary of the irrwell block
    supwell : dict {per: np.recarray}
        dictionary of the supwell block
    extension : str, optional
        default is .ag
    unitnumber : list, optional
        fortran unit number for modflow, default 69
    filenames : list, optional
        file name for ModflowAwu package to write input
    nper : int
        number of stress periods in the model

    Examples
    --------

    load a ModflowAg file

    >>> import gsflow
    >>> ml = gsflow.modflow.Modflow('awutest')
    >>> ag = gsflow.modflow.ModflowAg.load('test.awu', ml, nper=2)

    """

    _options = OrderedDict(
        [
            ("noprint", OptionBlock.simple_flag),
            (
                "irrigation_diversion",
                {
                    OptionBlock.dtype: np.bool_,
                    OptionBlock.nested: True,
                    OptionBlock.n_nested: 2,
                    OptionBlock.vars: OrderedDict(
                        [
                            ("numirrdiversions", OptionBlock.simple_int),
                            ("maxcellsdiversion", OptionBlock.simple_int),
                        ]
                    ),
                },
            ),
            (
                "irrigation_well",
                {
                    OptionBlock.dtype: np.bool_,
                    OptionBlock.nested: True,
                    OptionBlock.n_nested: 2,
                    OptionBlock.vars: OrderedDict(
                        [
                            ("numirrwells", OptionBlock.simple_int),
                            ("maxcellswell", OptionBlock.simple_int),
                        ]
                    ),
                },
            ),
            (
                "irrigation_pond",
                {
                    OptionBlock.dtype: np.bool_,
                    OptionBlock.nested: True,
                    OptionBlock.n_nested: 2,
                    OptionBlock.vars: OrderedDict(
                        [
                            ("numirrponds", OptionBlock.simple_int),
                            ("maxcellspond", OptionBlock.simple_int),
                        ]
                    ),
                },
            ),
            (
                "supplemental_well",
                {
                    OptionBlock.dtype: np.bool_,
                    OptionBlock.nested: True,
                    OptionBlock.n_nested: 2,
                    OptionBlock.vars: OrderedDict(
                        [
                            ("numsupwells", OptionBlock.simple_int),
                            ("maxdiversions", OptionBlock.simple_int),
                        ]
                    ),
                },
            ),
            (
                "maxwells",
                {
                    OptionBlock.dtype: np.bool_,
                    OptionBlock.nested: True,
                    OptionBlock.n_nested: 1,
                    OptionBlock.vars: OrderedDict(
                        [("nummaxwell", OptionBlock.simple_int)]
                    ),
                },
            ),
            (
                "maxponds",
                {
                    OptionBlock.dtype: np.bool_,
                    OptionBlock.nested: True,
                    OptionBlock.n_nested: 1,
                    OptionBlock.vars: OrderedDict(
                        [("nummaxpond", OptionBlock.simple_int)]
                    ),
                },
            ),
            ("tabfiles", OptionBlock.simple_tabfile),
            (
                "tabfileswell",
                {
                    OptionBlock.dtype: np.bool_,
                    OptionBlock.nested: True,
                    OptionBlock.n_nested: 2,
                    OptionBlock.vars: OrderedDict(
                        [
                            ("numtabwell", OptionBlock.simple_int),
                            ("maxvalwell", OptionBlock.simple_int),
                        ]
                    ),
                },
            ),
            (
                "tabfilespond",
                {
                    OptionBlock.dtype: np.bool_,
                    OptionBlock.nested: True,
                    OptionBlock.n_nested: 2,
                    OptionBlock.vars: OrderedDict(
                        [
                            ("numtabpond", OptionBlock.simple_int),
                            ("maxvalpond", OptionBlock.simple_int),
                        ]
                    ),
                },
            ),
            ("phiramp", OptionBlock.simple_flag),
            (
                "etdemand",
                {
                    OptionBlock.dtype: np.bool_,
                    OptionBlock.nested: True,
                    OptionBlock.n_nested: 1,
                    OptionBlock.vars: OrderedDict(
                        [("accel", OptionBlock.simple_float)]
                    ),
                },
            ),
            ("trigger", OptionBlock.simple_flag),
            ("timeseries_diversion", OptionBlock.simple_flag),
            ("timeseries_well", OptionBlock.simple_flag),
            ("timeseries_pond", OptionBlock.simple_flag),
            ("timeseries_diversionet", OptionBlock.simple_flag),
            ("timeseries_wellet", OptionBlock.simple_flag),
            ("timeseries_pondet", OptionBlock.simple_flag),
            (
                "diversionlist",
                {
                    OptionBlock.dtype: np.bool_,
                    OptionBlock.nested: True,
                    OptionBlock.n_nested: 1,
                    OptionBlock.vars: OrderedDict(
                        [("unit_diversionlist", OptionBlock.simple_int)]
                    ),
                },
            ),
            (
                "welllist",
                {
                    OptionBlock.dtype: np.bool_,
                    OptionBlock.nested: True,
                    OptionBlock.n_nested: 1,
                    OptionBlock.vars: OrderedDict(
                        [("unit_welllist", OptionBlock.simple_int)]
                    ),
                },
            ),
            (
                "wellirrlist",
                {
                    OptionBlock.dtype: np.bool_,
                    OptionBlock.nested: True,
                    OptionBlock.n_nested: 1,
                    OptionBlock.vars: OrderedDict(
                        [("unit_wellirrlist", OptionBlock.simple_int)]
                    ),
                },
            ),
            (
                "pondlist",
                {
                    OptionBlock.dtype: np.bool_,
                    OptionBlock.nested: True,
                    OptionBlock.n_nested: 1,
                    OptionBlock.vars: OrderedDict(
                        [("unit_pondlist", OptionBlock.simple_int)]
                    ),
                },
            ),
            (
                "pondirrlist",
                {
                    OptionBlock.dtype: np.bool_,
                    OptionBlock.nested: True,
                    OptionBlock.n_nested: 1,
                    OptionBlock.vars: OrderedDict(
                        [("unit_pondirrlist", OptionBlock.simple_int)]
                    ),
                },
            ),
            (
                "diversionirrlist",
                {
                    OptionBlock.dtype: np.bool_,
                    OptionBlock.nested: True,
                    OptionBlock.n_nested: 1,
                    OptionBlock.vars: OrderedDict(
                        [("unit_diversionirrlist", OptionBlock.simple_int)]
                    ),
                },
            ),
            (
                "wellcbc",
                {
                    OptionBlock.dtype: np.bool_,
                    OptionBlock.nested: True,
                    OptionBlock.n_nested: 1,
                    OptionBlock.vars: OrderedDict(
                        [("unitcbc", OptionBlock.simple_int)]
                    ),
                },
            ),
        ]
    )

    def __init__(
        self,
        model,
        options=None,
        time_series=None,
        well_list=None,
        pond_list=None,
        irrdiversion=None,
        irrwell=None,
        irrpond=None,
        supwell=None,
        extension="ag",
        unitnumber=None,
        filenames=None,
        nper=0,
    ):
        if unitnumber is None:
            unitnumber = 37

        self.irrigation_pond = False
        self.numirrponds = 0
        self.maxcellspond = 0
        self.maxponds = False
        self.nummaxpond = 0
        self.tabfileswell = False
        self.numtabwell = 0
        self.maxvalwell = 0
        self.tabfilespond = False
        self.numtabpond = 0
        self.maxvalpond = 0
        self.accel = 0
        self.timeseries_pond = False
        self.timeseries_pondet = False
        self.pondlist = False
        self.unit_pondlist = None
        self.pondirrlist = False
        self.unit_pondirrlist = None
        self.pond_list = pond_list
        self.irrpond = irrpond

        super(ModflowAg, self).__init__(
            model=model,
            options=options,
            time_series=time_series,
            well_list=well_list,
            irrdiversion=irrdiversion,
            irrwell=irrwell,
            supwell=supwell,
            extension=extension,
            unitnumber=unitnumber,
            filenames=filenames,
            nper=nper,
        )

    def write_file(self, check=False):
        """
        Write method for ModflowAwu

        Parameters
        ----------
        check: bool
            not implemented currently
        """
        if self.parent.version2 != "gsflow":
            super(ModflowAg, self).write_file(check=check)

        else:
            ws = self.parent.model_ws
            name = self.file_name[0]
            with open(os.path.join(ws, name), "w") as foo:
                foo.write(self.heading)

                # update options
                self.options.update_from_package(self)
                self.options.write_options(foo)

                # check if there is a timeseries block and write output
                if self.time_series is not None:
                    foo.write("# ag time series\n")
                    fmt = "{}   {:d}   {:d}\n"
                    foo.write("TIME SERIES \n")
                    for record in self.time_series:
                        if record["keyword"] in (
                            "welletall",
                            "wellall",
                            "pondetall",
                            "pondall",
                        ):
                            foo.write(
                                "{}   {:d}\n".format(
                                    record["keyword"], record["unit"]
                                ).upper()
                            )
                        else:
                            foo.write(fmt.format(*record).upper())

                    foo.write("END \n")

                # check if item 12 exists and write item 12 - 14
                if self.segment_list is not None:
                    foo.write("# segment list for irriagation diversions\n")
                    foo.write("SEGMENT LIST\n")
                    for iseg in self.segment_list:
                        foo.write("{:d}\n".format(iseg))

                    foo.write("END \n")

                # check if item 15 exists and write item 15 through 17
                if self.well_list is not None:
                    foo.write("# ag well list\n")
                    foo.write("WELL LIST \n")
                    if self.tabfiles or self.tabfileswell:
                        # item 16a
                        fmt16a = True
                        fmt16 = "{:d}   {:d}   {:d}   {:d}   {:d}\n"
                    else:
                        # item 16b
                        fmt16a = False
                        fmt16 = "{:d}   {:d}   {:d}   {:f}\n"

                    for record in self.well_list:
                        if fmt16a:
                            foo.write(
                                fmt16.format(
                                    record["unit"],
                                    record["tabval"],
                                    record["k"] + 1,
                                    record["i"] + 1,
                                    record["j"] + 1,
                                )
                            )
                        else:
                            foo.write(
                                fmt16.format(
                                    record["k"] + 1,
                                    record["i"] + 1,
                                    record["j"] + 1,
                                    record["flux"],
                                )
                            )

                    foo.write("END \n")

                # check if pond list exists and write items
                if self.pond_list is not None:
                    foo.write("# ag pond list\n")
                    foo.write("POND LIST \n")
                    if self.tabfilespond:
                        # item pond
                        fmt16 = "{:d}   {:d}   {:d}   {:d}\n"
                    else:
                        # item pond
                        fmt16 = "{:d}   {:f}   {:d}\n"

                    for record in self.pond_list:
                        if self.tabfilespond:
                            foo.write(
                                fmt16.format(
                                    record["unit"],
                                    record["tabval"],
                                    record["hru_id"] + 1,
                                    record["segid"],
                                )
                            )
                        else:
                            foo.write(
                                fmt16.format(
                                    record["hru_id"] + 1,
                                    record["q"],
                                    record["segid"],
                                )
                            )

                    foo.write("END \n")

                foo.write("# ag stress period data\n")
                for per in range(self._nper):
                    foo.write("STRESS PERIOD {}\n".format(per + 1))

                    # check for item 18 and write items 18 - 21
                    if self.irrdiversion is not None:
                        foo.write("IRRDIVERSION \n")

                        if self.trigger:
                            # item 20
                            fmt20 = "{:d}   {:d}   {:f}   {:f}\n"
                        else:
                            # item 20, if trigger option is false we need 0's
                            # for period and trigger fac.
                            fmt20 = "{:d}   {:d}   0   0\n"

                        if per in self.irrdiversion:
                            if np.isscalar(self.irrdiversion[per]):
                                # write item 19
                                foo.write("-1  \n")
                            else:
                                write = True
                                if per - 1 in self.irrdiversion and per > 0:
                                    if not np.isscalar(
                                        self.irrdiversion[per - 1]
                                    ):
                                        if (
                                            self.irrdiversion[per - 1].size
                                            == self.irrdiversion[per].size
                                            and self.irrdiversion[per].size > 0
                                            and self.irrdiversion[per - 1].size
                                            > 0
                                        ):
                                            if np.array_equal(
                                                self.irrdiversion[per],
                                                self.irrdiversion[per - 1],
                                            ):
                                                foo.write("-1 \n")
                                                write = False
                                if write:
                                    recarray = self.irrdiversion[per]

                                    # write item 19
                                    foo.write("{:d} \n".format(len(recarray)))

                                    # item 21b
                                    fmt21 = "{:d}   {:d}   {:f}   {:f}\n"

                                    for rec in recarray:
                                        num = rec["numcell"]
                                        if self.trigger:
                                            foo.write(
                                                fmt20.format(
                                                    rec["segid"],
                                                    rec["numcell"],
                                                    rec["period"],
                                                    rec["triggerfact"],
                                                )
                                            )
                                        else:
                                            foo.write(
                                                fmt20.format(
                                                    rec["segid"],
                                                    rec["numcell"],
                                                )
                                            )

                                        for i in range(num):
                                            foo.write(
                                                fmt21.format(
                                                    rec["hru_id{}".format(i)]
                                                    + 1,
                                                    rec["dum{}".format(i)] + 1,
                                                    rec[
                                                        "eff_fact{}".format(i)
                                                    ],
                                                    rec[
                                                        "field_fact{}".format(
                                                            i
                                                        )
                                                    ],
                                                )
                                            )

                        else:
                            # write item 19
                            foo.write("0  \n")

                    # check for item 22 and write 22 - 25
                    if self.irrwell is not None:
                        foo.write("IRRWELL \n")

                        if self.trigger:
                            # item 24
                            fmt24 = "{:d}   {:d}   {:f}   {:f}\n"
                        else:
                            # item 24
                            fmt24 = "{:d}   {:d}   0   0\n"

                        if per in self.irrwell:
                            if np.isscalar(self.irrwell[per]):
                                foo.write("-1  \n")

                            else:
                                write = True
                                if per - 1 in self.irrwell and per > 0:
                                    if not np.isscalar(self.irrwell[per - 1]):
                                        if (
                                            self.irrwell[per - 1].size
                                            == self.irrwell[per].size
                                            and self.irrwell[per - 1].size > 0
                                            and self.irrwell[per].size > 0
                                        ):
                                            if np.array_equal(
                                                self.irrwell[per],
                                                self.irrwell[per - 1],
                                            ):
                                                foo.write("-1 \n")
                                                write = False
                                if write:
                                    recarray = self.irrwell[per]

                                    # write item 23
                                    foo.write("{:d} \n".format(len(recarray)))

                                    fmt25 = "{:d}   {:d}   {:f}   {:f}\n"

                                    for rec in recarray:
                                        num = rec["numcell"]
                                        if self.trigger:
                                            foo.write(
                                                fmt24.format(
                                                    rec["wellid"] + 1,
                                                    rec["numcell"],
                                                    rec["period"],
                                                    rec["triggerfact"],
                                                )
                                            )
                                        else:
                                            foo.write(
                                                fmt24.format(
                                                    rec["wellid"] + 1,
                                                    rec["numcell"],
                                                )
                                            )

                                        for i in range(num):
                                            foo.write(
                                                fmt25.format(
                                                    rec["hru_id{}".format(i)]
                                                    + 1,
                                                    rec["dum{}".format(i)] + 1,
                                                    rec[
                                                        "eff_fact{}".format(i)
                                                    ],
                                                    rec[
                                                        "field_fact{}".format(
                                                            i
                                                        )
                                                    ],
                                                )
                                            )

                        else:
                            # write item 23
                            foo.write("0  \n")

                    # check if item 26 and write items 26 - 29
                    if self.supwell is not None:
                        foo.write("SUPWELL \n")

                        fmt28 = "{:d}   {:d}\n"

                        if per in self.supwell:
                            if np.isscalar(self.supwell[per]):
                                foo.write("-1  \n")
                            else:
                                write = True
                                if per - 1 in self.supwell and per > 0:
                                    if not np.isscalar(self.supwell[per - 1]):
                                        if (
                                            self.supwell[per - 1].size
                                            == self.supwell[per].size
                                            and self.supwell[per - 1].size > 0
                                            and self.supwell[per].size > 0
                                        ):
                                            if np.array_equal(
                                                self.supwell[per],
                                                self.supwell[per - 1],
                                            ):
                                                foo.write("-1 \n")
                                                write = False
                                if write:
                                    recarray = self.supwell[per]

                                    # write item 27
                                    foo.write("{:d} \n".format(len(recarray)))

                                    for rec in recarray:
                                        num = rec["numcell"]

                                        foo.write(
                                            fmt28.format(
                                                rec["wellid"] + 1,
                                                rec["numcell"],
                                            )
                                        )

                                        for i in range(num):
                                            if (
                                                rec["fracsupmax{}".format(i)]
                                                != -1e10
                                            ):
                                                foo.write(
                                                    "{:d}   {:f}   {:f}\n".format(
                                                        rec[
                                                            "segid{}".format(i)
                                                        ],
                                                        rec[
                                                            "fracsup{}".format(
                                                                i
                                                            )
                                                        ],
                                                        rec[
                                                            "fracsupmax{}".format(
                                                                i
                                                            )
                                                        ],
                                                    )
                                                )

                                            else:
                                                foo.write(
                                                    "{:d}   {:f}\n".format(
                                                        rec[
                                                            "segid{}".format(i)
                                                        ],
                                                        rec[
                                                            "fracsup{}".format(
                                                                i
                                                            )
                                                        ],
                                                    )
                                                )

                        else:
                            # write item 27
                            foo.write("0 \n")

                    if self.irrpond is not None and self.irrpond:
                        foo.write("IRRPOND \n")

                        if self.trigger:
                            # item 32
                            fmt32 = "{:d}   {:d}   {:f}   {:f}   {:d}\n"
                        else:
                            # item 32
                            fmt32 = "{:d}   {:d}   {:f}   {:d}\n"

                        fmt33 = "{:d}   {:d}   {:f}   {:f}\n"
                        if per in self.irrpond:
                            if np.isscalar(self.irrpond[per]):
                                foo.write("-1  \n")

                            else:
                                write = True
                                if per - 1 in self.irrpond:
                                    if not np.isscalar(self.irrpond[per - 1]):
                                        if (
                                            self.irrpond[per - 1].size
                                            == self.irrpond[per].size
                                            and self.irrpond[per - 1].size > 0
                                            and self.irrpond[per].size > 0
                                        ):
                                            if np.array_equal(
                                                self.irrpond[per],
                                                self.irrpond[per - 1],
                                            ):
                                                foo.write("-1 \n")
                                                write = False
                                if write:
                                    recarray = self.irrpond[per]
                                    # write item 31
                                    foo.write("{:d} \n".format(len(recarray)))
                                    for rec in recarray:
                                        num = rec["numcell"]
                                        if self.trigger:
                                            foo.write(
                                                fmt32.format(
                                                    rec["pond_id"] + 1,
                                                    rec["numcell"],
                                                    rec["period"],
                                                    rec["triggerfact"],
                                                    rec["flowthrough"],
                                                )
                                            )
                                        else:
                                            foo.write(
                                                fmt32.format(
                                                    rec["pond_id"] + 1,
                                                    rec["numcell"],
                                                    rec["period"],
                                                    rec["flowthrough"],
                                                )
                                            )

                                        for i in range(num):
                                            foo.write(
                                                fmt33.format(
                                                    rec["hru_id{}".format(i)]
                                                    + 1,
                                                    rec["dum{}".format(i)],
                                                    rec[
                                                        "eff_fact{}".format(i)
                                                    ],
                                                    rec[
                                                        "field_fact{}".format(
                                                            i
                                                        )
                                                    ],
                                                )
                                            )
                    foo.write("END\n")

    @staticmethod
    def get_empty(numrecords, maxells=0, block="well"):
        """
        Creates an empty record array corresponding to the block data type
        it is associated with.

        Parameters
        ----------
        numrecords : int
            number of records to create recarray with
        maxells : int, optional
            maximum number of irrigation links
        block : str
            str which indicates data set valid options are
            "well" ,
            "tabfile_well" ,
            "timeseries" ,
            "irrdiversion" ,
            "irrwell" ,
            "supwell"

        Returns:
            np.recarray

        """
        dtype = ModflowAg.get_default_dtype(maxells=maxells, block=block)
        return create_empty_recarray(numrecords, dtype, default_value=-1.0e10)

    @staticmethod
    def get_default_dtype(maxells=0, block="well"):
        """

        Parameters
        ----------
        maxells : int
             maximum number of irrigation links
        block : str
            str which indicates data set valid options are
            "well" ,
            "tabfile_well" ,
            "timeseries" ,
            "irrdiversion" ,
            "irrwell" ,
            "supwell"

        Returns
        -------
            dtype : (list, tuple)
        """
        if block == "well":
            dtype = [
                ("k", "int"),
                ("i", "int"),
                ("j", "int"),
                ("flux", "float"),
            ]

        elif block == "tabfile_well":
            dtype = [
                ("unit", "int"),
                ("tabval", "int"),
                ("k", "int"),
                ("i", "int"),
                ("j", "int"),
            ]

        elif block == "pond":
            dtype = [("hru_id", int), ("q", float), ("segid", int)]

        elif block == "tabfile_pond":
            dtype = [
                ("unit", int),
                ("tabval", int),
                ("hru_id", int),
                ("segid", int),
            ]

        elif block == "time series":
            dtype = [("keyword", "object"), ("id", "int"), ("unit", "int")]

        elif block == "irrdiversion":
            dtype = [
                ("segid", "int"),
                ("numcell", "int"),
                ("period", "float"),
                ("triggerfact", "float"),
            ]

            for i in range(maxells):
                dtype += [
                    ("hru_id{}".format(i), "int"),
                    ("dum{}".format(i), "int"),
                    ("eff_fact{}".format(i), "float"),
                    ("field_fact{}".format(i), "float"),
                ]

        elif block == "irrwell":
            dtype = [
                ("wellid", "int"),
                ("numcell", "int"),
                ("period", "float"),
                ("triggerfact", "float"),
            ]

            for i in range(maxells):
                dtype += [
                    ("hru_id{}".format(i), "int"),
                    ("dum{}".format(i), "int"),
                    ("eff_fact{}".format(i), "float"),
                    ("field_fact{}".format(i), "float"),
                ]
        elif block == "irrpond":
            dtype = [
                ("pond_id", int),
                ("numcell", int),
                ("period", float),
                ("triggerfact", float),
                ("flowthrough", int),
            ]

            for i in range(maxells):
                dtype += [
                    ("hru_id{}".format(i), int),
                    ("dum{}".format(i), int),
                    ("eff_fact{}".format(i), float),
                    ("field_fact{}".format(i), float),
                ]

        elif block == "supwell":
            dtype = [("wellid", "int"), ("numcell", "int")]

            for i in range(maxells):
                dtype += [
                    ("segid{}".format(i), "int"),
                    ("fracsup{}".format(i), "float"),
                    ("fracsupmax{}".format(i), "float"),
                ]

        else:
            raise NotImplementedError(
                "block type {}, not supported".format(block)
            )

        return np.dtype(dtype)

    @staticmethod
    def load(f, model, nper=0, ext_unit_dict=None):
        """
        Method to load the AG package from file

        Parameters
        ----------
        f : str
            filename
        model : gsflow.modflow.Modflow object
            model to attach the ag pacakge to
        nper : int
            number of stress periods in model
        ext_unit_dict : dict, optional

        Returns
        -------
            ModflowAwu object
        """
        if nper == 0:
            nper = model.nper

        with open(f) as mfag:

            # strip the file header if it exists
            while True:
                line = multi_line_strip(mfag)
                if line:
                    break

            # read the options block
            options = OptionBlock.load_options(mfag, ModflowAg)

            line = multi_line_strip(mfag)

            time_series = None
            if "time series" in line:
                # read time_series
                t = []
                while True:
                    line = multi_line_strip(mfag)
                    if line == "end":
                        line = multi_line_strip(mfag)
                        break

                    else:
                        t.append(line.split())

                if len(t) > 0:
                    nrec = len(t)
                    time_series = ModflowAg.get_empty(
                        nrec, block="time series"
                    )

                    for ix, rec in enumerate(t):
                        if rec[0] in (
                            "welletall",
                            "wellall",
                            "pondetall",
                            "pondall",
                        ):
                            time_series[ix] = (rec[0], -999, rec[-1])
                        else:
                            time_series[ix] = tuple(rec[:3])

            # read item 12-14
            if "segment list" in line:
                # read item 13, no need to store it's regenerated later
                while True:
                    line = multi_line_strip(mfag)
                    if line == "end":
                        line = multi_line_strip(mfag)
                        break

            # read item 15-17 well_list
            well = None
            if "well list" in line:
                # read item 16
                t = []
                while True:
                    line = multi_line_strip(mfag)
                    if line == "end":
                        line = multi_line_strip(mfag)
                        break

                    else:
                        t.append(line.split())

                if len(t) > 0:
                    nrec = len(t)

                    # check if this is block 16a
                    if isinstance(options.tabfileswell, np.recarray):
                        tf = True
                        well = ModflowAg.get_empty(nrec, block="tabfile_well")
                    else:
                        tf = False
                        well = ModflowAg.get_empty(nrec, block="well")

                    for ix, rec in enumerate(t):
                        if not tf:
                            k = int(rec[0]) - 1
                            i = int(rec[1]) - 1
                            j = int(rec[2]) - 1
                            well[ix] = (k, i, j, rec[3])
                        else:
                            k = int(rec[2]) - 1
                            i = int(rec[3]) - 1
                            j = int(rec[4]) - 1
                            well[ix] = (rec[0], rec[1], k, i, j)

            pond = None
            if "pond list" in line:
                t = []
                while True:
                    line = multi_line_strip(mfag)
                    if line == "end":
                        line = multi_line_strip(mfag)
                        break

                    else:
                        t.append(line.split())

                if len(t) > 0:
                    nrec = len(t)

                    # check if this is block 19a
                    if isinstance(options.tabfilespond, np.recarray):
                        tf = True
                        pond = ModflowAg.get_empty(nrec, block="tabfile_pond")
                    else:
                        tf = False
                        pond = ModflowAg.get_empty(nrec, block="pond")

                    for ix, rec in enumerate(t):
                        if not tf:
                            hru_id = int(rec[0]) - 1
                            pond[ix] = (hru_id, rec[1], rec[2])
                        else:
                            hru_id = int(rec[2]) - 1
                            pond[ix] = (rec[0], rec[1], hru_id, rec[3])

            maxcellsdiversion = 0
            if options.maxcellsdiversion is not None:
                maxcellsdiversion = options.maxcellsdiversion

            maxcellswell = 0
            if options.maxcellswell is not None:
                maxcellswell = options.maxcellswell

            maxdiversions = 0
            if options.maxdiversions is not None:
                maxdiversions = options.maxdiversions

            maxcellspond = 0
            if options.maxcellspond is not None:
                maxcellspond = options.maxcellspond

            irr_diversion = {}
            irr_well = {}
            irr_pond = {}
            sup_well = {}
            # get the stress period data from blocks 18 - 29
            for per in range(nper):
                while True:
                    if "stress period" in line:
                        line = multi_line_strip(mfag)

                    # block 18
                    elif "irrdiversion" in line:
                        # read block 19
                        nrec = int(multi_line_strip(mfag).split()[0])
                        if nrec == -1:
                            irr = np.copy(irr_diversion[per - 1])
                        else:
                            irr = ModflowAg.get_empty(
                                nrec,
                                maxells=maxcellsdiversion,
                                block="irrdiversion",
                            )

                            # read blocks 20 & 21
                            irr = _read_block_21_25_or_29(mfag, nrec, irr, 21)

                        irr_diversion[per] = irr
                        line = multi_line_strip(mfag)

                    # block 22
                    elif "irrwell" in line:
                        # read block 23
                        nrec = int(multi_line_strip(mfag).split()[0])
                        if nrec == -1:
                            irr = np.copy(irr_well[per - 1])
                        else:
                            irr = ModflowAg.get_empty(
                                nrec, maxells=maxcellswell, block="irrwell"
                            )

                            # read blocks 24 & 25
                            irr = _read_block_21_25_or_29(mfag, nrec, irr, 25)

                        irr_well[per] = irr
                        line = multi_line_strip(mfag)

                    elif "irrpond" in line:
                        nrec = int(multi_line_strip(mfag).split()[0])
                        if nrec == -1:
                            irrpond = np.copy(irr_pond[per - 1])
                        else:
                            irrpond = ModflowAg.get_empty(
                                nrec, maxells=maxcellspond, block="irrpond"
                            )

                            irrpond = _read_irrpond_block(
                                mfag, nrec, irrpond, options.trigger
                            )

                        irr_pond[per] = irrpond
                        line = multi_line_strip(mfag)

                    # block 26
                    elif "supwel" in line:
                        # read block 27
                        nrec = int(multi_line_strip(mfag).split()[0])
                        if nrec == -1:
                            sup = np.copy(sup_well[per - 1])

                        else:
                            sup = ModflowAg.get_empty(
                                nrec, maxells=maxdiversions, block="supwell"
                            )
                            # read blocks 28 & 29
                            sup = _read_block_21_25_or_29(mfag, nrec, sup, 29)

                        sup_well[per] = sup
                        line = multi_line_strip(mfag)

                    # block 30?
                    elif "end" in line:
                        if per == nper - 1:
                            break

                        line = multi_line_strip(mfag)
                        break

                    else:
                        raise ValueError(
                            "Something went wrong at: {}".format(line)
                        )

        return ModflowAg(
            model,
            options=options,
            time_series=time_series,
            well_list=well,
            pond_list=pond,
            irrwell=irr_well,
            irrdiversion=irr_diversion,
            irrpond=irr_pond,
            supwell=sup_well,
            nper=nper,
        )


def _read_irrpond_block(fobj, nrec, recarray, trigger):
    """
    Method to read irrigation pond block from AG package. GSFLOW only
    capability!

    Parameters
    ----------
    fobj : object
        open ag file object
    nrec : int
        number of records
    recarray : np.recarray
        empty recarray for pond block
    trigger : bool
        method to set trigger factor

    Returns
    -------
        np.recarray
    """
    t = []

    for _ in range(nrec):
        t1 = []
        ll = multi_line_strip(fobj).split()
        t1.append(int(ll[0]) - 1)
        t1 += ll[1:3]
        if trigger:
            t1 += ll[3:5]
        else:
            t1 += [0, ll[3]]

        for _ in range(int(t1[1])):
            tmp = multi_line_strip(fobj).split()[:4]
            tmp[0] = int(tmp[0]) - 1

            t1 += tmp

        t.append(t1)

    if len(t) > 0:
        for ix, rec in enumerate(t):
            for ix2, name in enumerate(recarray.dtype.names):
                if ix2 >= len(rec):
                    pass
                else:
                    recarray[name][ix] = rec[ix2]

    return recarray
