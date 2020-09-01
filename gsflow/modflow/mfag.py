import os
import numpy as np
import flopy
from flopy.modflow.mfag import _read_block_21_25_or_29
from flopy.utils.recarray_utils import create_empty_recarray
from gsflow.utils.io import multi_line_strip
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
    _options = OrderedDict([('noprint', OptionBlock.simple_flag),
                            ('irrigation_diversion',
                             {OptionBlock.dtype: np.bool_,
                              OptionBlock.nested: True,
                              OptionBlock.n_nested: 2,
                              OptionBlock.vars: OrderedDict(
                                  [(
                                   'numirrdiversions', OptionBlock.simple_int),
                                   ('maxcellsdiversion',
                                    OptionBlock.simple_int)]
                              )}),
                            ('irrigation_well', {OptionBlock.dtype: np.bool_,
                                                 OptionBlock.nested: True,
                                                 OptionBlock.n_nested: 2,
                                                 OptionBlock.vars: OrderedDict(
                                                     [("numirrwells",
                                                       OptionBlock.simple_int),
                                                      ('maxcellswell',
                                                       OptionBlock.simple_int)]
                                                 )}),
                            ('supplemental_well', {OptionBlock.dtype: np.bool_,
                                                   OptionBlock.nested: True,
                                                   OptionBlock.n_nested: 2,
                                                   OptionBlock.vars: OrderedDict(
                                                       [("numsupwells",
                                                         OptionBlock.simple_int),
                                                        ("maxdiversions",
                                                         OptionBlock.simple_int)]
                                                   )}),
                            ('maxwells', {OptionBlock.dtype: np.bool_,
                                          OptionBlock.nested: True,
                                          OptionBlock.n_nested: 1,
                                          OptionBlock.vars: OrderedDict(
                                              [('nummaxwell',
                                                OptionBlock.simple_int)]
                                          )}),
                            ('tabfiles', OptionBlock.simple_tabfile),
                            ('phiramp', OptionBlock.simple_flag),
                            ('etdemand', OptionBlock.simple_flag),
                            ('trigger', OptionBlock.simple_flag),
                            ('timeseries_diversion', OptionBlock.simple_flag),
                            ('timeseries_well', OptionBlock.simple_flag),
                            (
                            'timeseries_diversionet', OptionBlock.simple_flag),
                            ('timeseries_wellet', OptionBlock.simple_flag),
                            ('diversionlist', {OptionBlock.dtype: np.bool_,
                                               OptionBlock.nested: True,
                                               OptionBlock.n_nested: 1,
                                               OptionBlock.vars: OrderedDict(
                                                   [('unit_diversionlist',
                                                     OptionBlock.simple_int)]
                                               )}),
                            ('welllist', {OptionBlock.dtype: np.bool_,
                                          OptionBlock.nested: True,
                                          OptionBlock.n_nested: 1,
                                          OptionBlock.vars: OrderedDict(
                                              [('unit_welllist',
                                                OptionBlock.simple_int)]
                                          )}),
                            ('wellirrlist', {OptionBlock.dtype: np.bool_,
                                             OptionBlock.nested: True,
                                             OptionBlock.n_nested: 1,
                                             OptionBlock.vars: OrderedDict(
                                                 [('unit_wellirrlist',
                                                   OptionBlock.simple_int)]
                                             )}),
                            ('diversionirrlist', {OptionBlock.dtype: np.bool_,
                                                  OptionBlock.nested: True,
                                                  OptionBlock.n_nested: 1,
                                                  OptionBlock.vars: OrderedDict(
                                                      [(
                                                       'unit_diversionirrlist',
                                                       OptionBlock.simple_int)]
                                                  )}),
                            ('wellcbc', {OptionBlock.dtype: np.bool_,
                                         OptionBlock.nested: True,
                                         OptionBlock.n_nested: 1,
                                         OptionBlock.vars: OrderedDict(
                                             [('unitcbc',
                                               OptionBlock.simple_int)]
                                         )})
                            ])

    def __init__(self, model, options=None, time_series=None, well_list=None,
                 irrdiversion=None, irrwell=None, supwell=None,
                 extension="ag", unitnumber=None, filenames=None, nper=0):
        if unitnumber is None:
            unitnumber = 37

        super(ModflowAg, self).__init__(model=model, options=options,
                                        time_series=time_series,
                                        well_list=well_list,
                                        irrdiversion=irrdiversion,
                                        irrwell=irrwell,
                                        supwell=supwell,
                                        extension=extension,
                                        unitnumber=unitnumber,
                                        filenames=filenames, nper=nper)

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
                        if record["keyword"] in ('welletall', 'wellall'):
                            foo.write("{}   {:d}\n".format(
                                record['keyword'],
                                record['unit']).upper())
                        else:
                            foo.write(fmt.format(*record).upper())

                    foo.write("END \n")

                # check if item 12 exists and write item 12 - 14
                if self.segment_list is not None:
                    foo.write('# segment list for irriagation diversions\n')
                    foo.write("SEGMENT LIST\n")
                    for iseg in self.segment_list:
                        foo.write("{:d}\n".format(iseg))

                    foo.write("END \n")

                # check if item 15 exists and write item 15 through 17
                if self.well_list is not None:
                    foo.write("# ag well list\n")
                    foo.write("WELL LIST \n")
                    if self.tabfiles:
                        # item 16a
                        fmt16a = True
                        fmt16 = "{:d}   {:d}   {:d}   {:d}   {:d}\n"
                    else:
                        # item 16b
                        fmt16a = False
                        fmt16 = "{:d}   {:d}   {:d}   {:f}\n"

                    for record in self.well_list:
                        if fmt16a:
                            foo.write(fmt16.format(record["unit"],
                                                   record["tabval"],
                                                   record["k"] + 1,
                                                   record["i"] + 1,
                                                   record["j"] + 1))
                        else:
                            foo.write(fmt16.format(record["k"] + 1,
                                                   record["i"] + 1,
                                                   record["j"] + 1,
                                                   record["flux"]))

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
                                recarray = self.irrdiversion[per]

                                # write item 19
                                foo.write("{:d} \n".format(len(recarray)))

                                # item 21b
                                fmt21 = "{:d}   {:f}   {:f}\n"

                                for rec in recarray:
                                    num = rec['numcell']
                                    if self.trigger:
                                        foo.write(fmt20.format(
                                            rec['segid'],
                                            rec['numcell'],
                                            rec['period'],
                                            rec['triggerfact']))
                                    else:
                                        foo.write(fmt20.format(rec['segid'],
                                                               rec['numcell']))

                                    for i in range(num):
                                        foo.write(fmt21.format(
                                            rec['hru_id{}'.format(i)] + 1,
                                            rec["eff_fact{}".format(i)],
                                            rec['field_fact{}'.format(i)]))

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
                                recarray = self.irrwell[per]

                                # write item 23
                                foo.write("{:d} \n".format(len(recarray)))

                                fmt25 = "{:d}   {:d}   {:f}   {:f}\n"

                                for rec in recarray:
                                    num = rec['numcell']
                                    if self.trigger:
                                        foo.write(fmt24.format(
                                            rec['wellid'] + 1,
                                            rec['numcell'],
                                            rec['period'],
                                            rec['triggerfact']))
                                    else:
                                        foo.write(fmt24.format(
                                            rec['wellid'] + 1,
                                            rec['numcell']))

                                    for i in range(num):
                                        foo.write(fmt25.format(
                                            rec['hru_id{}'.format(i)] + 1,
                                            rec['dum{}'.format(i)] + 1,
                                            rec["eff_fact{}".format(i)],
                                            rec['field_fact{}'.format(i)]))

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
                                recarray = self.supwell[per]

                                # write item 27
                                foo.write("{:d} \n".format(len(recarray)))

                                for rec in recarray:
                                    num = rec['numcell']

                                    foo.write(fmt28.format(rec["wellid"] + 1,
                                                           rec["numcell"]))

                                    for i in range(num):
                                        if rec["fracsupmax{}".format(i)] != -1e+10:
                                            foo.write(
                                                "{:d}   {:f}   {:f}\n".format(
                                                    rec['segid{}'.format(i)],
                                                    rec['fracsup{}'.format(i)],
                                                    rec['fracsupmax{}'.format(i)]))

                                        else:
                                            foo.write("{:d}   {:f}\n".format(
                                                rec['segid{}'.format(i)],
                                                rec['fracsup{}'.format(i)]))

                        else:
                            # write item 27
                            foo.write("0 \n")

                    foo.write("END \n")

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
        return create_empty_recarray(numrecords, dtype, default_value=-1.0E+10)

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
            dtype = [('k', np.int), ('i', np.int),
                     ('j', np.int), ('flux', np.float)]

        elif block == "tabfile_well":
            dtype = [('unit', np.int), ('tabval', np.int),
                     ('k', np.int), ('i', np.int), ('j', np.int)]

        elif block == "time series":
            dtype = [('keyword', np.object), ('id', np.int),
                     ('unit', np.int)]

        elif block == "irrdiversion":
            dtype = [("segid", np.int), ("numcell", np.int),
                     ("period", np.float), ("triggerfact", np.float)]

            for i in range(maxells):
                dtype += [("hru_id{}".format(i), np.int),
                          ("eff_fact{}".format(i), np.float),
                          ("field_fact{}".format(i), np.float)]

        elif block == "irrwell":
            dtype = [("wellid", np.int), ("numcell", np.int),
                     ("period", np.float), ("triggerfact", np.float)]

            for i in range(maxells):
                dtype += [("hru_id{}".format(i), np.int),
                          ("dum{}".format(i), np.int),
                          ("eff_fact{}".format(i), np.float),
                          ("field_fact{}".format(i), np.float)]

        elif block == "supwell":
            dtype = [("wellid", np.int), ("numcell", np.int)]

            for i in range(maxells):
                dtype += [("segid{}".format(i), np.int),
                          ("fracsup{}".format(i), np.float),
                          ("fracsupmax{}".format(i), np.float)]

        else:
            raise NotImplementedError(
                "block type {}, not supported".format(block))

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
                    time_series = ModflowAg.get_empty(nrec,
                                                      block="time series")

                    for ix, rec in enumerate(t):
                        if rec[0] in ('welletall', 'wellall'):
                            time_series[ix] = (rec[0], -999, rec[-1])
                        else:
                            time_series[ix] = tuple(rec[:3])

            # read item 12-14
            segments = None
            if "segment list" in line:
                # read item 13
                t = []
                while True:
                    line = multi_line_strip(mfag)
                    if line == "end":
                        line = multi_line_strip(mfag)
                        break
                    else:
                        t.append(line.split())

                if len(t) > 0:
                    segments = []
                    for rec in t:
                        iseg = int(rec[0])
                        segments.append(iseg)

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
                    if isinstance(options.tabfiles, np.recarray):
                        tf = True
                        well = ModflowAg.get_empty(nrec, block='tabfile_well')
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

            maxcellsdiversion = 0
            if options.maxcellsdiversion is not None:
                maxcellsdiversion = options.maxcellsdiversion

            maxcellswell = 0
            if options.maxcellswell is not None:
                maxcellswell = options.maxcellswell

            maxdiversions = 0
            if options.maxdiversions is not None:
                maxdiversions = options.maxdiversions

            irr_diversion = {}
            irr_well = {}
            sup_well = {}
            # get the stress period data from blocks 18 - 29
            for per in range(nper):
                while True:
                    if 'stress period' in line:
                        line = multi_line_strip(mfag)

                    # block 18
                    elif 'irrdiversion' in line:
                        # read block 19
                        nrec = int(multi_line_strip(mfag).split()[0])
                        if nrec == -1:
                            irr = np.copy(irr_diversion[per - 1])
                        else:
                            if model.version2 == "gsflow":
                                irr = ModflowAg.get_empty(
                                    nrec,
                                    maxells=maxcellsdiversion,
                                    block="irrdiversion")
                            else:
                                irr = flopy.modflow.ModflowAg.get_empty(
                                    nrec,
                                    maxells=maxcellsdiversion,
                                    block="irrdiversion")

                            # read blocks 20 & 21
                            irr = _read_block_21_25_or_29(mfag, nrec, irr, 21)
                        irr_diversion[per] = irr
                        line = multi_line_strip(mfag)

                    # block 22
                    elif 'irrwell' in line:
                        # read block 23
                        nrec = int(multi_line_strip(mfag).split()[0])
                        if nrec == -1:
                            irr = np.copy(irr_well[per - 1])
                        else:
                            if model.version2 == "gsflow":
                                irr = ModflowAg.get_empty(
                                    nrec,
                                    maxells=maxcellswell,
                                    block="irrwell")
                            else:
                                irr = flopy.modflow.ModflowAg.get_empty(
                                    nrec,
                                    maxells=maxcellswell,
                                    block="irrwell")

                            # read blocks 24 & 25
                            irr = _read_block_21_25_or_29(mfag, nrec, irr, 25)

                        irr_well[per] = irr
                        line = multi_line_strip(mfag)

                    # block 26
                    elif 'supwel' in line:
                        # read block 27
                        nrec = int(multi_line_strip(mfag).split()[0])
                        if nrec == -1:
                            sup = np.copy(sup_well[per - 1])

                        else:
                            sup = ModflowAg.get_empty(
                                nrec,
                                maxells=maxdiversions,
                                block="supwell")
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
                            "Something went wrong at: {}".format(line))

        return ModflowAg(model, options=options, time_series=time_series,
                         well_list=well, irrwell=irr_well,
                         irrdiversion=irr_diversion,
                         supwell=sup_well, nper=nper)
