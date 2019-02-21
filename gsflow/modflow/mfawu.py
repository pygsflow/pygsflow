import flopy
import numpy as np
from flopy.pakbase import Package
from flopy.utils import MfList
from flopy.utils.recarray_utils import create_empty_recarray
from gsflow.utils.io import multi_line_strip
from flopy.utils.optionblock import OptionBlock
from collections import OrderedDict


class ModflowAg(Package):
    """
    The ModflowAg class is used to build read, write, and edit data
    from the MODFLOW-NWT AG package.

    Parameters
    ----------
    model
    extension
    options
    unitnumber
    filenames
    """
    _options = OrderedDict([('noprint', OptionBlock.simple_flag),
                            ('irrigation_diversion', {OptionBlock.dtype: np.bool_,
                                                      OptionBlock.nested: True,
                                                      OptionBlock.n_nested: 2,
                                                      OptionBlock.vars: OrderedDict(
                                                          [('numirrdiversions', OptionBlock.simple_int),
                                                           ('maxcellsdiversions', OptionBlock.simple_int)]
                                                      )}),
                            ('irrigation_well', {OptionBlock.dtype: np.bool_,
                                                 OptionBlock.nested: True,
                                                 OptionBlock.n_nested: 2,
                                                 OptionBlock.vars: OrderedDict(
                                                     [("numirrwells", OptionBlock.simple_int),
                                                      ('maxcellswell', OptionBlock.simple_int)]
                                                 )}),
                            ('supplemental_well', {OptionBlock.dtype: np.bool_,
                                                   OptionBlock.nested: True,
                                                   OptionBlock.n_nested: 2,
                                                   OptionBlock.vars: OrderedDict(
                                                       [("numsupwells", OptionBlock.simple_int),
                                                        ("maxdiversions", OptionBlock.simple_int)]
                                                   )}),
                            ('maxwell', {OptionBlock.dtype: np.bool_,
                                         OptionBlock.nested: True,
                                         OptionBlock.n_nested: 1,
                                         OptionBlock.vars: OrderedDict(
                                             [('nummaxwell', OptionBlock.simple_int)]
                                         )}),
                            ('tabfiles', OptionBlock.simple_tabfile),
                            ('phiramp', OptionBlock.simple_float),
                            ('etdemand', {OptionBlock.dtype: np.bool_,
                                          OptionBlock.nested: True,
                                          OptionBlock.n_nested: 1,
                                          OptionBlock.vars: OrderedDict(
                                              [('accel', OptionBlock.simple_float)]
                                          )}),
                            ('trigger', OptionBlock.simple_flag),
                            ('timeseries_diversion', OptionBlock.simple_flag),
                            ('timeseries_well', OptionBlock.simple_flag),
                            ('timeseries_diversionet', OptionBlock.simple_flag),
                            ('timeseries_wellet', OptionBlock.simple_flag),
                            ('diversionlist', {OptionBlock.dtype: np.bool_,
                                               OptionBlock.nested: True,
                                               OptionBlock.n_nested: 1,
                                               OptionBlock.vars: OrderedDict(
                                                   [('unit_diversionlist', OptionBlock.simple_int)]
                                               )}),
                            ('welllist', {OptionBlock.dtype: np.bool_,
                                          OptionBlock.nested: True,
                                          OptionBlock.n_nested: 1,
                                          OptionBlock.vars: OrderedDict(
                                              [('unit_diversionlist', OptionBlock.simple_int)]
                                          )}),
                            ('wellirrlist', {OptionBlock.dtype: np.bool_,
                                             OptionBlock.nested: True,
                                             OptionBlock.n_nested: 1,
                                             OptionBlock.vars: OrderedDict(
                                                 [('unit_wellirrlist', OptionBlock.simple_int)]
                                             )}),
                            ('diversionirrlist', {OptionBlock.dtype: np.bool_,
                                                  OptionBlock.nested: True,
                                                  OptionBlock.n_nested: 1,
                                                  OptionBlock.vars: OrderedDict(
                                                      [('unit_diversionirrlist', OptionBlock.simple_int)]
                                                  )}),
                            ('wellcbc', {OptionBlock.dtype: np.bool_,
                                         OptionBlock.nested: True,
                                         OptionBlock.n_nested: 1,
                                         OptionBlock.vars: OrderedDict(
                                             [('unitcbc', OptionBlock.simple_int)]
                                         )})
                            ])

    def __init__(self, model,
                 extension="ag", options=None, unitnumber=None, filenames=None):

        # setup the package parent class
        if unitnumber is None:
            unitnumber = ModflowAg.defaultunit()

        if filenames is None:
            filenames = [None]
        elif isinstance(filenames, str):
            filenames = [filenames]

        name = [ModflowAg.ftype()]
        units = [unitnumber]
        extra = [""]

        # set package name
        fname = [filenames[0]]

        super(ModflowAg, self).__init__(model,extension=extension,
                                        name=name, unit_number=units,
                                        extra=extra, filenames=fname)

        # set up class
        self.heading = "# {} package for {}, generated " \
                       "by gsflopy".format(self.name[0],
                                           model.version_types[model.version])
        self.url = "ag.htm"
        # todo: package specific code

    @staticmethod
    def defaultunit():
        return 69

    @staticmethod
    def ftype():
        return "AG"

    @staticmethod
    def get_empty(numrecords, maxells=0, block="well"):
        """
        Creates an empty record array corresponding to the block data type
        it is associated with.

        Parameters
        ----------
            numrecords :
            maxells :
            block :

        Returns:

        """
        dtype = ModflowAg.get_default_dtype(maxells=maxells, block=block)
        return create_empty_recarray(numrecords, dtype, default_value=-1.0E+10)

    @staticmethod
    def get_default_dtype(maxells=0, block="well"):
        """

        Parameters
        ----------
            maxells :
            block :

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

        elif block == "irrdiversion_modflow":
            dtype = [("segid", np.int), ("numcell", np.int),
                     ("period", np.float), ("triggerfact", np.float)]

            for i in range(maxells):
                dtype += [("i{}".format(i), np.int),
                          ("j{}".format(i), np.int),
                          ("eff_fact{}".format(i), np.float),
                          ("field_fact{}".format(i), np.float)]

        elif block == "irrdiversion_gsflow":
            dtype = [("segid", np.int), ("numcell", np.int),
                     ("period", np.float), ("triggerfact", np.float)]

            for i in range(maxells):
                dtype += [("hru_id{}".format(i), np.int),
                          ("dum{}".format(i), np.int),
                          ("hru_eff_fact{}".format(i), np.float),
                          ("hru_field_fact{}".format(i), np.float)]

        else:
            dtype = None

        return np.dtype(dtype)

    def write_file(self, check=False):
        pass


    @staticmethod
    def load(f, model, nper=0, method="gsflow", ext_unit_dict=None):
        """

        Parameters
        ----------
        f
        model
        ext_unit_dict

        Returns
        -------
            ModflowAg object
        """

        with open(f) as mfag:

            # strip the file header if it exists
            while True:
                line = multi_line_strip(mfag)
                if line:
                    break

            # read the options block
            options = OptionBlock.load_options(mfag, ModflowAg)

            nummaxwell = 0
            if options.nummaxwell is not None:
                nummaxwell = options.nummaxwell

            line = multi_line_strip(mfag)

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
                    # todo: set data to a timeseries object (recarray...)
                    nrec = len(t)
                    time_series = ModflowAg.get_empty(nrec, block="time series")

                    for ix, rec in enumerate(t):
                        if rec[0] in ('welletall', 'wellall'):
                            time_series[ix] = (rec[0], -999, rec[-1])
                        else:
                            time_series[ix] = tuple(rec[:3])

            # read item 1-2 well_list
            if "well list" in line:
                # read item 2
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

                    # check if this is block 2a
                    if isinstance(options.tabfiles, np.recarray):
                        well = ModflowAg.get_empty(nrec, block='tabfile_well')
                    else:
                        well = ModflowAg.get_empty(nrec, block="well")

                    for ix, rec in enumerate(t):
                        well[ix] = tuple(rec[:4])

            irr_diversion = {}
            irr_well = {}
            sup_well = {}
            for per in range(nper):
                if 'stress period' in line:
                    line = multi_line_strip(mfag)

                while True:
                    if 'irrdiversion' in line:
                        nrec = int(multi_line_strip(mfag).split()[0])
                        # model.version2 will need to be changed
                        # for pure flopy compatibility if migrated
                        if model.version2 == "gsflow":
                            irr = ModflowAg.get_empty(nrec, block="irrdiversion_gsflow")
                        else:
                            irr = ModflowAg.get_empty(nrec, block="irrdiversion_modflow")

                        irr = _read_block_6_or_10(mfag, nrec, irr)
                        irr_diversion[per] = irr
                        line = multi_line_strip(mfag)

                    if 'irrwell' in line:
                        pass

                    if 'supwel' in line:
                        pass

            print(line)
            print('break')



        return ModflowAg(model,)


def _read_block_6_or_10(fobj, nrec, recarray):
    t = []
    for _ in range(nrec):
        t1 = []
        ll = multi_line_strip(fobj).split()
        t1.append(ll[:4])

        for numcell in range(int(ll[1])):
            t1.append(multi_line_strip(fobj).split()[:4])

        t.append(t1)

    if len(t) > 0:
        for ix, rec in enumerate(t):
            for ix2, name in enumerate(recarray):
                recarray[ix][name] = rec[ix][ix2]

    return recarray
