import flopy
import numpy as np
from flopy.pakbase import Package
from flopy.utils import MfList
from flopy.utils.recarray_utils import create_empty_recarray
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
    def get_empty():
        pass

    @staticmethod
    def get_default_dtype():
        pass

    def write_file(self, check=False):
        pass

    def load(self, f, model, ext_unit_dict=None):

        return ModflowAg(model,)

