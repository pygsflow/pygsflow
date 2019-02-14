# -*- coding: utf-8 -*-
"""
@ comment JL

Author, email, and version information is supposed to
go into either the setup.py or a version.py, email.py, author.py

this is not an appropriate place for it
"""
__author__ = """Ayman Alzraiee"""
__email__ = 'aalzraiee@usgs.gov'
__version__ = '0.1.0'

from .control import Control, ControlFile, ControlRecord
from .prms.prms_model import Prms, PrmsModel
from .prms.prms_parameter import Parameters, PrmsParameters, ParameterRecord
from .prms.prms_data import Prms_data, PrmsData
from .gsflow import Gsflow, GsflowModel
from . import utils
from . import modflow
from . import output
from . import prms

"""
@ comment JL

what's up with all of the pass statements
they literally do nothing in your code. and I've counted
at least 20 of them that can be removed without changing 
functionality.
"""



