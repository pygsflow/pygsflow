# -*- coding: utf-8 -*-
"""
@ comment JL

Author, email, and version information is supposed to
go into either the setup.py or a version.py, email.py, author.py

this is not an appropriate place for it
"""
__author__ = """Ayman Alzraiee, Josh Larsen"""
__email__ = """aalzraiee@usgs.gov, jlarsen@usgs.gov"""
__version__ = '0.1.0'

from .control import Control, ControlFile, ControlRecord
from .prms.prms_model import PrmsModel
from .prms.prms_parameter import PrmsParameters, ParameterRecord
from .prms.prms_data import PrmsData
from .gsflow import GsflowModel
from . import utils
from . import modflow
from . import output
from . import prms
from . import modsim
