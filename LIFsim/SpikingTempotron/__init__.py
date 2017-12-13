from . import _ST
from . import pyFunctions as __pf
from . import __plotPatterns__ as __pp

_ST.SpikingTempotron.getVoltageTrace = __pf.__SpikingTempotron_getVoltageTrace
_ST.SpikingTempotron.getVoltageTraceFromInputLayer = __pf.__SpikingTempotron_getVoltageTrace_1
_ST.Tempotron.w = __pf.__tempotron__w
_ST.SpikeTrain.plot = __pp._st_plot

from ._ST import *
from .pyFunctions import times2SpikeTrain
