"""
data
====

GSadjust objects for non-Survey | Loop | Station objects (those are
represented as PyQt objects).
--------------------------------------------------------------------------------

The main data objects are:

ChannelList: Holds lists of observed data values (g, tilt, temp, etc.). The
lists are copied to ObsTreeStation objects.

Adjustment | AdjustmentOptions } Adjustment Results: The first contains
instances of the latter two. Holds the input data, options, and results of the
network adjustment.

Datum: Absolute-gravity observation or other reference for the relative-gravity
network. At least one datum is required for network adjustment.

Delta: Relative-gravity difference calculated from two station occupations.
May or may not include drift correction.

SimpleLoop | SimpleStation: This is a simplified representation of the main
ObsTreeSurvey | ObsTreeLoop | ObstreeStation PyQt datastore. The Simple...
objects exist because PyQt objects can't be serialized, which is mandatory for
saving the workspace using pickle (or json).

Tare: Represents an offset applied to the data

TimeSeries: Used only for tide correction.

This software is preliminary, provisional, and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition tha
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.
"""
from . import adjustment, analysis, channel, correction, datum, delta, tare, timeseries
from .adjustment import (
    AdjustedStation,
    Adjustment,
    AdjustmentOptions,
    AdjustmentResults,
)
from .channel import ChannelList
from .datum import Datum
from .delta import (
    Delta3Point,
    DeltaList,
    DeltaManual,
    DeltaNormal,
    create_delta_by_type,
)
from .tare import Tare
from .timeseries import TimeSeries
