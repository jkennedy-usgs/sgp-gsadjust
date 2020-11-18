"""
data/tare.py
===============

GSadjust objects for Tare: Represents an offset applied to the data

This software is preliminary, provisional, and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition tha
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.
"""
from PyQt5.QtCore import QDateTime


class Tare:
    """
    Object to store relative-gravity difference

    Deltas are calculated when one of the "populate delta table" menu commands
    is chosen. They may or may not be drift-corrected. Before a network
    adjustment is done, there is no residual and the default value of -999 is
    assigned.

    It's tempting to store just the Station objects that make up the delta,
    and calculate dg, sd, etc. as needed. But that doesn't work too well with
    the Roman method because the delta-g is calculated between (1) a station
    and (2) an interpolated value between two occupations at the other station.
    """

    def __init__(self, datetime, tare):
        self.datetime = datetime
        self.tare = float(tare)
        self.checked = 2

    def __str__(self):
        if self.checked == 2:
            in_use = 'x'
        else:
            in_use = 'o'
        return f"{in_use} {self.datetime} {self.tare}"

    def to_json(self):
        return {
            'checked': self.checked,
            'datetime': str(self.datetime),
            'tare': self.tare,
        }

