#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
data/datum.py
===============

GSadjust objects for Datum: Absolute-gravity observation or other reference for the relative-gravity
network. At least one datum is required for network adjustment.

This software is preliminary, provisional, and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition tha
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.
"""


class Datum:
    """
    Object to store data typically associated with an absolute
    gravity measurement.

    At least one datum observation is required for adjustment; the default datum
    value (50000) can be used at one of the relative-gravity stations to see
    relative adjusted-gravity values. The gradient and measurement height are
    considered when transferring between the absolute- and relative-gravimeter
    measuring heights. The date field is not used in the adjustment (datum
    observations must be manually assigned to a survey).
    """

    def __init__(
        self, station, g=50000, sd=5, date='1/1/2000', gradient=-3.0, meas_height=0, checked=2,
    ):
        self.station = station
        self.g = g
        self.sd = sd
        self.date = date
        self.gradient = gradient
        self.meas_height = meas_height
        self.residual = -999
        self.checked = checked
        self.n_sets = None

    def __str__(self):
        """
        Override the built-in method 'print' when applied to such object
        """
        if self.checked == 2:
            in_use = '1'
        else:
            in_use = '0'
        return f'{in_use} {self.station} {self.date} {self.g:.2f} {self.sd:.2f} {self.residual:.2f}'
