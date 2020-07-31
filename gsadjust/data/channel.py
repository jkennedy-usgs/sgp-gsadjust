#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
data_objects.py
===============

GSadjust objects for ChannelList: Holds lists of observed data values (g, tilt, temp, etc.). The
lists are copied to ObsTreeStation objects.

This software is preliminary, provisional, and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition tha
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.
"""


class ChannelList:
    """
    Stores relative-gravity data as lists. Inherited by station object.
    Each item in each list (line, station, etc.) is a sample (i.e., single reading) from
    a relative gravimeter.
    """

    def __init__(self):
        self.meter_type = None
        # Attributes that are populated when a file is read
        self.line = []  # OLine number (Scintrex only)
        self.station = []  # Station name
        self.elev = []  # Elevation
        self.lat = []  # Station latitude
        self.long = []  # Station longitude
        self.height = []  # meter height

        # A station gravity value (.grav or .gmean) is calculated "on the fly" from .raw_grav, .etc,
        # and .tare as a property in the ObsTreeFunction class
        self.raw_grav = []  # gravity reading WITHOUT Earth tide correction
        # (etc is removed when data are read)
        # self.corr_g = []           # Corrected gravity value, no longer used; still some legacy usage in OL code
        self.sd = []  # Standard deviation (Scintrex only)
        self.etc = []  # Earth tide correction - can be updated to another model
        self.meter_etc = []  # Earth tide correction from meter
        self.tare = []  # tare magnitude (included in corr_g)

        self.tiltx = []  # tilt x
        self.tilty = []  # tilt y
        self.temp = []  # temperature, deg. C
        self.dur = []  # Measurement duration (Scintrex only)
        self.rej = []  # Accept or reject flag (Scintrex only)
        self.t = []  # Measurement time
        self.oper = []  # Operator
        self.dial = []  # Dial setting
        self.feedback = []  # Feedback (Burris only)
        self.meter = []  # Meter serial number
        self.keepdata = []  # 0/1 flags indicating which samples are to be used

    def __iter__(self):
        """
        Function to iterate over ChannelList. Useful because Scintrex and Burris
        meters have different fields; this returns only non-blank fields.
        :return: tuple of a single non-empty field in the ChannelList object:
        (field, list of values)
        """
        for attr, value in self.__dict__.items():
            a = getattr(self, attr)
            if type(a) is list:
                if len(a) > 0:
                    yield attr, a

    def extract_subset_idx(self, start_idx, end_idx):
        """
        Function for extracting shorter time series based on location in list.
        :param start_idx: list index of first value to return.
        :param end_idx: list index of last value to return.
        :return: ChannelList subset
        """
        templist = ChannelList()
        # Burris and Scintrex have different fields. This only subsets non-empty fields.
        for field, value in self:
            temp = value[start_idx : end_idx + 1]
            setattr(templist, field, temp)
        return templist
