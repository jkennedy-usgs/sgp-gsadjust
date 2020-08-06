#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
data/delta.py
===============

GSadjust objects for Delta: Relative-gravity difference calculated from two station occupations.
May or may not include drift correction.

This software is preliminary, provisional, and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition tha
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.
"""
import numpy as np
from matplotlib.dates import num2date


def create_delta_by_type(delta_type, *args, **kwargs):
    """
    Create a new Delta subclass of the given type, passing in provided arguments.

    :returns: an instance of the Delta subclass.
    """
    cls = {
        'normal': DeltaNormal,
        'assigned': DeltaAssigned,
        'three_point': Delta3Point,
        'list': DeltaList,
    }.get(delta_type)
    if cls is None:
        raise TypeError
    return cls(*args, **kwargs)


class DeltaBase:
    """
    Object to store relative-gravity difference

    Deltas are calculated when one of the "populate delta table" menu commands
    is chosen. They may or may not be drift-corrected. Before a network
    adjustment is done, there is no residual and the default value of -999 is
    assigned.

    Four types of deltas implemented as subclasses.
    DeltaNormal ('normal'): calculated between two stations (g2 - g1)
    Delta3Point ('three_point'): calculated between one station, and an interpolated value
        betweeen two stations (Roman method).  self.station2 is a tuple with two
        stations, self.station1 is a single station.
    DeltaList ('list'): a delta is calculated from a list of deltas. Used with the Roman
        method to average three-point deltas. self.station2 is a list.
        self.station1 is None.
    DeltaAssigned ('assigned'): delta-g is manually assigned. Created automatically when a
        'normal' delta is edited on the network adjustment tab. Roman
        ('three_point' and 'list') deltas can't be converted to 'assigned'.

    ls_drift: records degree of least squares adjustment: they must be identical
        for all deltas to use Gravnet. Its a tuple with (loop name, degree of
        drift polynomial).
    adj_sd: stores the standard deviation used in the adjustment, which may
        differ from the one initially associated
        with the delta.
    assigned_sd: Manually-edited std. dev.
    """

    editable = False  # Only Normal deltas are editable.

    def __init__(
        self,
        sta1,
        sta2,
        driftcorr=0.0,
        ls_drift=None,
        checked=2,
        adj_sd=999,
        cal_coeff=1,
        loop=None,
    ):
        self.station1 = sta1
        self.station2 = sta2
        self.checked = checked
        self.adj_sd = adj_sd
        self.residual = -999.0
        self.driftcorr = driftcorr
        self.ls_drift = ls_drift  # degree of drift, if included in network adjustment
        self.cal_coeff = cal_coeff
        self.loop = loop
        self.assigned_dg = None

    # General methods, used on all subclasses.

    @property
    def sd_for_adjustment(self):
        """
        Standard deviation used in the network adjustment.

        If the standard deviation hasn't changed (e.g., it's the default
        value 999), return the standard deviation initially associated with
        the delta.

        Default value (i.e., if not set in adjustment options) is self.sd.
        :return: float
        """
        if self.adj_sd == 999:
            return self.sd
        else:
            return self.adj_sd

    def __str__(self):
        # if self.checked == 2:
        #     in_use = '1'
        # else:
        #     in_use = '0'
        # Rarely a station time will be -999 (if all samples are unchecked)
        # if self.sta1_t() != -999:
        #     delta_time = num2date(self.sta1_t()).strftime('%Y-%m-%d %H:%M:%S')
        # else:
        #     delta_time = '-999'
        return_str = '{} {} {} {} {:0.3f} {:0.3f} {:0.3f} {:0.3f}'.format(
            int(self.checked / 2),
            self.sta1,
            self.sta2,
            self.time_string(),
            self.dg,
            self.sd,
            self.sd_for_adjustment,
            self.residual,
        )
        return return_str

    @property
    def sta1(self):
        raise NotImplementedError

    @property
    def sta2(self):
        raise NotImplementedError

    @property
    def sta1_t(self):
        raise NotImplementedError

    @property
    def sta2_t(self):
        raise NotImplementedError

    @property
    def sd(self):
        raise NotImplementedError

    def time(self):
        if isinstance(self.station2, tuple):
            return self.sta1_t
        elif isinstance(self.station2, list):
            t = [delta.sta1_t for delta in self.station2]
            return np.mean(t)
        else:
            # FIXME: Default to ObsTreeStation so we don't need the import.
            # Should generally avoid type checks, may be able to remove
            # these remaining two.
            return (self.sta1_t + self.sta2_t) / 2

    def time_string(self):
        try:
            return num2date(self.time()).strftime('%Y-%m-%d %H:%M:%S')
        except Exception:
            return '-999'

    @property
    def dg(self):
        return -999

    @property
    def duration(self):
        if type(self.station2) is tuple:
            return self.station2[1].tmean() - self.station2[0].tmean()
        elif type(self.station2) is list:
            dgs = [delta.dg for delta in self.station2]
            return np.mean(dgs)
        else:
            return self.sta2_t - self.sta1_t


class DeltaNormal(DeltaBase):
    """
    'normal': calculated between two stations (g2 - g1)
    """

    type = 'normal'

    @property
    def meter(self):
        return self.sta2.meter[0]

    @property
    def sd(self):
        """
        Standard deviation determined from drift correction. Default value for
        network adjustment.
        :return: float
        """
        try:
            if hasattr(self.station1, 'asd'):
                if self.station1.assigned_sd:
                    s = np.sqrt(
                        self.station2.assigned_sd ** 2 + self.station1.assigned_sd ** 2
                    )
                else:
                    s = np.sqrt(self.station2.stdev() ** 2 + self.station1.stdev() ** 2)
            else:
                s = np.sqrt(self.station2.stdev() ** 2 + self.station1.stdev() ** 2)
            return s
        except:
            return -999

    @property
    def key(self):
        """
        Used to match up adjustment results with observations, and to store
        simple Deltas
        :return: str
        """
        return (
            f"{self.type}"
            f"{self.station1.station_name}"
            f"{self.station1.station_count}"
            f"{self.station2.station_name}"
            f"{self.station2.station_count}"
        )

    @property
    def sta1(self):
        return self.station1.station_name

    @property
    def sta2(self):
        return self.station2.station_name

    @property
    def dg(self):
        # Roman correction: dg requires three station-observations
        gm1, gm2, = self.station1.gmean(), self.station2.gmean()
        return gm2 - gm1 - self.driftcorr

    @property
    def sta1_t(self):
        return self.station1.tmean()

    @property
    def sta2_t(self):
        return self.station2.tmean()

    @property
    def duration(self):
        if type(self.station2) is tuple:
            return self.station2[1].tmean() - self.station2[0].tmean()
        elif type(self.station2) is list:
            dgs = [delta.dg for delta in self.station2]
            return np.mean(dgs)
        else:
            return self.sta2_t - self.sta1_t


class DeltaAssigned(DeltaNormal):
    """
        'assigned': delta-g is manually assigned. Created automatically when a
            'normal' delta is edited on the network adjustment tab. Roman
            ('three_point' and 'list') deltas can't be converted to 'assigned'.
    """

    type = 'assigned'

    @property
    def dg(self):
        return self.assigned_dg


class Delta3Point(DeltaBase):
    """
    'three_point': calculated between one station, and an interpolated value
        betweeen two stations (Roman method).  self.station2 is a tuple with two
        stations, self.station1 is a single station.
    """

    type = 'three_point'

    @property
    def meter(self):
        return self.sta2[1].meter[0]

    @property
    def sd(self):
        """
        Standard deviation determined from drift correction. Default value for
        network adjustment.
        :return: float
        """
        return 3.0

    @property
    def sta1(self):
        return self.station1.station_name

    @property
    def sta2(self):
        return self.station2[1].station_name

    @property
    def dg(self):
        # Roman correction: dg requires three station-observations
        gm1, gm2a, gm2b, = (
            self.station1.gmean(),
            self.station2[1].gmean(),
            self.station2[0].gmean(),
        )
        tm1, tm2a, tm2b = (
            self.station1.tmean(),
            self.station2[0].tmean(),
            self.station2[1].tmean(),
        )

        sta2_dg = gm2a - gm2b
        time_prorate = (tm1 - tm2a) / (tm2b - tm2a)
        return (gm2b + sta2_dg * time_prorate) - gm1

    @property
    def sta1_t(self):
        return self.station1.tmean()

    @property
    def sta2_t(self):
        return self.station2[0].tmean()


class DeltaList(DeltaBase):
    """
    'list': a delta is calculated from a list of deltas. Used with the Roman
        method to average three-point deltas. self.station2 is a list.
        self.station1 is None.
    """

    type = 'list'

    @property
    def meter(self):
        return self.sta2.meter[0]

    @property
    def sd(self):
        """
        Standard deviation determined from drift correction. Default value for
        network adjustment.
        :return: float
        """
        try:
            if len(self.station2) > 1:
                s = [np.abs(delta.dg) for delta in self.station2]
                return np.std(s)
            else:
                return 3.0
        except:
            return -999

    @property
    def key(self):
        """
        Used to match up adjustment results with observations, and to store
        simple Deltas
        :return: str
        """
        return f"{self.type}{self.sta1}{self.sta2}{self.dg}"

    @property
    def sta1(self):
        if self.station2[0].dg > 0:
            return self.station2[0].station1.station_name
        else:
            return self.station2[0].station2[1].station_name

    @property
    def sta2(self):
        if self.station2[0].dg > 0:
            return self.station2[0].station2[1].station_name
        else:
            return self.station2[0].station1.station_name

    @property
    def dg(self):
        # Roman correction: return average of individual deltas
        dg_all = [np.abs(delta.dg) for delta in self.station2]
        return np.mean(dg_all)

    @property
    def sta1_t(self):
        return self.station2[0].station1.tmean()

    @property
    def sta2_t(self):
        return self.station2[0].station1.tmean()

    @property
    def duration(self):
        if type(self.station2) is tuple:
            return self.station2[1].tmean() - self.station2[0].tmean()
        elif type(self.station2) is list:
            dgs = [delta.dg for delta in self.station2]
            return np.mean(dgs)
        else:
            return self.sta2_t - self.sta1_t
