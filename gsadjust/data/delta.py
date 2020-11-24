#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
data/delta.py
===============

Relative-gravity difference calculated from two station occupations. May or may not
include drift correction.

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
    """Create a new Delta subclass of the given type, passing in provided arguments.

    Parameters
    ----------
    delta_type : {'normal', 'three_point', 'list'}

    Returns
    -------
    Delta
        An instance of the Delta subclass.
    """
    cls = {
        "normal": DeltaNormal,
        "three_point": Delta3Point,
        "list": DeltaList,
    }.get(delta_type)
    if cls is None:
        raise TypeError
    # Filter kwargs; remove any None values:
    kwargs = {k: v for k, v in kwargs.items() if v is not None}
    return cls(*args, **kwargs)


class DeltaBase:
    """Object to store relative-gravity difference

    Deltas are calculated when one of the "populate delta table" menu commands
    is chosen. They may or may not be drift-corrected. Before a network
    adjustment is done, there is no residual and the default value of -999 is
    assigned.

    Three types of deltas implemented as subclasses.
    DeltaNormal ('normal'): calculated between two stations (g2 - g1)
    Delta3Point ('three_point'): calculated between one station, and an interpolated value
        betweeen two stations (Roman method).  self.station2 is a tuple with two
        stations, self.station1 is a single station.
    DeltaList ('list'): a delta is calculated from a list of deltas. Used with the Roman
        method to average three-point deltas. self.station2 is a list.
        self.station1 is None.

    3 standard deviation values are relevant:
    1) sd: from the drift correction method. This value is always shown on the Drift tab
    2) adj_sd: from a combination of the drift correction sd value and the options set under
        adjustment options (e.g., additive and multiplicative factors)
    3) assigned_sd: assigned by the user

    Either (2) or (3) is used in the network adjustment (i.e., shown on the Network
    Adjustment tab). If (3) is assigned, it overrides the value set by (2).

    (1) is a property; (2) and (3) are attributes, although 2 could also possibly be
    a property.

    Attributes
    ----------
    sta1 : ObsTreeStation or None
        Content varies depending on delta type.
    sta2 : ObsTreeStation, tuple, or list.
        Content varies depending on delta type.
    driftcorr : float
        Specified when a delta is created, based on the drift-correction type
    ls_drift : tuple
        Stores the degree of least squares adjustment.
        (loop name, degree of drift polynomial).
    adj_sd : float
        Stores the standard deviation used in the adjustment, which may
        differ from the one initially associated with the delta.
        (i.e., it can be user-specified)
    cal_coeff : int
        The user-specified calibration coefficient (default=1)
    loop : str
        Loop name
    assigned_sd : float
        User-specified delta-g
    checked : int
        For PyQt models

    """

    editable = False  # Only Normal deltas are editable.

    def __init__(
        self,
        sta1,
        sta2,
        driftcorr=0.0,
        ls_drift=None,
        checked=2,
        adj_sd=777,
        cal_coeff=1,
        loop=None,
        assigned_dg=None,
        **kwargs  # discard
    ):
        self.station1 = sta1
        self.station2 = sta2
        self.checked = checked
        self.adj_sd = adj_sd
        self.residual = -999.0
        self.driftcorr = driftcorr
        self.ls_drift = ls_drift
        self.cal_coeff = cal_coeff
        self.loop = loop
        self.assigned_dg = assigned_dg

    # General methods, used on all subclasses.

    def __str__(self):
        return_str = "{} {} {} {} {:0.3f} {:0.3f} {:0.3f} {:0.3f}".format(
            int(self.checked / 2),
            self.sta1,
            self.sta2,
            self.time_string(),
            self.dg,
            self.sd,
            self.adj_sd,
            self.residual,
        )
        return return_str

    def to_json(self):
        """Method for serializing.

        Deltas are re-created from station data and loop attributes when a workspace
        is loaded. Here, we only store the attributes that were potentially modified
        after the delta was created.

        Returns
        -------
        dict
            JSON-friently dict

        """
        return {
            "checked": self.checked,
            "type": self.type,
            "adj_sd": self.adj_sd,
            "assigned_dg": self.assigned_dg,
            "key": self.key,
            "loop": self.loop,
        }

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

    @property
    def key(self):
        """Used to match up saved delta attributes with the deltas that were newly
        created when a workspace is loaded.

        Returns
        -------
        str
        """
        return self.sta1 + self.sta2 + self.time_string()

    def time(self):
        if isinstance(self.station2, tuple):
            return self.sta1_t
        elif isinstance(self.station2, list):
            t = [delta.sta1_t for delta in self.station2]
            return np.mean(t)
        else:
            return (self.sta1_t + self.sta2_t) / 2

    def time_string(self):
        try:
            return num2date(self.time()).strftime("%Y-%m-%d %H:%M:%S")
        except Exception:
            return "-999"

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
    """DeltaNormal is a delta calculated between two stations (g2 - g1)

    Attributes
    ----------
    station1 : ObsTreeStation
    station2 : ObsTreeStation

    """

    type = "normal"

    @property
    def meter(self):
        return self.station2.meter[0]

    @property
    def sd(self):
        """
        Standard deviation, calculated by summing the station S.Ds in quadrature

        Returns
        -------
        float

        """
        try:
            if hasattr(self.station1, "assigned_sd"):
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
    def sta1(self):
        return self.station1.station_name

    @property
    def sta2(self):
        return self.station2.station_name

    @property
    def dg(self):
        # Roman correction: dg requires three station-observations
        gm1, gm2, = (
            self.station1.gmean(),
            self.station2.gmean(),
        )
        try:
            return gm2 - gm1 - self.driftcorr
        # If drift method == "netadj", self.driftcorr is a str
        except TypeError:
            return gm2 - gm2

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


class Delta3Point(DeltaBase):
    """Delta3Point provides a dg calculated between a station and an interpolated
    value at a second station.

    Dg is calculated by:
    -Interpolating the gravity value at a station, when that station (e.g., station
    A) has been visited twice, and some other station (e.g., station B) visited in
    between.
    -The value at station A is linearly interpolated at the time that station B was
    visited.
    -The interpolated value at A is differenced with the observed value at B.

    Attributes
    ----------
    station1 : ObsTreeStation
    station2 : tuple
        Tuple of length 2 with two observations (at different times) at a
        single station.
    """

    type = "three_point"
    default_SD = 3.0

    @property
    def meter(self):
        return self.station2[1].meter[0]

    @property
    def sd(self):
        """
        Default value. Ignored unless there is only one Delta3Point comprising a
        DeltaList.

        Delta3Points are not adjusted directly, only the derived dgs/sds stored in a
        DeltaList are used.

        Returns
        -------
        float
        """
        return self.default_SD

    @property
    def sta1(self):
        return self.station1.station_name

    @property
    def sta2(self):
        return self.station2[1].station_name

    @property
    def dg(self):
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
    """DeltaList provides a dg calculated from a list of deltas.

    Used with the Roman method to average three-point deltas.

    Attributes
    ----------
    station1 : None
    station2 : list
        List of Delta3Point deltas.

    """

    type = "list"

    @property
    def meter(self):
        return self.station2[0].meter

    @property
    def sd(self):
        """Standard deviation determined from the standard deviation of the 3Point
        dg values.

        If only one 3Point delta exists, returns the S.D for that delta (a default
        assigned value)

        Returns
        -------
        float

        """
        try:
            if len(self.station2) > 1:
                s = [np.abs(delta.dg) for delta in self.station2]
                return np.std(s)
            else:
                return self.station2[0].sd
        except:
            return -999

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
