#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
data_objects.py
===============

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
import datetime as dt

import numpy as np
from matplotlib.dates import date2num, num2date
from PyQt5 import QtCore
from scipy.interpolate import interp1d

import pyqt_models


###############################################################################
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


###############################################################################


def create_delta_by_type(delta_type, *args, **kwargs):
    """
    Create a new Delta subclass of the given type, passing in provided arguments.

    :returns: an instance of the Delta subclass.
    """
    cls = {
        'normal': DeltaNormal,
        'manual': DeltaManual,
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
    DeltaManual ('assigned'): delta-g is manually assigned. Created automatically when a
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
    def dg(self):
        raise NotImplementedError

    @property
    def sd(self):
        raise NotImplementedError

    def time(self):
        if type(self.station2) is pyqt_models.ObsTreeStation:
            return (self.sta1_t() + self.sta2_t()) / 2
        elif type(self.station2) is tuple:
            return self.sta1_t()
        elif type(self.station2) is list:
            t = [delta.sta1_t() for delta in self.station2]
            return np.mean(t)

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
            return self.sta2_t() - self.sta1_t()


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
            f"{self.station1.station_count"
            f"{self.station2.station_name}"
            f"{self.station2.station_count}"
        )

    @property
    def sta1(self):
        return self.station1.station_name

    @property
    def sta2(self):
        return self.station2.station_name

    def dg(self):
        # Roman correction: dg requires three station-observations
        gm1, gm2, = self.station1.gmean(), self.station2.gmean()
        dg = gm2 - gm1 - self.driftcorr

    def sta1_t(self):
        return self.station1.tmean()

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
            return self.sta2_t() - self.sta1_t()


class DeltaManual(DeltaNormal):
    """
        'assigned': delta-g is manually assigned. Created automatically when a
            'normal' delta is edited on the network adjustment tab. Roman
            ('three_point' and 'list') deltas can't be converted to 'assigned'.
    """

    type = 'manual'

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
        return sta2[1].meter[0]

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
        dg = (gm2b + sta2_dg * time_prorate) - gm1

    def sta1_t(self):
        return self.station1.tmean()

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

    def dg(self):
        # Roman correction: return average of individual deltas
        dg_all = [np.abs(delta.dg for delta in self.station2]
        dg = np.mean(dg_all)

    def sta1_t(self):
        return self.station2[0].station1.tmean()

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
            return self.sta2_t() - self.sta1_t()


###############################################################################
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

    def __init__(self, date, time, tare):
        self.date = date
        self.time = time
        self.tare = float(tare)
        self.checked = 2

    def __str__(self):
        if self.checked == 2:
            in_use = 'x'
        else:
            in_use = 'o'
        return f"{in_use} {self.date} {self.time} {self.tare}"

    @property
    def datetime(self):
        """
        :return: Python datetime object.
        """
        return QtCore.QDateTime(self.date, self.time).toPyDateTime()


###############################################################################
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
        self,
        station,
        g=50000,
        sd=5,
        date='1/1/2000',
        gradient=-3.0,
        meas_height=0,
        checked=2,
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
        return '{} {} {} {:.2f} {:.2f} {:.2f}'.format(
            in_use, self.station, self.date, self.g, self.sd, self.residual
        )


###############################################################################
class AdjustmentResults:
    """
    Object to store least-squares adjustment statistics.
    """

    def __init__(self):
        self.n_deltas, self.n_deltas_notused = 0, 0
        self.n_datums, self.n_datums_notused = 0, 0
        self.n_unknowns = 0
        self.max_dg_residual, self.min_dg_residual = 0, 0
        self.max_datum_residual, self.min_datum_residual = 0, 0
        self.avg_stdev = 0
        self.chi2, self.chi2c = 0, 0
        self.cal_dic, self.netadj_drift_dic = {}, {}
        self.text = []

    def __str__(self):
        return_str = ''
        for attr in vars(self):
            if attr != 'text':
                return_str += '{}: {}\n'.format(attr, getattr(self, attr))
        for line in self.text:
            return_str += line
        return return_str


###############################################################################
class AdjustedStation:
    """
    Object to store the network-adjusted gravity value at each station.
    """

    def __init__(self, station, g, sd):
        self.station = station
        self.g = g
        self.sd = sd

    def __str__(self):
        return '{} {:0.2f} {:0.2f}'.format(self.station, self.g, self.sd)


###############################################################################
# other functions and classes definitions
###############################################################################
class TimeSeries:
    """
    Time Series object used to store a simple time series, used in
    tide_correction.py could have been a ChannelList object
     Used for storing synthetic tides, or atmospheric time series for instance

    Properties:
    t:                      time vector: datetime format
    d:                      data

    Functions:
    - populateFromTsoftFile
    - populateFromEternaFile
    - interpolateOnGivenTimes
    """

    def __init__(self):
        """
        """
        self.t = []
        self.d = []

    def populate_from_tsoft_channel(self, filename, channel):
        """
        Load a .tsf file and populate the object's properties
        Assume data is in nm/s² and convert it to mgal!
        """
        i = 0
        try:
            # essaye d'ouvrir le fichier
            with open(filename, 'r') as fh:
                test = 0
                for line in fh:
                    i += 1
                    # Clean line
                    line = line.strip()
                    # Skip blank and comment lines
                    if (not line) or (line[0] == '/') or (line[0] == 'L'):
                        continue
                    # parse string line first with respect to '/' caracters (used in the date format),
                    # then with ':' (used for the time display), eventually with the classic ' '
                    vals = line.split()
                    if test == 1:
                        self.t.append(
                            dt.datetime(
                                int(vals[0]),
                                int(vals[1]),
                                int(vals[2]),
                                int(vals[3]),
                                int(vals[4]),
                                int(vals[5]),
                            )
                        )
                        # read nm/s² data and convert it to mgal:
                        self.d.append(float(vals[5 + channel]) / 10000.0)
                    if vals[0] == "[DATA]":
                        test = 1
        except IOError:
            # si ça ne marche pas, affiche ce message et continue le prog
            print('No file : {}'.format(filename))
        except ValueError:
            print('pb at line {:d} : Is it really .tsf (or .TSF) format? '.format(i))
        except IndexError:
            print(
                'pb at line {:d} : check raw data file: possibly last line?'.format(i)
            )

    def populate_from_eterna_file(self, filename, channel):
        """
        load an eterna file and populate the object's properties
        Assume data is in nm/s² and convert it to mgal!
        """
        i = 0
        try:
            with open(filename, 'r') as fh:
                test = 0
                for line in fh:
                    i += 1

                    # Clean line
                    line = line.strip()

                    # Skip blank and comment lines
                    if (not line) or (line[0] == '/') or (line[0] == 'L'):
                        continue

                    # parse string line first with respect to '/' caracters (used in the date format),
                    # then with ':' (used for the time display), eventually with the classic ' '
                    vals = line.split()
                    if vals[0] == "99999999":
                        test = 0

                    if test == 1:
                        # pad with 0
                        ttemp = vals[1].zfill(6)
                        self.t.append(
                            dt.datetime(
                                int(vals[0][0:4]),
                                int(vals[0][4:6]),
                                int(vals[0][6:8]),
                                int(ttemp[0:2]),
                                int(ttemp[2:4]),
                                int(ttemp[4:6]),
                            )
                        )
                        # read nm/s² data and convert it to mgal:
                        self.d.append(float(vals[1 + channel]) / 10000.0)
                    if vals[0] == "77777777":
                        test = 1
        except IOError:
            # si ça ne marche pas, affiche ce message et continue le prog
            print('No file : {}'.format(filename))
        except ValueError:
            print('pb at line {:d} : Is it really eterna format? '.format(i))
        except IndexError:
            print(
                'pb at line {:d} : check raw data file: possibly last line?'.format(i)
            )

    def interpolate_on_given_times(self, t):
        """
        Used in tide correction routines to interpolate the time series on the 
        user input time vector overlain the previous t and d fields 
        """
        tord = [date2num(tmp) for tmp in self.t]
        f = interp1d(tord, self.d, kind='linear', bounds_error=False)
        self.d = f([date2num(tmp) for tmp in t])
        self.t = t


class AdjustmentOptions:
    """
    Object that holds options for least-squares adjustment. Options are set 
    by gui_objects.AdjustOptions().
    """

    def __init__(self):
        # self.use_model_temp = False
        # self.model_temp_degree = 0
        self.use_sigma_prefactor = False
        self.sigma_prefactor = 1.0
        self.use_sigma_postfactor = False
        self.sigma_postfactor = 1.0
        self.use_sigma_add = False
        self.sigma_add = 0.0
        self.use_sigma_min = False
        self.sigma_min = 3.0
        self.alpha = 0.05
        self.cal_coeff = False
        self.adj_type = 'gravnet'
        self.specify_cal_coeff = False
        self.meter_cal_dict = None

    def __str__(self):
        return_str = ''
        for attr in vars(self):
            return_str += '{}: {}\n'.format(attr, getattr(self, attr))
        return return_str


###############################################################################
class Adjustment:
    """
    Object that holds least-squares adjustment input matrices and results. There
    is one Adjustment object per Survey().

    Properties:
    A                   model matrix (Nobs rel + NobsAbs)*Nunknowns
    P                   Weight matrix (Nobs rel + NobsAbs)*(Nobs rel + NobsAbs)
    Obs                 observation vector (Nobs rel + NobsAbs)
    S                   StX=0
    X                   Unknowns
    r                   residuals (V)
    var                 Diagonal of the a posteriori covariance matrix
    VtPV
    dof                 degree of freedom (Nobs rel + NobsAbs - Nunknowns)
    """

    def __init__(self):

        """
        Initializes

        """
        self.A = np.array([])
        self.P = np.array([])
        self.Obs = np.array([])
        self.S = np.array([])
        self.X = np.array([])
        self.r = np.array([])
        self.var = np.array([])
        self.VtPV = np.array([])
        self.SDaposteriori = 0
        self.dof = 0
        self.adjustmentoptions = AdjustmentOptions()
        self.adjustmentresults = AdjustmentResults()
        self.n_meters = 0
        self.meter_dic = dict()
        self.deltas = []
        self.datums = []
        self.g_dic = dict()
        self.sd_dic = dict()
        self.netadj_loop_keys = dict()
        self.sta_dic_ls = dict()

    def results_string(self):
        try:
            if self.adjustmentresults.text:  # Gravnet results
                return [x + '\n' for x in self.adjustmentresults.text]
            elif self.adjustmentresults.n_unknowns > 0:
                text_out = list()
                # text_out.append("Number of stations:                 {:d}\n".format(len(sta_dic_ls)))
                # text_out.append("Number of loops:                    {:d}\n".format(nloops))
                # text_out.append("Polynomial degree for time:         {:d}\n".format(drift_time))
                text_out.append(
                    "Number of unknowns:                  {:d}\n".format(
                        self.adjustmentresults.n_unknowns
                    )
                )
                text_out.append(
                    "Number of relative observations:     {:d}\n".format(
                        self.adjustmentresults.n_deltas
                    )
                )
                text_out.append(
                    "Number of absolute observations:     {:d}\n".format(
                        self.adjustmentresults.n_datums
                    )
                )
                text_out.append(
                    "Degrees of freedom (nobs-nunknowns): {:d}\n".format(self.dof)
                )
                text_out.append(
                    "SD a posteriori:                     {:4f}\n".format(
                        float(np.sqrt(self.VtPV / self.dof))
                    )
                )
                text_out.append(
                    "chi square value:                  {:6.2f}\n".format(
                        float(self.adjustmentresults.chi2)
                    )
                )
                text_out.append(
                    "critical chi square value:         {:6.2f}\n".format(
                        float(self.adjustmentresults.chi2c)
                    )
                )
                if float(self.adjustmentresults.chi2) < float(
                    self.adjustmentresults.chi2c
                ):
                    text_out.append("Chi-test accepted\n")
                else:
                    text_out.append("Chi-test rejected\n")

                text_out.append(
                    "Average standard deviation:          {:.2f}\n".format(
                        self.adjustmentresults.avg_stdev
                    )
                )
                text_out.append(
                    "Maximum delta residual:              {:.2f}\n".format(
                        self.adjustmentresults.max_dg_residual
                    )
                )
                text_out.append(
                    "Maximum absolute datum residual:     {:.2f}\n".format(
                        self.adjustmentresults.max_datum_residual
                    )
                )
                text_out.append(
                    "Minimum absolute datum residual:     {:.2f}\n".format(
                        self.adjustmentresults.min_datum_residual
                    )
                )
                if self.adjustmentoptions.cal_coeff:
                    for k, v in self.adjustmentresults.cal_dic.items():
                        text_out.append(
                            "Calibration coefficient for meter {}: {:.6f} +/- {:.6f}\n".format(
                                k, v[0], v[1]
                            )
                        )
                elif self.adjustmentoptions.specify_cal_coeff:
                    for k, v in self.adjustmentoptions.meter_cal_dict.items():
                        text_out.append(
                            "Specified calibration coefficient for meter {}: {:.6f}\n".format(
                                k, v
                            )
                        )
                if self.netadj_loop_keys:
                    text_out.append("Network adjustment drift coefficient(s):\n")
                    for (
                        loop
                    ) in (
                        self.netadj_loop_keys.items()
                    ):  # this dict only has loops with netadj option
                        text_out.append("Loop " + loop[0] + ": ")
                        for i in range(loop[1][1]):  # degree of polynomial
                            degree = self.X[
                                len(self.sta_dic_ls) + self.n_meters + loop[1][0] + i
                            ][0]
                            if degree == 1:
                                text_out.append(
                                    "Drift coefficient, degree {}: {:.3f}".format(
                                        i + 1, degree
                                    )
                                )
                            else:
                                text_out.append(
                                    "Drift coefficient, degree {}: {:.3f} (µGal/day)".format(
                                        i + 1, degree
                                    )
                                )
                return text_out
        except (AttributeError, TypeError):
            return ''

    @property
    def adjustment_gdic_str(self):
        return_str = ''
        for k, v in self.g_dic.items():
            return_str += '{} {:0.2f} {:0.2f}\n'.format(k, v, self.sd_dic[k])
        return return_str

    def python_lsq_inversion(self):
        """
        Solve system of equations using numpy linalg.inv. Similar to LSQ
        inversion from Hwang et al (2002)
        """
        At = np.transpose(self.A)
        # St = np.transpose(self.S)
        N = At.dot(self.P).dot(self.A)
        # original PyGrav solution:
        # self.X = np.linalg.inv(N+self.S.dot(St)).dot(At).dot(self.P).dot(self.Obs)
        self.X = np.linalg.inv(N).dot(At).dot(self.P).dot(self.Obs)
        self.r = self.A.dot(self.X) - self.Obs
        rt = np.transpose(self.r)
        self.VtPV = rt.dot(self.P).dot(self.r)
        var_post_norm = self.VtPV / self.dof
        self.SDaposteriori = np.sqrt(var_post_norm)
        cov_post = np.linalg.inv(N) * var_post_norm
        self.var = np.diagonal(cov_post)

    def lsq_statistics(self):
        """
        a priori variance of unit weight = 1
        """
        alpha = self.adjustmentoptions.alpha
        self.adjustmentresults.chi2 = self.VtPV

        t = np.sqrt(2 * np.log(1 / alpha))
        chi_1_alpha = t - (2.515517 + 0.802853 * t + 0.010328 * t ** 2) / (
            1 + 1.432788 * t + 0.189269 * t ** 2 + 0.001308 * t ** 3
        )
        dof = float(self.dof)
        self.adjustmentresults.chi2c = (
            dof * (chi_1_alpha * np.sqrt(2 / (9 * dof)) + 1 - 2.0 / (9 * dof)) ** 3
        )
