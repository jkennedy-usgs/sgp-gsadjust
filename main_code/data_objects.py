#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
data_objects.py
===============

GSadjust objects for non-Survey | Loop | Station objects (those are represented as PyQt objects).
--------------------------------------------------------------------------------------------------------------------

The main data objects are:

ChannelList: Holds lists of observed data values (g, tilt, temp, etc.). The lists are copied to ObsTreeStation objects.

Adjustment | AdjustmentOptions } Adjustment Results: The first contains instances of the latter two. Holds the input
  data, options, and results of the network adjustment.
 
Datum: Absolute-gravity observation or other reference for the relative-gravity network. At
  least one datum is required for network adjustment.
 
Delta: Relative-gravity difference calculated from two station occupations. May or may not include drift correction.

SimpleSurvey | SimpleLoop | SimpleStation: This is a simplified representation of the main ObsTreeSurvey | ObsTreeLoop |
  ObstreeStation PyQt datastore. The Simple... objects exist because PyQt objects can't be serialized, which is
  mandatory for saving the workspace using pickle (or json).

Tare: Represents an offset applied to the data

TimeSeries: Used only for tide correction.

This software is preliminary, provisional, and is subject to revision. It is being provided to meet the need for
timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related
material nor shall the fact of release constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or
unauthorized use of the software.
"""
import numpy as np
import datetime as dt
import pickle
import logging
from PyQt5 import QtCore, QtWidgets

from scipy.interpolate import interp1d
from collections import OrderedDict
from matplotlib.dates import date2num
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
        self.line = []
        self.station = []
        self.elev = []
        self.raw_grav = []
        self.sd = []
        self.tiltx = []
        self.tilty = []
        self.temp = []
        self.etc = []
        self.meter_etc = []
        self.dur = []
        self.rej = []
        self.t = []
        self.oper = []
        self.dial = []
        self.feedback = []
        self.tide = []
        self.meter = []
        self.lat = []
        self.long = []
        self.corr_g = []
        self.tare = []
        self.keepdata = []  # List of 0's and 1's indicating which samples are to be used

    def __iter__(self):
        """
        Function to iterate over ChannelList. Useful because Scintrex and Burris meters have different fields; this
        returns the fields relevant to the particular meter and skips blank fields.
        :return: tuple of a single non-empty field in the ChannelList object: (field, list of values)
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
            temp = value[start_idx:end_idx + 1]
            setattr(templist, field, temp)
        return templist

    def extract_subset_time(self, start_date, end_date):
        """
        Function for extracting shorter time series based on dates.
        Provide a starting and ending dates and the function will return a ChannelList object with a subset of
        the original.
        :param start_date: list index of first value to return.
        :param end_date: list index of last value to return.
        :return: ChannelList subset
        """
        templist = ChannelList()
        tsub = [t for t in self.t if date2num(start_date) <= date2num(t) <= date2num(end_date)]
        indexes = [self.t.index(t) for t in tsub]
        indexes.sort()
        # Burris and Scintrex have different fields. This only subsets non-empty fields.
        for field, value in self:
            temp = [value[ind] for ind in indexes]
            setattr(templist, field, temp)
        return templist


###############################################################################
class Delta:
    """
    Object to store relative-gravity difference
    
    Deltas are calculated when one of the "populate delta table" menu commands is chosen. They may or may not be
    drift-corrected. Before a network adjustment is done, there is no residual and the default value of -999 is 
    assigned.

    Types of deltas:
    'normal': calculated between two stations
    'three_point': calculated between one station, and an interpolated value betweeen two stations (Roman method).
      self.station2 is a tuple with two stations, self.station1 is a single station.
    'list': a delta is calculated from a list of deltas. Used with the Roman method to average three-point deltas.
      self.station2 is a list. self.station1 is None.

    ls_drift: records degree of least squares adjustment: they must be identical for all deltas to use Gravnet.
      Its a tuple with (loop name, degree of drift polynomial).
    """
    def __init__(self, sta1, sta2, driftcorr=0., sd=3., ls_drift=None, delta_type='normal', checked=2,
                 adj_stdev=999):
        self.station1 = sta1
        self.station2 = sta2
        self.stdev = float(sd)
        self.adj_stdev = adj_stdev
        self.checked = checked
        self.type = delta_type
        if type(sta2) is pyqt_models.ObsTreeStation:
            self.meter = sta2.meter[0]
        # Roman correction: list of deltas
        elif type(sta2) is list:
            self.meter = sta2[0].meter
        # Roman correction: single dg
        elif type(sta2) is tuple:
            self.meter = sta2[1].meter[0]
        self.residual = -999.
        self.driftcorr = driftcorr
        self.ls_drift = ls_drift  # degree of drift, if included in network adjustment

    @classmethod
    def from_list(cls, list_of_deltas):
        return cls(None, list_of_deltas, delta_type='list')

    # Used to match up adjustment results with observations, and to store simple Deltas
    @property
    def key(self):
        if self.type is 'list':
            return self.type + self.sta1 + self.sta2 + str(self.sta1_t)# self.sta1 + self.sta2 + str(self.dg) + str(self.sd)
        elif self.type is 'normal':
            return self.type + \
                   self.station1.station_name + \
                   self.station1.station_count + \
                   self.station2.station_name + \
                   self.station2.station_count

    @property
    def adj_sd(self):
        if self.adj_stdev == 999:
            return self.sd
        else:
            return self.adj_stdev

    @property
    def sd(self):
        if type(self.station2) is list:
            s = [np.abs(delta.dg) for delta in self.station2]
            return np.std(s)
        elif type(self.station2) is pyqt_models.ObsTreeStation:
            s = np.sqrt(self.station2.stdev**2 + self.station1.stdev**2)
            return s
        else:
            return self.stdev

    @property
    def sta1(self):
        if self.type == 'list':
            if self.station2[0].dg > 0:
                return self.station2[0].station1.station_name
            else:
                return self.station2[0].station2[1].station_name
        elif self.type == 'three_point':
            return self.station1.station_name
        else:
            return self.station1.station_name

    @property
    def sta2(self):
        # Roman correction
        if self.type == 'list':
            if self.station2[0].dg > 0:
                return self.station2[0].station2[1].station_name
            else:
                return self.station2[0].station1.station_name
        elif self.type == 'three_point':
            return self.station2[1].station_name
        else:
            return self.station2.station_name

    @property
    def dg(self):
        # Roman correction: dg requires three station-observations
        if self.type == 'three_point':
            sta2_dg = self.station2[1].gmean - self.station2[0].gmean
            time_prorate = (self.station1.tmean - self.station2[0].tmean)/(self.station2[1].tmean -
                                                                           self.station2[0].tmean)
            dg = (self.station2[0].gmean + sta2_dg * time_prorate) - self.station1.gmean
            return dg
        # Roman correction: return average of individual deltas
        elif self.type == 'list':
            dg = [np.abs(delta.dg) for delta in self.station2]
            return np.mean(dg)
        # Normal delta
        else:
            return self.station2.gmean - self.station1.gmean - self.driftcorr

    @property
    def sta1_t(self):
        if self.type == 'three_point':
            return self.station1.tmean
        elif self.type == 'list':
            return self.station2[0].station1.tmean
        else:
            return self.station1.tmean

    @property
    def sta2_t(self):
        if self.type == 'three_point':
            return self.station2[0].tmean
        elif self.type == 'list':
            return self.station2[0].station1.tmean
        else:
            return self.station2.tmean

    @property
    def duration(self):
        if type(self.station2) is tuple:
            return self.station2[1].tmean - self.station2[0].tmean
        elif type(self.station2) is list:
            dgs = [delta.dg for delta in self.station2]
            return np.mean(dgs)
        else:
            return self.sta2_t - self.sta1_t

    def __str__(self):
        if self.checked == 2:
            in_use = 'x'
        else:
            in_use = 'o'
        return ' '.join([in_use, self.sta1, self.sta2, str(self.sta1_t), str(self.sta2_t),
                         str(self.dg), str(self.sd)])

    def time(self):
        if type(self.station2) is pyqt_models.ObsTreeStation:
            return (self.sta1_t + self.sta2_t) / 2
        elif type(self.station2) is tuple:
            return self.sta1_t
        elif type(self.station2) is list:
            t = [delta.sta1_t for delta in self.station2]
            return np.mean(t)


###############################################################################
class Tare:
    """
    Object to store relative-gravity difference

    Deltas are calculated when one of the "populate delta table" menu commands is chosen. They may or may not be
    drift-corrected. Before a network adjustment is done, there is no residual and the default value of -999 is
    assigned.

    It's tempting to store just the Station objects that make up the delta, and calculate dg, sd, etc. as needed. But
    that doesn't work too well with the Roman method because the delta-g is calculated between (1) a station and (2) an
    interpolated value between two occupations at the other station.
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
        return ' '.join([in_use, self.date, self.time, str(self.tare)])

    @property
    def datetime(self):
        """
        :return: Python datetime object.
        """
        return QtCore.QDateTime(self.date, self.time).toPyDateTime()


###############################################################################
class Datum:
    """
    Object to store data typically associated with an absolute gravity measurement. 
    
    At least one datum observation is required for adjustment; the default datum value (50000) can be used at one of 
    the relative-gravity stations to see relative adjusted-gravity values. The gradient and measurement height are 
    considered when transferring between the absolute- and relative-gravimeter measuring heights. The date field is 
    not used in the adjustment (datum observations must be manually assigned to a survey). 
    """

    def __init__(self, station, g=50000, sd=5, date='1/1/2000', gradient=-3., meas_height=0, checked=2):
        self.station = station
        self.g = g
        self.sd = sd
        self.date = date
        self.gradient = gradient
        self.meas_height = meas_height
        self.residual = -999
        self.checked = checked

    def __str__(self):
        """
        Override the built-in method 'print' when applied to such object
        """
        if self.checked == 2:
            in_use = 'x '
        else:
            in_use = 'o '
        return in_use + self.station + ' ' + str(self.g) + ' ' + str(self.sd) + ' ' + self.date


###############################################################################
class AdjustmentResults:
    """
    Object to store least-squares adjustment statistics.
    """
    def __init__(self):
        self.n_deltas, self.n_deltas_notused = 0, 0
        self.n_datums, self.n_datums_notused = 0, 0
        self.max_dg_residual, self.min_dg_residual = 0, 0
        self.max_datum_residual, self.min_datum_residual = 0, 0
        self.avg_stdev = 0
        self.text = None


###############################################################################
class AdjustedStation:
    """
    Object to store the network-adjusted gravity value at each station.
    """

    def __init__(self, station, g, sd):
        self.station = station
        self.g = g
        self.sd = sd


###############################################################################
# other functions and classes definitions
###############################################################################
class TimeSeries:
    """Time Series object
    used to store a simple time series, used in tide_correction.py
    could have been a ChannelList object
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
                        self.t.append(dt.datetime(int(vals[0]), int(vals[1]),
                                                  int(vals[2]), int(vals[3]),
                                                  int(vals[4]), int(vals[5])))
                        # read nm/s² data and convert it to mgal:
                        self.d.append(float(vals[5 + channel]) / 10000.)
                    if vals[0] == "[DATA]":
                        test = 1
        except IOError:
            # si ça ne marche pas, affiche ce message et continue le prog
            print('No file : {}'.format(filename))
        except ValueError:
            print('pb at line {:d} : Is it really .tsf (or .TSF) format? '.format(i))
        except IndexError:
            print('pb at line {:d} : check raw data file: possibly last line?'.format(i))

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
                        self.t.append(dt.datetime(int(vals[0][0:4]), int(vals[0][4:6]),
                                                  int(vals[0][6:8]), int(ttemp[0:2]),
                                                  int(ttemp[2:4]), int(ttemp[4:6])))
                        # read nm/s² data and convert it to mgal:
                        self.d.append(float(vals[1 + channel]) / 10000.)
                    if vals[0] == "77777777":
                        test = 1
        except IOError:
            # si ça ne marche pas, affiche ce message et continue le prog
            print('No file : {}'.format(filename))
        except ValueError:
            print('pb at line {:d} : Is it really eterna format? '.format(i))
        except IndexError:
            print('pb at line {:d} : check raw data file: possibly last line?'.format(i))

    def interpolate_on_given_times(self, t):
        """
        Used in tide correction routines to interpolate the time series on the user input time vector overlain the 
        previous t and d fields 
        """
        tord = [date2num(tmp) for tmp in self.t]
        f = interp1d(tord, self.d, kind='linear', bounds_error=False)
        self.d = f([date2num(tmp) for tmp in t])
        self.t = t


class AdjustmentOptions:
    """
    Object that holds options for least-squares adjustment. Options are set by gui_objects.AdjustOptions().
    """

    def __init__(self):
        self.use_model_temp = False
        self.model_temp_degree = 0
        self.use_sigma_factor = False
        self.sigma_factor = 1.0
        self.use_sigma_add = False
        self.sigma_add = 0.0
        self.use_sigma_min = False
        self.sigma_min = 3.0
        self.alpha = 0.05
        self.woutfiles = False
        self.cal_coeff = False
        self.adj_type = 'gravnet'


###############################################################################
class Adjustment:
    """
    Object that holds least-squares adjustment input matrices and results. There is one Adjustment object per Survey().

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
        self.chi2 = 0
        self.chi2c = 0
        self.adjustmentoptions = AdjustmentOptions()
        self.adjustmentresults = AdjustmentResults()
        self.meter_dic = dict()
        self.deltas = []
        self.datums = []

    def python_lsq_inversion(self):
        """
        Solve system of equations using numpy linalg.inv. Similar to LSQ inversion from Hwang et al (2002)
        """
        At = np.transpose(self.A)
        St = np.transpose(self.S)
        N = At.dot(self.P).dot(self.A)

        # original solution:
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
        self.chi2 = self.VtPV

        t = np.sqrt(2 * np.log(1 / alpha))
        chi_1_alpha = t - (2.515517 + 0.802853 * t + 0.010328 * t ** 2) / (
            1 + 1.432788 * t + 0.189269 * t ** 2 + 0.001308 * t ** 3)
        dof = float(self.dof)
        self.chi2c = dof * (chi_1_alpha * np.sqrt(2 / (9 * dof)) + 1 - 2. / (9 * dof)) ** 3

        print("SD a posteriori: {:6.2f}".format(float(self.SDaposteriori)))
        print("chi square value: {:6.2f}".format(float(self.chi2)))
        print("critical chi square value: {:6.2f}".format(float(self.chi2c)))

        if self.chi2 < self.chi2c:
            print("Chi-test accepted")
        else:
            print("Chi-test rejected")


###############################################################################
# Simple.... objects have no PyQT elements, so they can be pickled
###############################################################################
class SimpleSurvey:
    """
    Qt objects can't be pickled. For the delta and datum tables, data are only stored in the respective tables,
    so we need to rewrite them as python objects. For the loop_dic and station_dics, they duplicate the data
    in the Qt models, so we can't just write the dics and ignore the Qt objects - they'll be recreated when the
    file loads.

    This also clears adjustment results, which we don't want to save.
    """

    def __init__(self, survey):
        self.loops = []
        self.deltas = []
        self.datums = []
        # Remove ObsTreeStation objects from deltas in the survey delta_model (which is different than the individual
        # loop delta_models; those are recreated when the workspace is loaded.
        for i in range(survey.delta_model.rowCount()):
            ind = survey.delta_model.createIndex(i, 0)
            delta = survey.delta_model.data(ind, QtCore.Qt.UserRole)
            simpledelta = SimpleDelta(delta)
            self.deltas.append(simpledelta)
        for i in range(survey.datum_model.rowCount()):
            ind = survey.datum_model.createIndex(i, 0)
            self.datums.append(survey.datum_model.data(ind, QtCore.Qt.UserRole))
        for i in range(survey.rowCount()):
            obstreeloop = survey.child(i)
            simpleloop = SimpleLoop(obstreeloop)
            self.loops.append(simpleloop)

        self.name = survey.name
        self.adjoptions = survey.adjustment.adjustmentoptions


class SimpleDelta:
    """
    Here we remove the ObsTreeStation objects from the Delta object (because they can't be pickled). Store a reference
    to the station name and station count so that when a workspace is loaded, the delta is recreated with the
    appropriate stations.
    """
    def __init__(self, delta):
        # Normal delta
        sta1 = None
        sta2 = None
        self.key = delta.key
        self.adj_sd = delta.adj_sd
        self.type = delta.type
        self.ls_drift = delta.ls_drift
        self.driftcorr = delta.driftcorr
        self.checked = delta.checked

class SimpleLoop:
    """
    Qt objects can't be pickled.

    Loops have a tare_model PyQt object that must be removed. We don't need to worry about the delta_model because it
    is recreated when the workspace is loaded.
    """

    def __init__(self, loop):
        self.__dict__ = loop.__dict__
        self.stations = []
        self.tares = []
        if loop.tare_model is not None:
            for i in range(loop.tare_model.rowCount()):
                ind = loop.tare_model.createIndex(i, 0)
                self.tares.append(loop.tare_model.data(ind, QtCore.Qt.UserRole))
        for i in range(loop.rowCount()):
            station = SimpleStation(loop.child(i))
            self.stations.append(station)
        self.delta_model = None
        self.tare_model = None


class SimpleStation:
    def __init__(self, obstreestation):
        self.__dict__ = obstreestation.__dict__
        self.cellcolor = None
