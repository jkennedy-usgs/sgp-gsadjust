#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
data_correction.py
==================
"""
import os
import shutil
import subprocess
from copy import deepcopy
from PyQt5 import QtGui

from .data_objects import TimeSeries
from .synthetic_tides import *






def load_tide_time_series(self):
    """
    Load time series, and apply the correction to all the data. Time series should be either a .TSF (Tsoft) or an
    eterna formatted file. It is assumed that the synthetic tide is stored in the first data column.
    """
    self.statusBar().showMessage("Loading synthetic tide data")
    fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file',
                                              'home')

    Tides = TimeSeries()
    if fname[len(fname) - 4:len(fname)] == ".tsf" or fname[len(fname) - 4:len(fname)] == ".TSF":
        # tsoftfile
        Tides.populate_from_tsoft_channel(fname, 1)
    else:
        # assume it's eterna: JPBoy loading calculations are usually on channel 3
        Tides.populate_from_eterna_file(fname, 3)

    self.applyTideCorrection(Tides)




def ocean_loading_correction(self):
    """
    Get the ocean Loading tidal parameters, compute the correction, and apply to the dataset

    read amplitude and phases of ocean loading from standard BLQ file as obtained from
    http://holt.oso.chalmers.se/loading/ for a single station.

    amplitudes & phases: lists with amplitudes and phases. List order is:
    M2,S2,K1,O1,N2,P1,K2,Q1,Mf,Mm,Ssa
    (amp,phases)
    """

    # open file
    fname = QtGui.QFileDialog.getOpenFileName(self, 'Open tidal parameters file (BLQ)', self.data_path)

    # read amplitudes and phases:
    fh = open(fname, 'r')
    test = 0
    i = 0
    amp = []
    phases = []
    for line in fh:
        i += 1
        # Clean line
        line = line.strip()
        # Skip blank and comment lines
        if (not line) or (line == '$$ END TABLE'): continue
        vals = line.split()
        if len(vals) >= 3:
            if vals[2] == 'GRAV':
                test = 1
                lon = float(vals[5])
                i = 0
        if test == 1 and i == 1:
            for j in range(11):
                amp.append(float(vals[j]))
        if test == 1 and i == 4:
            for j in range(11):
                phases.append(float(vals[j]))
            break

    self.campaigndata.oceanCorrectionAgnew(amp, phases, lon)


def atmospheric_correction(self):
    """
    Load a time series and apply the correction

    Time series should be either a .TSF (Tsoft) or an eterna formatted file. It is assumed that the synthetic tide is
    stored in the first data column. This function populates/updats the corrg field of all gravity objects within
    self.campaigndata and same for grav field.

    IMPORTANT: atmospheric correction is currently not stored in the data structure,
    nor saved in an output file, this could be updated
    """
    self.statusBar().showMessage("Loading Atmospheric loading time series")
    fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file',
                                              'home')

    Atm = TimeSeries()
    if fname[len(fname) - 4:len(fname)] == ".tsf" or fname[len(fname) - 4:len(fname)] == ".TSF":
        # tsoftfile
        Atm.populate_from_tsoft_channel(fname, 1)
    else:
        # assume it's eterna
        Atm.populate_from_eterna_file(fname, 1)

    # if data is loaded from already processed data, self.campaigndata time
    # series are not populated, and correction should not be applied at such
    # levels
    if self.campaigndata.t:
        Atm.interpolate_on_given_times(self.campaigndata.t)
        print(len(self.campaigndata.corr_g))
        print(len(Atm.d))

    self.campaigndata.corr_g = self.campaigndata.grav - Atm.d
    self.campaigndata.grav = self.campaigndata.corr_g

    i = 1
    for keysurv, surv in self.campaigndata.survey_dic.iteritems():
        i += 1
        Atmtemp = deepcopy(Atm)
        if self.campaigndata.survey_dic[keysurv].t:
            Atmtemp.interpolate_on_given_times(self.campaigndata.survey_dic[keysurv].t)
            self.campaigndata.survey_dic[keysurv].corr_g = \
                self.campaigndata.survey_dic[keysurv].grav - Atmtemp.d
            self.campaigndata.survey_dic[keysurv].grav = self.campaigndata.survey_dic[keysurv].corr_g
        j = 1
        for keyloop, loop in self.campaigndata.survey_dic[keysurv].loop_dic.iteritems():
            j += 1
            Atmtemp2 = deepcopy(Atmtemp)
            if self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].t:
                Atmtemp2.interpolate_on_given_times(self.campaigndata.
                                                    survey_dic[keysurv].loop_dic[keyloop].t)
                self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].corr_g = \
                    self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].grav - Atmtemp2.d
                self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].grav = \
                    self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].corr_g
            k = 1
            for keysta, sta in self.campaigndata.survey_dic[keysurv]. \
                    loop_dic[keyloop].station_dic.iteritems():
                k += 1
                Atmtemp3 = deepcopy(Atmtemp2)
                Atmtemp3.interpolate_on_given_times(self.campaigndata.survey_dic[keysurv].
                                                    loop_dic[keyloop].station_dic[keysta].t)
                self.campaigndata.survey_dic[keysurv].loop_dic[keyloop]. \
                    station_dic[keysta].corr_g = self.campaigndata.survey_dic[keysurv]. \
                                                    loop_dic[keyloop].station_dic[keysta].grav - Atmtemp3.d
                self.campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta]. \
                    grav = self.campaigndata.survey_dic[keysurv].loop_dic[keyloop]. \
                    station_dic[keysta].corr_g
