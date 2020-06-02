#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
tide_corretion.py
=================

Module to calculate tide corrections for GSadjust, a graphical user interface for interactive network adjustmnent of
relative-gravity surveys.
---------------------------------------------------------------------------------------------------------------------

Actual tide computations are carried out in synthetic_tides.py.

Except for minor changes, this software is provided as written by B. Hector. It should be considered
preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science.
The software has not received final approval by the U.S. Geological Survey (USGS). No warranty, expressed or implied,
is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the
fact of release constitute any such warranty. The software is provided on the condition that neither the USGS nor the
U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the
software.
"""
import logging
import os
import shutil
import subprocess

import numpy as np
from PyQt5 import QtWidgets

from data_objects import TimeSeries
from synthetic_tides import earth_tide, ocean_loading


# TODO: Only Agnew is tested


def use_meter_tide_correction(MainProg):
    """
    Function for using internal CG5 tide correction option.
    update the Campaign().corr_g list
    """
    tide_correction_meter(MainProg.campaigndata)


def use_tide_time_series(MainProg):
    """
    Load time series, and apply the correction to all the data. time series should be either a .TSF (Tsoft) or an
    eterna formatted file It is assumed that the synthetic tide is stored in the first data column.
    """
    fname, _ = QtWidgets.QFileDialog.getOpenFileName(None, 'Open file',
                                                     'home')

    tides = TimeSeries()
    if fname[len(fname) - 4:len(fname)] == ".tsf" or fname[len(fname) - 4:len(fname)] == ".TSF":
        # tsoftfile
        tides.populate_from_tsoft_channel(fname, 1)
    else:
        # assume it's eterna: JPBoy loading calculations are usually on channel 3
        tides.populate_from_eterna_file(fname, 3)
    tide_correction_ts(MainProg.campaigndata, tides)


def use_predict_tides():
    """
    Generate synthetic tides from a list of gravimetric factors for tide corrections
    """

    enterCoordinates = QtWidgets.QWidget()
    self.statusBar().showMessage("Please enter survey coordinates")

    enterCoordinates.lat = QtWidgets.QLabel('latitude')
    enterCoordinates.lon = QtWidgets.QLabel('longitude')
    enterCoordinates.alt = QtWidgets.QLabel('altitude')
    enterCoordinates.latEdit = QtWidgets.QLineEdit()
    enterCoordinates.lonEdit = QtWidgets.QLineEdit()
    enterCoordinates.altEdit = QtWidgets.QLineEdit()
    # create buttons and actions
    enterCoordinates.btn1 = QtWidgets.QPushButton('Use default tidal parameters', self)
    enterCoordinates.btn1.clicked.connect(lambda: self.launchPredict(enterCoordinates, option=1))
    enterCoordinates.btn2 = QtWidgets.QPushButton('Load tidal parameters file', self)
    enterCoordinates.btn2.clicked.connect(lambda: self.launchPredict(enterCoordinates, option=2))

    # locations
    grid = QtWidgets.QGridLayout()
    grid.addWidget(enterCoordinates.lat, 1, 0)
    grid.addWidget(enterCoordinates.latEdit, 1, 1)
    grid.addWidget(enterCoordinates.lon, 2, 0)
    grid.addWidget(enterCoordinates.lonEdit, 2, 1)
    grid.addWidget(enterCoordinates.alt, 3, 0)
    grid.addWidget(enterCoordinates.altEdit, 3, 1)
    grid.addWidget(enterCoordinates.btn1, 4, 0)
    grid.addWidget(enterCoordinates.btn2, 4, 1)
    enterCoordinates.setLayout(grid)
    enterCoordinates.setWindowTitle('Survey coordinates')
    self.popup = enterCoordinates
    enterCoordinates.show()


def launch_predict(self, cwin, option):
    """
    Launch predict program: write a standard project file.
    IMPORTANT: the format of the tide groups definition in the
    .ini is very important; specifications can be found in the Eterna doc if problems are encountered here.
    """
    self.lat = float(cwin.latEdit.text())
    self.lon = float(cwin.lonEdit.text())
    self.alt = float(cwin.altEdit.text())

    t = self.campaigndata.t
    dur = t[len(t) - 1] - t[0]  # dur is a timedelta object
    dur = int(dur.days * 24 + dur.seconds / 60 / 60 + 24)  #
    projname = "PR000000.ini"
    predict_input_file = open(self.output_root_dir + os.sep + projname, 'w')
    predict_input_file.write("TEXTHEADER= MAREES THEORIQUES\n")
    predict_input_file.write("TEXTHEADER= %d (1 minute)\n" % t[0].year)
    predict_input_file.write("TEXTHEADER= FILTRE GGP2 (DLAG = 17.18 sec)\n\n")
    predict_input_file.write("SENSORNAME=   CG5         #earth tide sensor name\n")
    predict_input_file.write("SAMPLERATE=     60         #sampling interval in seconds \n")
    predict_input_file.write("STATLATITU=   %f       #stations latitude  in degree\n" % self.lat)
    predict_input_file.write("STATLONITU=    %f       #stations longitude in degree\n" % self.lon)
    predict_input_file.write("STATELEVAT=  %f         #stations elevation in meter\n" % self.alt)
    predict_input_file.write("STATGRAVIT=    0.0\n")
    predict_input_file.write("TIDALCOMPO=      0         #tidal component, see manual\n")
    predict_input_file.write("TIDALPOTEN=      7         #Earth Tide Potential\n")
    predict_input_file.write("AMTRUNCATE=1.D-10\n")
    predict_input_file.write("INITIALEPO=%d  %02d  %02d\n" % (t[0].year, t[0].month, t[0].day))
    predict_input_file.write("PREDICSPAN=%d          #time span in hours\n " % dur)
    predict_input_file.write("PRINTDEVEL=      1        #ANALYZE print param. for tidal development (1=yes)\n")
    predict_input_file.write("SEARDATLIM=     -1        #ANALYZE search for data error threshold \n")
    predict_input_file.write("NUMHIGPASS=      0        #ANALYZE highpass filtering = 1 \n")
    predict_input_file.write("PRINTOBSER=      0        #ANALYZE print parameter for observations (1=yes)\n")
    predict_input_file.write("RIGIDEARTH=      0        #ANALYZE parameter for rigid earth model (1=yes)\n")
    predict_input_file.write("HANNWINDOW=      0        #ANALYZE parameter for Hann-window (1=yes)\n")
    predict_input_file.write("QUICKLOOKA=      0        #ANALYZE parameter for quick look analysis (1=yes)\n")
    predict_input_file.write("POLETIDCOR=      0        #ANALYZE parameter for pole corrections (1=yes)\n")
    predict_input_file.write("LODTIDECOR=      0        #ANALYZE parameter for LOD corrections (1=yes)\n")
    predict_input_file.write("STORENEQSY=      1\n\n")

    if option == 1:
        # use default tidal parameters: load predict_input_file
        tides = open("./eterna_files/200D.INI", 'r')
        for line in tides:
            predict_input_file.write(line)

        tides.close()
    if option == 2:
        # Ask the user to provide tide files
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Load tidal parameters file',
                                                         self.settings.value('tide_path'))

        tides = open(fname, 'r')
        for line in tides:
            # Clean line
            line = line.strip()
            vals = line.split()
            predict_input_file.write("TIDALPARAM=")
            predict_input_file.write("%10.6f%10.6f%10.6f%10.6f %4s " % (float(vals[0]),
                                                                        float(vals[1]),
                                                                        float(vals[2]),
                                                                        float(vals[3]),
                                                                        vals[4]))
            predict_input_file.write(" #ANALYZE wave group\n")
        tides.close()
    predict_input_file.close()
    file2 = open(self.output_root_dir + os.sep + "project", 'w')
    file2.write(projname)
    file2.close()

    # cp predict.exe to output directory
    shutil.copyfile("./eterna_files/predict.exe", self.output_root_dir + os.sep + "predict.exe")

    # run predict
    curr_dir = os.getcwd()
    os.chdir(self.output_root_dir)
    subprocess.call(["predict.exe"])
    os.chdir(curr_dir)

    # load results and apply tide corrections
    tides = TimeSeries()
    tides.populate_from_eterna_file(self.output_root_dir + os.sep + "PR000000.prd", 1)

    self.campaigndata.tide_correction_predict(tides)

    self.popup.close()
    self.statusBar().showMessage("Process some data?")


def launch_agnew(cwin, MainProg):
    """
    Launch the Agnew tide correction
    """
    MainProg.campaigndata.tide_lat = float(cwin.latEdit.text())
    MainProg.campaigndata.tide_lon = float(cwin.lonEdit.text())
    MainProg.campaigndata.tide_alt = float(cwin.elevEdit.text())

    # apply tide correction
    MainProg.tide_popup.close()
    tide_correction_agnew(MainProg,
                          MainProg.campaigndata.tide_lat,
                          MainProg.campaigndata.tide_lon,
                          MainProg.campaigndata.tide_alt)
    MainProg.update_data_tab()


def tide_correction_ts(campaigndata, tide_in):
    """
    Apply tide correction from previously-loaded Tsoft file. Assume that units are mGal (already converted from nanogal
    when the file was loaded). The sign of the correction is opposite that of Agnew.
    :param campaigndata: data
    :param tides: [t,s]
    """
    logging.info('New tide correction from time series')
    for keysurv, surv in campaigndata.survey_dic.items():
        for keyloop, loop in campaigndata.survey_dic[keysurv].loop_dic.items():
            for keysta, sta in campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic.items():
                tides = np.round(np.array([ts_tide_corr(tide_in, t)
                                           for t in campaigndata.survey_dic[keysurv].
                                          loop_dic[keyloop].
                                          station_dic[keysta].t]) * 10000) / 10000.
                tides *= 1000  # convert milligal to microgal
                campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].grav = \
                    [campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].grav[i]
                     - campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].etc[i]  # remove previous
                     - tides[i]  # add new correction
                     for i in range(len(campaigndata.survey_dic[keysurv].
                                        loop_dic[keyloop].
                                        station_dic[keysta].t))]
                campaigndata.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].etc = \
                    [tides[i] * -1 for i in range(len(campaigndata.survey_dic[keysurv].
                                                      loop_dic[keyloop].
                                                      station_dic[keysta].t))]


def ts_tide_corr(tides, time):
    idx = tides.t.index(min(tides.t, key=lambda x: abs(x - time)))
    return tides.d[idx]


def tide_correction_meter(MainProg):
    """
    Function for using internal CG5 tide correction option.
    update the Campaign().corr_g list
    """
    for i in range(MainProg.obsTreeModel.invisibleRootItem().rowCount()):
        survey = MainProg.obsTreeModel.invisibleRootItem().child(i)
        for ii in range(survey.rowCount()):
            loop = survey.child(ii)
            for iii in range(loop.rowCount()):
                station = loop.child(iii)
                station.etc = station.meter_etc


def tide_correction_agnew(MainProg, lat, lon, alt):
    """
    Apply Agnew tide correction.

    - Agnew, D.C., 2007, 3.06 - Earth Tides, in Schubert, G. ed.,
    Treatise on Geophysics, Amsterdam, Elsevier, p. 163–195.
    - Agnew, D.C., 2012, SPOTL: Some Programs for Ocean-Tide Loading

    """
    logging.info('New tide correction, Lat: %f Long: %f Elevation: %f ' % (lat, lon, alt))
    for i in range(MainProg.obsTreeModel.invisibleRootItem().rowCount()):
        survey = MainProg.obsTreeModel.invisibleRootItem().child(i)
        for ii in range(survey.rowCount()):
            loop = survey.child(ii)
            for iii in range(loop.rowCount()):
                station = loop.child(iii)
                tides = np.round(np.array([earth_tide(lat, lon, t)
                                           for t in station.t]) * 10000) / 10000.
                tides *= 1000  # convert milligal to microgal
                station.etc = tides


def ocean_correction_agnew(self, amp, phases, lon):
    """
    compute and apply ocean loading correction to the dataset:

    update the gravity field but not the etc field for avoiding possible later conflict with earth tide correction


    Arguments:
    amplitudes & phases: lists with amplitudes and phases. List order is:
    M2,S2,K1,O1,N2,P1,K2,Q1,Mf,Mm,Ssa
    (amp,phases)
    lon:     site longitude
    """

    # if any(self.t):
    #     # get tides and round to µgal level
    #     tides = np.round(np.array([ocean_loading(t, amp, phases, lon) for t in self.t]) * 1000) / 1000.
    #     g = self.grav()
    #     self.corr_g = [g[i] + tides[i] for i in range(len(self.t))]
        # self.grav = self.corr_g
    i = 1
    for keysurv, surv in self.survey_dic.items():
        i += 1
        if any(self.survey_dic[keysurv].t):
            # get tides and round to µgal level
            tides = np.round(np.array([ocean_loading(t, amp, phases, lon)
                                       for t in self.survey_dic[keysurv].t]) * 1000) / 1000.
            self.survey_dic[keysurv].corr_g = [self.survey_dic[keysurv].grav[i] + tides[i]
                                               for i in range(len(self.survey_dic[keysurv].t))]
            self.survey_dic[keysurv].grav = self.survey_dic[keysurv].corr_g
        j = 1
        for keyloop, loop in self.survey_dic[keysurv].loop_dic.items():
            j += 1
            if any(self.survey_dic[keysurv].loop_dic[keyloop].t):
                # get tides and round to µgal level
                tides = np.round(np.array([ocean_loading(t, amp, phases, lon)
                                           for t in self.survey_dic[keysurv].loop_dic[keyloop].t]) * 1000) / 1000.
                self.survey_dic[keysurv].loop_dic[keyloop].corr_g = \
                    [self.survey_dic[keysurv].loop_dic[keyloop].grav[i] + tides[i]
                     for i in range(len(self.survey_dic[keysurv].loop_dic[keyloop].t))]
                self.survey_dic[keysurv].loop_dic[keyloop].grav = \
                    self.survey_dic[keysurv].loop_dic[keyloop].corr_g
            k = 1
            for keysta, sta in self.survey_dic[keysurv].loop_dic[keyloop].station_dic.items():
                k += 1
                # get tides and round to µgal level
                tides = np.round(np.array([ocean_loading(t, amp, phases, lon)
                                           for t in self.survey_dic[keysurv].
                                          loop_dic[keyloop].
                                          station_dic[keysta].t]) * 1000) / 1000.
                self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].corr_g = \
                    [self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].grav[i] + tides[i]
                     for i in range(len(self.survey_dic[keysurv].
                                        loop_dic[keyloop].
                                        station_dic[keysta].t))]
                self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].grav = \
                    self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].corr_g


def ocean_loading_correction(self):
    """
    get the ocean Loading tidal parameters, compute the correction and
    apply to the dataset

    read amplitude and phases of ocean loading from standard BLQ file as
    obtained from http://holt.oso.chalmers.se/loading/ for a single station.

    amplitudes & phases: lists with amplitudes and phases. List order is:
    M2,S2,K1,O1,N2,P1,K2,Q1,Mf,Mm,Ssa
    (amp,phases)
    lon:     site longitude
    """

    # open file
    fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open tidal parameters file (BLQ)', self.data_path)

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
        if (not line) or (line == '$$ END TABLE'):
            continue
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

    self.campaigndata.ocean_correction_agnew(amp, phases, lon)
