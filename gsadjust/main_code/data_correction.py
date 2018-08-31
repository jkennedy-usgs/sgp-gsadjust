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


def tide_correction(self):
    """
    - set the Window for tide Correction options
    - link selection options to appropriate tide correction functions
    Three correction functions are availabel:
        - use_meter_tide_correction: uses the Burris/Scintrex internal tide correction. Only valid
          if appropriate survey coordinates have been input in the meter prior to the survey
        - use_predict: calculate a synthetic tide based on tidal parameters and survey coordinates.
          It is currently not possible to apply different tide corrections to different stations or loops.
        - load_tide_time_series: load a synthetic tide from file to use as a correction.
    """
    tidecorrwin = QtGui.QWidget()
    self.statusBar().showMessage("Please choose tidal correction method")

    tidecorrwin.btn1 = QtGui.QPushButton('Use CG5 tide correction', self)
    tidecorrwin.btn1.clicked.connect(self.use_meter_tide_correction)
    tidecorrwin.btn2 = QtGui.QPushButton('Use synthetic tides from predict', self)
    tidecorrwin.btn2.clicked.connect(self.use_predict)
    tidecorrwin.btn3 = QtGui.QPushButton('Use synthetic tides from Agnew', self)
    tidecorrwin.btn3.clicked.connect(self.use_agnew)
    tidecorrwin.btn4 = QtGui.QPushButton('Load time series', self)
    tidecorrwin.btn4.clicked.connect(self.load_tide_time_series)

    # button locations
    grid = QtGui.QGridLayout()
    grid.addWidget(tidecorrwin.btn1, 0, 0, 1, 1)
    grid.addWidget(tidecorrwin.btn2, 1, 0, 1, 1)
    grid.addWidget(tidecorrwin.btn3, 2, 0, 1, 1)
    grid.addWidget(tidecorrwin.btn4, 3, 0, 1, 1)
    tidecorrwin.setLayout(grid)

    tidecorrwin.setWindowTitle('Tide correction method')
    tidecorrwin.setGeometry(50, 50, 350, 300)
    self.popup = tidecorrwin
    self.popup.show()


def use_meter_tide_correction(self):
    """
    Function for using internal CG5 tide correction option.
    update the Campaign().corr_g list
    """
    self.campaigndata.corr_g = self.campaigndata.grav
    Tides = TimeSeries()
    Tides.d = deepcopy(self.campaigndata.etc)
    # should be negative to match the further application of applyTideCorrection(Tides): native CG5 etc value is not
    # the synthetic tide that is subtracted but already the correction (hence - the synt tide....)
    Tides.d = [-d for d in Tides.d]
    Tides.t = deepcopy(self.campaigndata.t)
    self.applyTideCorrection(Tides)


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


def use_predict(self):
    """
    Generate synthetic tides from a list of gravimetric factors for tide corrections
    """
    enterCoordinates = QtGui.QWidget()
    self.statusBar().showMessage("Please enter survey coordinates")

    enterCoordinates.lat = QtGui.QLabel('latitude')
    enterCoordinates.lon = QtGui.QLabel('longitude')
    enterCoordinates.alt = QtGui.QLabel('altitude')
    enterCoordinates.latEdit = QtGui.QLineEdit()
    enterCoordinates.lonEdit = QtGui.QLineEdit()
    enterCoordinates.altEdit = QtGui.QLineEdit()
    # create buttons and actions
    enterCoordinates.btn1 = QtGui.QPushButton('Use default tidal parameters', self)
    enterCoordinates.btn1.clicked.connect(lambda: self.launchPredict(enterCoordinates, option=1))
    enterCoordinates.btn2 = QtGui.QPushButton('Load tidal parameters file', self)
    enterCoordinates.btn2.clicked.connect(lambda: self.launchPredict(enterCoordinates, option=2))

    # locations
    grid = QtGui.QGridLayout()
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

    IMPORTANT the format of the tide groups definition in the .ini is very important; specifications can be found in
    the Eterna doc if problems are encountered here.
    """
    self.lat = float(cwin.latEdit.text())
    self.lon = float(cwin.lonEdit.text())
    self.alt = float(cwin.altEdit.text())

    t = self.campaigndata.t
    dur = t[len(t) - 1] - t[0]  # dur is a timedelta object
    dur = int(dur.days * 24 + dur.seconds / 60 / 60 + 24)  #
    projname = "PR000000.ini"
    predict_file = open(self.output_root_dir + os.sep + projname, 'w')
    predict_file.write("TEXTHEADER= MAREES THEORIQUES\n")
    predict_file.write("TEXTHEADER= %d (1 minute)\n" % t[0].year)
    predict_file.write("TEXTHEADER= FILTRE GGP2 (DLAG = 17.18 sec)\n\n")
    predict_file.write("SENSORNAME=   CG5         #earth tide sensor name\n")
    predict_file.write("SAMPLERATE=     60         #sampling interval in seconds \n")
    predict_file.write("STATLATITU=   %f       #stations latitude  in degree\n" % self.lat)
    predict_file.write("STATLONITU=    %f       #stations longitude in degree\n" % self.lon)
    predict_file.write("STATELEVAT=  %f         #stations elevation in meter\n" % self.alt)
    predict_file.write("STATGRAVIT=    0.0\n")
    predict_file.write("TIDALCOMPO=      0         #tidal component, see manual\n")
    predict_file.write("TIDALPOTEN=      7         #Earth Tide Potential\n")
    predict_file.write("AMTRUNCATE=1.D-10\n")
    predict_file.write("INITIALEPO=%d  %02d  %02d\n" % (t[0].year, t[0].month, t[0].day))
    predict_file.write("PREDICSPAN=%d          #time span in hours\n " % dur)
    predict_file.write("PRINTDEVEL=      1        #ANALYZE print param. for tidal development (1=yes)\n")
    predict_file.write("SEARDATLIM=     -1        #ANALYZE search for data error threshold \n")
    predict_file.write("NUMHIGPASS=      0        #ANALYZE highpass filtering = 1 \n")
    predict_file.write("PRINTOBSER=      0        #ANALYZE print parameter for observations (1=yes)\n")
    predict_file.write("RIGIDEARTH=      0        #ANALYZE parameter for rigid earth model (1=yes)\n")
    predict_file.write("HANNWINDOW=      0        #ANALYZE parameter for Hann-window (1=yes)\n")
    predict_file.write("QUICKLOOKA=      0        #ANALYZE parameter for quick look analysis (1=yes)\n")
    predict_file.write("POLETIDCOR=      0        #ANALYZE parameter for pole corrections (1=yes)\n")
    predict_file.write("LODTIDECOR=      0        #ANALYZE parameter for LOD corrections (1=yes)\n")
    predict_file.write("STORENEQSY=      1\n\n")

    if option == 1:
        # use default tidal parameters: load file
        tides = open("./eterna_files/200D.INI", 'r')
        for line in tides:
            predict_file.write(line)
        tides.close()

    if option == 2:
        # Ask the user to provide tide files
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Load tidal parameters file',
                                                  self.data_path)
        tides = open(fname, 'r')
        for line in tides:
            # Clean line
            line = line.strip()
            vals = line.split()
            predict_file.write("TIDALPARAM=")
            predict_file.write("%10.6f%10.6f%10.6f%10.6f %4s " % (float(vals[0]),
                                                          float(vals[1]), float(vals[2]), float(vals[3]), vals[4]))
            predict_file.write(" #ANALYZE wave group\n")
        tides.close()
    predict_file.close()
    predict_file2 = open(self.output_root_dir + os.sep + "project", 'w')
    predict_file2.write(projname)
    predict_file2.close()

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
    self.campaigndata.tideCorrectionPredict(tides)

    self.popup.close()
    self.statusBar().showMessage("Process some data?")


def use_agnew(self):
    """
    generate synthetic tides from the Agnew approach
    see
    - Agnew, D.C., 2007, 3.06 - Earth Tides, in Schubert, G. ed.,
    Treatise on Geophysics, Amsterdam, Elsevier, p. 163–195.
    - Agnew, D.C., 2012, SPOTL: Some Programs for Ocean-Tide Loading
    """

    enterCoordinates = QtGui.QWidget()

    self.statusBar().showMessage("Please enter survey coordinates")

    enterCoordinates.lat = QtGui.QLabel('Latitude')
    enterCoordinates.lon = QtGui.QLabel('Longitude')
    enterCoordinates.alt = QtGui.QLabel('Altitude')
    enterCoordinates.latEdit = QtGui.QLineEdit()
    enterCoordinates.lonEdit = QtGui.QLineEdit()
    enterCoordinates.altEdit = QtGui.QLineEdit()
    # create buttons and actions
    enterCoordinates.btn1 = QtGui.QPushButton('OK', self)
    enterCoordinates.btn1.clicked.connect(lambda: self.launchAgnew(enterCoordinates, option=1))

    # locations
    grid = QtGui.QGridLayout()
    grid.addWidget(enterCoordinates.lat, 1, 0)
    grid.addWidget(enterCoordinates.latEdit, 1, 1)
    grid.addWidget(enterCoordinates.lon, 2, 0)
    grid.addWidget(enterCoordinates.lonEdit, 2, 1)
    grid.addWidget(enterCoordinates.alt, 3, 0)
    grid.addWidget(enterCoordinates.altEdit, 3, 1)
    grid.addWidget(enterCoordinates.btn1, 4, 0)
    # grid.addWidget(enterCoordinates.btn2,4,1)
    enterCoordinates.setLayout(grid)
    enterCoordinates.setWindowTitle('Survey coordinates')
    enterCoordinates.show()
    self.popup = enterCoordinates
    self.popup.show()


def launch_agnew(self, cwin):
    """
    Launch the Agnew tide correction
    """
    self.lat = float(cwin.latEdit.text())
    self.lon = float(cwin.lonEdit.text())
    self.alt = float(cwin.altEdit.text())

    # apply tide correction
    self.campaigndata.tideCorrectionAgnew(self.lat, self.lon, self.alt)

    self.popup.close()
    self.statusBar().showMessage("Process some data?")


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
