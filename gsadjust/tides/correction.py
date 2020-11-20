"""
tides/correction.py
=================

Handles tide corrections for GSadjust.

Only Agnew tide correction is enabled. Ocean Loading correction is not implemented or
verified.

Actual tide computations are carried out in synthetic_tides.py.

Except for minor changes, this software is provided as written by B. Hector. It
should be considered preliminary or provisional and is subject to revision. It
is being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.
"""
import logging

import numpy as np
from PyQt5 import QtWidgets

from .synthetic import earth_tide, ocean_loading


def use_meter_tide_correction(MainProg):
    """
    Function for using internal CG5 tide correction option.
    update the Campaign().corr_g list
    """
    tide_correction_meter(MainProg.campaigndata)


def launch_agnew(cwin, MainProg):
    """
    Launch the Agnew tide correction
    """
    MainProg.campaigndata.tide_lat = float(cwin.latEdit.text())
    MainProg.campaigndata.tide_lon = float(cwin.lonEdit.text())
    MainProg.campaigndata.tide_alt = float(cwin.elevEdit.text())

    # apply tide correction
    MainProg.tide_popup.close()
    tide_correction_agnew(
        MainProg,
        MainProg.campaigndata.tide_lat,
        MainProg.campaigndata.tide_lon,
        MainProg.campaigndata.tide_alt,
    )
    MainProg.update_data_tab()


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
    logging.info("New tide correction, Lat: %f Long: %f Elevation: %f ", lat, lon, alt)
    for i in range(MainProg.obsTreeModel.invisibleRootItem().rowCount()):
        survey = MainProg.obsTreeModel.invisibleRootItem().child(i)
        for ii in range(survey.rowCount()):
            loop = survey.child(ii)
            for iii in range(loop.rowCount()):
                station = loop.child(iii)
                tides = (
                    np.round(
                        np.array([earth_tide(lat, lon, t) for t in station.t]) * 10000
                    )
                    / 10000.0
                )
                tides *= 1000  # convert milligal to microgal
                station.etc = tides.tolist()


def ocean_correction_agnew(self, amp, phases, lon):
    """
    compute and apply ocean loading correction to the dataset:

    NOT IMPLEMENTED.

    update the gravity field but not the etc field for avoiding possible later
    conflict with earth tide correction


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
            tides = (
                np.round(
                    np.array(
                        [
                            ocean_loading(t, amp, phases, lon)
                            for t in self.survey_dic[keysurv].t
                        ]
                    )
                    * 1000
                )
                / 1000.0
            )
            self.survey_dic[keysurv].corr_g = [
                self.survey_dic[keysurv].grav[i] + tides[i]
                for i in range(len(self.survey_dic[keysurv].t))
            ]
            self.survey_dic[keysurv].grav = self.survey_dic[keysurv].corr_g
        j = 1
        for keyloop, loop in self.survey_dic[keysurv].loop_dic.items():
            j += 1
            if any(self.survey_dic[keysurv].loop_dic[keyloop].t):
                # get tides and round to µgal level
                tides = (
                    np.round(
                        np.array(
                            [
                                ocean_loading(t, amp, phases, lon)
                                for t in self.survey_dic[keysurv].loop_dic[keyloop].t
                            ]
                        )
                        * 1000
                    )
                    / 1000.0
                )
                self.survey_dic[keysurv].loop_dic[keyloop].corr_g = [
                    self.survey_dic[keysurv].loop_dic[keyloop].grav[i] + tides[i]
                    for i in range(len(self.survey_dic[keysurv].loop_dic[keyloop].t))
                ]
                self.survey_dic[keysurv].loop_dic[keyloop].grav = (
                    self.survey_dic[keysurv].loop_dic[keyloop].corr_g
                )
            k = 1
            for keysta, sta in (
                self.survey_dic[keysurv].loop_dic[keyloop].station_dic.items()
            ):
                k += 1
                # get tides and round to µgal level
                tides = (
                    np.round(
                        np.array(
                            [
                                ocean_loading(t, amp, phases, lon)
                                for t in self.survey_dic[keysurv]
                                .loop_dic[keyloop]
                                .station_dic[keysta]
                                .t
                            ]
                        )
                        * 1000
                    )
                    / 1000.0
                )
                self.survey_dic[keysurv].loop_dic[keyloop].station_dic[
                    keysta
                ].corr_g = [
                    self.survey_dic[keysurv]
                    .loop_dic[keyloop]
                    .station_dic[keysta]
                    .grav[i]
                    + tides[i]
                    for i in range(
                        len(
                            self.survey_dic[keysurv]
                            .loop_dic[keyloop]
                            .station_dic[keysta]
                            .t
                        )
                    )
                ]
                self.survey_dic[keysurv].loop_dic[keyloop].station_dic[keysta].grav = (
                    self.survey_dic[keysurv]
                    .loop_dic[keyloop]
                    .station_dic[keysta]
                    .corr_g
                )


def ocean_loading_correction(self):
    """
    NOT IMPLEMENTED.

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
    fname, _ = QtWidgets.QFileDialog.getOpenFileName(
        self, "Open tidal parameters file (BLQ)", self.data_path
    )

    # read amplitudes and phases:
    fh = open(fname, "r")
    test = 0
    i = 0
    amp = []
    phases = []
    for line in fh:
        i += 1
        # Clean line
        line = line.strip()
        # Skip blank and comment lines
        if (not line) or (line == "$$ END TABLE"):
            continue
        vals = line.split()
        if len(vals) >= 3:
            if vals[2] == "GRAV":
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
