#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""
pyqt_modules.py
===============

PyQt models for GSadjust. Handles assembling input matrices for network adjustment.
--------------------------------------------------------------------------------------------------------------------

NB: PyQt models follow the PyQt CamelCase naming convention. All other methods/functions in GSadjust use PEP-8
lowercase_underscore convention.

This software is preliminary, provisional, and is subject to revision. It is being provided to meet the need for
timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related
material nor shall the fact of release constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or
unauthorized use of the software.
"""
import copy
import datetime as dt
import logging
import os
import json, jsons
import numpy as np
from PyQt5 import QtGui, QtCore, QtWidgets
from PyQt5.QtCore import QVariant
from matplotlib.dates import num2date, date2num

import data_objects
from data_analysis import numpy_inversion

# Constants for column headers
STATION_NAME, STATION_DATETIME, STATION_MEAN = range(3)
LOOP_NAME = 0
SURVEY_NAME = 0

DATUM_STATION, DATUM_DATE, DATUM_G, DATUM_SD, N_SETS, MEAS_HEIGHT, GRADIENT, DATUM_RESIDUAL = range(8)
DELTA_STATION1, DELTA_STATION2, LOOP, DELTA_TIME, DELTA_G, DELTA_DRIFTCORR, DELTA_SD, DELTA_ADJ_SD, DELTA_RESIDUAL = range(
    9)
ADJSTA_STATION, ADJSTA_G, ADJSTA_SD = range(3)
SCINTREX_STATION, SCINTREX_DATE, SCINTREX_G, SCINTREX_SD, SCINTREX_X_TILT, SCINTREX_Y_TILT, SCINTREX_TEMP, SCINTREX_DUR, SCINTREX_REJ = range(
    9)
BURRIS_STATION, BURRIS_OPER, BURRIS_METER, BURRIS_DATE, BURRIS_G, BURRIS_DIAL, BURRIS_FEEDBACK, BURRIS_TIDE, BURRIS_ELEVATION, BURRIS_LAT, BURRIS_LONG = range(
    11)
TARE_DATE, TARE_TIME, TARE_TARE = range(3)
ROMAN_FROM, ROMAN_TO, ROMAN_DELTA = range(3)

def format_numeric_column(column):
    """
    Format fn for simple numeric columns.
    """
    return column + 1


class ObsTreeItem(QtGui.QStandardItem):
    """
    Basic tree-view item used to populate data tree, used for Surveys, Loops, and Stations.
    Not used directly but inherited by ObsTreeStation, ...Loop, and ...Survey
    """

    def __init__(self):
        super(ObsTreeItem, self).__init__()
        self.setFlags(self.flags() |
                      QtCore.Qt.ItemIsEnabled |
                      QtCore.Qt.ItemIsEditable |
                      QtCore.Qt.ItemIsUserCheckable)

        self.setCheckState(QtCore.Qt.Checked)
        self.fontweight = QtGui.QFont.Normal
        self.cellcolor = QtCore.Qt.white


class ObsTreeStation(ObsTreeItem):
    """
    PyQt model for stations. Station data is stored as a Station object in .station
    """

    def __init__(self, k, station_name, station_count):
        """
        Create a new station from ChannelList object. The fields (lists) of the ChannelList are copies to the
        ObsTreeStation.
        :param k: ChannelList object.
        :param station_name: String
        :param station_count: Int
        """
        super(ObsTreeStation, self).__init__()  # call properties from the baseclass
        self.__dict__ = copy.deepcopy(k.__dict__)
        self.station_name = station_name
        self.station_count = station_count
        if not hasattr(k, 'asd'):
            asd = None  # Records the number of times the station is occupied in a loop
        if hasattr(k, 'checked'):
            self.setCheckState(k.checked)
        # For legacy .p files
        # if len(self.corr_g) == 0:
        #     self.corr_g = self.raw_grav

    def __str__(self):
        return self.station_name

    def _weights_(self):
        """
        weights for Scintrex data are calculated based on the sample standard deviation divided by the duration.
        Neither of these are reported by the Burris meter.  Instead, duration is set to a constant (5) when the file
        is read. Parameter sd is calculated from the multiple samples comprising a station. Therefore w is a constant
        array if Burris data, variable if Scintrex.

        TODO: test effect of weighting if Burris and Scintrex data are combined in a survey
        :return:
        """
        if self.meter_type == 'Burris' or self.meter_type == 'CG6Tsoft':
            return [1 for i in self.keepdata if i == 1]
        else:
            sdtmp = [self.sd[i] / np.sqrt(self.dur[i]) for i in range(len(self.t))]
            w = np.array([1. / (sdtmp[i] * sdtmp[i]) for i in range(len(self.t)) if self.keepdata[i] == 1])
            return w

    @classmethod
    def from_station(cls, station):
        """
        Create a station from an existing station. Used to clone stations.
        :param station: ObsTreeStation object to be copied
        :return: ObsTreeStation
        """
        temp = cls(station, station.station_name, station.station_count)
        temp.__dict__ = station.__dict__
        return temp

    @property
    def key(self):
        return (self.station_name, self.tmean)

    @property
    def grav(self):
        """
        Applies tares and earth tide correction to raw_grav
        :return: List
        """
        # if len(self.corr_g) == 0:
        #     self.corr_g = self.raw_grav
        data = np.array(self.raw_grav) - np.array(self.tare) + np.array(self.etc)
        return data.tolist()

    def display_name(self):
        return self.station_name + ' (' + self.station_count + ')'

    def display_datetime_or_tmean(self):
        return num2date(self.tmean).strftime('%Y-%m-%d %H:%M:%S') if self.tmean != -999 else self.tmean

    def display_mean_stddev(self):
        return '{:.1f} ± {:.1f}'.format(self.gmean, self.stdev)

    @property
    def display_column_map(self):
        return {
            STATION_NAME: (self.display_name, ),
            STATION_DATETIME: (self.display_datetime_or_tmean, ),
            STATION_MEAN: (self.display_mean_stddev, )
        }

    @property
    def tooltip(self):
        return None

    @property
    def gmean(self):
        """
        The try-except block handles errors when all keepdata == 0.
        """
        try:
            gtmp = np.array([self.grav[i] for i in range(len(self.t)) if self.keepdata[i] == 1])
            gmean = sum(gtmp * self._weights_()) / sum(self._weights_())
            return gmean
        except:
            return -999

    @property
    def tmean(self):
        """
        The try-except block handles errors when all keepdata == 0.
        """
        try:
            ttmp = np.array([self.t[i] for i in range(len(self.t)) if self.keepdata[i] == 1])
            tmean = sum(ttmp * self._weights_()) / sum(self._weights_())
            return tmean
        except:
            return -999

    @property
    def t_stdev(self):
        try:
            if self.sd[0] == -999:
                return 1
            else:
                ttmp = np.array([self.t[i] for i in range(len(self.t)) if self.keepdata[i] == 1])
                num = np.zeros(len(ttmp))
                for i in range(len(ttmp)):
                    num[i] = self._weights_()[i] * (ttmp[i] - np.mean(ttmp)) ** 2
                sd = np.sqrt(
                    np.sum(num) / ((len(self._weights_()) - 1) * np.sum(self._weights_())))  # np.sqrt(1. / sum(
                # self._weights_()))
                return sd
        except:
            return -999

    @property
    def original_sd(self):
        sdtmp = np.array([self.sd[i] for i in range(len(self.t)) if self.keepdata[i] == 1])
        gtmp = np.array([self.grav[i] for i in range(len(self.t)) if self.keepdata[i] == 1])
        num = np.zeros(len(sdtmp))
        for i in range(len(sdtmp)):
            num[i] = self._weights_()[i] * (gtmp[i] - np.mean(gtmp)) ** 2
        sd = np.sqrt(np.sum(num) / ((len(self._weights_()) - 1) * np.sum(self._weights_())))
        return sd

    @property
    def stdev(self):
        """
        The try-except block handles errors when all keepdata == 0.

        The Scintrex meter reports an S.D. for each sample; the Burris meter does not report S.D. (with the Palm PDA,
        it shows S.D. on the display but does not record it).
        """
        try:
            if hasattr(self, 'asd'):
                if self.asd:
                    return self.asd
            if self.sd[0] == -999:  # Burris meter: return the S.D. calculated from all the samples at a station
                gtmp = np.array([self.grav[i] for i in range(len(self.t)) if self.keepdata[i] == 1])
                if len(gtmp) == 1:  # Can't take S.D. if there's only one sample
                    return 3.0
                else:
                    return float(np.std(gtmp))
            else:  # Scintrex: return the average SD of the samples
                sdtmp = np.array([self.sd[i] for i in range(len(self.t)) if self.keepdata[i] == 1])
                sd = np.mean(sdtmp)  #
                gtmp = np.array([self.grav[i] for i in range(len(self.t)) if self.keepdata[i] == 1])
                num = np.zeros(len(sdtmp))
                for i in range(len(sdtmp)):
                    num[i] = self._weights_()[i] * (gtmp[i] - np.mean(gtmp)) ** 2
                sd = np.sqrt(
                    np.sum(num) / ((len(self._weights_()) - 1) * np.sum(self._weights_())))  # np.sqrt(1. / sum(
                # self._weights_()))
                return sd
        except:
            return -999

    @property
    def summary_str(self):
        if self.tmean == -999:
            tm = self.tmean
        else:
            tm = num2date(self.tmean).strftime('%Y-%m-%d %H:%M:%S')
        summary_str = '{} {} {} {} {} {} {} {:0.3f} {:0.3f}\n'.format(self.display_name(),
                                                                      tm,
                                                                      self.oper[0],
                                                                      self.meter_type,
                                                                      self.lat[0],
                                                                      self.long[0],
                                                                      self.elev[0],
                                                                      self.gmean,
                                                                      self.stdev)
        return summary_str

    def iter_samples(self):
        """
        Iterator that returns print statements for each sample, used when writing adjustment summary
        :return: All of the stations in a campaign
        """
        for i in range(len(self.raw_grav)):
            if self.meter_type == 'Burris':
                return_str = '{} {} {:0.2f} {:0.2f} {:0.2f} {:0.2f}\n'.format(self.keepdata[i],
                                                                              self.station[i],
                                                                              self.raw_grav[i],
                                                                              self.grav[i],
                                                                              self.sd[i],
                                                                              self.etc[i])
            else:
                return_str = '{} {} {:0.2f} {} {}\n'.format(self.keepdata[i],
                                                            self.station[i],
                                                            self.raw_grav[i],
                                                            self.sd[i],
                                                            self.etc[i])
            yield return_str

    def to_json(self):
        jsonalizable_dict = self.__dict__
        jsonalizable_dict['checked'] = self.checkState()
        return self.__dict__


class ObsTreeLoop(ObsTreeItem):
    """
    PyQt model for Loops, parent item for stations. Loop attributes stored in .loop
    """

    def __init__(self, name):
        super(ObsTreeLoop, self).__init__()
        self.delta_model = DeltaTableModel()
        self.tare_model = TareTableModel()
        self.name = name
        self.drift_method = 'none'
        self.drift_cont_method = 0  # If continuous model, also need to keep track of which type of model
        self.drift_cont_startend = 0  # behavior at start/end. 0: extrapolate, 1: constant
        self.drift_netadj_method = 2  # If netadj method, keep track of polynomial degree
        self.meter = ''  # Meter S/N, for the case where multiple meters are calibrated
        self.comment = ''
        self.oper = ''

    def __str__(self):
        return 'Loop: {}, ' \
               'Drift method: {}, ' \
               'Meter type: {}, ' \
               'Continuous drift method: {}, ' \
               'Continuous drift start/end method: {}, ' \
               'Netadj method: {}\n'.format(self.name,
                                            self.drift_method,
                                            self.meter_type,
                                            self.drift_cont_method,
                                            self.drift_cont_startend,
                                            self.drift_netadj_method)

    @property
    # To accomodate legacy files, which might have meter and oper set to None:
    def tooltip(self):
        return 'Loop: {}\n' \
               'Drift method: {}\n' \
               'Meter: {}\n' \
               'Operator: {}\n' \
               'Comment: {}'.format(
            self.name,
            self.__dict__.get('drift_method', ''),
            self.__dict__.get('meter', ''),
            self.__dict__.get('oper', ''),
            self.__dict__.get('comment', ''),
        )

    @property
    def display_column_map(self):
        return {
            # Column name map
            # index: name
            LOOP_NAME: (str, self.name),
            1: (str,''),
            2: (str,'')
        }

    @property
    def meter_type(self):
        station = self.child(0)
        return station.meter_type

    @classmethod
    def from_json(cls, data):
        temp = cls(data['name'])
        temp.__dict__ = data
        temp.delta_model = DeltaTableModel()
        temp.tare_model = TareTableModel()
        for tare_dict in data['tares']:
            tare_object = data_objects.Tare(tare_dict['datetime'].date(), tare_dict['datetime'].time(),
                                            tare_dict['tare'])
            temp.tare_model.insertRows(tare_object, 0)
        if 'checked' in data:
            temp.setCheckState(data['checked'])
        return temp

    def rename(self, from_name, to_name):
        for i in range(self.rowCount()):
            item = self.child(i, 0)
            if item.station_name == from_name:
                item.station_name = to_name
                for iiii in range(len(item.station)):
                    item.station[iiii] = to_name

    def populate(self, data):
        """
        Populate loop dictionary using data passed as an option. For now, only the 'single' option works.

        :param data: Passed to Loop.populate_station_dic()
        """
        self.oper = data.oper[0]
        self.meter = data.meter[0]
        prev_sta = data.station[0]
        ind_start = 0
        station_count_dic = dict()

        # loop over samples:
        for i in range(len(data.station)):
            curr_sta = data.station[i]
            if curr_sta != prev_sta:
                if prev_sta not in station_count_dic.keys():
                    station_count_dic[prev_sta] = 1
                else:
                    station_count_dic[prev_sta] += 1
                ind_end = i - 1
                # create a new ChannelList object:
                temp_sta = data.extract_subset_idx(ind_start, ind_end)
                obstreestation = ObsTreeStation(temp_sta, prev_sta, str(station_count_dic[prev_sta]))
                obstreestation.meter_type = data.meter_type
                logging.info('Station added: ' + prev_sta + str(station_count_dic[prev_sta]))
                ind_start = i
                self.appendRow([obstreestation, QtGui.QStandardItem('0'), QtGui.QStandardItem('0')])
            if i == len(data.station) - 1:
                # enter last data line
                ind_end = i
                if curr_sta not in station_count_dic.keys():
                    station_count_dic[curr_sta] = 1
                else:
                    station_count_dic[curr_sta] += 1

                # create a new Loop-type object:
                temp_sta = data.extract_subset_idx(ind_start, ind_end)
                obstreestation = ObsTreeStation(temp_sta, curr_sta, str(station_count_dic[curr_sta]))
                obstreestation.meter_type = data.meter_type
                logging.info('Station added: ' + curr_sta + '(' + str(station_count_dic[curr_sta]) + ')')
                self.appendRow([obstreestation, QtGui.QStandardItem('0'), QtGui.QStandardItem('0')])
            prev_sta = curr_sta

    def get_station_reoccupation(self, station):
        """
        Get index of station reoccupation. To be used iteratively within the populate_station_dic function
        :param station: Station for which a repeat occupation is sought.
        :return: index of station reoccupation
        """
        station_reoc = [sta[0] for sta, stakey in self.station_dic.items() if sta[0] == station]
        return len(station_reoc)

    def get_data_for_plot(self):
        """
        Retrieves data from loop that is used for plotting: occupations must be checked, and there must be
        more than one occupation at a station.
        :return: list of plot_data lists
        """
        unique_stations = set()
        plot_data = []
        for i in range(self.rowCount()):
            station = self.child(i)
            if self.child(i).data(role=QtCore.Qt.CheckStateRole) == 2:
                if (station.station_name, station.meter[0]) not in unique_stations:
                    unique_stations.add(station.station_name)

        for unique_station in unique_stations:
            x, y, y_sd, t_sd = [], [], [], []
            for i in range(self.rowCount()):
                station = self.child(i)
                if self.child(i).data(role=QtCore.Qt.CheckStateRole) == 2:
                    if station.station_name == unique_station:
                        x.append(station.tmean)
                        y.append(station.gmean)
                        y_sd.append(station.original_sd)
                        t_sd.append(station.t_stdev)
            new_x, new_y, new_y_sd, new_t_sd = [], [], [], []
            # -999's can occur when all samples at a station are unchecked
            for idx, i in enumerate(x):
                if i != -999:
                    new_x.append(i)
                    new_y.append(y[idx])
                    new_y_sd.append(y_sd[idx])
                    new_t_sd.append(t_sd[idx])
            plot_data.append([new_x, new_y, unique_station, new_y_sd, new_t_sd])

        # sort plot_data by initial x in each line
        plot_data_sorted = []
        exes = [x[0][0] for x in plot_data]

        # get index of sorted list
        idx = sorted(range(len(exes)), key=lambda m: exes[m])
        for i in idx:
            plot_data_sorted.append(plot_data[i])
        return plot_data_sorted

    def checked_stations(self):
        """
        Retrieves checked data from loop; used for calculating delta-gs in the "None" and "Continuous" drift options.
        :return: list of checked stations
        """
        data = []
        for i in range(self.rowCount()):
            station = self.child(i)
            if self.child(i).data(role=QtCore.Qt.CheckStateRole) == 2:
                data.append(station)
        return data

    def deltas(self):
        """
        Called when 'Populate delta table... menu item is selected. Applies SD settings in options.
        """
        deltas = []
        table = self.delta_model
        for i in range(table.rowCount()):
            delta = table.data(table.index(i, 0), role=QtCore.Qt.UserRole)
            deltas.append(delta)
            if self.parent().adjustment.adjustmentoptions.use_sigma_min:
                delta.adj_sd = max(self.parent().adjustment.adjustmentoptions.sigma_min, delta.sd)
        return deltas

    def to_json(self):
        # Copy tares and stations from PyQt models to lists
        stations, tares = [], []
        if self.tare_model is not None:
            for i in range(self.tare_model.rowCount()):
                ind = self.tare_model.createIndex(i, 0)
                tares.append(self.tare_model.data(ind, QtCore.Qt.UserRole))
        for i in range(self.rowCount()):
            station = self.child(i).to_json()
            stations.append(station)

        self.delta_model = None
        self.tare_model = None
        return {
            'checked': self.checkState(),
            'delta_model': None,
            'tare_model': None,
            'stations': stations,
            'tares': tares,
            'name': self.name,
            'drift_method': self.drift_method,
            'drift_cont_method': self.drift_cont_method,
            # If continuous model, also need to keep track of which type of model
            'drift_cont_startend': self.drift_cont_startend,  # behavior at start/end. 0: extrapolate, 1: constant
            'drift_netadj_method': self.drift_netadj_method,  # If netadj method, keep track of polynomial degree
            'meter': self.meter,  # Meter S/N, for the case where multiple meters are calibrated
            'comment': self.comment,
            'oper': self.oper
        }


class ObsTreeSurvey(ObsTreeItem):
    """
    A survey object is a little different than Loop or Station objects in that it doesn't have a corresponding
    data object. All relevant info is in this object.
    """

    def __init__(self, name):
        super(ObsTreeSurvey, self).__init__()
        # Both surveys and loops have delta models. The survey delta_model holds the observations used in an
        # adjustment, and may have deltas from one or more loops.
        self.name = name
        self.delta_model = DeltaTableModel()
        self.datum_model = DatumTableModel()
        self.results_model = ResultsTableModel()
        self.adjustment = data_objects.Adjustment()

    def __str__(self):
        return self.name

    @property
    def display_column_map(self):
        return {
            # Column name map
            # index: name
            SURVEY_NAME: (str, self.name),
            1: (str, ''),
            2: (str, '')
        }

    @classmethod
    def from_json(cls, data):
        """
        When loading a workspace, repopulate PyQt models
        """
        temp = cls(data['name'])

        deltas = []
        for delta in data['deltas']:
            # Converting list of dicts to list of SimpleNamespace is to accommodate the loading routine,
            # which expects a delta-like object, not a dict. Could probably keep change it to "temp.deltas = data[
            # 'deltas'] if not for that.
            from types import SimpleNamespace
            sd = SimpleNamespace()
            sd.adj_sd = delta['adj_sd']
            sd.checked = delta['checked']
            sd.driftcorr = delta['driftcorr']
            sd.loop = delta['loop']
            sd.ls_drift = delta['ls_drift']
            sd.type = delta['type']
            try:
                sd.assigned_dg = delta['assigned_dg']
            except:
                pass
            try:
                sd.sta1 = delta['sta1']
                sd.sta2 = delta['sta2']
            except KeyError as e:
                # Raised if delta type is 'assigned'
                pass
            deltas.append(sd)

        temp.deltas = deltas

        for datum in data['datums']:
            d = data_objects.Datum(datum['station'])
            d.__dict__ = datum
            temp.datum_model.insertRows(d, 0)

        ao = data_objects.AdjustmentOptions()
        ao.__dict__ = data['adjoptions']
        temp.adjustment.adjustmentoptions = ao

        if 'checked' in data:
            temp.setCheckState(data['checked'])
        return temp

    def to_json(self):
        loops, datums = [], []
        # Remove ObsTreeStation objects from deltas in the survey delta_model (which is different than the individual
        # loop delta_models; those are ignored and not save (they're recreated when the workspace is loaded).
        for i in range(self.datum_model.rowCount()):
            ind = self.datum_model.createIndex(i, 0)
            datums.append(self.datum_model.data(ind, QtCore.Qt.UserRole))
        for i in range(self.rowCount()):
            loops.append(self.child(i).to_json())
        return {
            'loops': loops,
            'deltas': jsons.dump(self.delta_model),
            'datums': jsons.dump(datums),
            'checked': self.checkState(),
            'name': self.name,
            'adjoptions': jsons.dump(self.adjustment.adjustmentoptions)
        }

    @property
    def tooltip(self):
        return 'Survey: {}\n' \
               'Meters: {}\n' \
               'Number of loops: {}'.format(self.name,
                                            self.unique_meters,
                                            self.loop_count)

    @property
    def unique_meters(self):
        """
        :return: List of unique gravity-meter IDs (serial numbers) used in the survey
        """
        meters = []
        for station in self.iter_stations():
            meters.append(station.meter[0])  # Get the first entry; Assume meter number can't change at a station
        return list(set(meters))

    @property
    def loops(self):
        return [self.child(i) for i in range(self.rowCount())]

    def deltas(self):
        deltas = []
        for i in range(self.delta_model.rowCount()):
            ind = self.delta_model.createIndex(i, 0)
            delta = self.delta_model.data(ind, QtCore.Qt.UserRole)
            deltas.append(delta)
        return deltas

    @property
    def datums(self):
        datums = []
        for i in range(self.datum_model.rowCount()):
            ind = self.datum_model.createIndex(i, 0)
            datums.append(self.datum_model.data(ind, QtCore.Qt.UserRole))
        return datums

    def _handler_edit_ObsTreeSurvey(self, item, value):
        """
        Handle editing of ObsTreeLoop objects (renaming).
        """
        new_name = str(value)
        old_name = item.name
        logging.info('Survey renamed from {} to {}'.format(old_name, new_name))
        item.name = new_name
        return True

    @property
    def loop_count(self):
        a = self.loop_names
        return len(a)

    @property
    def loop_names(self):
        loops = []
        for i in range(self.rowCount()):
            loop = self.child(i)
            loops.append(loop.name)
        return loops

    def get_loop_by_name(self, name):
        for i in range(self.rowCount()):
            loop = self.child(i, 0)
            if loop.name == name:
                return loop
        raise KeyError

    def rename(self, from_name, to_name):
        for i in range(self.rowCount()):
            loop = self.child(i, 0)
            loop.rename(from_name, to_name)

    def populate(self, survey_data, name='0'):
        """
        Called from open_raw_data. Loads all survey_data into a single loop.
        :param survey_data: data returned from read_raw_data_file
        :param name: Survey name
        :return: None
        """
        obstreeloop = ObsTreeLoop(name)
        obstreeloop.populate(survey_data)
        self.appendRow([obstreeloop, QtGui.QStandardItem('0'), QtGui.QStandardItem('0')])
        logging.info('Survey added')

    def iter_stations(self):
        """
        Iterator that returns stations
        :return: All of the stations in a campaign
        """
        for i in range(self.rowCount()):
            obstreeloop = self.child(i)
            for ii in range(obstreeloop.rowCount()):
                obstreestation = obstreeloop.child(ii)
                yield obstreestation

    def run_inversion(self, adj_type='PyLSQ', write_out_files='n', output_root_dir='./', ):
        """
        Prepares inversion data and calls appropriate method.

        Adds an AdjustmentResults object to each Adjustment object of each Survey.
            -Counts n_deltas, n_deltas_notused, n_datums, n_datums_notused
            -Corrects datum observations for gradient
            -Runs inversion
        :param adj_type: 'PyLSQ' (numpy) or 'Gravnet'
        :param write_out_files: Write output files, 'y' or 'n'
        :param output_root_dir: Where to write output files
        """
        from gui_objects import show_message
        if self.data(role=QtCore.Qt.CheckStateRole) == 2:
            deltas = []
            datums = []
            adjustmentresults = data_objects.AdjustmentResults()

            # Form datasets for adjustment from selected table items, count number of deltas and datums (used when
            # writing metadata about adjustment).
            specify_cal_coeff = False
            cal_dic = None
            # To deal with old files
            if hasattr(self.adjustment.adjustmentoptions, 'specify_cal_coeff'):
                if self.adjustment.adjustmentoptions.specify_cal_coeff:
                    specify_cal_coeff = True
                    cal_dic = self.adjustment.adjustmentoptions.meter_cal_dict
            for ii in range(self.delta_model.rowCount()):
                ind = self.delta_model.createIndex(ii, 0)
                chk = self.delta_model.data(ind, QtCore.Qt.CheckStateRole)
                if chk == 2:  # Checkbox checked
                    delta = self.delta_model.data(ind, QtCore.Qt.UserRole)
                    if specify_cal_coeff:
                        delta.cal_coeff = cal_dic[delta.meter]
                    else:
                        delta.cal_coeff = 1
                    if delta.type == 'normal':
                        try:
                            # Don't include deltas in which one of the stations has been unchecked
                            if delta.station1.data(role=QtCore.Qt.CheckStateRole) == 2 and \
                                    delta.station2.data(role=QtCore.Qt.CheckStateRole) == 2:
                                deltas.append(delta)
                                adjustmentresults.n_deltas += 1
                            else:
                                adjustmentresults.n_deltas_notused += 1
                        except:
                            self.msg = show_message("Delta station not found. Was it deleted?", "GSadjust error")
                    else:
                        deltas.append(delta)
                        adjustmentresults.n_deltas += 1
                else:
                    adjustmentresults.n_deltas_notused += 1

            for ii in range(self.datum_model.rowCount()):
                ind = self.datum_model.createIndex(ii, 0)
                chk = self.datum_model.data(ind, QtCore.Qt.CheckStateRole)
                if chk == 2:  # Checkbox checked
                    # account for vertical gradient
                    datum = copy.copy(self.datum_model.data(ind, QtCore.Qt.UserRole))
                    if hasattr(datum, 'gradient'):
                        datum.g = datum.g - datum.gradient * datum.meas_height
                    datums.append(datum)
                    adjustmentresults.n_datums += 1
                else:
                    adjustmentresults.n_datums_notused += 1

            self.adjustment.deltas = deltas
            self.adjustment.datums = datums
            self.adjustment.adjustmentresults = adjustmentresults

            try:
                self.results_model.clearResults()
                if len(self.adjustment.datums) == 0:
                    self.msg = show_message(
                        "Survey {}: At least one datum must be specified".format(self.name),
                        "Inversion error")
                    return
                if len(self.adjustment.deltas) == 0:
                    self.msg = show_message(
                        "Survey {}: At least one relative-gravity difference must be specified".format(self.name),
                        "Inversion error")
                    return
                if adj_type == 'PyLSQ':
                    logging.info('Numpy inversion, Survey: {}'.format(self.name))
                    self.start_inversion(output_root_dir, write_out_files)
                elif adj_type == 'Gravnet':
                    logging.info('Gravnet inversion, Survey: {}'.format(self.name))
                    self.gravnet_inversion()
            except ZeroDivisionError:
                self.msg = show_message('Are there standard deviations that are zero or very small?',
                                        'ZeroDivisionError')
            except KeyError as e:
                self.msg = show_message('Key error: {}'.format(e.args[0]),
                                        'KeyError')
            except Exception as e:
                logging.exception(e, exc_info=True)
                self.msg = show_message("Inversion error", "unknown error")

    def gravnet_inversion(self):
        """
        Writes input files for Gravnet.exe, runs the executable, and reads in the results
        """
        # TODO: method could probably be made static
        from gui_objects import show_message
        from data_objects import AdjustedStation
        dir_changed1, dir_changed2 = False, False
        # Check that executable exists in current directory (it should, if run as compiled .exe
        if not os.path.exists('.\\gravnet.exe'):
            # if not found, it might be in the dist directory (i.e., running GSadjust.py)
            if os.path.exists('..\\dist\\gravnet.exe'):
                os.chdir('..\\dist')
                dir_changed1 = True
            if os.path.exists('.\\dist\\gravnet.exe'):
                os.chdir('.\\dist')
                dir_changed2 = True
            # not found at all: error
            else:
                self.msg = show_message('Gravnet.exe not found, aborting', 'Inversion error')
                logging.error('Inversion error, Gravnet.exe not found')
                return

        # If using gravnet with netadj drift option (including drift term in the network adjustment), a single
        # drift model is used for all observations in the adjustment. Check that the method is set to netadj
        # for all loops, otherwise show an error message and return.
        ls_degree = []
        drift_term = ''
        for delta in self.adjustment.deltas:
            if delta.ls_drift is not None:
                ls_degree.append(delta.ls_drift[1])
        unique_ls = list(set(ls_degree))
        if len(unique_ls) > 1:
            self.msg = show_message(
                'It appears that more than one polynomial degree was specified for different loops for the '
                'network, or that some loops are not using the ' +
                'adjustment drift option. When using Gravnet, all loops must have the same degree drift ' +
                'model. Aborting.',
                'Inversion error')
            return
        if len(unique_ls) == 1:
            if unique_ls[0] is not None:
                drift_term = '-T' + str(unique_ls[0])

        # Remove old gravnet files
        for ext in ['gra', 'err', 'his', 'met', 'res', 'sta']:
            try:
                os.remove('{}.{}'.format(self.name, ext))
            except OSError:
                pass

        # Warn if station names will be truncated
        truncate_warning = False
        for delta in self.adjustment.deltas:
            if len(delta.sta1) > 6 or len(delta.sta2) > 6:
                truncate_warning = True
        if truncate_warning:
            self.msg = show_message(
                'One or more station names is longer than 6 characters. Names will be truncated to 6 '
                'characters in the Gravnet input file. Please verify that names will still be unique '
                'after truncating.',
                'Inversion warning')
        # Write delta-g observation file
        dg_file = self.name + '_dg.obs'
        with open(dg_file, 'w') as fid:
            fid.write('start, end, difference (mgal), mjd (from), mjd (to),' +
                      ' reading (from, CU), reading (to, CU), standard deviation (mgal)\n')
            # Gravnet station names are limited to 6 characters; units are mGal
            for delta in self.adjustment.deltas:
                fid.write('{} {} {:0.6f} {} {} {:0.6f} {} {:0.6f}\n'.format(delta.sta1[:6],
                                                                            delta.sta2[:6],
                                                                            delta.dg / 1000. * delta.cal_coeff,
                                                                            delta.sta1_t,
                                                                            delta.sta2_t,
                                                                            delta.dg / 1000, '0',
                                                                            delta.sd_for_adjustment / 1000.))
        # Write absolute-g (aka datum, aka fix) file
        fix_file = self.name + '_fix.txt'
        with open(fix_file, 'w') as fid:
            for datum in self.adjustment.datums:
                fid.write(
                    '{} {:0.6f} {:0.3f}\n'.format(datum.station[:6], float(datum.g) / 1000., float(datum.sd) / 1000.))

        # Check if calibration coefficient is calculated
        cal_dic = {}
        if self.adjustment.adjustmentoptions.cal_coeff:
            if len(self.unique_meters) > 1:
                self.msg = show_message("It appears more than one meter was used on the survey. Gravnet calculates a " +
                                        "calibration coefficient for a single meter only. Use Numpy inversion " +
                                        "to calculate meter-specific calibration coefficients",
                                        "Inversion warning")
            if len(self.adjustment.datums) == 1:
                self.msg = show_message(
                    "Two or more datum observations are required to calculate a calibration coefficient." +
                    " Aborting.", "Inversion warning")
                return
            # Run gravnet with calibration coefficient
            os.system('gravnet -D' + dg_file + ' -N' + self.name + ' -M2 -C1 ' + drift_term + ' -F' + fix_file)
            with open(self.name + '.met', 'r') as fid:
                symbol = ''
                while symbol != 'Y_l':
                    line = fid.readline().split()
                    symbol = line[0]
                meter_calib_params = fid.readline().split()
                cal_dic[self.unique_meters[0]] = (1 + float(meter_calib_params[0]), float(meter_calib_params[1]))
        else:
            # Run gravnet without calibration coefficient
            os.system('gravnet -D' + dg_file + ' -N' + self.name + ' -M2 ' + drift_term + ' -F' + fix_file)

        # Read drift coefficients
        if drift_term != '':
            meter_drift_params = []
            with open(self.name + '.met', 'r') as fid:
                symbol = ''
                while symbol != 'Coefficients':  # Indicates 'Coefficients of drift' in .met file
                    line = fid.readline().split()
                    symbol = line[0]
                order = int(line[-1])
                for i in range(order):
                    meter_drift_params.append(fid.readline().split())

        self.results_model.setData(0, 0, QtCore.Qt.UserRole)  # clears results table

        # Read gravnet results
        g_dic, sd_dic = {}, {}
        with open(self.name + '.gra', 'r') as fid:
            _ = fid.readline()  # Header line
            all_sd = []
            while True:
                fh = fid.readline()
                if len(fh) == 0:
                    break
                parts = fh.split()
                sta = parts[1]
                g = float(parts[2]) * 1000
                sd = float(parts[3]) * 1000
                all_sd.append(sd)
                if len(parts) == 4:
                    g_dic[sta] = g
                    sd_dic[sta] = sd
                    self.results_model.insertRows(AdjustedStation(sta, g, sd), 0)
            self.adjustment.adjustmentresults.avg_stdev = np.mean(all_sd)

        # Match up residuals with input data
        self.adjustment.g_dic = g_dic
        self.adjustment.sd_dic = sd_dic
        self.match_inversion_results(inversion_type='gravnet', cal_dic=cal_dic)

        # Add calibration coefficient and/or drift coefficients to output text
        self.adjustment.adjustmentresults.text = []
        with open(self.name + '.sta', 'r') as fid:
            while True:
                fh = fid.readline()
                if len(fh) == 0:
                    break
                self.adjustment.adjustmentresults.text.append(fh)
            if self.adjustment.adjustmentoptions.cal_coeff:
                self.adjustment.adjustmentresults.text.append("\n\nGravimeter calibration coefficient: " +
                                                              str(1 + float(meter_calib_params[0])) + ' ± ' +
                                                              meter_calib_params[1])
            if drift_term != '':
                self.adjustment.adjustmentresults.text.append("\n\nGravimeter drift coefficient(s):\n ")
                for coeffs in meter_drift_params:
                    self.adjustment.adjustmentresults.text.append(coeffs[0] + ' ± ' + coeffs[1] + '\n')

        if dir_changed1:
            os.chdir('..\\gsadjust')
        elif dir_changed2:
            os.chdir('..')

    def start_inversion(self, output_root_dir, write_out_files='n'):
        """
        Least-squares network adjustment using numpy

        Parameters
        write_out_files:    (y/n): write output files for drift adjustment
                            (similar to MCGravi output files)
        output_root_dir     directory for writing output files
        """
        from gui_objects import show_message
        from data_objects import AdjustedStation
        self.adjustment.adjustmentresults.text = []

        # sta_dic_LS is a dictionary, key: station name, value: column for A matrix
        self.adjustment.sta_dic_ls = self.get_station_indices()
        self.adjustment.nloops = self.loop_count

        # number of unknowns:
        # TODO: temperature drift function
        drift_time = 0  # self.adjustment.adjustmentoptions.drift_time
        if self.adjustment.adjustmentoptions.use_model_temp:
            self.adjustment.drift_temp = self.adjustment.adjustmentoptions.model_temp_degree
        else:
            self.adjustment.drift_temp = 0

        if self.adjustment.adjustmentoptions.cal_coeff:
            if len(self.adjustment.datums) == 1:
                self.msg = show_message(
                    "Two or more datum observations are required to calculate a calibration coefficient." +
                    " Aborting.", "Inversion warning")
                return
            n_meters = len(self.unique_meters)
            self.adjustment.meter_dic = dict(zip(self.unique_meters, range(n_meters + 1)))
        else:
            n_meters = 0
        self.adjustment.n_meters = n_meters

        # Drift may be included in the network adjustment for any or all loops. This section of code
        # builds some dicts used later to assemble the A matrix.
        ndrift = 0

        # Temporary var to compile all delta.ls_drift tuples: (loop.name, degree of drift model)
        ls_drift_list = []

        # dict of tuples, used to identify column of drift observation in A matrix:
        # (loop.name, (column relative to end of A matrix, drift degree)
        netadj_loop_keys = dict()

        # dict of unique delta.ls_drift values, used to cross-reference loop name with poly. degree
        # (the loop object only knows if the drift method is netadj, or not. The delta object stores
        # the degree of the drift polynomial, but doesn't know if the drift correction is netadj (a
        # delta could have a non (0,0) ls_drift value, but use some other drift corretion method.
        loop_ls_dict = dict()
        for delta in self.adjustment.deltas:
            # If a workspace was loaded from json, need to convert list to tuple
            if type(delta.ls_drift) == list:
                delta.ls_drift = (delta.ls_drift[0], delta.ls_drift[1])
            ls_drift_list.append(delta.ls_drift)
        if not all(v is None for v in ls_drift_list):
            ls_drift_list = [x for x in ls_drift_list if x is not None]
            ls_drift_dict = dict(set(ls_drift_list))
            active_loops = set([x[0] for x in ls_drift_list])
            for lp in active_loops:
                obstreeloop = self.get_loop_by_name(lp)
                loop_ls_dict[obstreeloop.name] = obstreeloop.drift_method
                if obstreeloop.drift_method == 'netadj':
                    ls_degree = ls_drift_dict[obstreeloop.name]
                    netadj_loop_keys[obstreeloop.name] = (ndrift, ls_degree)
                    ndrift += ls_degree
        self.adjustment.loop_ls_dict = loop_ls_dict
        self.adjustment.ndrift = ndrift
        self.adjustment.netadj_loop_keys = netadj_loop_keys
        cal_dic = numpy_inversion(self.adjustment, self.results_model, output_root_dir, write_out_files='n')
        self.match_inversion_results('numpy', cal_dic)
        # self.results_model.

    def match_inversion_results(self, inversion_type, cal_dic=None):
        """
        Populates delta and datum table residuals from inversion results
        :param inversion_type: 'gravnet' or 'numpy'
        cal_dic: has an entry for each meter if a calibration coefficient was calculated. Different than
                 delta.cal_coeff, which is an a priori specified coefficient.
        """
        from gui_objects import show_message
        datum_residuals, dg_residuals = [], []
        g_dic = self.adjustment.g_dic
        # Reset residuals to all -999's
        for i in range(self.delta_model.rowCount()):
            ind = self.delta_model.createIndex(i, 0)
            tempdelta = self.delta_model.data(ind, QtCore.Qt.UserRole)
            tempdelta.residual = -999.
            self.delta_model.setData(ind, tempdelta, QtCore.Qt.UserRole)

        # Matchup adjustment residuals with observations
        for i in range(self.delta_model.rowCount()):
            ind = self.delta_model.createIndex(i, 0)
            chk = self.delta_model.data(ind, QtCore.Qt.CheckStateRole)
            tempdelta = self.delta_model.data(ind, QtCore.Qt.UserRole)
            if chk == 2:
                try:
                    if inversion_type == 'gravnet':
                        station1_name = tempdelta.sta1[:6]
                        station2_name = tempdelta.sta2[:6]
                    elif inversion_type == 'numpy':
                        station1_name = tempdelta.sta1
                        station2_name = tempdelta.sta2
                    if tempdelta.meter in cal_dic:
                        cal_adj_dg = tempdelta.dg * cal_dic[tempdelta.meter][0] * tempdelta.cal_coeff
                    else:
                        cal_adj_dg = tempdelta.dg * tempdelta.cal_coeff
                    adj_g1 = g_dic[station1_name]
                    adj_g2 = g_dic[station2_name]
                    adj_dg = adj_g2 - adj_g1
                    tempdelta.residual = adj_dg - cal_adj_dg
                    dg_residuals.append(tempdelta.residual)
                except KeyError:
                    self.msg = show_message("Key error", "Key error")
                    return
            else:
                tempdelta.residual = -999.
            self.delta_model.setData(ind, tempdelta, QtCore.Qt.UserRole)

        for i in range(self.datum_model.rowCount()):
            ind = self.datum_model.createIndex(i, 0)
            tempdatum = self.datum_model.data(ind, QtCore.Qt.UserRole)
            if inversion_type == 'numpy':
                station_name = tempdatum.station
            elif inversion_type == 'gravnet':
                station_name = tempdatum.station[0:6]
            if station_name in g_dic.keys():
                adj_g1 = g_dic[station_name]
                residual = adj_g1 - (tempdatum.g - tempdatum.meas_height * tempdatum.gradient)
                datum_residuals.append(residual)
            else:
                residual = -999
            tempdatum.residual = residual
            self.datum_model.setData(ind, tempdatum, QtCore.Qt.UserRole)

        self.adjustment.adjustmentresults.min_dg_residual = np.min(dg_residuals)
        self.adjustment.adjustmentresults.max_dg_residual = np.max(dg_residuals)
        self.adjustment.adjustmentresults.min_datum_residual = np.min(datum_residuals)
        self.adjustment.adjustmentresults.max_datum_residual = np.max(datum_residuals)

    def get_station_indices(self):
        """
        Creates a dictionary with unique (key=stationname, value=integer) pairs.
        :return: Dict with A-matrix index for each station
        """
        station_list = []
        sta_dic_ls = {}
        for delta in self.adjustment.deltas:
            station_list.append(delta.sta1)
            station_list.append(delta.sta2)
        station_list = list(set(station_list))
        for i in range(len(station_list)):
            sta_dic_ls[station_list[i]] = i
        return sta_dic_ls

    def return_obstreestation(self, station_id):
        """
        returns ObsTreeStation object corresponding to delta_id. Used when recreating deltas from a saved workspace.
        :param station_id: tuple, (station_name, time)
        :return: ObsTreeStation
        """
        for station in self.iter_stations():
            if station.station_name == station_id[0]:
                if abs(station.tmean - station_id[1]) < 0.0001:
                    return station
        return None

    def populate_delta_model(self, loop=None, clear=True):
        """
        Copy deltas from the delta_model shown on the drift tab to the model shown on the adjustment tab.
        :param loop:
        :return:
        """
        if clear:
            self.delta_model.clearDeltas()
        # If just a single loop
        if type(loop) is ObsTreeLoop:
            for ii in range(loop.delta_model.rowCount()):
                if loop.delta_model.data(loop.delta_model.index(ii, 0), role=QtCore.Qt.CheckStateRole) == 2:
                    delta = loop.delta_model.data(loop.delta_model.index(ii, 0), role=QtCore.Qt.UserRole)
                    self.delta_model.insertRows(delta, 0)
        # Populate all loops
        elif loop is None:
            for i in range(self.rowCount()):
                loop = self.child(i)
                if loop.checkState() == QtCore.Qt.Checked:
                    try:
                        for ii in range(loop.delta_model.rowCount()):
                            if loop.delta_model.data(loop.delta_model.index(ii, 0), role=QtCore.Qt.CheckStateRole) == 2:
                                delta = loop.delta_model.data(loop.delta_model.index(ii, 0), role=QtCore.Qt.UserRole)
                                # Need to create a new delta here instead of just putting the original one, from the
                                # drift tab, in the net adj. tab. Otherwise checking/unchecking on the net adj tab
                                # overrides repopulating the delta table.
                                if type(delta.station2) == list:
                                    newdelta = data_objects.Delta.from_list(delta.station2)
                                else:
                                    newdelta = data_objects.Delta(delta.station1, delta.station2)
                                newdelta.ls_drift = delta.ls_drift
                                newdelta.driftcorr = delta.driftcorr
                                newdelta.type = delta.type
                                newdelta.loop = delta.loop
                                self.delta_model.insertRows(newdelta, 0)
                    except Exception as e:
                        logging.exception(e, exc_info=True)
                        # Sometimes the delta table isn't created when a workspace is loaded
                        from gui_objects import show_message
                        self.msg = show_message("Error populating delta table. Please check the drift correction " +
                                                "for survey " + self.name + ", loop " + loop.name,
                                                "GSadjust error")
        return True

class ObsTreeModel(QtGui.QStandardItemModel):
    """
    Tree model that shows station name, date, and average g value.

    The model is populated by appending stations to loops, and loops to surveys. The parent-child relationship is not
    explicitly stored.
    """

    # TODO: Implement drag and drop.
    def __init__(self):
        super(ObsTreeModel, self).__init__()
        self.setColumnCount(3)
        self.setHorizontalHeaderLabels(['Name', 'Date', 'g (\u00b5Gal)'])
        self.station_coords = None

    signal_refresh_view = QtCore.pyqtSignal()
    signal_name_changed = QtCore.pyqtSignal()
    signal_delta_update_required = QtCore.pyqtSignal()

    def columnCount(self, QModelIndex_parent=None, *args, **kwargs):
        return 3

    def flags(self, QModelIndex):
        if not QModelIndex.isValid():
            return QtCore.Qt.NoItemFlags
        return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable | \
               QtCore.Qt.ItemIsEditable | \
               QtCore.Qt.ItemIsUserCheckable

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if index.model() is not None:
            column = index.column()
            if role == QtCore.Qt.DisplayRole:
                
                if column > 0:
                    m = index.model().itemFromIndex(index.sibling(index.row(), 0))
                else:
                    m = index.model().itemFromIndex(index)

                fn, *args = m.display_column_map.get(column, (format_numeric_column, column))
                return fn(*args)

            elif role == QtCore.Qt.CheckStateRole:
                if column == 0:
                    m = index.model().itemFromIndex(index)
                    return m.checkState()
            elif role == QtCore.Qt.ToolTipRole:
                m = index.model().itemFromIndex(index)
                try:
                    return m.tooltip
                except AttributeError:
                    return ''

    def setData(self, index, value, role):
        if role == QtCore.Qt.CheckStateRole and index.column() == 0:
            m = index.model().itemFromIndex(index)
            m.setCheckState(value)
            if type(m) is ObsTreeLoop:
                self.signal_delta_update_required.emit()

            self.dataChanged.emit(index, index)
            # self.layoutChanged.emit()
            self.signal_refresh_view.emit()
            return True

        if role == QtCore.Qt.EditRole:
            if index.isValid() and index.column() == 0:
                item = index.model().itemFromIndex(index)
                return {
                    ObsTreeStation: self._handler_edit_ObsTreeStation,
                    ObsTreeLoop: self._handler_edit_ObsTreeLoop,
                    ObsTreeSurvey: self._handler_edit_ObsTreeSurvey
                }[type(item)](item, value)

    def _handler_edit_ObsTreeStation(self, item, value):
        """
        Handle editing of ObsTreeStation objects (renaming).
        """
        from gui_objects import rename_dialog
        old_name = item.station_name
        new_name = str(value)
        if new_name is not item.station_name:
            rename_type = rename_dialog(old_name, new_name)
            if rename_type == 'Loop':
                loop = item.parent()
                loop.rename(old_name, new_name)

            if rename_type == 'Survey':
                loop = item.parent()
                survey = loop.parent()
                survey.rename(old_name, new_name)

            if rename_type == 'Campaign':
                campaign = item.model().invisibleRootItem()
                for i in range(campaign.rowCount()):
                    survey = campaign.child(i, 0)
                    survey.rename(old_name, new_name)

            if rename_type == 'Station':
                item.station_name = new_name
                for i in range(len(item.station)):
                    item.station[i] = new_name

            logging.info('Stations renamed from {} to {} in {}'.format(old_name, new_name,
                                                                       rename_type))
            self.signal_name_changed.emit()
        return True

    def _handler_edit_ObsTreeLoop(self, item, value):
        """
        Handle editing of ObsTreeLoop objects (renaming).
        """
        new_name = str(value)
        old_name = item.name
        logging.info('Loop renamed from {} to {}'.format(old_name, new_name))
        item.name = new_name
        return True

    def _handler_edit_ObsTreeSurvey(self, item, value):
        new_name = str(value)
        try:
            item.name = dt.datetime.strptime(new_name, 'yyyy-MM-dd')
        except Exception as e:
            pass
        old_name = item.name
        logging.info('Loop renamed from {} to {}'.format(old_name, new_name))
        item.name = new_name
        return True

    def checkState(self, index):
        if index.column > 0:
            m = index.model().itemFromIndex(index.sibling(index.row(), 0))
        else:
            m = index.model().itemFromIndex(index)
        return m.checkState()

    def insertRows(self, station, position, rows=1, index=QtCore.QModelIndex()):
        self.beginInsertRows(index, position, position + rows - 1)
        self.appendRow(station)
        self.endInsertRows()

    def load_workspace(self, fname):
        """
        Load previously-save workspace. Need to recreate PyQt models.

        Importantly, there are two types of deltas: those on the Drift tab and those on the network adjustment tab.
        deltas are COPIED from one to another by discrete menu commands, it's not automatic.

        The workflow:
        1) Crete obstreesurvey object so we have somewhere to store loop and station objects.
        2) deltas on the drift tab are created strictly from the station objects and specified options (drift
        correction method, etc.). We don't store a corresponding delta object in the saved workspace/json.
        3) deltas on the network adjustment tab have additional information that must be stored in the saved
        workspace/json (checked state, std. dev. for adj., etc.). When loading a workspace, this delta is INDEPENDENT of
        the delta stored in the delta table on the Drift tab. The Net Adj. tab delta table is re-created later,
        not here, because Roman-method deltas depend on the drift-tab deltas.
        :param fname:
        :return: (ObsTreeSurvey, delta_table, coords)
        """
        logging.info("Workspace loaded: " + fname)
        delta_tables, obstreesurveys = [], []
        coords, surveys = None, None
        with open(fname, "r") as f:
            data = jsons.load(json.load(f))
            if all(isinstance(x, ObsTreeSurvey) for x in data):
                surveys = data
            elif len(data) > 1:
                coords = data[1]
                surveys = data[0]
        # Populate PyQt objects. First, create obstreesurvey
        for survey in surveys:
            obstreesurvey = ObsTreeSurvey.from_json(survey)
            for loop in survey['loops']:
                obstreeloop = ObsTreeLoop.from_json(loop)
                for station in loop['stations']:
                    if 'station_name' in station:  # Sometimes blank stations are generated, not sure why?
                        temp_station = tempStation(station)
                        obstreestation = ObsTreeStation(temp_station, temp_station.station_name,
                                                        temp_station.station_count)
                        if type(obstreestation.t[0]) == dt.datetime:
                            obstreestation.t = [date2num(i) for i in obstreestation.t]
                        obstreeloop.appendRow([obstreestation,
                                               QtGui.QStandardItem('0'),
                                               QtGui.QStandardItem('0')])
                obstreesurvey.appendRow([obstreeloop,
                                         QtGui.QStandardItem('0'),
                                         QtGui.QStandardItem('0')])
            obstreesurveys.append(obstreesurvey)
            delta_tables.append(obstreesurvey.deltas)

        return (obstreesurveys, delta_tables, coords)

    def datums(self):
        datum_list = []
        for i in range(self.rowCount()):
            obstreesurvey = self.itemFromIndex(self.index(i, 0))
            for ii in range(obstreesurvey.datum_model.rowCount()):
                idx = obstreesurvey.datum_model.index(ii, 0)
                datum = obstreesurvey.datum_model.data(idx, role=QtCore.Qt.UserRole)
                datum_list.append(datum.station)
        return list(set(datum_list))

    def deltas(self):
        delta_list = []
        for i in range(self.rowCount()):
            obstreesurvey = self.itemFromIndex(self.index(i, 0))
            for ii in range(obstreesurvey.datum_model.rowCount()):
                idx = obstreesurvey.datum_model.index(ii, 0)
                delta = obstreesurvey.datum_model.data(idx, role=QtCore.Qt.UserRole)
                delta_list.append(delta.station)
        return delta_list

    def save_workspace(self, fname):
        try:
            surveys = []
            for i in range(self.rowCount()):
                obstreesurvey = self.itemFromIndex(self.index(i, 0))
                surveys.append(obstreesurvey)
            workspace_data = [surveys, self.station_coords]
            if fname[-4:] != '.gsa':
                fname += '.gsa'

            with open(fname, "w") as f:
                json.dump(jsons.dump(workspace_data), f)
            logging.info('Saving JSON workspace to {}'.format(fname))
        except Exception:
            return False
        return fname


class tempStation():
    def __init__(self, station):
        self.__dict__ = station


class RomanTableModel(QtCore.QAbstractTableModel):
    """
    Model to store the individual delta-g's calculated using the Roman-method. Shown on the drift tab.

    Shown on the drift tab.
    """
    _headers = {
        ROMAN_FROM: 'From',
        ROMAN_TO: 'To',
        ROMAN_DELTA: 'Delta g'
    }

    def __init__(self):
        super(RomanTableModel, self).__init__()
        self.dgs = []

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        if role == QtCore.Qt.TextAlignmentRole:
            if orientation == QtCore.Qt.Horizontal:
                return QVariant(int(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter))
            return QVariant(int(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter))

        if role == QtCore.Qt.DisplayRole and orientation == QtCore.Qt.Horizontal:
            return self._headers.get(section, section + 1)

    def insertRows(self, delta, position, rows=1, index=QtCore.QModelIndex()):
        self.beginInsertRows(QtCore.QModelIndex(), position, position + rows - 1)
        self.dgs.append(delta)
        self.endInsertRows()

    def rowCount(self, parent=None):
        return len(self.dgs)

    def columnCount(self, parent=None):
        return 3

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if index.isValid():
            delta = self.dgs[index.row()]
            column = index.column()
            if role == QtCore.Qt.DisplayRole:
                fn, *args = {
                    ROMAN_FROM: (str, delta.sta1),
                    ROMAN_TO: (str, delta.sta2),
                    ROMAN_DELTA: (format, delta.dg, "0.1f")
                }.get(column)
                return fn(*args)

            elif role == QtCore.Qt.UserRole:
                # check status definition
                return delta

    def setData(self, index, value, role):
        # type: (object, object, object) -> object
        """
        If a row is unchecked, update the keepdata value to 0 setData launched when role is acting value is
        QtCore.Qt.Checked or QtCore.Qt.Unchecked
        """
        if role == QtCore.Qt.CheckStateRole and index.column() == 0:
            datum = self.datums[index.row()]
            datum.checked = value
            return True

        if role == QtCore.Qt.EditRole:
            if index.isValid() and index.row() >= 0 and value:
                datum = self.datums[index.row()]
                column = index.column()

                if column == DATUM_STATION:
                    datum.station = value
                elif column == DATUM_G:
                    datum.g = float(value)
                elif column == DATUM_SD:
                    datum.sd = float(value)
                elif column == DATUM_DATE:
                    datum.date = value
                self.dataChanged.emit(index, index)

            return True

    def flags(self, index):
        return (QtCore.Qt.ItemIsEnabled |
                QtCore.Qt.ItemIsSelectable |
                QtCore.Qt.ItemIsEditable)


# noinspection PyUnresolvedReferences
class DatumTableModel(QtCore.QAbstractTableModel):
    """
    Model to store Datums, shown on the adjust tab.
    """

    _headers = {  # As map, so do not need to be kept in order with the above.
        DATUM_STATION: 'Station',
        DATUM_G: 'g',
        DATUM_SD: 'Std. dev.',
        DATUM_DATE: 'Date',
        MEAS_HEIGHT: 'Meas. height',
        GRADIENT: 'Gradient',
        DATUM_RESIDUAL: 'Residual',
        N_SETS: '# sets',
    }

    _attrs = {  # From column constants to object attributes, for setting.
        DATUM_STATION: ('station', str),
        DATUM_G: ('g', float),
        DATUM_SD: ('sd', float),
        DATUM_DATE: ('date', lambda x: x),  # pass through
        MEAS_HEIGHT: ('meas_height', float),
        GRADIENT: ('gradient', float),
    }

    signal_adjust_update_required = QtCore.pyqtSignal()

    def __init__(self):
        super(DatumTableModel, self).__init__()
        self.datums = []

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        if role == QtCore.Qt.TextAlignmentRole:
            if orientation == QtCore.Qt.Horizontal:
                return QVariant(int(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter))
            return QVariant(int(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter))

        if role == QtCore.Qt.DisplayRole and orientation == QtCore.Qt.Horizontal:
            return self._headers.get(section, section + 1)

    def insertRows(self, datum, position, rows=1, index=QtCore.QModelIndex()):
        self.beginInsertRows(QtCore.QModelIndex(), position, position + rows - 1)
        self.datums.append(datum)
        self.endInsertRows()

    def removeRow(self, index):
        datum = self.data(index, role=QtCore.Qt.UserRole)
        self.beginRemoveRows(index, index.row(), 1)
        self.datums.remove(datum)
        self.endRemoveRows()
        self.beginResetModel()
        self.endResetModel()

    def rowCount(self, parent=None):
        return len(self.datums)

    def columnCount(self, parent=None):
        return 8

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if index.isValid():
            datum = self.datums[index.row()]
            column = index.column()
            if role == QtCore.Qt.DisplayRole:
                # To accommodate old save files
                def get_nsets():
                    try:
                        return datum.n_sets
                    except:
                        return 'NA'
                fn, *args = {
                    DATUM_SD: (format, datum.sd, '0.2f'),
                    DATUM_G: (format, datum.g, '8.1f'),
                    DATUM_STATION: (str, datum.station),
                    DATUM_DATE: (QtCore.QDate.fromString, datum.date, 'yyyy-MM-dd'),
                    MEAS_HEIGHT: (format, datum.meas_height, '0.2f'),
                    GRADIENT: (format, datum.gradient, '0.2f'),
                    DATUM_RESIDUAL: (format, datum.residual, '0.1f'),
                    N_SETS: (get_nsets,)
                }.get(column, (format_numeric_column, column))

                return fn(*args)

            elif role == QtCore.Qt.CheckStateRole:
                # check status definition
                if index.column() == 0:
                    return self.checkState(datum)

            elif role == QtCore.Qt.UserRole:
                # check status definition
                return datum

    def checkState(self, datum):
        """
        By default, everything is checked. If keepdata property from the ChannelList object is 0, it is unchecked
        """
        if datum.checked == 0:
            # self.unchecked[index]= QtCore.Qt.Unchecked
            # return self.unchecked[index]
            return QtCore.Qt.Unchecked
        else:
            return QtCore.Qt.Checked

    def setData(self, index, value, role):
        # type: (object, object, object) -> object
        """
        If a row is unchecked, update the keepdata value to 0 setData launched when role is acting value is
        QtCore.Qt.Checked or QtCore.Qt.Unchecked
        """
        if role == QtCore.Qt.CheckStateRole and index.column() == 0:
            datum = self.datums[index.row()]
            if value == QtCore.Qt.Checked:
                datum.checked = 2
            elif value == QtCore.Qt.Unchecked:
                datum.checked = 0
            self.dataChanged.emit(index, index)
            self.signal_adjust_update_required.emit()
            return True

        if role == QtCore.Qt.EditRole:
            if index.isValid() and 0 <= index.row():
                if value:
                    datum = self.datums[index.row()]
                    column = index.column()
                    attr, vartype = self._attrs.get(column, (None, None))
                    if attr:
                        setattr(datum, attr, vartype(value))
                    self.dataChanged.emit(index, index)

            return True

        if role == QtCore.Qt.UserRole:
            self.datums[index.row()] = value
            self.dataChanged.emit(index, index)

    def flags(self, index):
        return (QtCore.Qt.ItemIsUserCheckable |
                QtCore.Qt.ItemIsEnabled |
                QtCore.Qt.ItemIsSelectable |
                QtCore.Qt.ItemIsEditable)

    def clearDatums(self):
        self.beginRemoveRows(self.index(0, 0), 0, self.rowCount())
        self.datums = []
        self.endRemoveRows()
        # The ResetModel calls is necessary to remove blank rows from the table view.
        self.beginResetModel()
        self.endResetModel()
        return QVariant()

    def datum_names(self):
        dn = []
        for datum in self.datums:
            dn.append(datum.station)
        return dn


class TareTableModel(QtCore.QAbstractTableModel):
    """
    Model to store tares (offsets)
    """

    _headers = {
        TARE_DATE: 'Date',
        TARE_TIME: 'Time',
        TARE_TARE: 'Tare (\u00b5Gal)'
    }

    def __init__(self):
        super(TareTableModel, self).__init__()
        self.tares = []

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        if role == QtCore.Qt.TextAlignmentRole:
            if orientation == QtCore.Qt.Horizontal:
                return QVariant(int(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter))
            return QVariant(int(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter))

        if role == QtCore.Qt.DisplayRole and orientation == QtCore.Qt.Horizontal:
            return self._headers.get(section, section + 1)

    def insertRows(self, tare, position=0, rows=1, index=QtCore.QModelIndex()):
        self.beginInsertRows(QtCore.QModelIndex(), position, position + rows - 1)
        self.tares.append(tare)
        self.endInsertRows()

    def removeRow(self, index):
        tare = self.data(index, role=QtCore.Qt.UserRole)
        self.beginRemoveRows(index, index.row(), 1)
        self.tares.remove(tare)
        self.endRemoveRows()
        self.beginResetModel()
        self.endResetModel()

    def rowCount(self, parent=None):
        return len(self.tares)

    def columnCount(self, parent=None):
        return 3

    def data(self, index, role):
        if index.isValid():
            tare = self.tares[index.row()]
            column = index.column()
            # print 'c'+str(column)
            if role == QtCore.Qt.DisplayRole:
                fn, *args = {
                    TARE_DATE: (str, tare.date),
                    TARE_TIME: (str, tare.time),
                    TARE_TARE: (format, tare.tare, "0.1f")
                }.get(column)
                return fn(*args)

            elif role == QtCore.Qt.CheckStateRole:
                if index.column() == 0:
                    return self.checkState(tare)

            elif role == QtCore.Qt.UserRole:
                return tare

    def setData(self, index, value, role):
        """
        If a row is unchecked, update the keepdata value to 0 setData launched when role is acting value is
        QtCore.Qt.Checked or QtCore.Qt.Unchecked
        """
        if role == QtCore.Qt.CheckStateRole and index.column() == 0:
            tare = self.tares[index.row()]
            if value == QtCore.Qt.Checked:
                tare.checked = 2
            elif value == QtCore.Qt.Unchecked:
                tare.checked = 0
            return True
        if role == QtCore.Qt.EditRole:
            if index.isValid() and 0 <= index.row():
                tare = self.tares[index.row()]
                column = index.column()
                if len(str(value)) > 0:
                    if column == 0:
                        tare.date = float(value)
                    elif column == 1:
                        tare.time = float(value)
                    elif column == 2:
                        tare.tare = float(value)
                    self.dataChanged.emit(index, index)
                return True
        if role == QtCore.Qt.UserRole:
            self.tares[index.row()] = value

    def deleteTare(self, idx):
        self.beginRemoveRows(idx, 0, 1)
        self.endRemoveRows

    def clearTares(self):
        self.beginRemoveRows(self.index(0, 0), 0, self.rowCount())
        self.tares = []
        self.endRemoveRows()
        # The ResetModel calls is necessary to remove blank rows from the table view.
        self.beginResetModel()
        self.endResetModel()
        return QVariant()

    def flags(self, index):
        return QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled | \
               QtCore.Qt.ItemIsEditable | QtCore.Qt.ItemIsSelectable

    def checkState(self, tare):
        """
        By default, everything is checked. If keepdata property from the ChannelList object is 0, it is unchecked
        """
        if tare.checked == 0:
            return QtCore.Qt.Unchecked
        else:
            return QtCore.Qt.Checked


class ResultsTableModel(QtCore.QAbstractTableModel):
    """
    Model to store network-adjusted gravity values.

    There is one ResultsTableModel per survey.
    """

    _headers = {
        ADJSTA_STATION: 'Station',
        ADJSTA_G: 'g',
        ADJSTA_SD: 'Std. dev'
    }

    def __init__(self):
        super(ResultsTableModel, self).__init__()
        self.adjusted_stations = []

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        if role == QtCore.Qt.TextAlignmentRole:
            if orientation == QtCore.Qt.Horizontal:
                return QVariant(int(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter))
            return QVariant(int(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter))

        if role == QtCore.Qt.DisplayRole and orientation == QtCore.Qt.Horizontal:
            return self._headers.get(section, section + 1)

    def insertRows(self, adjusted_station, position, rows=1, index=QtCore.QModelIndex()):
        self.beginInsertRows(QtCore.QModelIndex(), position, position + rows - 1)
        self.adjusted_stations.append(adjusted_station)
        self.endInsertRows()

    def rowCount(self, parent=None):
        return len(self.adjusted_stations)

    def columnCount(self, parent=None):
        return 3

    def data(self, index, role):
        if index.isValid():
            sta = self.adjusted_stations[index.row()]
            column = index.column()
            if role == QtCore.Qt.DisplayRole:
                fn, *args ={
                    ADJSTA_STATION: (str, sta.station),
                    ADJSTA_G: (format, sta.g, '8.1f'),
                    ADJSTA_SD: (format, sta.sd, '1.1f')
                }.get(column, (format_numeric_column, column))
                return fn(*args)

            if role == QtCore.Qt.UserRole:
                return sta

    def setData(self, index, value, role):
        # type: (object, object, object) -> object
        """
        If a row is unchecked, update the keepdata value to 0 setData launched when role is acting value is
        QtCore.Qt.Checked or QtCore.Qt.Unchecked
        """
        if role == QtCore.Qt.UserRole:
            self.adjusted_stations = []
            return QVariant()
        # return QtCore.QAbstractTableModel.setData(self, index, value, role)

    def flags(self, index):
        return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable

    def copyToClipboard(self):
        """
        Copies results tble to clipboard. Table needs to be selected first by clicking in the upper left corner.
        """
        clipboard = ""

        for r in range(self.rowCount()):
            for c in range(self.columnCount()):
                idx = self.index(r, c)
                clipboard += str(self.data(idx, role=QtCore.Qt.DisplayRole))
                if c != (self.columnCount() - 1):
                    clipboard += '\t'
            clipboard += '\n'

        # copy to the system clipboard
        sys_clip = QtWidgets.QApplication.clipboard()
        sys_clip.setText(clipboard)

    def clearResults(self):

        self.beginRemoveRows(self.index(0, 0), 0, self.rowCount())
        self.adjusted_stations = []
        self.endRemoveRows()
        # The ResetModel calls is necessary to remove blank rows from the table view.
        self.beginResetModel()
        self.endResetModel()
        return QVariant()


class DeltaTableModel(QtCore.QAbstractTableModel):
    """
    Model that stores delta-g's used in network adjustment. Used on the network adjustment tab and as the Roman
    average table (bottom right table when using Roman method drift correction).

    There is one DeltaTableModel per survey. The delta-g's that populate the model depend on the drift method used;
    if the Roman method is used, the delta-g's are the average of the individual delta-g's in the RomanTableModel.
    """

    _headers = {
        DELTA_STATION1: 'From',
        DELTA_STATION2: 'To',
        LOOP: 'Loop',
        DELTA_TIME: 'Time',
        DELTA_G: 'Delta-g',
        DELTA_DRIFTCORR: 'Drift corr.',
        DELTA_SD: 'Std. dev.',
        DELTA_ADJ_SD: 'SD for adj.',
        DELTA_RESIDUAL: 'Residual',
    }

    def __init__(self):
        super(DeltaTableModel, self).__init__(parent=None)
        self._deltas = []

    signal_adjust_update_required = QtCore.pyqtSignal()

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        if role == QtCore.Qt.TextAlignmentRole:
            if orientation == QtCore.Qt.Horizontal:
                return QVariant(int(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter))
            return QVariant(int(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter))

        if role == QtCore.Qt.DisplayRole and orientation == QtCore.Qt.Horizontal:
            return self._headers.get(section, section + 1)

    def insertRows(self, delta, position, rows=1, index=QtCore.QModelIndex()):
        self.beginInsertRows(QtCore.QModelIndex(), position, position + rows - 1)
        self._deltas.append(delta)
        self.endInsertRows()

    def rowCount(self, parent=None):
        if parent is None:
            return len(self._deltas)
        elif parent.isValid():
            return 0
        else:
            return len(self._deltas)

    def columnCount(self, parent=None):
        if parent is None:
            return 9
        elif parent.isValid():
            return 0
        else:
            return 9

    def data(self, index, role):
        if index.isValid():
            delta = self._deltas[index.row()]
            column = index.column()
            if role == QtCore.Qt.ForegroundRole:
                brush = QtGui.QBrush(QtCore.Qt.black)
                if delta.type == 'normal':
                    try:
                        if delta.station1.data(role=QtCore.Qt.CheckStateRole) == 2 and \
                                delta.station2.data(role=QtCore.Qt.CheckStateRole) == 2:
                            brush = QtGui.QBrush(QtCore.Qt.black)
                        else:
                            if delta.station1.data(role=QtCore.Qt.CheckStateRole) == 2:
                                if column == 0:
                                    brush = QtGui.QBrush(QtCore.Qt.darkGray)
                                else:
                                    brush = QtGui.QBrush(QtCore.Qt.lightGray)
                            elif delta.station2.data(role=QtCore.Qt.CheckStateRole) == 2:
                                if column == 1:
                                    brush = QtGui.QBrush(QtCore.Qt.darkGray)
                                else:
                                    brush = QtGui.QBrush(QtCore.Qt.lightGray)
                            else:
                                brush = QtGui.QBrush(QtCore.Qt.lightGray)
                    except:
                        catch = 1
                elif delta.type == 'assigned':
                    if column == DELTA_G:
                        brush = QtGui.QBrush(QtCore.Qt.red)
                return brush
            if role == QtCore.Qt.DisplayRole:
                def delta_station_loop():
                    if delta.loop is None:
                        if type(delta.station2) == list:
                            return delta.station2[0].loop
                        else:
                            return "NA"
                    else:
                        return delta.loop

                def format_datetime(date):
                    return dt.datetime.utcfromtimestamp((date - 719163) * 86400.).strftime('%Y-%m-%d %H:%M:%S')

                fn, *args = {
                    DELTA_STATION1: (str, delta.sta1),
                    DELTA_STATION2: (str, delta.sta2),
                    LOOP: (str, delta_station_loop()),
                    DELTA_TIME: (format_datetime, delta.time()),
                    DELTA_G: (format, delta.dg, "0.1f"),
                    DELTA_DRIFTCORR: (format, delta.driftcorr, "0.1f"),
                    DELTA_SD: (format, delta.sd, "0.1f"),
                    DELTA_ADJ_SD: (format, delta.sd_for_adjustment, "0.1f"),
                    DELTA_RESIDUAL: (format, delta.residual, "0.1f")
                }.get(column, (str, "NA"))

                return fn(*args)

            elif role == QtCore.Qt.CheckStateRole:
                if index.column() == 0:
                    return self.checkState(delta)

            elif role == QtCore.Qt.UserRole:
                return delta

    def setData(self, index, value, role):
        """
        If a row is unchecked, update the keepdata value to 0 setData launched when role is acting value is
        QtCore.Qt.Checked or QtCore.Qt.Unchecked
        """
        if role == QtCore.Qt.CheckStateRole and index.column() == 0:
            delta = self._deltas[index.row()]
            if value == QtCore.Qt.Checked:
                delta.checked = 2
            elif value == QtCore.Qt.Unchecked:
                delta.checked = 0
            self.signal_adjust_update_required.emit()
            return True
        if role == QtCore.Qt.EditRole:
            if index.isValid() and 0 <= index.row():
                delta = self._deltas[index.row()]
                column = index.column()
                if len(str(value)) > 0:
                    if column == DELTA_ADJ_SD:
                        delta.adj_sd = float(value)
                    if column == DELTA_G:
                        delta.type = 'assigned'
                        delta.assigned_dg = float(value)
                    self.dataChanged.emit(index, index)
                    self.signal_adjust_update_required.emit()
                return True
        if role == QtCore.Qt.UserRole:
            self._deltas[index.row()] = value

    def clearDeltas(self):
        self.beginRemoveRows(self.index(0, 0), 0, self.rowCount())
        for delta in self._deltas:
            delta.residual = -999
        self._deltas = []
        self.endRemoveRows()
        # The ResetModel calls is necessary to remove blank rows from the table view.
        self.beginResetModel()
        self.endResetModel()
        return QVariant()

    def updateTable(self):
        self.beginResetModel()
        self.endResetModel()
        return QVariant()

    def flags(self, QModelIndex):
        if not QModelIndex.isValid():
            return QtCore.Qt.NoItemFlags
        return QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled | \
               QtCore.Qt.ItemIsEditable | QtCore.Qt.ItemIsSelectable

    def checkState(self, delta):
        """
        By default, everything is checked. If keepdata property from the ChannelList object is 0, it is unchecked
        """
        if delta.checked == 0:
            return QtCore.Qt.Unchecked
        else:
            return QtCore.Qt.Checked

    @classmethod
    def from_json(cls, data):
        """
        When loading a workspace, repopulate PyQt models
        """
        temp = cls(data['name'])

        deltas = []
        for delta in data['deltas']:
            from types import SimpleNamespace
            sd = SimpleNamespace()
            sd.adj_sd = delta['adj_sd']
            sd.checked = delta['checked']
            sd.driftcorr = delta['driftcorr']
            sd.loop = delta['loop']
            sd.ls_drift = delta['ls_drift']
            sd.type = delta['type']
            try:
                sd.assigned_dg = delta['assigned_dg']
            except:
                pass
            # d = data_objects.SimpleDelta(sd)
            try:
                sd.sta1 = delta['sta1']
                sd.sta2 = delta['sta2']
            except KeyError as e:
                # Raised if delta type is 'assigned'
                pass
            deltas.append(sd)

        temp.deltas = deltas

        for datum in data['datums']:
            d = data_objects.Datum(datum['station'])
            d.__dict__ = datum
            temp.datum_model.insertRows(d, 0)

        ao = data_objects.AdjustmentOptions()
        ao.__dict__ = data['adjoptions']
        temp.adjustment.adjustmentoptions = ao

        if 'checked' in data:
            temp.setCheckState(data['checked'])
        return temp

    def to_json(self):
        # Normal delta
        json_deltas = []
        for delta in self._deltas:
            if delta.type == 'normal' or delta.type == 'assigned':
                sta1 = delta.station1.key
                sta2 = delta.station2.key
            elif delta.type == 'list':
                sta1 = None
                sta2 = []
                for threepoint in delta.station2:
                    stations = [threepoint.station1.key, threepoint.station2[0].key, threepoint.station2[1].key]
                    sta2.append(stations)
            adj_sd = delta.adj_sd
            delta_type = delta.type
            ls_drift = delta.ls_drift
            if delta.loop is None:
                if isinstance(delta.station2, list):
                    loop = delta.station2[0].loop
                else:
                    loop = "NA"
            else:
                loop = delta.loop
            driftcorr = delta.driftcorr
            checked = delta.checked
            temp_delta = {
                'sta1': sta1,
                'sta2': sta2,
                'adj_sd': adj_sd,
                'ls_drift': ls_drift,
                'driftcorr': driftcorr,
                'checked': checked,
                'loop': loop,
                'type': delta_type
            }
            json_deltas.append(temp_delta)
        return json_deltas


class ScintrexTableModel(QtCore.QAbstractTableModel):
    """
    Model to store Scintrex data.

    There is one ScintrexTableModel (or BurrisTableModel) per station occupation. The model is created dynamically
    each time a station is selected in the tree view on the data tab, rather than stored in memory.


    The station position in the data hierarchy are stored, so that if a modification is triggered, the original data
    can be accessed and changed accordingly (keysurv,keyloop,keysta)

    by default, all table entries are checked (this can be modified to allow pre-check based on user criteria (tiltx,
    tilty,...)). Then, if one is unchecked, the keepdata property of the ChannelList object at the table row position
    is set to 0

    properties:
    _headers:             table header
    unchecked:            a dictionnary of unchecked items. Keys are item
                          indexes, entries are states
    ChannelList_obj:      an object of ChannelList-type: used to store the
                          table data as the structured data
    arraydata:            an array representation of the data from the
                          ChannelList_obj
    keysurv:              the key of a survey object within the grav_obj
                          hierarchy
    keyloop:              the key of a loop object within the grav_obj
                          hierarchy
    keysta:               the key of a station object within the grav_obj
                          hierarchy
    """

    _headers = {
        SCINTREX_STATION: "Station",
        SCINTREX_DATE: "Date",
        SCINTREX_G: "g (\u00b5gal)",
        SCINTREX_SD: "sd (\u00b5gal)",
        SCINTREX_X_TILT: "X Tilt",
        SCINTREX_Y_TILT: "Y Tilt",
        SCINTREX_TEMP: "Temp (K)",
        SCINTREX_DUR: "dur (s)",
        SCINTREX_REJ: "rej"
    }

    signal_update_coordinates = QtCore.pyqtSignal()
    signal_adjust_update_required = QtCore.pyqtSignal()
    signal_uncheck_station = QtCore.pyqtSignal()
    signal_check_station = QtCore.pyqtSignal()

    def __init__(self, ChannelList_obj, parent=None):
        QtCore.QAbstractTableModel.__init__(self, parent)
        self.unchecked = {}
        self.createArrayData(ChannelList_obj)

    def createArrayData(self, ChannelList_obj):
        """
        Create the np array data for table display, and update the ChannelList_obj. This function can be called from
        outside to update the table display
        """
        self.ChannelList_obj = ChannelList_obj
        self.arraydata = np.concatenate((ChannelList_obj.station,
                                         np.array(ChannelList_obj.t),
                                         np.array(ChannelList_obj.grav),
                                         np.array(ChannelList_obj.sd), ChannelList_obj.tiltx,
                                         ChannelList_obj.tilty, ChannelList_obj.temp, ChannelList_obj.dur,
                                         ChannelList_obj.rej)).reshape(len(ChannelList_obj.t), 9, order='F')

    def rowCount(self, parent=None):
        return len(self.ChannelList_obj.t)

    def columnCount(self, parent):
        return len(self._headers)

    def flags(self, index):
        return QtCore.Qt.ItemIsUserCheckable | \
               QtCore.Qt.ItemIsEnabled | \
               QtCore.Qt.ItemIsSelectable

    def data(self, index, role):
        if not index.isValid():
            return None
        if role == QtCore.Qt.DisplayRole:
            # view definition
            row = index.row()
            column = index.column()
            try:
                value = float(self.arraydata[row][column])
            except ValueError:
                value = self.arraydata[row][column]

            def format_datetime(dt):
                return num2date(float(dt)).strftime('%Y-%m-%d %H:%M:%S')

            fn, *args = {
                SCINTREX_DATE: (format_datetime, value),
                SCINTREX_REJ: (format, value, "2.0f"),
                SCINTREX_DUR: (format, value, "3.0f"),
                SCINTREX_G: (format, value, "8.1f"),
                SCINTREX_SD: (format, value, "8.1f")
            }.get(column, (str, value))

            return fn(*args)

        if role == QtCore.Qt.CheckStateRole:
            # check status definition
            if index.column() == 0:
                return self.checkState(index)

    def checkState(self, index):
        """
        By default, everything is checked. If keepdata property from the ChannelList object is 0, it is unchecked
        """
        if self.ChannelList_obj.keepdata[index.row()] == 0:
            self.unchecked[index] = QtCore.Qt.Unchecked
            return self.unchecked[index]
        else:
            return QtCore.Qt.Checked

    def setData(self, index, value, role):
        # type: (object, object, object) -> object
        """
        if a row is unchecked, update the keepdata value to 0 setData launched when role is acting value is
        QtCore.Qt.Checked or QtCore.Qt.Unchecked
        """
        if role == QtCore.Qt.CheckStateRole and index.column() == 0:
            if value == QtCore.Qt.Checked:
                self.ChannelList_obj.keepdata[index.row()] = 1
                self.signal_check_station.emit()
            elif value == QtCore.Qt.Unchecked:
                self.unchecked[index] = value
                self.ChannelList_obj.keepdata[index.row()] = 0
                if not any(self.ChannelList_obj.keepdata):
                    self.signal_uncheck_station.emit()
            self.signal_adjust_update_required.emit()
            self.dataChanged.emit(index, index)
            return True

    def headerData(self, section, orientation, role):
        if role == QtCore.Qt.DisplayRole and orientation == QtCore.Qt.Horizontal:
            return self._headers.get(section, section + 1)


# noinspection PyUnresolvedReferences
class BurrisTableModel(QtCore.QAbstractTableModel):
    """
    Model to store Burris data.

    There is one BurrisTableModel (or ScintrexTableModel) per station occupation. The model is created dynamically
    each time a station is selected in the tree view on the data tab, rather than stored in memory.

    The station position in the data hierarchy are stored, so that if a modification is triggered, the original data
    can be accessed and changed accordingly (keysurv,keyloop,keysta)

    by default, all table entries are checked (this can be modified to allow pre-check based on user criteria (tiltx,
    tilty,...)). Then, if one is unchecked, the keepdata property of the ChannelList object at the table row position
    is set to 0

    properties:
    _headers:             table header
    unchecked:            a dictionary of unchecked items. Keys are item
                          indexes, entries are states
    ChannelList_obj:      an object of ChannelList-type: used to store the
                          table data as the structured data
    arraydata:            an array representation of the data from the
                          ChannelList_obj
    keysurv:              the key of a survey object within the grav_obj
                          hierarchy
    keyloop:              the key of a loop object within the grav_obj
                          hierarchy
    keysta:               the key of a station object within the grav_obj
                          hierarchy
    """

    _headers = {
        BURRIS_STATION: "Station",
        BURRIS_OPER: "Oper",
        BURRIS_METER: "Meter",
        BURRIS_DATE: "Date",
        BURRIS_G: "g (\u00b5gal)",
        BURRIS_DIAL: "Dial setting",
        BURRIS_FEEDBACK: "Feedback",
        BURRIS_TIDE: "Tide corr.",
        BURRIS_ELEVATION: "Elev.",
        BURRIS_LAT: "Lat",
        BURRIS_LONG: "Long"
    }

    signal_update_coordinates = QtCore.pyqtSignal()
    signal_adjust_update_required = QtCore.pyqtSignal()
    signal_uncheck_station = QtCore.pyqtSignal()
    signal_check_station = QtCore.pyqtSignal()

    def __init__(self, ChannelList_obj, parent=None):
        QtCore.QAbstractTableModel.__init__(self, parent)
        self.unchecked = {}
        self.ChannelList_obj = None
        self.createArrayData(ChannelList_obj)
        self.arraydata = None

    def createArrayData(self, ChannelList_obj):
        """
        Create the np array data for table display, and update the ChannelList_obj. This function can be called from
        outside to update the table display
        """
        # channel list attributes aren't specified in advance, so don't worry whether it's Burris or CG5 data
        self.ChannelList_obj = ChannelList_obj
        try:
            self.arraydata = np.concatenate((ChannelList_obj.station,
                                             ChannelList_obj.oper,
                                             ChannelList_obj.meter,
                                             ChannelList_obj.t,
                                             np.array(ChannelList_obj.grav),
                                             np.array(ChannelList_obj.dial),
                                             np.array(ChannelList_obj.feedback),
                                             np.array(ChannelList_obj.etc),
                                             np.array(ChannelList_obj.elev),
                                             np.array(ChannelList_obj.lat),
                                             np.array(ChannelList_obj.long))). \
                reshape(len(ChannelList_obj.t), 11, order='F')
        except:
            return

    def rowCount(self, parent=None):
        return len(self.ChannelList_obj.t)

    def columnCount(self, parent=None):
        return len(self._headers)

    def flags(self, index):
        return QtCore.Qt.ItemIsUserCheckable | \
               QtCore.Qt.ItemIsEnabled | \
               QtCore.Qt.ItemIsSelectable | \
               QtCore.Qt.ItemIsEditable

    def data(self, index, role):
        if not index.isValid():
            return None
        if role == QtCore.Qt.DisplayRole:
            row = index.row()
            column = index.column()
            try:
                value = float(self.arraydata[row][column])
            except ValueError:
                value = self.arraydata[row][column]

            def format_datetime(dt):
                return num2date(float(dt)).strftime('%Y-%m-%d %H:%M:%S')

            fn, *args = {
                BURRIS_DATE: (format_datetime, value),
                BURRIS_TIDE: (format, value, "0.1f"),
            }.get(column, (str, value))

            return fn(*args)

        if role == QtCore.Qt.CheckStateRole:
            # check status definition
            if index.column() == 0:
                return self.checkState(index)

    def checkState(self, index):
        """
        By default, everything is checked. If keepdata property from the ChannelList object is 0, it is unchecked
        """
        if self.ChannelList_obj.keepdata[index.row()] == 0:
            self.unchecked[index] = QtCore.Qt.Unchecked
            return self.unchecked[index]
        else:
            return QtCore.Qt.Checked

    def setData(self, index, value, role, silent=False):
        """
        If a row is unchecked, update the keepdata value to 0 setData launched when role is acting value is
        QtCore.Qt.Checked or QtCore.Qt.Unchecked
        """
        if role == QtCore.Qt.CheckStateRole and index.column() == 0:

            if value == QtCore.Qt.Checked:
                self.ChannelList_obj.keepdata[index.row()] = 1
                self.signal_check_station.emit()
            elif value == QtCore.Qt.Unchecked:
                self.unchecked[index] = value
                self.ChannelList_obj.keepdata[index.row()] = 0
                if not any(self.ChannelList_obj.keepdata):
                    self.signal_uncheck_station.emit()
            if not silent:
                self.signal_adjust_update_required.emit()
                self.dataChanged.emit(index, index)
            return True

        elif role == QtCore.Qt.EditRole:
            if index.isValid() and 0 <= index.row():
                if len(str(value)) > 0:
                    attr = None
                    if index.column() == 1:  # Oper
                        self.arraydata[index.row()][index.column()] = value
                        attr = 'oper'
                    if index.column() == 2:  # Meter
                        self.arraydata[index.row()][index.column()] = value
                        attr = 'meter'
                    if index.column() == 9:  # Lat
                        self.arraydata[index.row()][index.column()] = float(value)
                        attr = 'lat'
                        self.signal_update_coordinates.emit()
                    if index.column() == 10:  # Long
                        self.arraydata[index.row()][index.column()] = float(value)
                        attr = 'long'
                        self.signal_update_coordinates.emit()
                    if attr is not None:
                        self.dataChanged.emit(index, index)
                        return True

    def headerData(self, section, orientation, role):
        if role == QtCore.Qt.DisplayRole and orientation == QtCore.Qt.Horizontal:
            return self._headers.get(section, section + 1)


class GravityChangeModel(QtCore.QAbstractTableModel):
    """
    Model to store gravity change between surveys.

    There is only one such model per campaign. Gravity change is calculated when the respective menu item is chosen.
    """

    def __init__(self, header, table, table_type='simple', parent=None):
        QtCore.QAbstractTableModel.__init__(self, parent)
        self._headers = {n: col for n, col in enumerate(header)}
        self.createArrayData(table, table_type)

    def createArrayData(self, table, table_type):
        """
        Create the np array data for table display, and update the ChannelList_obj. This function can be called from
        outside to update the table display
        """
        if table_type == 'simple' or table_type == 'list':
            array = np.array(table).transpose()
        elif table_type == 'full':
            array = np.array(table)
        self.arraydata = array

    def rowCount(self, parent=None):
        return self.arraydata.shape[0]

    def columnCount(self, parent=None):
        return len(self._headers)

    def flags(self, index):
        return QtCore.Qt.ItemIsEnabled | \
               QtCore.Qt.ItemIsSelectable

    def data(self, index, role):
        if not index.isValid():
            return None
        if role == QtCore.Qt.DisplayRole:
            row = index.row()
            column = index.column()
            return str(self.arraydata[row][column])

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        if role == QtCore.Qt.DisplayRole and orientation == QtCore.Qt.Horizontal:
            return self._headers.get(section, section + 1)


class CustomSortingModel(QtCore.QSortFilterProxyModel):
    """
    Used to sort by date in importAbsG dialog
    """

    def lessThan(self, left, right):
        col = left.column()
        dataleft = left.data()
        dataright = right.data()

        if col == 3:
            dataleft = QtCore.QDate.fromString(dataleft, "MM/dd/yy").addYears(100)
            dataright = QtCore.QDate.fromString(dataright, "MM/dd/yy").addYears(100)
        elif not col == 0:
            dataleft = float(dataleft)
            dataright = float(dataright)

        return dataleft < dataright


class MeterCalibrationModel(QtGui.QStandardItemModel):
    def __init__(self):
        super(MeterCalibrationModel, self).__init__()
        self.setColumnCount(2)
        self.setHorizontalHeaderLabels(['Meter', 'Calibration factor'])


def survey_serializer(obj, cls, **kwargs):
    """
    Handle serialization of ObsTreeSurvey via .to_json() method.
    """
    return obj.to_json()


jsons.set_serializer(survey_serializer, ObsTreeSurvey)
jsons.set_serializer(survey_serializer, DeltaTableModel)
