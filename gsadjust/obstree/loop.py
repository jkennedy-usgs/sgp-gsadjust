"""
obstree/loop.pu
===============

PyQt models for loops in GSadjust tree view.
--------------------------------------------------------------------------------

NB: PyQt models follow the PyQt CamelCase naming convention. All other
methods/functions in GSadjust use PEP-8 lowercase_underscore convention.

This software is preliminary, provisional, and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.
"""

import logging
from collections import defaultdict

from PyQt5 import QtGui
from PyQt5.QtCore import Qt

from .base import ObsTreeItemBase
from .station import ObsTreeStation
from ..data import Tare

# Constants for column headers
STATION_NAME, STATION_DATETIME, STATION_MEAN = range(3)
LOOP_NAME = 0
SURVEY_NAME = 0


class ObsTreeLoop(ObsTreeItemBase):
    """
    PyQt model for Loops, parent item for stations.

    Loops are groups of station data that have a common drift method (and
    accompanying parameters), observed in order.
    """

    def __init__(self, name):
        super(ObsTreeLoop, self).__init__()

        self.deltas = []
        self.tares = []

        self.name = name
        self.drift_method = "none"  # 'none', 'netadj', 'roman', or 'continuous'

        # If continuous model, also need to keep track of which type of model
        self.drift_cont_method = 0

        # behavior at start/end. 0: extrapolate, 1: constant
        self.drift_cont_startend = 0

        # check box check state
        self.drift_cont_weighting = 0

        # If netadj method, keep track of polynomial degree
        self.drift_netadj_method = 2

        # Ignore repeat observations with delta-t longer than
        self.time_extent_check = 0
        self.time_extent_time = (1, 0)

        self.comment = ""  # String that can be specified from GUI

        # TODO: Import comment from Burris file?
        self.source = ""  # Filename of raw input file

    def __str__(self):
        if self.drift_method == "roman" or self.drift_method == "none":
            return "Loop: {}, Drift method: {}, Meter type: {}\n".format(
                self.name, self.drift_method, self.meter_type
            )

        elif self.drift_method == "continuous":
            return (
                "Loop: {}, "
                "Drift method: {}, "
                "Meter type: {}, "
                "Continuous drift method: {}, "
                "Continuous drift start/end method: {}, "
                "Weighting: {}\n".format(
                    self.name,
                    self.drift_method,
                    self.meter_type,
                    self.drift_cont_method,
                    self.drift_cont_startend,
                    self.drift_cont_weighting,
                )
            )
        elif self.drift_method == "netadj":
            return (
                "Loop: {}, "
                "Drift method: {}, "
                "Meter type: {}, "
                "Netadj method: {}, "
                "Time extent filter: {}, {}\n".format(
                    self.name,
                    self.drift_method,
                    self.meter_type,
                    self.drift_netadj_method,
                    self.time_extent_check,
                    self.time_extent_time,
                )
            )

    @property
    def oper(self):
        try:
            return self.child(0).oper[0]
        except AttributeError:
            return ""

    @property
    def meter(self):
        try:
            return self.child(0).meter[0]
        except AttributeError:
            return ""

    @property
    def tooltip(self):
        return (
            "Loop: {}\n"
            "Drift method: {}\n"
            "Meter: {}\n"
            "Operator: {}\n"
            "Comment: {}\n"
            "Source: {}".format(
                self.name,
                self.__dict__.get("drift_method", ""),
                self.meter,
                self.oper,
                self.__dict__.get("comment", ""),
                self.__dict__.get("source", ""),
            )
        )

    @property
    def display_column_map(self):
        return {
            # Column name map
            # index: name
            LOOP_NAME: (str, self.name),
            1: (str, ""),
            2: (str, ""),
        }

    @property
    def meter_type(self):
        station = self.child(0)
        # To deal with old-format files (2020-06). Can probably be deleted someday.
        if station.meter_type == "Scintrex":
            return "CG5"
        else:
            return station.meter_type

    @classmethod
    def from_json(cls, data):
        temp = cls(data["name"])
        temp.__dict__.update(data)
        potential_missing_fields = {
            "comment": "",
            "oper": "",
            "source": "",
            "drift_cont_weighting": 0,
        }
        tare_objects = []
        for tare_dict in data["tares"]:
            tare_object = Tare(tare_dict["datetime"], tare_dict["tare"])
            tare_object.checked = tare_dict["checked"]
            tare_objects.append(tare_object)
        temp.tares = tare_objects
        if "checked" in data:
            temp.setCheckState(data["checked"])
        for k, v in potential_missing_fields.items():
            if not hasattr(temp, k):
                setattr(temp, k, v)
        return temp

    def rename(self, from_name, to_name):
        for i in range(self.rowCount()):
            item = self.child(i, 0)
            if item.station_name == from_name:
                item.station_name = to_name
                for iiii in range(len(item.station)):
                    item.station[iiii] = to_name

    def update_tide(self, *args, **kwargs):
        for i in range(self.rowCount()):
            item = self.child(i, 0)
            item.update_tide(*args, **kwargs)

    def populate(self, data):
        """
        Populate loop dictionary using data passed as an option. For now, only
        the 'single' option works.

        Parameters
        ----------
        data: ChannelList object
        """
        prev_sta = data.station[0]
        ind_start = 0

        # Default dict, will default entries to zero when accessed.
        station_count_dic = defaultdict(int)

        # loop over samples:
        for i in range(len(data.station)):
            curr_sta = data.station[i]
            if curr_sta != prev_sta:
                station_count_dic[prev_sta] += 1
                ind_end = i - 1
                # create a new ChannelList object:
                temp_sta = data.extract_subset_idx(ind_start, ind_end)
                obstreestation = ObsTreeStation(
                    temp_sta, prev_sta, str(station_count_dic[prev_sta])
                )
                obstreestation.meter_type = data.meter_type
                logging.info(
                    "Station added: " + prev_sta + str(station_count_dic[prev_sta])
                )
                ind_start = i
                self.appendRow(
                    [obstreestation, QtGui.QStandardItem("0"), QtGui.QStandardItem("0")]
                )
            if i == len(data.station) - 1:
                # enter last data line
                ind_end = i
                station_count_dic[prev_sta] += 1

                # create a new Loop-type object:
                temp_sta = data.extract_subset_idx(ind_start, ind_end)
                obstreestation = ObsTreeStation(
                    temp_sta, curr_sta, str(station_count_dic[curr_sta])
                )
                obstreestation.meter_type = data.meter_type
                logging.info(
                    "Station added: "
                    + curr_sta
                    + "("
                    + str(station_count_dic[curr_sta])
                    + ")"
                )
                self.appendRow(
                    [obstreestation, QtGui.QStandardItem("0"), QtGui.QStandardItem("0")]
                )
            prev_sta = curr_sta

    def get_data_for_plot(self):
        """
        Retrieves data from loop that is used for plotting: occupations must be
        checked, and there must be more than one occupation at a station.

        Returns
        -------
        list
            One entry per station: [[plot x],[plot y], station name, [y s.d.]]
        """
        unique_stations = set()
        plot_data = []
        for i in range(self.rowCount()):
            station = self.child(i)
            if self.child(i).data(role=Qt.CheckStateRole) == 2:
                if (station.station_name, station.meter[0]) not in unique_stations:
                    unique_stations.add(station.station_name)

        for unique_station in unique_stations:
            x, y, y_sd, t_sd = [], [], [], []
            for i in range(self.rowCount()):
                station = self.child(i)
                if self.child(i).data(role=Qt.CheckStateRole) == 2:
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
        Retrieves checked stations from loop; used for calculating delta-gs in
        the "None" and "Continuous" drift options.

        Returns
        -------
        list
            list of checked stations
        """
        stations = []
        for i in range(self.rowCount()):
            station = self.child(i)
            if self.child(i).data(role=Qt.CheckStateRole) == 2:
                stations.append(station)
        return stations

    def stations(self):
        """
        Convenience method to get all stations in a loop.

        Returns
        -------
        list
            list of station from PyQt model
        """
        stations = []
        for i in range(self.rowCount()):
            stations.append(self.child(i))
        return stations

    def n_unique_stations(self):
        sn = set(s.station_name for s in self.stations())
        return len(sn)

    def to_json(self):
        """
        For serializing (PyQt models can't be serialized directly)

        Returns
        -------
        dict
            json-friendly dict
        """

        # Copy stations from PyQt models to lists
        stations_json = [s.to_json() for s in self.stations()]
        tares_json = [t.to_json() for t in self.tares]

        return {
            "checked": self.checkState(),
            "delta_model": None,
            "tare_model": None,
            "stations": stations_json,
            "tares": tares_json,
            "name": self.name,
            "drift_method": self.drift_method,
            # If continuous model, also need to keep track of which type of model
            "drift_cont_method": self.drift_cont_method,
            # behavior at start/end. 0: extrapolate, 1: constant
            "drift_cont_startend": self.drift_cont_startend,
            # If netadj method, keep track of polynomial degree
            "drift_netadj_method": self.drift_netadj_method,
            # Meter S/N, for the case where multiple meters are calibrated
            "meter": self.meter,
            "time_extent": self.time_extent_time,
            "comment": self.comment,
            "oper": self.oper,
            "source": self.source,
        }
