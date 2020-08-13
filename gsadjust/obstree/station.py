"""
pyqt_modules.py
===============

PyQt models for GSadjust. Handles assembling input matrices for
network adjustment.
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
import copy

import numpy as np
from matplotlib.dates import num2date

from .base import ObsTreeItemBase

# Constants for column headers
STATION_NAME, STATION_DATETIME, STATION_MEAN = range(3)


class ObsTreeStation(ObsTreeItemBase):
    """
    PyQt model for stations. Station data is stored as a Station object in .station
    """

    def __init__(self, k, station_name, station_count):
        """
        Create a new station from ChannelList object. The fields (lists) of the
        ChannelList are copies to the ObsTreeStation.
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
        weights for Scintrex data are calculated based on the sample standard
        deviation divided by the duration. Neither of these are reported by the
        Burris meter.  Instead, duration is set to a constant (5) when the file
        is read. Parameter sd is calculated from the multiple samples comprising
        a station. Therefore w is a constant array if Burris data, variable
        if Scintrex.

        TODO: test effect of weighting if Burris and Scintrex data are combined
        in a survey
        :return:
        """
        if self.meter_type == 'Burris' or self.meter_type == 'CG6Tsoft':
            return [1 for i in self.keepdata if i == 1]
        else:
            # sdtmp = [self.sd[i] / np.sqrt(self.dur[i]) for i in range(len(self.t))]
            w = [
                1.0 / (self.sd[i] * self.sd[i])
                for i in range(len(self.t))
                if self.keepdata[i] == 1
            ]
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
        return (self.station_name, self.tmean())

    # @property
    def grav(self):
        """
        Applies tares and earth tide correction to raw_grav
        :return: List
        """
        # if len(self.corr_g) == 0:
        #     self.corr_g = self.raw_grav
        # data = np.array(self.raw_grav) - np.array(self.tare) + np.array(self.etc)
        data = [a - b + c for (a, b, c) in zip(self.raw_grav, self.tare, self.etc)]
        return data

    def display_name(self):
        return self.station_name + ' (' + self.station_count + ')'

    def display_datetime_or_tmean(self):
        return (
            num2date(self.tmean()).strftime('%Y-%m-%d %H:%M:%S')
            if self.tmean() != -999
            else self.tmean()
        )

    def display_mean_stddev(self):
        return '{:.1f} ± {:.1f}'.format(self.gmean(), self.stdev())

    @property
    def display_column_map(self):
        return {
            STATION_NAME: (self.display_name,),
            STATION_DATETIME: (self.display_datetime_or_tmean,),
            STATION_MEAN: (self.display_mean_stddev,),
        }

    @property
    def tooltip(self):
        return None

    def _filter(self, d):
        """
        Filter data by the .keepdata list.
        """
        return [v for i, v in enumerate(d) if self.keepdata[i] == 1]

    # @property
    def gmean(self):
        """
        The try-except block handles errors when all keepdata == 0.
        """
        g = self.grav()
        try:
            gtmp = self._filter(g)
            d = sum(self._weights_())
            w = self._weights_()
            wg = [g * w for (g, w) in zip(gtmp, w)]
            gmean = sum(wg) / d
            return gmean
        except:
            return -999

    # @property
    def tmean(self):
        """
        The try-except block handles errors when all keepdata == 0.
        """
        try:
            ttmp = self._filter(self.t)
            d = sum(self._weights_())
            w = self._weights_()
            wt = [g * w for (g, w) in zip(ttmp, w)]
            tmean = sum(wt) / d
            return tmean
        except:
            return -999

    @property
    def t_stdev(self):
        try:
            if self.sd[0] == -999:
                return 1
            else:
                ttmp = np.array(self._filter(self.t))
                num = np.zeros(len(ttmp))
                for i in range(len(ttmp)):
                    num[i] = self._weights_()[i] * ((ttmp[i] - np.mean(ttmp)) * 24) ** 2
                sd = np.sqrt(
                    np.sum(num)
                    / ((len(self._weights_()) - 1) * np.sum(self._weights_()))
                )  # np.sqrt(1. / sum(
                # self._weights_()))
                return sd
        except:
            return -999

    @property
    def original_sd(self):
        g = self.grav()
        sdtmp = np.array(self._filter(self.sd))
        gtmp = np.array(self._filter(g))
        num = np.zeros(len(sdtmp))
        for i in range(len(sdtmp)):
            num[i] = self._weights_()[i] * (gtmp[i] - np.mean(gtmp)) ** 2
        sd = np.sqrt(
            np.sum(num) / ((len(self._weights_()) - 1) * np.sum(self._weights_()))
        )
        return sd

    # @property
    def stdev(self):
        """
        The try-except block handles errors when all keepdata == 0.

        The Scintrex meter reports an S.D. for each sample; the Burris meter
        does not report S.D. (with the Palm PDA, it shows S.D. on the display
        but does not record it).
        """
        try:
            g = self.grav()
            if hasattr(self, 'assigned_sd'):
                if self.assigned_sd:
                    return self.assigned_sd
            if self.sd[0] == -999:
                # Burris meter: return the S.D. calculated from all the samples at a station
                gtmp = np.array(self._filter(g))
                if len(gtmp) == 1:  # Can't take S.D. if there's only one sample
                    return 3.0
                else:
                    return float(np.std(gtmp))
            else:  # Scintrex: return the average SD of the samples
                sdtmp = np.array(self._filter(self.sd))
                # sd = np.mean(sdtmp)  #
                gtmp = np.array(self._filter(g))
                num = np.zeros(len(sdtmp))

                for i in range(len(sdtmp)):
                    num[i] = self._weights_()[i] * (gtmp[i] - np.mean(gtmp)) ** 2
                sd = np.sqrt(
                    np.sum(num)
                    / ((len(self._weights_()) - 1) * np.sum(self._weights_()))
                )  # np.sqrt(1. / sum(
                return sd
        except:
            return -999

    @property
    def summary_str(self):
        if self.tmean() == -999:
            tm = self.tmean()
        else:
            tm = num2date(self.tmean()).strftime('%Y-%m-%d %H:%M:%S')
        if len(self.lat) == 0:
            self.lat, self.long, self.elev = [0], [0], [0]
        summary_str = '{} {} {} {} {} {} {} {:0.3f} {:0.3f}\n'.format(
            self.display_name(),
            tm,
            self.oper[0],
            self.meter_type,
            self.lat[0],
            self.long[0],
            self.elev[0],
            self.gmean(),
            self.stdev(),
        )
        return summary_str

    def iter_samples(self):
        """
        Iterator that returns print statements for each sample, used when
        writing adjustment summary
        :return: All of the stations in a campaign
        """
        g = self.grav()
        for i in range(len(self.raw_grav)):
            yield '{} {} {:0.2f} {:0.2f} {:0.2f} {:0.2f}\n'.format(
                self.keepdata[i],
                self.station[i],
                self.raw_grav[i],
                self.etc[i],
                g[i],
                self.sd[i],
            )

    def to_json(self):
        jsonalizable_dict = self.__dict__
        jsonalizable_dict['checked'] = self.checkState()
        return self.__dict__
