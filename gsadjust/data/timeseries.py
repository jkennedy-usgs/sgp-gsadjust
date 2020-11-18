"""
data/timeseries.py
===============

GSadjust objects for TimeSeries: Used only for tide correction.

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

from matplotlib.dates import date2num
from scipy.interpolate import interp1d


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
            print('pb at line {:d} : check raw data file: possibly last line?'.format(i))

    def interpolate_on_given_times(self, t):
        """
        Used in tide correction routines to interpolate the time series on the
        user input time vector overlain the previous t and d fields
        """
        tord = [date2num(tmp) for tmp in self.t]
        f = interp1d(tord, self.d, kind='linear', bounds_error=False)
        self.d = f([date2num(tmp) for tmp in t])
        self.t = t
