"""
file/read.py
============

GSadjust code for importing relative-gravity meter data.
--------------------------------------------------------------------------------

This software is preliminary, provisional, and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.
"""

import datetime as dt
import logging

from matplotlib.dates import date2num

from ..data import ChannelList, Datum
from ..utils import index_or_none


class InvalidMeterException(Exception):
    """
    A loader for an unknown meter type was requested.
    """


def file_reader(meter_type, fh):
    """
    Generic file reader, which given a meter type and a file handle will
    pass loading off to an appropriate reader, returning the result.

    Parameters
    ----------
    meter_type : {'csv', 'CG5', 'Burris', 'CG6', CG6Tsoft'}
        The type of data tp load
    fh : TextIOWrapper
        Open file handle to the file to load.

    Returns
    -------
    ChannelList

    """
    read_fn = {
        "csv": read_csv,
        "CG5": read_cg5,
        "Burris": read_burris,
        "CG6": read_cg6,
        "CG6Tsoft": read_cg6tsoft,
    }.get(meter_type)
    if read_fn is None:
        raise InvalidMeterException

    # We have a valid reader, return the contents of the file.
    data = read_fn(fh)
    return data


def read_csv(fh):
    """Read arbitraty data in a csv file.

    TODO: Make more general

    CSV columns (must be in this order:
      Station name
      Latitude
      Longitude
      ALtitude
      Gravity
      Standard deviation
      Tilt X
      Tilt Y
      Temperature
      Tide
      Duration
      Rejected
      Time
      Decimal Date.time (ignored)
      Terrain correction (ignored)
      Date

    Representative line:
    TEST,26.3724,-112.4124,110.0,6134.455,0.016,-1,0.2,-1.49,-0.014,60,0,...
    14:11:20,42572.59026,0,8/22/2016
    """
    i = 0
    meter, oper = None, None
    all_survey_data = ChannelList()

    # Skip 1 header line
    _ = fh.readline()
    for i, orig_line in enumerate(fh, 1):
        try:
            line = orig_line.strip()
            # Skip blank and comment lines
            if (not line) or (line[0] == "%"):
                continue

            # parse string line first with respect to '/' characters (used in the
            # date format), then with ':' (used for the time display), eventually
            # with the classic ' '
            vals_temp1 = line.split("/")
            vals_temp2 = vals_temp1[0].split(":")
            vals_temp3 = vals_temp2[0].split(",")
            vals_temp4 = vals_temp2[2].split(",")

            # fill object properties:
            all_survey_data.station.append(vals_temp3[0].strip())
            all_survey_data.lat.append(float(vals_temp3[1]))
            all_survey_data.long.append(float(vals_temp3[2]))
            all_survey_data.elev.append(float(vals_temp3[3]))
            all_survey_data.raw_grav.append(
                float(vals_temp3[4]) * 1000.0 - float(vals_temp3[9]) * 1000.0
            )  # convert to microGal; remove etc
            all_survey_data.tare.append(0)
            all_survey_data.sd.append(float(vals_temp3[5]) * 1000.0)
            all_survey_data.tiltx.append(float(vals_temp3[6]))
            all_survey_data.tilty.append(float(vals_temp3[7]))
            all_survey_data.temp.append(float(vals_temp3[8]))
            all_survey_data.etc.append(float(vals_temp3[9]) * 1000.0)
            all_survey_data.meter_etc.append(float(vals_temp3[9]) * 1000.0)
            all_survey_data.dur.append(int(vals_temp3[10]))
            all_survey_data.rej.append(int(vals_temp3[11]))
            all_survey_data.t.append(
                date2num(
                    dt.datetime(
                        int(vals_temp1[-1]),
                        int(vals_temp4[-1]),
                        int(vals_temp1[-2]),
                        int(vals_temp3[-1]),
                        int(vals_temp2[-2]),
                        int(vals_temp4[0]),
                    )
                )
            )

            all_survey_data.meter.append(meter or "-999")
            all_survey_data.oper.append(oper or "-999")

            all_survey_data.keepdata.append(1)

        except (IndexError, ValueError) as e:
            logging.exception("Error loading CSV file at line %d", i)
            logging.info("LINE: %s", line)

            e.i = i
            e.line = orig_line
            raise e
    all_survey_data.meter_type = "csv"
    return all_survey_data


def read_cg5(fh):
    """
    Read CG5 formatted file, from given open file handle.

    Parameters
    ----------
    fh : TextIOWrapper
       open file handle, of type CG5

    Returns
    -------
    ChannelList

    """
    meter, oper = None, None
    all_survey_data = ChannelList()

    for i, orig_line in enumerate(fh, 1):
        try:
            # Clean line
            line = orig_line.strip()

            # Skip blank and comment lines
            if (not line) or (line[0] == "L"):
                continue

            # Header line; look for useful information
            if line[0] == "/":
                vals_temp = line.split()
                if len(vals_temp) > 1:
                    if vals_temp[1] == "Instrument":
                        meter = vals_temp[-1]
                    if vals_temp[1] == "Operator:":
                        oper = vals_temp[-1]
                continue

            # parse string line first with respect to '/' characters (used in the date format),
            # then with ':' (used for the time display), eventually with the classic ' '
            vals_temp1 = line.split("/")
            vals_temp2 = vals_temp1[0].split(":")
            vals_temp3 = vals_temp2[0].split()
            vals_temp4 = vals_temp2[2].split()

            # fill object properties:
            all_survey_data.line.append(float(vals_temp3[0]))
            s = vals_temp3[1].replace(".0000000", "")
            all_survey_data.station.append(s.strip())
            all_survey_data.elev.append(float(vals_temp3[2]))
            all_survey_data.raw_grav.append(
                float(vals_temp3[3]) * 1000.0 - float(vals_temp3[8]) * 1000.0
            )  # convert to microGal; remove tide correction
            all_survey_data.tare.append(0)
            all_survey_data.sd.append(float(vals_temp3[4]) * 1000.0)
            all_survey_data.tiltx.append(float(vals_temp3[5]))
            all_survey_data.tilty.append(float(vals_temp3[6]))
            all_survey_data.temp.append(float(vals_temp3[7]))
            all_survey_data.etc.append(float(vals_temp3[8]) * 1000.0)
            all_survey_data.meter_etc.append(float(vals_temp3[8]) * 1000.0)
            all_survey_data.dur.append(int(vals_temp3[9]))
            all_survey_data.rej.append(int(vals_temp3[10]))
            all_survey_data.t.append(
                date2num(
                    dt.datetime(
                        int(vals_temp4[3]),
                        int(vals_temp1[1]),
                        int(vals_temp1[2]),
                        int(vals_temp3[11]),
                        int(vals_temp2[1]),
                        int(vals_temp4[0]),
                    )
                )
            )

            all_survey_data.meter.append(meter or "-999")
            all_survey_data.oper.append(oper or "-999")

            all_survey_data.keepdata.append(1)
        except (IndexError, ValueError) as e:
            logging.exception("Error loading CG5 file at line %d", i)
            logging.info("LINE: %s", line)

            e.i = i
            e.line = orig_line
            raise e
        except ValueError as e:
            e.i = i
            e.line = orig_line
            raise e
    all_survey_data.meter_type = "CG5"
    return all_survey_data


def read_burris(fh):
    """
    Read Burris formatted file, from given open file handle.

    Accepts comma or tab-separated files.

    Parameters
    ----------
    fh : TextIOWrapper
       open file handle

    Returns
    -------
    ChannelList

    """

    all_survey_data = ChannelList()

    for i, orig_line in enumerate(fh, 1):
        try:
            line = orig_line.strip()
            if line.find(",") != -1:
                vals_temp = line.split(",")
                if vals_temp[0] == "Station ID" or vals_temp[0] == "Station":
                    continue
            elif line.find("\t") != -1:
                vals_temp = line.split("\t")
            else:
                vals_temp = line.split()
            if vals_temp[0] == "Station ID" or vals_temp[0] == "Station":
                continue

            if len(vals_temp) == 15:  # no meter operator specified
                (
                    c_station,
                    c_meter,
                    c_date,
                    c_time,
                    c_grav,
                    c_dial,
                    c_feedback,
                    c_tide,
                    c_tilt,
                    _,
                    _,
                    c_height,
                    c_elev,
                    c_lat,
                    c_long,
                ) = range(
                    15
                )  # 0 - 15
                all_survey_data.oper.append("None")

            else:  # 16 values, includes meter operator.
                # Numbers are columns in the imported file
                (
                    c_station,
                    c_oper,
                    c_meter,
                    c_date,
                    c_time,
                    c_grav,
                    c_dial,
                    c_feedback,
                    c_tide,
                    c_tilt,
                    _,
                    _,
                    c_height,
                    c_elev,
                    c_lat,
                    c_long,
                ) = range(
                    16
                )  # 0 - 14

                all_survey_data.oper.append(vals_temp[c_oper])

            if line.find("/") != -1:
                date_temp = vals_temp[c_date].split("/")
            elif line.find("-") != -1:
                date_temp = vals_temp[c_date].split("-")
            else:
                date_temp = []
            if int(date_temp[2]) > 999:
                date_temp = [date_temp[2], date_temp[0], date_temp[1]]
            elif int(date_temp[0]) > 999:
                date_temp = [date_temp[0], date_temp[1], date_temp[2]]
            # Else raise date error

            time_temp = vals_temp[c_time].split(":")

            # fill object properties:
            all_survey_data.station.append(vals_temp[c_station].strip())
            all_survey_data.elev.append(float(vals_temp[c_elev]))
            all_survey_data.height.append(float(vals_temp[c_height]))
            all_survey_data.lat.append(float(vals_temp[c_lat]))
            all_survey_data.long.append(float(vals_temp[c_long]))
            # remove Earth tide correction; it's added in using the @grav property
            all_survey_data.raw_grav.append(
                float(vals_temp[c_grav]) * 1000.0 - float(vals_temp[c_tide]) * 1000.0
            )
            all_survey_data.tare.append(0)
            all_survey_data.etc.append(float(vals_temp[c_tide]) * 1000.0)
            all_survey_data.meter_etc.append(float(vals_temp[c_tide]) * 1000.0)
            all_survey_data.dial.append(float(vals_temp[c_dial]))
            all_survey_data.feedback.append(float(vals_temp[c_feedback]))
            all_survey_data.sd.append(-999)  # Burris doesn't ouput SD, tiltx, tilty
            all_survey_data.meter.append(vals_temp[c_meter])
            all_survey_data.tiltx.append(float(vals_temp[c_tilt]) * 1000.0)
            all_survey_data.tilty.append(0.0)
            all_survey_data.temp.append(0.0)
            all_survey_data.dur.append(5)
            all_survey_data.rej.append(5)
            all_survey_data.t.append(
                date2num(
                    dt.datetime(
                        int(date_temp[0]),
                        int(date_temp[1]),
                        int(date_temp[2]),
                        int(time_temp[0]),
                        int(time_temp[1]),
                        int(time_temp[2]),
                    )
                )
            )
            all_survey_data.keepdata.append(1)

        except (IndexError, ValueError) as e:
            logging.exception("Error loading Burris file at line %d", i)
            logging.info("LINE: %s", line)

            e.i = i
            e.line = orig_line
            raise e

    all_survey_data.meter_type = "Burris"
    return all_survey_data


def read_cg6(fh):
    """
    Read CG6 formatted file, from given open file handle.

    Parameters
    ----------
    fh : TextIOWrapper
       open file handle

    Returns
    -------
    ChannelList

    """
    meter, oper = None, None
    all_survey_data = ChannelList()

    for i, orig_line in enumerate(fh, 1):
        try:
            line = orig_line.strip()
            vals_temp = line.split("\t")
            if line[0] == "/":
                vals_temp = line.split()
                if len(vals_temp) > 1:
                    if vals_temp[1] == "Instrument":
                        meter = vals_temp[-1]
                    if vals_temp[1] == "Operator:":
                        oper = vals_temp[-1]
                continue
            # Numbers are columns in the imported file
            c_station, c_date, c_time, c_sd = 0, 1, 2, 5
            c_tiltx, c_tilty = 8, 9
            c_tide, c_tilt, c_temp = 11, 12, 13
            c_dur = 15
            c_grav, c_elev, c_lat, c_long = 3, 19, 17, 18

            date_temp = vals_temp[c_date].split("-")
            time_temp = vals_temp[c_time].split(":")

            # fill object properties:
            all_survey_data.line.append(0.0)
            all_survey_data.station.append(vals_temp[c_station].strip())
            all_survey_data.elev.append(float(vals_temp[c_elev]))
            all_survey_data.lat.append(float(vals_temp[c_lat]))
            all_survey_data.long.append(float(vals_temp[c_long]))
            all_survey_data.raw_grav.append(
                float(vals_temp[c_grav]) * 1000.0 - float(vals_temp[c_tide]) * 1000.0
            )
            all_survey_data.tare.append(0)
            all_survey_data.etc.append(float(vals_temp[c_tide]) * 1000.0)
            all_survey_data.meter_etc.append(float(vals_temp[c_tide]) * 1000.0)
            all_survey_data.sd.append(float(vals_temp[c_sd]) * 1000.0)
            all_survey_data.meter.append(meter)
            all_survey_data.tiltx.append(float(vals_temp[c_tiltx]) * 1000.0)
            all_survey_data.tilty.append(float(vals_temp[c_tilty]) * 1000.0)
            all_survey_data.temp.append(float(vals_temp[c_temp]) * 1000.0)
            all_survey_data.dur.append(int(vals_temp[c_dur]))
            all_survey_data.rej.append(5)
            all_survey_data.t.append(
                date2num(
                    dt.datetime(
                        int(date_temp[0]),
                        int(date_temp[1]),
                        int(date_temp[2]),
                        int(time_temp[0]),
                        int(time_temp[1]),
                        int(time_temp[2]),
                    )
                )
            )

            all_survey_data.meter.append(meter or "-999")
            all_survey_data.oper.append(oper or "-999")

            all_survey_data.keepdata.append(1)

        except (IndexError, ValueError) as e:
            logging.exception("Error loading CG6 file at line %d", i)
            logging.info("LINE: %s", line)
            e.i = i
            e.line = orig_line
            raise e
    all_survey_data.meter_type = "CG6"
    return all_survey_data


def read_cg6tsoft(fh):
    """
    Read CG6 TSoft formatted file, from given open file handle.

    Parameters
    ----------
    fh : TextIOWrapper
       open file handle, of type CG5

    Returns
    -------
    ChannelList

    """

    meter, oper = None, None
    all_survey_data = ChannelList()
    station_name = None
    for i, orig_line in enumerate(fh, 1):
        try:
            line = orig_line.strip()
            vals_temp = line.split()

            if line[0] == "/":
                if len(vals_temp) > 1:
                    if vals_temp[1] == "Instrument":
                        meter = vals_temp[-1]
                    if vals_temp[1] == "Operator:":
                        oper = vals_temp[-1]
                    if vals_temp[1] == "Station:":
                        station_name = vals_temp[-1]
                continue

            # Numbers are columns in the imported file
            c_tiltx, c_tilty = 12, 13
            c_tide, c_tilt, c_temp = 19, 18, 14
            c_lat, c_long, c_elev = 7, 8, 9
            c_grav = 11  # CorrGravity channel

            # fill object properties:
            if station_name:
                all_survey_data.station.append(station_name)
                all_survey_data.elev.append(float(vals_temp[c_elev]))
                all_survey_data.lat.append(float(vals_temp[c_lat]))
                all_survey_data.long.append(float(vals_temp[c_long]))
                all_survey_data.raw_grav.append(float(vals_temp[c_grav]) * 1000.0)
                all_survey_data.tare.append(0)
                all_survey_data.etc.append(float(vals_temp[c_tide]) * 1000.0)
                all_survey_data.meter_etc.append(float(vals_temp[c_tide]) * 1000.0)
                all_survey_data.sd.append(
                    -999
                )  # SD not exported in Tsoft format?? It is in regular format
                all_survey_data.meter.append(meter)
                all_survey_data.tiltx.append(float(vals_temp[c_tiltx]) * 1000.0)
                all_survey_data.tilty.append(float(vals_temp[c_tilty]) * 1000.0)
                all_survey_data.temp.append(float(vals_temp[c_temp]) * 1000.0)
                # cols 0, 1, 2, 3, 4, 5 = year, month, day, hour, minute, second
                temp_date_ints = (int(i) for i in vals_temp[:6])
                all_survey_data.t.append(date2num(dt.datetime(*temp_date_ints)))

                all_survey_data.meter.append(meter or "-999")
                all_survey_data.oper.append(oper or "-999")

                all_survey_data.keepdata.append(1)
                all_survey_data.dur.append(-999)
                all_survey_data.rej.append(-999)

        except (IndexError, ValueError) as e:
            logging.exception("Error loading CG6TSoft file %s, at line %d", fname, i)
            logging.info("LINE: %s", line)
            e.i = i
            e.line = orig_line
            raise e

    if all_survey_data.raw_grav:
        all_survey_data.meter_type = "CG6Tsoft"
        return all_survey_data

    else:
        raise ValueError


def import_abs_g_complete(fname):
    """
    Imports absolute gravity data from comma-separated file output by fg5_parse.py.
    Output will be displayed in datum model on Adjust tab.

    Parameters
    ----------
    fname : str
        filename to open (complete path)

    Returns
    -------
    list of Datum objects

    """

    datums = list()

    with open(fname, "r") as fh:
        # Read header line
        line = fh.readline()
        parts = [p.strip() for p in line.split("\t")]

        g_idx = index_or_none(parts, "Gravity")
        n_idx = index_or_none(parts, "Station Name")
        s_idx = index_or_none(parts, "Set Scatter")
        d_idx = index_or_none(parts, "Date")
        th_idx = index_or_none(parts, "Transfer Height")
        gr_idx = index_or_none(parts, "Gradient")

        for line in fh:
            try:
                if all([g_idx, n_idx, s_idx, d_idx, th_idx]):
                    parts = line.split("\t")
                    # Gradient can be optional
                    if gr_idx:
                        gr = float(parts[gr_idx])
                    else:
                        gr = -3.0

                    datum = Datum(
                        parts[n_idx],
                        g=float(parts[g_idx]),
                        sd=float(parts[s_idx]),
                        date=parts[d_idx],
                        meas_height=float(parts[th_idx]),
                        gradient=gr
                    )
                    datums.append(datum)
            except ValueError:
                logging.exception(
                    "Error loading absolute gravity data from %s", fname
                )
                return []
    return datums


def import_abs_g_simple(fname):
    """
    Imports absolute data from a three-column space-separated file:
      station g stdev.

    Output will be displayed in datum model on Adjust tab.

    Parameters
    ----------
    fname : str
        filename to open (complete path)

    Returns
    -------
    List of Datum objects
    """
    datums = list()
    with open(fname, "r") as fh:
        line = fh.readline()
        while True:
            if not line:
                break
            parts = line.split(" ")
            datum = Datum(parts[0], float(parts[1]), float(parts[2]))
            datums.append(datum)
            line = fh.readline()
    return datums
