"""
utils.py
========

GSadjust utility functions
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


def index_or_none(lst, i):
    if i not in lst:
        return None
    return lst.index(i)


def init_cal_coeff_dict(obstreemodel):
    """
    Initiate dict for storing meter calibration coefficients.

    Parameters
    ----------
    obstreemodel : ObsTreeModel

    Returns
    -------
    dict
        key: Meter (str), value: float

    """
    try:
        meter_list = {}
        for i in range(obstreemodel.invisibleRootItem().rowCount()):
            survey = obstreemodel.invisibleRootItem().child(i)
            for ii in range(survey.rowCount()):
                loop = survey.child(ii)
                if loop.meter not in meter_list:
                    meter_list[loop.meter] = 1.000
        return meter_list
    except Exception:
        return None


def init_station_coords_dict(obstreemodel):
    """
    Stores a single set of coordinates for each station with the obsTreeModel
    object. The coordinates of the last
    Station in the Survey > Loop > Station hierarchy will be used.
    """
    station_coords = dict()
    for i in range(obstreemodel.invisibleRootItem().rowCount()):
        survey = obstreemodel.invisibleRootItem().child(i)
        for ii in range(survey.rowCount()):
            loop = survey.child(ii)
            for iii in range(loop.rowCount()):
                station = loop.child(iii)
                try:
                    station_coords[station.station_name] = (
                        station.long[0],
                        station.lat[0],
                        station.elev[0],
                    )
                except Exception:
                    station_coords[station.station_name] = (0, 0, 0)
    return station_coords
