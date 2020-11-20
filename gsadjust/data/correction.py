"""
data/correction.py
==================

GSadjust module for applying a time correction to gravity data.
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
import logging


def time_correction(
    obstreemodel, correction_type, t_offset, survey_index, loop_index, station_indexes
):
    """Perform the specified time correction.

    Parameters
    ----------
    obstreemodel : ObsTreeModel
        Data structure from MainProg.
    correction_type : {"all", "survey", "loop", "station}
        Str indicating where to apply time correction.
    t_offset : float
        Time offset, in minutes. Can be a decimal value.
    survey_index : QModelIndex
        index of currently selected survey (only used if applying to all surveys)
    loop_index : QModelIndex
        index of currently selected loop (only used if applying to all surveys)
    station_indexes : list
        List of selected (highlighted) station indexes.

    Returns
    -------
    None

    """
    if correction_type == "all":
        _time_correction_all(obstreemodel, t_offset)

    elif correction_type == "survey":
        obstreesurvey = obstreemodel.itemFromIndex(survey_index)
        _time_correction_survey(obstreesurvey, t_offset)

    elif correction_type == "loop":
        obstreeloop = obstreemodel.itemFromIndex(loop_index)
        _time_correction_loop(obstreeloop, t_offset)

    elif correction_type == "station":
        # Because each tree item has three columns, len(indexes) equals the number
        # of items selected * 3. The next line takes every 3rd index.
        indexes = station_indexes[0::3]
        for index in indexes:
            obstreestation = obstreemodel.itemFromIndex(index)
            obstreestation.t = _time_correction_station(obstreestation, t_offset)


def _time_correction_all(obstreemodel, t_offset):
    for i in range(obstreemodel.invisibleRootItem().rowCount()):
        obstreesurvey = obstreemodel.invisibleRootItem().child(i)
        _time_correction_survey(obstreesurvey, t_offset)


def _time_correction_survey(obstreesurvey, t_offset):
    for ii in range(obstreesurvey.rowCount()):
        obstreeloop = obstreesurvey.child(ii)
        _time_correction_loop(obstreeloop, t_offset)


def _time_correction_loop(obstreeloop, t_offset):
    for iii in range(obstreeloop.rowCount()):
        obstreestation = obstreeloop.child(iii)
        obstreestation.t = _time_correction_station(obstreestation, t_offset)


def _time_correction_station(obstreestation, t_offset):
    logging.info(
        "%s minute offset added to station %s",
        t_offset,
        obstreestation.display_name(),
    )
    return [t + t_offset / 1440 for t in obstreestation.t]
