#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
data_correction.py
==================

GSadjust module for applying a time correction to gravity data.
--------------------------------------------------------------------------------------------------------------------


This software is preliminary, provisional, and is subject to revision. It is being provided to meet the need for
timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related
material nor shall the fact of release constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or
unauthorized use of the software.
"""

import logging
import datetime as dt


def time_correction(obstreemodel, correction_type, t_offset, survey_index, loop_index, station_indexes):
    if correction_type == 'all':
        for i in range(obstreemodel.invisibleRootItem().rowCount()):
            obstreesurvey = obstreemodel.invisibleRootItem().child(i)
            for ii in range(obstreesurvey.rowCount()):
                obstreeloop = obstreesurvey.child(ii)
                for iii in range(obstreeloop.rowCount()):
                    obstreestation = obstreeloop.child(iii)
                    obstreestation.t = [t + dt.timedelta(t_offset / 1440, 0, 0) for t in obstreestation.t]
                    logging.info("{} minute offset added to station {}".format(t_offset,
                                                                               obstreestation.display_name()))
    elif correction_type == 'survey':
        obstreesurvey = obstreemodel.itemFromIndex(survey_index)
        for ii in range(obstreesurvey.rowCount()):
            obstreeloop = obstreesurvey.child(ii)
            for iii in range(obstreeloop.rowCount()):
                obstreestation = obstreeloop.child(iii)
                obstreestation.t = [t + dt.timedelta(t_offset / 1440, 0, 0) for t in obstreestation.t]
                logging.info(
                    "{} minute offset added to station {}".format(t_offset, obstreestation.display_name()))
    elif correction_type == 'loop':
        obstreeloop = obstreemodel.itemFromIndex(loop_index)
        for iii in range(obstreeloop.rowCount()):
            obstreestation = obstreeloop.child(iii)
            obstreestation.t = [t + dt.timedelta(t_offset / 1440, 0, 0) for t in obstreestation.t]
            logging.info(
                "{} minute offset added to station {}".format(t_offset, obstreestation.display_name()))
    elif correction_type == 'station':
        # Because each tree item has three columns, len(indexes) equals the number of items selected * 3.
        # The next line takes every 3rd index.
        indexes = station_indexes[0::3]
        for index in indexes:
            obstreestation = obstreemodel.itemFromIndex(index)
            obstreestation.t = [t + dt.timedelta(t_offset / 1440, 0, 0) for t in obstreestation.t]
            logging.info("{} minute offset added to station {}".format(t_offset,
                                                                       obstreestation.display_name()))
