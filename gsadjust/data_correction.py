#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
data_correction.py
==================
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
                                                                               obstreestation.display_name))
    elif correction_type == 'survey':
        obstreesurvey = obstreemodel.itemFromIndex(survey_index)
        for ii in range(obstreesurvey.rowCount()):
            obstreeloop = obstreesurvey.child(ii)
            for iii in range(obstreeloop.rowCount()):
                obstreestation = obstreeloop.child(iii)
                obstreestation.t = [t + dt.timedelta(t_offset / 1440, 0, 0) for t in obstreestation.t]
                logging.info(
                    "{} minute offset added to station {}".format(t_offset, obstreestation.display_name))
    elif correction_type == 'loop':
        obstreeloop = obstreemodel.itemFromIndex(loop_index)
        for iii in range(obstreeloop.rowCount()):
            obstreestation = obstreeloop.child(iii)
            obstreestation.t = [t + dt.timedelta(t_offset / 1440, 0, 0) for t in obstreestation.t]
            logging.info(
                "{} minute offset added to station {}".format(t_offset, obstreestation.display_name))
    elif correction_type == 'station':
        # Because each tree item has three columns, len(indexes) equals the number of items selected * 3.
        # The next line takes every 3rd index.
        indexes = station_indexes[0::3]
        for index in indexes:
            obstreestation = obstreemodel.itemFromIndex(index)
            obstreestation.t = [t + dt.timedelta(t_offset / 1440, 0, 0) for t in obstreestation.t]
            logging.info("{} minute offset added to station {}".format(t_offset,
                                                                       obstreestation.display_name))
