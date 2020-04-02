#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""
data_export.py
===============

Data export code for GSadjust.
--------------------------------------------------------------------------------------------------------------------

This software is preliminary, provisional, and is subject to revision. It is being provided to meet the need for
timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related
material nor shall the fact of release constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or
unauthorized use of the software.
"""
import logging
import time
import os
from PyQt5 import QtCore
import csv

from data_analysis import compute_gravity_change

def export_metadata(obsTreeModel, data_path):
    """
    Write metadata text to file. Useful for USGS data releases.
    """
    fn = os.path.join(data_path, 'GSadjust_MetadataText_' + time.strftime("%Y%m%d-%H%M") + '.txt')
    results_written, first = False, False
    fmt = '%Y-%m-%d'
    with open(fn, 'w') as fid:
        for i in range(obsTreeModel.rowCount()):
            survey = obsTreeModel.invisibleRootItem().child(i)
            if survey.adjustment.adjustmentresults.text:  # check that there are results
                if first:
                    fid.write('Attribute accuracy is evaluated from the least-squares network adjustment results. ')
                    first = False
                results_written = True
                fid.write('For the {} survey, the minimum and maximum gravity-difference residuals were {:0.1f} '
                          'and {:0.1f} '.format(survey.name,
                                                survey.adjustment.adjustmentresults.min_dg_residual,
                                                survey.adjustment.adjustmentresults.max_dg_residual))
                fid.write(
                    'microGal, respectively. The minimum and maximum datum (absolute-gravity station) residuals ')
                fid.write(
                    'were {:0.1f} and {:0.1f} microGal, respectively. '.format(
                        survey.adjustment.adjustmentresults.min_datum_residual,
                        survey.adjustment.adjustmentresults.max_datum_residual))
                fid.write('The average standard deviation of the adjusted gravity values at each station')
                fid.write(
                    ' (derived from the network adjustment) was {:0.1f} microGal. '.format(
                        survey.adjustment.adjustmentresults.avg_stdev))
                # TODO: account for instance of 1 outlier ('1 was removed')
                fid.write('{} out of {} possible gravity differences were used in the adjustment ({} were removed '
                          .format(survey.adjustment.adjustmentresults.n_deltas -
                                  survey.adjustment.adjustmentresults.n_deltas_notused,
                                  survey.adjustment.adjustmentresults.n_deltas,
                                  survey.adjustment.adjustmentresults.n_deltas_notused))
                fid.write('as outliers). {} out of {} possible datum observations were used. '.format(
                    survey.adjustment.adjustmentresults.n_datums -
                    survey.adjustment.adjustmentresults.n_datums_notused,
                    survey.adjustment.adjustmentresults.n_datums))
        logging.info('Metadata text written to file')

    if not results_written:
        from gui_objects import show_message
        msg = show_message('No network adjustment results', 'Write error')


def export_summary(obsTreeModel, data_path):
    """
    Write summary of procesing to text file. Can be used to reproduce results.

    Parameters
    ----------
    obsTreeModel

    Returns
    -------

    """
    fn = os.path.join(data_path, 'GSadjust_Summary_' + time.strftime("%Y%m%d-%H%M") + '.txt')
    # Write header info
    with open(fn, 'w') as fid:
        fid.write('# GSadjust processing summary, {}\n#\n'.format(time.strftime("%Y%m%d-%H%M")))
        fid.write('# Station data\n')
        for i in range(obsTreeModel.rowCount()):
            survey = obsTreeModel.invisibleRootItem().child(i)
            for ii in range(survey.rowCount()):
                loop = survey.child(ii)
                fid.write('Survey {}, Loop {}\n'.format(survey.name, loop.name))
                for iii in range(loop.rowCount()):
                    station = loop.child(iii)
                    fid.write(station.summary_str)
                    for sample_str in station.iter_samples():
                        fid.write(sample_str)
        fid.write('# Loop data\n')
        for i in range(obsTreeModel.rowCount()):
            survey = obsTreeModel.invisibleRootItem().child(i)
            for ii in range(survey.rowCount()):
                loop = survey.child(ii)
                fid.write('Survey {}, Loop {}\n'.format(survey.name, loop.name))
                fid.write(str(loop))
                for delta in loop.delta_model.deltas:
                    fid.write('{}\n'.format(delta))
        fid.write('# Adjustment data\n')
        for i in range(obsTreeModel.rowCount()):
            survey = obsTreeModel.invisibleRootItem().child(i)
            fid.write('# Adjustment options, survey: {}\n'.format(survey.name))
            fid.write(str(survey.adjustment.adjustmentoptions))
            fid.write('# Deltas in the adjustment, survey: {}\n'.format(survey.name))
            for delta in survey.adjustment.deltas:
                fid.write('{}\n'.format(delta))
            fid.write('# Datums in the adjustment, survey: {}\n'.format(survey.name))
            for datum in survey.adjustment.datums:
                fid.write('{}\n'.format(datum))
            fid.write('# Adjustment results, survey: {}\n'.format(survey.name))
            fid.write(str(survey.adjustment.adjustmentresults))
            fid.write('# Adjusted station values and std. dev., survey: {}\n'.format(survey.name))
            for ii in range(survey.results_model.rowCount()):
                adj_sta = survey.results_model.data(survey.results_model.index(ii, 0), role=QtCore.Qt.UserRole)
                fid.write('{}\n'.format(str(adj_sta)))

def export_data(obstreemodel, data_path):
    """
    Export gravity change table to csv file
    """
    fn = os.path.join(data_path, 'GSadjust_TabularData_' + time.strftime("%Y%m%d-%H%M") + '.csv')
    table = compute_gravity_change(obstreemodel, table_type='full')

    with open(fn, 'w', newline='\n') as fid:
        wr = csv.writer(fid)
        wr.writerow(table[0])
        for row in table[1]:
            wr.writerow(row)