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
import datetime as dt
import json
import logging
import os

import jsons
import numpy as np
from matplotlib.dates import date2num, num2date
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt, QVariant

from ..data import (AdjustedStation, Adjustment, AdjustmentOptions,
                    AdjustmentResults, Datum, DeltaList, DeltaNormal, Tare)
from ..data.analysis import InversionError, numpy_inversion
from ..gui.messages import show_message
from .base import ObsTreeItemBase
from .loop import ObsTreeLoop

# Constants for column headers
STATION_NAME, STATION_DATETIME, STATION_MEAN = range(3)
LOOP_NAME = 0
SURVEY_NAME = 0


class ObsTreeSurvey(ObsTreeItemBase):
    """
    A survey object is a little different than Loop or Station objects in that
    it doesn't have a corresponding data object. All relevant info is in
    this object.
    """

    def __init__(self, name):
        super(ObsTreeSurvey, self).__init__()
        # Both surveys and loops have delta models. The survey delta_model holds the observations used in an
        # adjustment, and may have deltas from one or more loops.
        self.name = name

        # Store related data.
        self.deltas = []
        self.datums = []
        self.results = []

        self.adjustment = Adjustment()

    def __str__(self):
        return self.name

    @property
    def display_column_map(self):
        return {
            # Column name map
            # index: name
            SURVEY_NAME: (str, self.name),
            1: (str, ''),
            2: (str, ''),
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
            except Exception:
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
            d = Datum(datum['station'])
            d.__dict__.update(datum)
            d.residual = -999
            temp.datums.insert(0, d)

        ao = AdjustmentOptions()
        ao.__dict__.update(data['adjoptions'])
        if not hasattr(ao, 'use_sigma_prefactor') and hasattr(ao, 'sigma_factor'):
            ao.sigma_prefactor = ao.sigma_factor
            ao.use_sigma_prefactor = ao.use_sigma_factor
            ao.sigma_postfactor = 1.0
            ao.use_sigma_postfactor = False
        temp.adjustment.adjustmentoptions = ao

        if 'checked' in data:
            temp.setCheckState(data['checked'])
        return temp

    def to_json(self):
        loops, datums = [], []
        # Remove ObsTreeStation objects from deltas in the survey delta_model (which is different than the individual
        # loop delta_models; those are ignored and not save (they're recreated when the workspace is loaded).
        for i in range(self.rowCount()):
            loops.append(self.child(i).to_json())
        return {
            'loops': loops,
            'deltas': jsons.dump(self.deltas),
            'datums': jsons.dump(self.datums),
            'checked': self.checkState(),
            'name': self.name,
            'adjoptions': jsons.dump(self.adjustment.adjustmentoptions),
        }

    @property
    def tooltip(self):
        return (
            'Survey: {}\n'
            'Meters: {}\n'
            'Number of loops: {}'.format(self.name, self.unique_meters, self.loop_count)
        )

    @property
    def unique_meters(self):
        """
        :return: List of unique gravity-meter IDs (serial numbers) used in the survey
        """
        meters = []
        for station in self.iter_stations():
            meters.append(
                station.meter[0]
            )  # Get the first entry; Assume meter number can't change at a station
        return list(set(meters))

    def loops(self):
        return [self.child(i) for i in range(self.rowCount())]

    def loops_with_deltas(self):
        return list(set([l.loop for l in self.deltas]))

    # @property
    # def datums(self):
    #     return self.datums[::]

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
        return len(self.loop_names)

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

    def populate(self, survey_data, name='0', source='NA'):
        """
        Called from open_raw_data. Loads all survey_data into a single loop.
        :param survey_data: data returned from read_raw_data_file
        :param name: Survey name
        :return: None
        """
        obstreeloop = ObsTreeLoop(name)
        obstreeloop.populate(survey_data)
        obstreeloop.source = source
        self.appendRow(
            [obstreeloop, QtGui.QStandardItem('0'), QtGui.QStandardItem('0')]
        )
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

    def run_inversion(
        self, adj_type='PyLSQ', write_out_files='n', output_root_dir='./',
    ):
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

        if self.data(role=Qt.CheckStateRole) == 2:
            deltas = []
            datums = []
            adjustmentresults = AdjustmentResults()

            self.adjustment.adjustmentoptions.adj_type = adj_type

            # Form datasets for adjustment from selected table items, count number of deltas and datums (used when
            # writing metadata about adjustment).
            specify_cal_coeff = False
            cal_dic = None
            # To deal with old files
            if hasattr(self.adjustment.adjustmentoptions, 'specify_cal_coeff'):
                if self.adjustment.adjustmentoptions.specify_cal_coeff:
                    specify_cal_coeff = True
                    cal_dic = self.adjustment.adjustmentoptions.meter_cal_dict
            for delta in self.deltas:
                if delta.checked:  # Checkbox checked
                    if specify_cal_coeff:
                        delta.cal_coeff = cal_dic[delta.meter]
                    else:
                        delta.cal_coeff = 1
                    if delta.type == 'normal':
                        try:
                            # Don't include deltas in which one of the stations has been unchecked
                            if (
                                delta.station1.data(role=Qt.CheckStateRole) == 2
                                and delta.station2.data(role=Qt.CheckStateRole) == 2
                            ):
                                deltas.append(delta)
                                adjustmentresults.n_deltas += 1
                            else:
                                adjustmentresults.n_deltas_notused += 1
                        except:
                            self.msg = show_message(
                                "Delta station not found. Was it deleted?",
                                "GSadjust error",
                            )
                    else:
                        deltas.append(delta)
                        adjustmentresults.n_deltas += 1
                else:
                    adjustmentresults.n_deltas_notused += 1

            for datum in self.datums:
                if datum.checked:  # Checkbox checked
                    # account for vertical gradient
                    datum = copy.copy(datum)
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
                self.results = []
                if not self.adjustment.datums:  # empty.
                    self.msg = show_message(
                        "Survey {}: At least one datum must be specified".format(
                            self.name
                        ),
                        "Inversion error",
                    )
                    return
                if not self.adjustment.deltas:  # empty.
                    self.msg = show_message(
                        "Survey {}: At least one relative-gravity difference must be specified".format(
                            self.name
                        ),
                        "Inversion error",
                    )
                    return
                if adj_type == 'PyLSQ':
                    logging.info('Numpy inversion, Survey: {}'.format(self.name))
                    self.start_numpy_inversion(output_root_dir, write_out_files)
                elif adj_type == 'Gravnet':
                    logging.info('Gravnet inversion, Survey: {}'.format(self.name))
                    self.gravnet_inversion()
            except ZeroDivisionError:
                self.msg = show_message(
                    'Unable to adjust. Are there standard deviations that are zero or very small?',
                    'ZeroDivisionError',
                )
            except KeyError as e:
                self.msg = show_message(
                    'Station error\nSurvey: {}\nStation: {}'.format(
                        self.name, e.args[0]
                    ),
                    'KeyError',
                )
            except Exception as e:
                logging.exception(e, exc_info=True)
                self.msg = show_message("Inversion error", "GSadjust")

    def gravnet_inversion(self):
        """
        Writes input files for Gravnet.exe, runs the executable, and reads in the results
        """
        # TODO: method could probably be made static

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
                self.msg = show_message(
                    'Gravnet.exe not found, aborting', 'Inversion error'
                )
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
                'It appears that more than one polynomial degree was specified for different loops for the network, '
                'or that some loops are not using the adjustment drift option. When using Gravnet, all loops must '
                'have the same degree drift model. Aborting.',
                'Inversion error',
            )
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
                'Inversion warning',
            )
        # Write delta-g observation file
        dg_file = self.name + '_dg.obs'
        with open(dg_file, 'w') as fid:
            fid.write(
                'start, end, difference (mgal), mjd (from), mjd (to),'
                + ' reading (from, CU), reading (to, CU), standard deviation (mgal)\n'
            )
            # Gravnet station names are limited to 6 characters; units are mGal
            for delta in self.adjustment.deltas:
                fid.write(
                    '{} {} {:0.6f} {} {} {:0.6f} {} {:0.6f}\n'.format(
                        delta.sta1[:6],
                        delta.sta2[:6],
                        delta.dg / 1000.0 * delta.cal_coeff,
                        delta.sta1_t,
                        delta.sta2_t,
                        delta.dg / 1000,
                        '0',
                        delta.sd_for_adjustment / 1000.0,
                    )
                )
        # Write absolute-g (aka datum, aka fix) file
        fix_file = self.name + '_fix.txt'
        with open(fix_file, 'w') as fid:
            for datum in self.adjustment.datums:
                fid.write(
                    '{} {:0.6f} {:0.3f}\n'.format(
                        datum.station[:6],
                        float(datum.g) / 1000.0,
                        float(datum.sd) / 1000.0,
                    )
                )

        # Check if calibration coefficient is calculated
        cal_dic = {}
        if self.adjustment.adjustmentoptions.cal_coeff:
            if len(self.unique_meters) > 1:
                self.msg = show_message(
                    "It appears more than one meter was used on the survey. Gravnet calculates a "
                    + "calibration coefficient for a single meter only. Use Numpy inversion "
                    + "to calculate meter-specific calibration coefficients",
                    "Inversion warning",
                )
            if len(self.adjustment.datums) == 1:
                self.msg = show_message(
                    "Two or more datum observations are required to calculate a calibration coefficient."
                    + " Aborting.",
                    "Inversion warning",
                )
                return
            # Run gravnet with calibration coefficient
            os.system(
                'gravnet -D{} -N{} -M2 -C1 {} -F{}'.format(
                    dg_file, self.name, drift_term, fix_file
                )
            )
            with open(self.name + '.met', 'r') as fid:
                symbol = ''
                while symbol != 'Y_l':
                    line = fid.readline().split()
                    symbol = line[0]
                meter_calib_params = fid.readline().split()
                cal_dic[self.unique_meters[0]] = (
                    1 + float(meter_calib_params[0]),
                    float(meter_calib_params[1]),
                )
        else:
            # Run gravnet without calibration coefficient
            os.system(
                'gravnet -D{} -N{} -M2 {} -F{}'.format(
                    dg_file, self.name, drift_term, fix_file
                )
            )

        # Read drift coefficients
        if drift_term != '':
            meter_drift_params = []
            with open(self.name + '.met', 'r') as fid:
                symbol = ''
                while (
                    symbol != 'Coefficients'
                ):  # Indicates 'Coefficients of drift' in .met file
                    line = fid.readline().split()
                    symbol = line[0]
                order = int(line[-1])
                for i in range(order):
                    meter_drift_params.append(fid.readline().split())

        self.results = []

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
                    self.results.insert(0, AdjustedStation(sta, g, sd))
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
                self.adjustment.adjustmentresults.text.append(fh.strip())
            if self.adjustment.adjustmentoptions.cal_coeff:
                self.adjustment.adjustmentresults.text.append(
                    "Gravimeter calibration coefficient: {:.6} ± {:.6}".format(
                        1 + float(meter_calib_params[0]), meter_calib_params[1]
                    )
                )
            elif self.adjustment.adjustmentoptions.specify_cal_coeff:
                for k, v in self.adjustment.adjustmentoptions.meter_cal_dict.items():
                    self.adjustment.adjustmentresults.text.append(
                        "Specified calibration coefficient for meter {}: {:.6f}".format(
                            k, v
                        )
                    )
            if drift_term != '':
                self.adjustment.adjustmentresults.text.append(
                    "Gravimeter drift coefficient(s):"
                )
                for coeffs in meter_drift_params:
                    self.adjustment.adjustmentresults.text.append(
                        "{:.3} ± {:.3}".format(coeffs[0], coeffs[1])
                    )

        if dir_changed1:
            os.chdir('..\\gsadjust')
        elif dir_changed2:
            os.chdir('..')

    def start_numpy_inversion(self, output_root_dir, write_out_files='n'):
        """
        Least-squares network adjustment using numpy

        Parameters
        write_out_files:    (y/n): write output files for drift adjustment
                            (similar to MCGravi output files)
        output_root_dir     directory for writing output files
        """

        self.adjustment.adjustmentresults.text = []

        # sta_dic_LS is a dictionary, key: station name, value: column for A matrix
        self.adjustment.sta_dic_ls = self.get_station_indices()
        self.adjustment.nloops = self.loop_count

        # number of unknowns:
        # TODO: temperature drift function
        drift_time = 0  # self.adjustment.adjustmentoptions.drift_time
        # if self.adjustment.adjustmentoptions.use_model_temp:
        #     self.adjustment.drift_temp = self.adjustment.adjustmentoptions.model_temp_degree
        # else:
        #     self.adjustment.drift_temp = 0

        if self.adjustment.adjustmentoptions.cal_coeff:
            if len(self.adjustment.datums) == 1:
                self.msg = show_message(
                    "Two or more datum observations are required to calculate a calibration coefficient."
                    + " Aborting.",
                    "Inversion warning",
                )
                return
            n_meters = len(self.unique_meters)
            self.adjustment.meter_dic = dict(
                zip(self.unique_meters, range(n_meters + 1))
            )
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
        self.adjustment.adjustmentresults.text = ''

        try:
            cal_dic = numpy_inversion(
                self.adjustment, self.results, output_root_dir, write_out_files='n',
            )
        except InversionError:
            show_message(InversionError, "Inversion Error")
        self.match_inversion_results('numpy', cal_dic)

    def match_inversion_results(self, inversion_type, cal_dic=None):
        """
        Populates delta and datum table residuals from inversion results
        :param inversion_type: 'gravnet' or 'numpy'
        cal_dic: has an entry for each meter if a calibration coefficient was
            calculated. Different than
            delta.cal_coeff, which is an a priori specified coefficient.
        """
        datum_residuals, dg_residuals = [], []
        g_dic = self.adjustment.g_dic
        # Reset residuals to all -999's
        for delta in self.deltas:
            delta.residual = -999.0

        # Matchup adjustment residuals with observations
        for delta in self.deltas:
            if delta.checked:
                try:
                    if inversion_type == 'gravnet':
                        station1_name = delta.sta1[:6]
                        station2_name = delta.sta2[:6]
                    elif inversion_type == 'numpy':
                        station1_name = delta.sta1
                        station2_name = delta.sta2
                    if delta.meter in cal_dic:
                        cal_adj_dg = (
                            delta.dg * cal_dic[delta.meter][0] * delta.cal_coeff
                        )
                    else:
                        cal_adj_dg = delta.dg * delta.cal_coeff
                    adj_g1 = g_dic[station1_name]
                    adj_g2 = g_dic[station2_name]
                    adj_dg = adj_g2 - adj_g1
                    delta.residual = adj_dg - cal_adj_dg
                    dg_residuals.append(delta.residual)
                except KeyError:
                    self.msg = show_message("Key error", "Key error")
                    return
            else:
                delta.residual = -999.0

        for datum in self.datums:
            if inversion_type == 'numpy':
                station_name = datum.station
            elif inversion_type == 'gravnet':
                station_name = datum.station[0:6]
            if station_name in g_dic.keys():
                adj_g1 = g_dic[station_name]
                residual = adj_g1 - (datum.g - datum.meas_height * datum.gradient)
                datum_residuals.append(residual)
            else:
                residual = -999
            datum.residual = residual

        self.adjustment.adjustmentresults.min_dg_residual = np.min(np.abs(dg_residuals))
        self.adjustment.adjustmentresults.max_dg_residual = np.max(np.abs(dg_residuals))
        self.adjustment.adjustmentresults.min_datum_residual = np.min(
            np.abs(datum_residuals)
        )
        self.adjustment.adjustmentresults.max_datum_residual = np.max(
            np.abs(datum_residuals)
        )

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

        sta_dic_ls = {station: i for i, station in enumerate(station_list)}
        return sta_dic_ls

    def return_obstreestation(self, station_id):
        """
        returns ObsTreeStation object corresponding to delta_id. Used when recreating deltas from a saved workspace.
        :param station_id: tuple, (station_name, time)
        :return: ObsTreeStation
        """
        for station in self.iter_stations():
            if station.station_name == station_id[0]:
                if abs(station.tmean() - station_id[1]) < 0.0001:
                    return station
        return None

    def populate_delta_model(self, loop=None, clear=True):
        """
        Copy deltas from the delta_model shown on the drift tab to the model
        shown on the adjustment tab.
        :param loop: if ObsTreeLoop, only populate deltas from the selected
            loop; otherwise use all loops
        :param clear: if True, clear delta table first; if not, append deltas
        :return:
        """
        if clear:
            # Clear existing list object: do not reassign [], as this breaks the
            # link to the model object.
            self.deltas.clear()

        # If just a single loop
        if isinstance(loop, ObsTreeLoop):
            try:
                for delta in loop.deltas:
                    # if loop.delta_model.data(loop.delta_model.index(ii, 0), role=Qt.CheckStateRole) == 2:
                    self.deltas.insert(0, delta)
            except Exception as e:
                self.msg = show_message(
                    "Error populating delta table. Please check the drift correction "
                    + "for survey "
                    + self.name
                    + ", loop "
                    + loop.name,
                    "GSadjust error",
                )
        # Populate all loops
        elif loop is None:
            for loop in self.loops():
                if loop.checkState() == Qt.Checked:
                    try:
                        self.deltas += loop.deltas


                            # # Need to create a new delta here instead of just putting the original one, from the
                            # # drift tab, in the net adj. tab. Otherwise checking/unchecking on the net adj tab
                            # # overrides repopulating the delta table.
                            # if type(delta.station2) == list:
                            #     newdelta = DeltaList(None, delta.station2)
                            # else:
                            #     newdelta = DeltaNormal(delta.station1, delta.station2)
                            # newdelta.ls_drift = delta.ls_drift
                            # newdelta.driftcorr = delta.driftcorr
                            # newdelta.type = delta.type
                            # newdelta.loop = delta.loop
                            # self.delta.insert(0, newdelta)
                    except Exception as e:
                        logging.exception(e, exc_info=True)
                        # Sometimes the delta table isn't created when a workspace is loaded

                        self.msg = show_message(
                            "Error populating delta table. Please check the drift correction "
                            + "for survey "
                            + self.name
                            + ", loop "
                            + loop.name,
                            "GSadjust error",
                        )
        return True
