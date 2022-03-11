"""
data/analysis.py
===============

Analysis methods for gravity surveys.

For use by GSadjust, software for the interactive network adjustment of
relative-gravity networks
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
import numpy as np
from .adjustment import AdjustedStation

class InversionError(Exception):
    pass


def _equal_or_truncated(a, b):
    """For matching Gravnet output (Gravnet truncates to 6 characters)"""
    return a == b or a[:6] == b


def numpy_inversion(adjustment):
    """
    Pre- and post-processing data for network adjustment.

    The actual adjustment is carried out in adjustment.python_lsq_inversion()

    Parameters
    ----------
    adjustment : Adjustment object
        Stores information needed for the

    """

    adjustment.adjustmentresults.text = []

    # sta_dic_LS is a dictionary, key: station name, value: column for A matrix
    sta_dic_ls = adjustment.sta_dic_ls
    loop_ls_dict = adjustment.loop_ls_dict

    # get number of observations:
    n_rel_obs = len(adjustment.deltas)
    n_abs_obs = len(adjustment.datums)
    n_meters = adjustment.n_meters
    ndrift = adjustment.ndrift

    # dict of tuples, used to identify column of drift observation in A matrix:
    # (loop.name, (column relative to end of A matrix, drift degree)
    netadj_loop_keys = adjustment.netadj_loop_keys

    # Initialize least squares matrices
    # number of unknowns
    nb_x = len(sta_dic_ls) + ndrift + n_meters
    adjustment.adjustmentresults.n_unknowns = nb_x
    # model matrix:
    A = np.zeros((n_rel_obs + n_abs_obs, nb_x))
    # weight matrix:
    P = np.zeros((n_rel_obs + n_abs_obs, n_rel_obs + n_abs_obs))  # pas sur
    # observation matrix:
    Obs = np.zeros((n_rel_obs + n_abs_obs, 1))
    # datum-free constraint vector:
    S = np.zeros((nb_x, 1))

    row = 0
    delta_keys = []

    # Populate least squares matrices
    for delta in adjustment.deltas:
        delta_keys.append(delta.__hash__())
        dg = delta.dg if delta.type != "assigned" else delta.assigned_dg
        Obs[row] = dg * delta.cal_coeff
        P[row, row] = 1.0 / (delta.adj_sd ** 2)
        A[row, sta_dic_ls[delta.sta1]] = -1
        A[row, sta_dic_ls[delta.sta2]] = 1

        # Populate 1 column per gravimeter for calibration coefficient
        if adjustment.adjustmentoptions.cal_coeff:
            meter = delta.meter
            A[row, adjustment.meter_dic[meter] + len(sta_dic_ls)] = delta.dg

        # Populate column(s) for drift, if included in network adjustment
        if delta.ls_drift is not None:
            loop_name = delta.ls_drift[0]

            # It's possible for ls_drift to have been set, but the loop method to be
            # something other than netadj
            if loop_name:
                if loop_ls_dict[loop_name] == "netadj":
                    for i in range(delta.ls_drift[1]):  # degree of polynomial
                        A[
                            row,
                            len(sta_dic_ls)
                            + n_meters
                            + netadj_loop_keys[loop_name][0]
                            + i,
                        ] = (delta.sta2_t - delta.sta1_t) ** (i + 1)

        S[sta_dic_ls[delta.sta1]] = 1
        S[sta_dic_ls[delta.sta2]] = 1
        row += 1

    # add datum observation (absolute station(s) or station(s) with fixed values)
    # Key errors handled by calling routine
    i = 0
    for datum in adjustment.datums:
        A[n_rel_obs + i, sta_dic_ls[datum.station]] = 1
        P[n_rel_obs + i, n_rel_obs + i] = 1.0 / datum.sd ** 2
        Obs[n_rel_obs + i] = datum.g
        i += 1

    # Do the inversion
    adjustment.A = A
    adjustment.P = P
    adjustment.Obs = Obs
    adjustment.S = S
    adjustment.dof = n_rel_obs + n_abs_obs - nb_x
    adjustment.g_dic = dict()
    # zero-division errors are caught by caller
    adjustment.python_lsq_inversion()

    # Populate results
    results, sd_all = [], []
    for i in range(len(sta_dic_ls)):
        for key, val in sta_dic_ls.items():
            if val == i:
                try:
                    g = float(adjustment.X[i])
                    sd = float(np.sqrt(adjustment.var[i]))
                    t = AdjustedStation(key, g, sd)
                    results.append(t)
                    adjustment.g_dic[key] = g
                    adjustment.sd_dic[key] = sd
                    sd_all.append(sd)
                except Exception:
                    raise InversionError("Bad variance in results.")
                    return
    adjustment.adjustmentresults.avg_stdev = np.mean(sd_all)

    # Retrieve calibration coefficient(s)
    cal_dic = dict()
    if adjustment.adjustmentoptions.cal_coeff:
        for k, v in adjustment.meter_dic.items():
            cal_dic[k] = (
                float(1 - adjustment.X[len(sta_dic_ls) + v]),
                float(np.sqrt(adjustment.var[len(sta_dic_ls) + v])),
            )
    else:
        for k, v in adjustment.meter_dic.items():
            cal_dic[k] = (1.0, 0.0)
    adjustment.adjustmentresults.cal_dic = cal_dic

    # calculate and display statistics:
    adjustment.lsq_statistics()

    return results


def compute_gravity_change(obstreemodel, table_type="simple"):
    """Calculates gravity change between surveys.

    Parameters
    ----------
    obstreemodel : ObsTreeModel
        Tree structure with the data.
    table_type : {"simple", "full", "list"}
        Controls what is shown in the table.
        Simple: One column per time interval. Dg calculated between each interval, and
        between first and each subsequent interval.
        Full: Shows g and sd at each station, and gravity changes converted to meters
        of water.
        List: Shows station, date, g, sd for each adjusted station.

    Returns
    -------
    tuple
        header, table, dates

    """

    # Check that all values are positive (it should work either way, but it avoids confusion)
    # for i in range(obstreemodel.invisibleRootItem().rowCount()):
    #     survey = obstreemodel.invisibleRootItem().child(i)
    #     for ii in range(survey.results_model.rowCount()):
    #         adj_station = survey.results_model.data(survey.results_model.index(ii, 0),
    #                                                 role=256)  # 256=Qt.UserRole
    #         if adj_station.g < 0:
    #             return False

    compare_station, initial_station, iteration_station, iteration_name = (
        None,
        None,
        None,
        None,
    )
    logging.info("Calculating gravity change")
    first = True
    # unique_station_names = set()
    # for survey in obstreemodel.checked_surveys():
    #     for ii in range(survey.rowCount()):
    #         loop = survey.child(ii)
    #         for iii in range(loop.rowCount()):
    #             station = loop.child(iii)
    #             unique_station_names.add(station.station_name)
    unique_station_names = obstreemodel.results_stations()
    out_table_iteration, out_table_cumulative = [], []
    header1, header2 = [], []
    lat, lon, elev, all_g = [], [], [], []
    dates, header = [], []

    # Simple list of Station, Date, g, sd
    # Copied from results table on NA tab
    if table_type == "list":
        date_col, station_col, sd_col = [], [], []
        for survey in obstreemodel.checked_surveys():
            dates.append(dt.datetime.strptime(survey.name, "%Y-%m-%d"))
            header = ["Station", "Date", "g", "Std. dev."]
            for adj_station in survey.results:
                station_col.append(adj_station.station)
                date_col.append(survey.name)
                all_g.append(adj_station.g)
                sd_col.append(adj_station.sd)
        table = [station_col, date_col, all_g, sd_col]
        return header, table, dates

    # get station coordinate and g/sd data to later append to dg table
    if table_type == "full":
        for station in unique_station_names:
            station_g = []
            g_header = []
            # station_coords can not exist during testing, maybe other times also?
            try:
                lonc, latc, elevc = obstreemodel.station_coords[station]
                lon.append(f"{lonc:.5f}")
                lat.append(f"{latc:.5f}")
                elev.append(f"{elevc:.5f}")
            except TypeError:
                lat.append(-999)
                lon.append(-999)
                elev.append(-999)
            except KeyError:
                lat.append(-999)
                lon.append(-999)
                elev.append(-999)
            for survey in obstreemodel.checked_surveys():
                g_header.append(survey.name + "_g")
                g_header.append(survey.name + "_sd")
                for adj_station in survey.results:
                    # This might break on gravnet output where station names are
                    # truncated to 6 characters
                    if adj_station.station == station:
                        station_g.append("{:0.1f}".format(adj_station.g))
                        station_g.append("{:0.1f}".format(adj_station.sd))
                        break
                else:
                    station_g.append("")
                    station_g.append("")
            all_g.append(station_g)

    for survey in obstreemodel.checked_surveys():
        dates.append(dt.datetime.strptime(survey.name, "%Y-%m-%d"))
        diff_cumulative = []
        diff_iteration = []
        if table_type == "full":
            diff_cumulative_sd, diff_iteration_sd = [], []
        if first:
            # Calculate both the between-survey change and the change from the initial survey
            initial_survey = survey.results
            iteration_reference = initial_survey
            reference_name = survey.name
            iteration_name = reference_name
            first = False
        else:
            if table_type == "simple":
                header1.append(f"{iteration_name}_to_{survey.name}")
                header2.append(f"{reference_name}_to_{survey.name}")

            elif table_type == "full":
                header1.append(f"dH2O_{iteration_name}_to_{survey.name}")
                header1.append(f"dH2O_sd_{iteration_name}_to_{survey.name}")
                header2.append(f"dH2O_{reference_name}_to_{survey.name}")
                header2.append(f"dH2O_sd_{reference_name}_to_{survey.name}")

            compare_survey = survey.results
            for station_name in unique_station_names:

                for initial_station in initial_survey:
                    # Iterate through, look for matching station. 'if' statements deal
                    # with Gravnet, which truncates station names to 6 characters
                    if _equal_or_truncated(initial_station.station, station_name):
                        break
                else:
                    # If we get to the end without breaking, set it to None.
                    initial_station = None

                for iteration_station in iteration_reference:
                    if _equal_or_truncated(iteration_station.station, station_name):
                        break

                else:
                    # If we get to the end without breaking, set it to None.
                    iteration_station = None

                for compare_station in compare_survey:
                    if _equal_or_truncated(compare_station.station, station_name):
                        break
                else:
                    # If we get to the end without breaking, set it to None.
                    compare_station = None

                if initial_station is not None and compare_station is not None:
                    if table_type == "simple":
                        diff_cumulative.append(
                            "{:0.1f}".format(compare_station.g - initial_station.g)
                        )
                    elif table_type == "full":
                        diff_cumulative.append(
                            "{:0.2f}".format(
                                (compare_station.g - initial_station.g) / 41.9
                            )
                        )
                        var = (
                            np.sqrt(compare_station.sd ** 2 + initial_station.sd ** 2)
                            / 41.9
                        )
                        if np.isnan(var):
                            diff_cumulative_sd.append("")
                        else:
                            diff_cumulative_sd.append("{:0.2f}".format(var))
                else:
                    diff_cumulative.append("")
                    if table_type == "full":
                        diff_cumulative_sd.append("")  # for sd column
                if iteration_station is not None and compare_station is not None:
                    if table_type == "simple":
                        diff_iteration.append(
                            "{:0.1f}".format(compare_station.g - iteration_station.g)
                        )
                    elif table_type == "full":
                        diff_iteration.append(
                            "{:0.2f}".format(
                                (compare_station.g - iteration_station.g) / 41.9
                            )
                        )
                        var = (
                            np.sqrt(compare_station.sd ** 2 + iteration_station.sd ** 2)
                            / 41.9
                        )
                        if np.isnan(var):
                            diff_iteration_sd.append("")
                        else:
                            diff_iteration_sd.append("{:0.2f}".format(var))
                else:
                    diff_iteration.append("")
                    if table_type == "full":
                        diff_iteration_sd.append("")  # for sd column
            out_table_iteration.append(diff_iteration)
            out_table_cumulative.append(diff_cumulative)
            if table_type == "full":
                out_table_iteration.append(diff_iteration_sd)
                out_table_cumulative.append(diff_cumulative_sd)
            iteration_reference = compare_survey
            iteration_name = survey.name
    out_table = (
        [list(unique_station_names)] + out_table_iteration + out_table_cumulative
    )

    if table_type == "simple":
        # deal with 2-survey case
        if header1 == header2:
            header = ["station"] + header1
            table = out_table[:-1]
        else:
            header = ["station"] + header1 + header2
            table = out_table
        return header, table, dates
    elif table_type == "full":
        # transpose table
        g = [list(i) for i in zip(*all_g)]
        table = [unique_station_names, lat, lon, elev]
        table += g
        # deal with 2-survey case
        if header1 == header2:
            header = (
                ["Station", "Latitude", "Longitude", "Elevation"] + g_header + header1
            )
            table += out_table_iteration
        else:
            header = (
                ["Station", "Latitude", "Longitude", "Elevation"]
                + g_header
                + header1
                + header2
            )
            table += out_table_iteration
            table += out_table_cumulative
        # transpose back
        table = [list(i) for i in zip(*table)]
        return header, table, dates
