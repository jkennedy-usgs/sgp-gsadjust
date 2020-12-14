"""
drift/continuous.py
===================

GSadjust code for calculating continuous-model drift correction.
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

import numpy as np
from scipy.interpolate import UnivariateSpline

from ..data import DeltaNormal

N_PTS_IN_INTERPOLATION = 300
N_PTS_IN_EXTRAPOLATION = 200


def drift_continuous(
    data,
    plot_data,
    drift_x,
    drift_rate,
    method_key,
    tension_slider_value,
    extrapolation_type,
    weight_obs,
    min_time,
    max_time,
    loop_name,
):
    """Interpolate drift model: polynomial, spline, etc. at N_PTS_IN_INTERPOLATION
    points, plus N_PTS_IN_EXTRAPOLATION on either side.

    These values need to be relatively small for decent performance.

    Parameters
    ----------
    data : list
        list of ObsTreeStations
    plot_data : list
        One entry per station. Each entry is a list, in the order:
        [[plot time values],
         [plot g values],
         station name,
         [g standard deviation],
         [time standard deviation]]
    drift_x : list
        Drift time observations (x location of points on bottom plot)
    drift_rate : list
        Drift rate observations (y location of points on bottom plot)
    method_key : {0, 1, 2, 3, 4}
        Indicates type of drift correction.
        0: Constant
        1: Spline
        2-4: Polynomial (degree is one less than the value)
    tension_slider_value : int
        Controls tension on the interpolated spline
    extrapolation_type : {1, anything else}
        Controls how interpolation is extended from the outermost data.
        1: Constant
        not 1: linearly extend the fitted curve at the same slope as the first/last
               2 data points
    weight_obs : int
        Controls if observations are weighted when fitting a constant drift rate.
        Only used if drift is set to constant, not for other methods.
        0: no weighting
        not 0: weighted
    min_time : float
        Time to extrapolate at the beginning. Should be the time of the first station
        occupation of the loop.
    max_time : float
        Time to extrapolate at the end. Should be the time of the last station
        occupation of the loop.
    loop_name : str
        loop name, for creating deltas

    Returns
    -------
    delta_list : list
        List of deltas
    xp : ndarray
        For plotting the bottom plot
    yp : ndarray
        For plotting the bottom plot
    z_main : (mean_drift, sigma)
        These are displayed on the plot.

    """

    xp = np.linspace(min(drift_x), max(drift_x), N_PTS_IN_INTERPOLATION)  # constant
    drift_stats = None
    z_main = []
    if method_key == 0:  # constant drift correction
        if weight_obs == 0:
            mean_drift = sum(drift_rate) / len(drift_rate)
            sigma = np.std(drift_rate) / np.sqrt(len(drift_rate))
            yp = np.zeros(xp.size) + mean_drift
            z_main = [(mean_drift, sigma)]
        # Weight observations according to NGA method
        else:
            drifts, drift_w = [], []
            for station_data in plot_data:
                t, R, Rsd, tsd = (
                    station_data[0],
                    station_data[1],
                    station_data[3],
                    station_data[4],
                )
                if len(t) > 1:
                    for i in range(1, len(t)):
                        dr = R[i] - R[0]
                        dt = (t[i] - t[0]) * 24
                        sdr = np.sqrt(Rsd[i] ** 2 + Rsd[0] ** 2)
                        sdt = np.sqrt(tsd[i] ** 2 + tsd[0] ** 2)
                        drifts.append(dr / dt)
                        drift_sd = np.sqrt(
                            sdr ** 2 / dt ** 2 + dr ** 2 * sdt ** 2 / dt ** 4
                        )
                        drift_w.append(1 / drift_sd ** 2)
            num = []
            for idx, w in enumerate(drift_w):
                num.append(w * drifts[idx])
            mean_drift = np.sum(num) / np.sum(drift_w)
            num = []
            for idx, w in enumerate(drift_w):
                num.append(w * (drifts[idx] - mean_drift) ** 2)
            sigma_d = np.sqrt(np.sum(num) / ((len(drift_w) - 1) * np.sum(drift_w)))
            drift_stats = dict()
            drift_stats["t0"] = plot_data[0][0][0]
            drift_stats["sigma_d"] = sigma_d
            drift_stats["mean_drift"] = mean_drift
            yp = np.zeros(xp.size) + mean_drift
            z_main = [(mean_drift, sigma_d)]
    else:
        x0 = [f - np.min(drift_x) for f in drift_x]
        xp0 = [f - np.min(xp) for f in xp]
        idx = sorted(range(len(x0)), key=lambda xpt: x0[xpt])
        x_sorted, drift_rate_sorted = [], []
        for i in idx:
            x_sorted.append(x0[i])
            drift_rate_sorted.append(drift_rate[i])
        x0 = x_sorted
        drift_rate = drift_rate_sorted
        if method_key == 9:
            pass
        if method_key == 1:  # spline
            try:
                s = UnivariateSpline(x0, drift_rate, k=3, s=tension_slider_value)
                xs = np.linspace(x0[0], x0[-1], N_PTS_IN_INTERPOLATION)
                yp = s(xs)
                logging.info(
                    "Spline drift correction, tension={}".format(tension_slider_value)
                )
            except Exception:
                raise IndexError
        else:
            # Polynomial interpolation. Degree is one less than the method key, e.g.,
            #     method_key == 2 is 1st order polynomial, etc.
            try:
                z_main = np.polyfit(x0, drift_rate, method_key - 1)
                p = np.poly1d(z_main)
                yp = p(xp0)
                logging.info(
                    "Polynomial drift correction degree {}".format(method_key - 1)
                )
            except np.linalg.LinAlgError as e:
                return np.linalg.LinAlgError

    # Method for extrapolating beyond fitted drift curve extent
    if extrapolation_type == 1:  # constant
        new_xp = np.linspace(min_time, min(drift_x), N_PTS_IN_EXTRAPOLATION)
        new_xp = np.append(new_xp, xp)
        new_xp = np.append(new_xp, np.linspace(max(drift_x), max_time, 200))
        xp = new_xp
        new_yp = np.ones(200) * yp[0]
        new_yp = np.append(new_yp, yp)
        new_yp = np.append(new_yp, np.ones(200) * yp[-1])
        yp = new_yp
    else:  # linear extrapolation from first two (and last two) points
        # get first two points
        x = xp[:2]
        y = yp[:2]
        z = np.polyfit(x, y, 1)
        p = np.poly1d(z)
        new_xp1 = np.linspace(min_time, min(drift_x), 200)
        yp1 = p(new_xp1)
        x = xp[-2:]
        y = yp[-2:]
        z = np.polyfit(x, y, 1)
        p = np.poly1d(z)
        new_xp2 = np.linspace(max(drift_x), max_time, 200)
        yp2 = p(new_xp2)
        xp_temp = np.append(new_xp1, xp)
        xp = np.append(xp_temp, new_xp2)
        yp_temp = np.append(yp1, yp)
        yp = np.append(yp_temp, yp2)
    delta_list = calc_cont_dg(xp, yp, data, loop_name, drift_stats)
    return delta_list, xp, yp, z_main


def calc_cont_dg(xp, yp, data, loop_name, drift_stats):
    """
    Calculates delta-g's while removing drift using the input drift model

    Parameters
    ----------
    xp : ndarray
        times of continuous drift model
    yp : ndarray
        continuous drift model
    data : list
        list of ObsTreeStations
    loop_name : str
    drift_stats : dict or None
        If dict, observations are weighted; None if not

     Returns
     -------
    list
        list of deltas
    """
    first = True
    ypsum = [0]
    delta_list = []
    for x, drift_rate in zip(xp, yp):
        if first:
            first = False
            prev_x = x
        else:
            prev_sum = ypsum[-1]
            interval = (x - prev_x) * 24
            prev_x = x
            ypsum.append(prev_sum + drift_rate * interval)

    xp = xp.tolist()
    yp = ypsum  # yp = yp.tolist()
    prev_station = data.pop(0)
    for station in data:
        # If using weighted dg
        if drift_stats:
            station.assigned_sd = np.sqrt(
                station.original_sd ** 2
                + ((station.tmean() - drift_stats["t0"]) * 24) ** 2
                * drift_stats["sigma_d"] ** 2
                + np.sqrt(station.t_stdev ** 2 + data[0].t_stdev ** 2)
                * drift_stats["mean_drift"] ** 2
            )
        else:
            station.assigned_sd = None
        drift1_idx = min(
            range(len(xp)), key=lambda i: abs(xp[i] - prev_station.tmean())
        )
        drift1 = yp[drift1_idx]
        drift2_idx = min(range(len(xp)), key=lambda i: abs(xp[i] - station.tmean()))
        drift2 = yp[drift2_idx]
        delta = DeltaNormal(
            prev_station, station, driftcorr=drift2 - drift1, loop=loop_name
        )
        delta_list.append(delta)
        prev_station = station
    return delta_list
