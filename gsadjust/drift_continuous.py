import numpy as np
import logging
from scipy.interpolate import UnivariateSpline
from data_objects import Delta


def drift_continuous(plot_data, drift_x, drift_rate, method_key, tension_slider_value, extrapolation_type, min_time,
                     max_time, loop_name):
    # Interpolate drift model: polynomial, spline, etc. at xp number of points. xp needs to remain
    # relatively low to maintain performance.
    xp = np.linspace(min(drift_x), max(drift_x), 50)  # constant
    # method_key = self.drift_polydegree_combobox.currentIndex()
    if method_key == 0:  # constant drift correction
        mean_drift = sum(drift_rate) / len(drift_rate)
        yp = np.zeros(xp.size) + mean_drift
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
                xs = np.linspace(x0[0], x0[-1], 50)
                yp = s(xs)
                logging.info('Spline drift correction, tension={}'.format(tension_slider_value))
            except Exception as e:
                raise IndexError
        else:
            # Polynomial interpolation. Degree is one less than the method key, e.g.,
            #     method_key == 2 is 1st orderpolynomial, etc.
            try:
                z = np.polyfit(x0, drift_rate, method_key - 1)
                p = np.poly1d(z)
                yp = p(xp0)
                logging.info('Polynomial drift correction degree {}'.format(method_key - 1))
                # obstreeloop.drift_cont_method = 0
            except np.linalg.LinAlgError as e:
                return np.linalg.LinAlgError

    # Method for extrapolating beyond fitted drift curve extene
    if extrapolation_type == 1:  # constant
        new_xp = np.linspace(min_time, min(drift_x), 20)
        new_xp = np.append(new_xp, xp)
        new_xp = np.append(new_xp, np.linspace(max(drift_x), max_time, 20))
        xp = new_xp
        new_yp = np.ones(20) * yp[0]
        new_yp = np.append(new_yp, yp)
        new_yp = np.append(new_yp, np.ones(20) * yp[-1])
        yp = new_yp
    else:  # linear extrapolation from first two (and last two) points
        # get first two points
        x = xp[:2]
        y = yp[:2]
        z = np.polyfit(x, y, 1)
        p = np.poly1d(z)
        new_xp1 = np.linspace(min_time, min(drift_x), 20)
        yp1 = p(new_xp1)
        x = xp[-2:]
        y = yp[-2:]
        z = np.polyfit(x, y, 1)
        p = np.poly1d(z)
        new_xp2 = np.linspace(max(drift_x), max_time, 20)
        yp2 = p(new_xp2)
        xp_temp = np.append(new_xp1, xp)
        xp = np.append(xp_temp, new_xp2)
        yp_temp = np.append(yp1, yp)
        yp = np.append(yp_temp, yp2)
    delta_list = calc_cont_dg(xp, yp, plot_data, loop_name)
    return delta_list, xp, yp

def calc_cont_dg(xp, yp, data, loop_name):
    """
    Calculates delta-g's while removing drift using the input drift model
    :param xp: times of continuous drift model
    :param yp: continuous drift model
    :param data: plot_data list
    :return: PyQt DeltaTableModel
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
    first = True
    for station in data:
        if first:
            first = False
            prev_station = station
            continue
        drift1_idx = min(range(len(xp)), key=lambda i: abs(xp[i] - prev_station.tmean))
        drift1 = yp[drift1_idx]
        drift2_idx = min(range(len(xp)), key=lambda i: abs(xp[i] - station.tmean))
        drift2 = yp[drift2_idx]
        delta = Delta(prev_station,
                      station,
                      driftcorr=drift2 - drift1,
                      loop=loop_name)
        delta_list.append(delta)
        # delta_model.insertRows(delta, 0)
        prev_station = station
    return delta_list