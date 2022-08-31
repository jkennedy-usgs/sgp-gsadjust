"""
get_nwis_data Retrieve groundwater-level data from the USGS National Water Information System.

"""
from dateutil import parser
import requests, json
import numpy as np
import matplotlib.pyplot as plt
import datetime
from collections import OrderedDict, defaultdict
from math import radians, sin, cos, atan2, sqrt
from scipy.interpolate import interp1d


def plot_hydrograph_with_gravity(site_id, dates):
    fig, ax = plt.subplots()
    well_data = nwis_get_data(site_id)

    if well_data['continuous_x']:
        ydata_meters = [x * 0.3048 for x in well_data['continuous_y']]
        ax.plot(well_data['continuous_x'], ydata_meters, label=site_id, linewidth=0.5)
    if well_data['discrete_x']:
        ydata_meters = [x * 0.3048 for x in well_data['discrete_y']]
        ax.plot(well_data['discrete_x'], ydata_meters, 'k.', label=site_id)
    ax.invert_yaxis()
    ax.set_ylabel("Depth to groundwater (meters)")
    ylim = ax.get_ylim()
    d1 = datetime.datetime.strptime(dates[0], '%Y-%m-%d')
    d2 = datetime.datetime.strptime(dates[1], '%Y-%m-%d')
    ax.plot([d1, d1], ylim, 'k')
    ax.plot([d2, d2], ylim, 'k')
    ax.set_ylim(ylim)
    plt.title(site_id)
    # for tick in ax.get_xticklabels():
    #     tick.set_fontname("Times New Roman")
    # for tick in ax.get_yticklabels():
    #     tick.set_fontname("Times New Roman")
    # L = plt.legend()
    # plt.setp(L.texts, family='Times New Roman')
    # plt.xlim([datetime.datetime(2006,1,1,0,0,0), datetime.datetime(2019,1,1,0,0,0)])
    plt.show()


def plot_wells(cross_ref_file, site_IDs):
    fig, ax = plt.subplots()
    for well in site_IDs:
        well_data = nwis_get_data(cross_ref_file, well)

        if well_data['continuous_x']:
            ydata_meters = [x * 0.3048 for x in well_data['continuous_y']]
            ax.plot(well_data['continuous_x'], ydata_meters, label=well, linewidth=0.5)
        else:
            ydata_meters = [x * 0.3048 for x in well_data['discrete_y']]
            ax.plot(well_data['discrete_x'], ydata_meters, label=well, linewidth=0.5)
    ax.invert_yaxis()
    ax.set_ylabel("Depth to groundwater (meters)", fontname="Times New Roman", fontsize=12)
    for tick in ax.get_xticklabels():
        tick.set_fontname("Times New Roman")
    for tick in ax.get_yticklabels():
        tick.set_fontname("Times New Roman")
    L = plt.legend()
    plt.setp(L.texts, family='Times New Roman')
    plt.xlim([datetime.datetime(2006,1,1,0,0,0), datetime.datetime(2019,1,1,0,0,0)])
    plt.show()
    input()


def get_NWIS_ID(cross_ref_file, gravity_station_ID):
    if gravity_station_ID.isnumeric() and len(gravity_station_ID) == 15:
        nwis_ID = gravity_station_ID
        grav_ID = gravity_station_ID
    else:
        with open(cross_ref_file, 'r') as fid:
            for line in fid:

                grav_ID, nwis_ID = line.strip().split(',')
                if grav_ID.upper() == gravity_station_ID.upper():
                    if len(nwis_ID) != 15:
                        print('Invalid 15-digit ID for site {}'.format(gravity_station_ID))
                        return 0
                    else:
                        break
            else:
                print('Gravity Site-ID {} not found in cross-ref file'.format(gravity_station_ID))
                return 0


def nwis_get_data(nwis_ID):
    """Gets NWIS groundwater-level data via REST API
    
    :param nwis_ID: USGS 15-digit ID
    """
    out_dic = dict()
    name_dic = dict()
    # rdb_meas retrieval is preferred, it returns both discrete and continuous measurements.
    today = datetime.datetime.today().strftime('%Y-%m-%d')
    nwis_URL = 'http://nwis.waterdata.usgs.gov/nwis/dv?cb_72019=on&format=rdb_meas' + \
               f'&site_no={nwis_ID}' + \
               '&referred_module=gw&period=&begin_date=2010-10-01&end_date=%s' % today
    # print('Retrieving rdb data for {} from {}'.format(grav_ID, nwis_URL))
    r = requests.get(nwis_URL)
    # If there is continuous data, it will start with '# ----... WARNING ---...'
    nwis_data = r.text.split('\n')
    if nwis_data[-2].split()[0] != "agency_cd":
        c_x, c_y, d_x, d_y = parse_continuous_data(nwis_data)
    if len(d_x) == 0:
        # if no continuous data, retrieve discrete data
        nwis_URL: str = f'https://nwis.waterdata.usgs.gov/nwis/gwlevels/?site_no={nwis_ID}' + \
                   '&format=rdb_meas'
        # print('No continuous data for {}. Retrieving discrete data from {}'.format(grav_ID, nwis_URL))
        r = requests.get(nwis_URL)
        nwis_data = r.text.split('\n')
        c_x, c_y, d_x, d_y = parse_discrete_data(nwis_data)


    out_dic['continuous_x'] = c_x
    out_dic['continuous_y'] = c_y
    out_dic['discrete_x'] = d_x
    out_dic['discrete_y'] = d_y
    # if not out_dic:
    #     print('No NWIS data found for site {}'.format(gravity_station_ID))
    return out_dic


def parse_discrete_data(nwis_data):
    discrete_x, discrete_y = [], []
    continuous_x, continuous_y = [], []
    name_dic = dict()
    for nwis_line in nwis_data:
        line_elems = nwis_line.split('\t')
        # Need to test for null strings because it's possible for there to be a date without a measurement.
        try:  # the rdb format has changed; parsing by '\t' barely works with the fixed-width fields
            # The parameter we want, 72019_00003 (mean), isn't always in the same column, it
            # depends on what other parameters are output. Search for it here
            if line_elems[0] == u'agency_cd':
                gwl_disc_index = \
                [idx for idx, s in enumerate(line_elems) if s == 'lev_va'][0]
                continue
            elif line_elems[0] == "#" and line_elems[1] == "USGS":
                if len(line_elems) > 3:
                    site_name = " ".join(line_elems[3:])
                else:
                    site_name = ""
                name_dic[line_elems[2]] = site_name
                continue
            elif line_elems[0] == u'USGS':  # Doesn't work for other agency codes
                if line_elems[gwl_disc_index] != u'':
                    discrete_x.append(parser.parse(line_elems[3]))
                    discrete_y.append(float(line_elems[gwl_disc_index]))
                # elif line_elems[gwl_disc_index] != u'':  # Get depth to water
                #     discrete_x.append(parser.parse(line_elems[2]))
                #     discrete_y.append(float(line_elems[gwl_disc_index]))
        except Exception as e:
            continue
    return continuous_x, continuous_y, discrete_x, discrete_y


def parse_continuous_data(nwis_data):
    discrete_x, discrete_y = [], []
    continuous_x, continuous_y = [], []
    name_dic = dict()
    for nwis_line in nwis_data:
        line_elems = nwis_line.split('\t')
        # Need to test for null strings because it's possible for there to be a date without a measurement.
        try:  # the rdb format has changed; parsing by '\t' barely works with the fixed-width fields
            # The parameter we want, 72019_00003 (mean), isn't always in the same column, it
            # depends on what other parameters are output. Search for it here
            if line_elems[0] == u'agency_cd':
                gwl_cont_index = \
                [idx for idx, s in enumerate(line_elems) if '72019_00002' in s][0]
                gwl_disc_index = \
                [idx for idx, s in enumerate(line_elems) if '72019_00011' in s][0]
            elif line_elems[0] == "#" and line_elems[1] == "USGS":
                if len(line_elems) > 3:
                    site_name = str.join(" ", line_elems[3:])
                else:
                    site_name = ""
                name_dic[line_elems[2]] = site_name
            elif line_elems[0] == u'USGS':  # Doesn't work for other agency codes
                if line_elems[2] == u'GW':
                    if line_elems[10] != u'':
                        discrete_x.append(parser.parse(line_elems[3]))
                        discrete_y.append(float(line_elems[10]))
                elif line_elems[gwl_disc_index] != u'':  # Get depth to water
                    discrete_x.append(parser.parse(line_elems[2]))
                    discrete_y.append(float(line_elems[gwl_disc_index]))
                elif line_elems[gwl_cont_index] != u'':
                    continuous_x.append(parser.parse(line_elems[2]))
                    continuous_y.append(float(line_elems[gwl_cont_index]))
        except Exception as e:
            continue
    return continuous_x, continuous_y, discrete_x, discrete_y


def search_nwis(coords):
    degree_buffer = 0.01
    ID_dict = defaultdict(list)
    try:
        lon_min = coords[0] - degree_buffer
        lon_max = coords[0] + degree_buffer
        lat_min = coords[1] - degree_buffer
        lat_max = coords[1] + degree_buffer
        nwis_URL = "https://waterservices.usgs.gov/nwis/site/?format=rdb&bBox=" \
              f"{lon_min:.5f},{lat_min:.5f},{lon_max:.5f},{lat_max:.5f}" \
              "&outputDataTypeCd=gw&siteType" \
              "=GW&siteStatus=all&hasDataTypeCd=gw"
        r = requests.get(nwis_URL)

        nwis_data = r.text.split('\n')
        for nwis_line in nwis_data:
            line_elems = nwis_line.split('\t')

            # Need to test for null strings because it's possible for there to be a date without a measurement.
            try:  # the rdb format has changed; parsing by '\t' barely works with the fixed-width fields
                if line_elems[0] == u'USGS' and line_elems[13] == '72019':  # Doesn't work for other agency codes
                    lat = float(line_elems[4])
                    lon = float(line_elems[5])
                    distance = calc_distance(lat,
                                                 lon,
                                                 coords[1],
                                                 coords[0])
                    depth = get_site_property(int(line_elems[1]), prop='wellDepth')
                    ID_dict[line_elems[1]] = [line_elems[2],
                                              f"{distance:.0f}",
                                              line_elems[21],
                                              line_elems[22],
                                              depth]
            except:
                continue
    except:
        return {}
    ID_dict = OrderedDict(sorted(ID_dict.items(), key=lambda x: int(x[1][1])))
    return ID_dict


def calc_distance(lat1, lon1, lat2, lon2):
    """

    Returns
    -------
    object
    """
    R = 6373000
    dlon = radians(lon2) - radians(lon1)
    dlat = radians(lat2) - radians(lat1)

    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))

    distance_haversine_formula = R * c
    return distance_haversine_formula


def get_gwls_for_date_and_coords(date, coords, buffer=365):
    dt = datetime.datetime.strptime(date, '%Y-%m-%d')
    window_end = (dt + datetime.timedelta(days=buffer)).strftime('%Y-%m-%d')
    window_start = (dt - datetime.timedelta(days=buffer)).strftime('%Y-%m-%d')

    # Get DVs
    url = 'https://waterservices.usgs.gov/nwis/dv/?format=rdb&bBox=' + \
          f'{coords[0]:0.3f},{coords[2]:.3f},{coords[1]:.3f},{coords[3]:.3f}' + \
          f'&startDT={window_start}&endDT={window_end}' + \
          '&parameterCd=72019&statisticCd=00002&siteType=GW&siteStatus=all'
          # Anza project only publishes statisticCd=00002 (minimum)
          # Albuquerque project publishes 00001,00002,00003 (max, min, mean)

    s = requests.get(url, verify=False).text.split('\n')
    data = []
    for line in s:
        if line[:4] == 'USGS':
            data.append(line.split('\t'))

    # get gwlevels
    url = r'https://waterservices.usgs.gov/nwis/gwlevels/?format=rdb&bBox=' + \
          f'{coords[0]:0.3f},{coords[2]:.3f},{coords[1]:.3f},{coords[3]:.3f}' + \
          f'&startDT={window_start}&endDT={window_end}&parameterCd=72019'
    s = requests.get(url, verify=False).text.split('\n')
    for line in s:
        if line[:4] == 'USGS':
            line_elems = line.split('\t')
            data.append([line_elems[0],
                         line_elems[1],
                         line_elems[3],
                         line_elems[6],
                         'A'])

    df1 = np.array(data)

    return df1


def get_gwl_change(dates, coords, convert=False, **kwargs):
    df1 = get_gwls_for_date_and_coords(dates[0], coords, **kwargs)
    df2 = get_gwls_for_date_and_coords(dates[1], coords, **kwargs)

    out, lat, lon = dict(), dict(), dict(),
    date1, date2 = dict(), dict()
    stations = np.unique(df1[:,1])
    stations2 = np.unique(df2[:,1])
    for station in stations:
        if station not in stations2:
            continue
        try:
            wl1, d1 = get_closest_wl(df1[df1[:, 1] == station, :], dates[0])
            wl2, d2 = get_closest_wl(df2[df2[:, 1] == station, :], dates[1])
            if d1 == d2:
                continue
            # Subtracting depths to groundwater, hence wl1 - wl2 (not wl2 - wl1)
            if convert:
                out[station] = (wl1 - wl2) * .3048
            else:
                out[station] = wl1 - wl2
            lat[station] = float(station[:2]) + float(station[2:4]) / 60 + float(
                station[4:6]) / 3600
            lon[station] = -1 * (float(station[6:9]) + float(station[9:11]) / 60 + float(
                station[11:13]) / 3600)
            date1[station] = d1
            date2[station] = d2
        except:
            continue

    return out.keys(), out, lat, lon, date1, date2


def get_closest_wl(df, date):
    a = np.array([datetime.datetime.strptime(x, '%Y-%m-%d') for x in df[:, 2]])
    d = datetime.datetime.strptime(date, '%Y-%m-%d')
    b = np.absolute(a - d)
    index_of_min = np.argmin(b)
    closest_wl = df[index_of_min, 3]
    closest_date = df[index_of_min, 2]
    return float(closest_wl), closest_date


def filter_ts(nd, window):
    try:
        ts = smooth(np.array(nd['continuous_y']), window_len=window)
        nd['continuous_y'] = ts
        # plt.plot(v['continuous_x'], v['continuous_y'])
        # plt.plot(v['continuous_x'], ts)
        # plt.title(k)
        # plt.show()
    except:
        pass
    try:
        x = np.array([toTimestamp(x) for x in (nd['discrete_x'])])
        F = interp1d(x, nd['discrete_y'])
        interpx = np.arange(x[0], x[-1],
                            1)
        interpy = F(interpx)
        # plt.plot(interpx, interpy)
        ts = smooth(interpy, window_len=window, window='flat')
        nd['discrete_x'] = [datetime.datetime.fromtimestamp(x*86400) for x in interpx]
        nd['discrete_y'] = ts
        # plt.plot(v['discrete_x'], v['discrete_y'])
        # plt.plot(interpx, ts)
        # plt.title(k)
        # plt.show()
    except:
        pass
    return nd


def toTimestamp(d):
    return ((d - datetime.datetime(1970, 1, 1)) / datetime.timedelta(seconds=1)) / 86400


def smooth(x, window_len=11, window='hanning'):
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len < 3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError(
            "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
    s = np.r_[2 * x[0] - x[window_len - 1::-1], x, 2 * x[-1] - x[-1:-window_len:-1]]
    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')
    y = np.convolve(w / w.sum(), s, mode='same')
    return y[window_len:-window_len + 1]


def get_site_property(ID, prop='monitoringLocationAltitudeLandSurface'):
    usgs_endpoint = r'https://labs.waterdata.usgs.gov/api/observations/collections'
    try:
        geo_json = requests.get('/'.join([usgs_endpoint,
                                          'monitoring-locations',
                                          'items',
                                          f'USGS-{ID}']))
        station = json.loads(geo_json.text)
        station_prop = station['properties'][prop]
        if station_prop:
            return station_prop
        else:
            return ""
    except json.decoder.JSONDecodeError as e:
        return ""


# Wrote this and didn't need it, to get station name from header.
#
# if line_elems[0][0] == "#":
#     line_elems = nwis_line.split()
#     if len(line_elems) > 3:
#         if line_elems[0] == "#" and line_elems[1] == "USGS":
#             site_name = str.join(" ", line_elems[3:])
#         else:
#             site_name = ""
#         name_dic[line_elems[2]] = site_name