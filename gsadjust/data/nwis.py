"""
get_nwis_data Retrieve groundwater-level data from the USGS National Water Information System.

"""
from dateutil import parser
import requests
import matplotlib.pyplot as plt
import datetime
from collections import OrderedDict, defaultdict
from math import radians, sin, cos, atan2, sqrt
# import string


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
    
    :param cross_ref_file: csv-separated file with [gravitystation name], [USGS 15-digit ID]
    :param gravity_station_ID: gravity station name (e.g., RM109)
    :return: Dictionary with fields 'continuous_x', 'continuous_y', 'discrete_x', and 'discrete 'y'
    """
    out_dic = dict()
    name_dic = dict()
    # rdb_meas retrieval is preferred, it returns both discrete and continuous measurements.
    today = datetime.datetime.today().strftime('%Y-%m-%d')
    nwis_URL = 'http://nwis.waterdata.usgs.gov/nwis/dv?cb_72019=on&format=rdb_meas' + \
               f'&site_no={nwis_ID}' + \
               '&referred_module=gw&period=&begin_date=1999-10-01&end_date=%s' % today
    # print('Retrieving rdb data for {} from {}'.format(grav_ID, nwis_URL))
    r = requests.get(nwis_URL)
    # If there is continuous data, it will start with '# ----... WARNING ---...'
    nwis_data = r.text.split('\n')
    if nwis_data[-2].split()[0] != "agency_cd":
        c_x, c_y, d_x, d_y = parse_continuous_data(nwis_data)
    else:
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
    degree_buffer = 0.1
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
                if line_elems[0] == u'USGS' and line_elems[13] != '72019':  # Doesn't work for other agency codes
                    lat = float(line_elems[4])
                    lon = float(line_elems[5])
                    distance = calc_distance(lat,
                                                 lon,
                                                 coords[1],
                                                 coords[0])
                    ID_dict[line_elems[1]] = [line_elems[2],
                                              f"{distance:.0f}",
                                              line_elems[21],
                                              line_elems[22]]
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


if __name__ == "__main__":
    plot_wells('SiteIDcrossref.csv', ['T2-S2','AAC-17','PK-1','T1-S5','324421114482101'])

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