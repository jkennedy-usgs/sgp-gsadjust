import datetime as dt

from matplotlib.dates import date2num

from data_objects import ChannelList


def read_csv(fh):
    i = 0
    meter, oper = None, None
    all_survey_data = ChannelList()

    _ = fh.readline()
    for line in fh:
        i += 1
        line = line.strip()
        # Skip blank and comment lines
        if (not line) or (line[0] == '%'):
            continue

        # parse string line first with respect to '/' caracters (used in the date format),
        # then with ':' (used for the time display), eventually with the classic ' '
        vals_temp1 = line.split('/')
        vals_temp2 = vals_temp1[0].split(':')
        vals_temp3 = vals_temp2[0].split(',')
        vals_temp4 = vals_temp2[2].split(',')

        # fill object properties:
        all_survey_data.station.append(vals_temp3[0].strip())
        all_survey_data.lat.append(float(vals_temp3[1]))
        all_survey_data.long.append(float(vals_temp3[2]))
        all_survey_data.elev.append(float(vals_temp3[3]))
        all_survey_data.raw_grav.append(float(vals_temp3[4]) * 1000. -
                                        float(vals_temp3[9]) * 1000.)  # convert to microGal; remove etc
        all_survey_data.tare.append(0)
        all_survey_data.sd.append(float(vals_temp3[5]) * 1000.)
        all_survey_data.tiltx.append(float(vals_temp3[6]))
        all_survey_data.tilty.append(float(vals_temp3[7]))
        all_survey_data.temp.append(float(vals_temp3[8]))
        all_survey_data.etc.append(float(vals_temp3[9]) * 1000.)
        all_survey_data.meter_etc.append(float(vals_temp3[9]) * 1000.)
        all_survey_data.dur.append(int(vals_temp3[10]))
        all_survey_data.rej.append(int(vals_temp3[11]))
        all_survey_data.t.append(date2num(dt.datetime(int(vals_temp1[-1]),
                                                      int(vals_temp4[-1]),
                                                      int(vals_temp1[-2]),
                                                      int(vals_temp3[-1]),
                                                      int(vals_temp2[-2]),
                                                      int(vals_temp4[0]))))
        if meter:
            all_survey_data.meter.append(meter)
        else:
            all_survey_data.meter.append('-999')
        if oper:
            all_survey_data.oper.append(oper)
        else:
            all_survey_data.oper.append('-999')
        all_survey_data.keepdata.append(1)
    all_survey_data.meter_type = 'csv'
    return all_survey_data


def read_scintrex(fh):
    i = 0
    meter, oper = None, None
    all_survey_data = ChannelList()

    for line in fh:
        i += 1

        # Clean line
        line = line.strip()

        # Skip blank and comment lines
        if (not line) or (line[0] == 'L'):
            continue

        # Header line; look for useful information
        if line[0] == '/':
            vals_temp = line.split()
            if len(vals_temp) > 1:
                if vals_temp[1] == 'Instrument':
                    meter = vals_temp[-1]
                if vals_temp[1] == 'Operator:':
                    oper = vals_temp[-1]
            continue

        # parse string line first with respect to '/' caracters (used in the date format),
        # then with ':' (used for the time display), eventually with the classic ' '
        vals_temp1 = line.split('/')
        vals_temp2 = vals_temp1[0].split(':')
        vals_temp3 = vals_temp2[0].split()
        vals_temp4 = vals_temp2[2].split()

        # fill object properties:
        all_survey_data.line.append(float(vals_temp3[0]))
        s = vals_temp3[1].replace('.0000000', '')
        all_survey_data.station.append(s.strip())
        all_survey_data.elev.append(float(vals_temp3[2]))
        all_survey_data.raw_grav.append(float(vals_temp3[3]) * 1000. -
                                        float(vals_temp3[8]) * 1000.)  # convert to microGal; remove tide correction
        all_survey_data.tare.append(0)
        all_survey_data.sd.append(float(vals_temp3[4]) * 1000.)
        all_survey_data.tiltx.append(float(vals_temp3[5]))
        all_survey_data.tilty.append(float(vals_temp3[6]))
        all_survey_data.temp.append(float(vals_temp3[7]))
        all_survey_data.etc.append(float(vals_temp3[8]) * 1000.)
        all_survey_data.meter_etc.append(float(vals_temp3[8]) * 1000.)
        all_survey_data.dur.append(int(vals_temp3[9]))
        all_survey_data.rej.append(int(vals_temp3[10]))
        all_survey_data.t.append(date2num(dt.datetime(int(vals_temp4[3]),
                                                      int(vals_temp1[1]),
                                                      int(vals_temp1[2]),
                                                      int(vals_temp3[11]),
                                                      int(vals_temp2[1]),
                                                      int(vals_temp4[0]))))
        if meter:
            all_survey_data.meter.append(meter)
        else:
            all_survey_data.meter.append('-999')
        if oper:
            all_survey_data.oper.append(oper)
        else:
            all_survey_data.oper.append('-999')
        all_survey_data.keepdata.append(1)
    all_survey_data.meter_type = 'Scintrex'
    return all_survey_data


def read_burris(fh):
    i = 0
    all_survey_data = ChannelList()

    for line in fh:
        i += 1
        line = line.strip()
        vals_temp = line.split()
        # Numbers are columns in the imported file
        c_station, c_oper, c_meter, c_date, c_time = 0, 1, 2, 3, 4
        c_dial, c_feedback, c_tide, c_tilt = 6, 7, 8, 9
        c_grav, c_elev, c_lat, c_long = 5, 13, 14, 15
        if len(vals_temp) == 15:  # no meter operator specified
            c_dial -= 1
            c_feedback -= 1
            c_tide -= 1
            c_tilt -= 1
            c_meter -= 1
            c_date -= 1
            c_time -= 1
            c_grav -= 1
            c_elev -= 1
            c_lat -= 1
            c_long -= 1
            all_survey_data.oper.append('None')
        else:
            all_survey_data.oper.append(vals_temp[c_oper])

        date_temp = vals_temp[c_date].split('/')
        if int(date_temp[2]) > 999:
            date_temp = [date_temp[2], date_temp[0], date_temp[1]]
        time_temp = vals_temp[c_time].split(':')

        # fill object properties:
        all_survey_data.station.append(vals_temp[0].strip())
        all_survey_data.elev.append(float(vals_temp[c_elev]))
        all_survey_data.lat.append(float(vals_temp[c_lat]))
        all_survey_data.long.append(float(vals_temp[c_long]))
        # remove Earth tide correction; it's added in using the @grav property
        all_survey_data.raw_grav.append(float(vals_temp[c_grav]) * 1000. -
                                        float(vals_temp[c_tide]) * 1000.)
        all_survey_data.tare.append(0)
        all_survey_data.etc.append(float(vals_temp[c_tide]) * 1000.)
        all_survey_data.meter_etc.append(float(vals_temp[c_tide]) * 1000.)
        all_survey_data.dial.append(float(vals_temp[c_dial]))
        all_survey_data.feedback.append(float(vals_temp[c_feedback]))
        all_survey_data.sd.append(-999)  # Burris doesn't ouput SD, tiltx, tilty
        all_survey_data.meter.append(vals_temp[c_meter])
        all_survey_data.tiltx.append(float(vals_temp[c_tilt]) * 1000.)
        all_survey_data.tilty.append(0.)
        all_survey_data.temp.append(0.)
        all_survey_data.dur.append(5)
        all_survey_data.rej.append(5)
        all_survey_data.t.append(date2num(dt.datetime(int(date_temp[0]),
                                                      int(date_temp[1]),
                                                      int(date_temp[2]),
                                                      int(time_temp[0]),
                                                      int(time_temp[1]),
                                                      int(time_temp[2]))))
        all_survey_data.keepdata.append(1)
    all_survey_data.meter_type = 'Burris'
    return all_survey_data


def read_cg6(fh):
    i = 0
    meter, oper = None, None
    all_survey_data = ChannelList()

    for line in fh:
        i += 1
        line = line.strip()
        vals_temp = line.split()
        if line[0] == '/':
            vals_temp = line.split()
            if len(vals_temp) > 1:
                if vals_temp[1] == 'Instrument':
                    meter = vals_temp[-1]
                if vals_temp[1] == 'Operator:':
                    oper = vals_temp[-1]
            continue
        # Numbers are columns in the imported file
        c_station, c_date, c_time, c_sd = 0, 1, 2, 5
        c_tiltx, c_tilty = 8, 9
        c_tide, c_tilt, c_temp = 11, 12, 13
        c_grav, c_elev, c_lat, c_long = 3, 19, 17, 18

        date_temp = vals_temp[c_date].split('-')
        time_temp = vals_temp[c_time].split(':')

        # fill object properties:
        all_survey_data.line.append(0.)
        all_survey_data.station.append(vals_temp[0].strip())
        all_survey_data.elev.append(float(vals_temp[c_elev]))
        all_survey_data.lat.append(float(vals_temp[c_lat]))
        all_survey_data.long.append(float(vals_temp[c_long]))
        all_survey_data.raw_grav.append(float(vals_temp[c_grav]) * 1000.)
        all_survey_data.tare.append(0)
        all_survey_data.etc.append(float(vals_temp[c_tide]) * 1000.)
        all_survey_data.meter_etc.append(float(vals_temp[c_tide]) * 1000.)
        all_survey_data.sd.append(float(vals_temp[c_sd]) * 1000.)
        all_survey_data.meter.append(meter)
        all_survey_data.tiltx.append(float(vals_temp[c_tiltx]) * 1000.)
        all_survey_data.tilty.append(float(vals_temp[c_tilty]) * 1000.)
        all_survey_data.temp.append(float(vals_temp[c_temp]) * 1000.)
        all_survey_data.dur.append(5)
        all_survey_data.rej.append(5)
        all_survey_data.t.append(date2num(dt.datetime(int(date_temp[0]),
                                                      int(date_temp[1]),
                                                      int(date_temp[2]),
                                                      int(time_temp[0]),
                                                      int(time_temp[1]),
                                                      int(time_temp[2]))))
        if meter:
            all_survey_data.meter.append(meter)
        else:
            all_survey_data.meter.append('-999')
        if oper:
            all_survey_data.oper.append(oper)
        else:
            all_survey_data.oper.append('-999')
        all_survey_data.keepdata.append(1)
    all_survey_data.meter_type = 'CG6'
    return all_survey_data


def read_cg6tsoft(fh):
    i = 0
    meter, oper = None, None
    all_survey_data = ChannelList()

    for line in fh:
        i += 1
        line = line.strip()
        vals_temp = line.split()
        station_name = None
        if line[0] == '/':
            vals_temp = line.split()
            if len(vals_temp) > 1:
                if vals_temp[1] == 'Instrument':
                    meter = vals_temp[-1]
                if vals_temp[1] == 'Operator:':
                    oper = vals_temp[-1]
                if vals_temp[1] == 'Station:':
                    station_name = vals_temp[-1]
            continue

        # Numbers are columns in the imported file
        c_year, c_month, c_day, c_hour, c_minute, c_second = 0, 1, 2, 3, 4, 5
        c_tiltx, c_tilty = 12, 13
        c_tide, c_tilt, c_temp = 19, 18, 14
        c_lat, c_long, c_elev = 7, 8, 9
        c_grav = 11  # CorrGravity channel

        # fill object properties:
        if station_name:
            all_survey_data.station.append(station_name)
            all_survey_data.elev.append(float(vals_temp[c_elev]))
            all_survey_data.lat.append(float(vals_temp[c_lat]))
            all_survey_data.long.append(float(vals_temp[c_long]))
            all_survey_data.raw_grav.append(float(vals_temp[c_grav]) * 1000.)
            all_survey_data.tare.append(0)
            all_survey_data.etc.append(float(vals_temp[c_tide]) * 1000.)
            all_survey_data.meter_etc.append(float(vals_temp[c_tide]) * 1000.)
            all_survey_data.sd.append(-999)  # SD not exported in Tsoft format?? It is in regular format
            all_survey_data.meter.append(meter)
            all_survey_data.tiltx.append(float(vals_temp[c_tiltx]) * 1000.)
            all_survey_data.tilty.append(float(vals_temp[c_tilty]) * 1000.)
            all_survey_data.temp.append(float(vals_temp[c_temp]) * 1000.)
            all_survey_data.t.append(date2num(dt.datetime(int(vals_temp[c_year]),
                                                          int(vals_temp[c_month]),
                                                          int(vals_temp[c_day]),
                                                          int(vals_temp[c_hour]),
                                                          int(vals_temp[c_minute]),
                                                          int(vals_temp[c_second]))))
            if meter:
                all_survey_data.meter.append(meter)
            else:
                all_survey_data.meter.append('-999')
            if oper:
                all_survey_data.oper.append(oper)
            else:
                all_survey_data.oper.append('-999')
            all_survey_data.keepdata.append(1)
            all_survey_data.dur.append(-999)
            all_survey_data.rej.append(-999)
    all_survey_data.meter_type = 'CG6Tsoft'

    return all_survey_data
