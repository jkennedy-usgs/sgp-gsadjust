from data_objects import Delta

def drift_roman(data, loop_name, time_threshold=None):
    # roman_dg_model = RomanTableModel()
    deltas, vert_lines = [], []
    station_list = [i.station_name for i in data]
    unique_stations = list(set(station_list))

    # store initial value at each station
    initial_g = dict()

    # Easiest case: all g values are relative to the initial g at that station
    if time_threshold is None:
        for station_name in unique_stations:
            for station in data:
                if station.station_name == station_name:
                    initial_g[station_name] = station.gmean()
                    break
    else:
        # If time_threshold is specified (checked in the GUI) we need to build a list of possible initial g values
        # for each station. Possible values are those occurring after a gap >= time_threshold (i.e., if there is a
        # gap, reset the initial g that's subtracted from the measurements. This makes the lines start at y = 0 on
        # the plots.
        #
        # Builds the dictionary:
        # initial_g{key:station_name value:(time, g)}
        for station_name in unique_stations:
            stations = []
            for station in data:
                if station.station_name == station_name:
                    stations.append(station)
            iter_stations = iter(stations)
            first_station = next(iter_stations)
            initial_xy = [(first_station.tmean(), first_station.gmean())]
            for station in iter_stations:
                if (station.tmean() - first_station.tmean()) * 1440 > time_threshold:
                    initial_xy.append((station.tmean(), station.gmean()))
                first_station = station
            initial_g[station_name] = initial_xy

    # For each station in data, loop over all the other stations looking for two observations that bracket the
    # first station
    for station in data:
        for other_station in unique_stations:
            # Ignore it if its the same station
            if other_station == station.station_name:
                continue
            else:
                # get all occurrences of the other station
                other_stations = [i for i in data if i.station_name == other_station]
                if len(other_stations) > 1:
                    iter_stations = iter(other_stations)
                    other1 = next(iter_stations)
                    for other2 in iter_stations:
                        # Check for 3-point configuration (2 observations at other station bracket the initial obs)
                        if other1.tmean() < station.tmean() < other2.tmean():
                            # Check that time_threshold is met, or not set
                            if time_threshold is None or \
                                    (other2.tmean() - other1.tmean()) * 1440 < time_threshold:
                                delta = Delta(station,
                                              (other1, other2),
                                              delta_type='three_point',
                                              loop=loop_name)
                                sta2_dg = other2.gmean() - other1.gmean()
                                # this is the drift correction
                                time_prorate = (station.tmean() - other1.tmean()) / (
                                        other2.tmean() - other1.tmean())
                                # Look for previous occupation at same station. If there is a break > time_threshold
                                # between the previous and current occupation, we need to account for the shift in
                                # initial g. Each station has a unique initial g (that might change,
                                # depending on the time_threshold).
                                if time_threshold is not None:
                                    initial_gees = initial_g[other1.station_name]
                                    other_initial_g = initial_gees[0][1]
                                    if len(initial_gees) > 1:
                                        for initial_xy in initial_gees[1:]:
                                            if other1.tmean() >= initial_xy[0]:
                                                other_initial_g = initial_xy[1]
                                    initial_gees = initial_g[station.station_name]
                                    station_initial_g = initial_gees[0][1]
                                    if len(initial_gees) > 1:
                                        for initial_xy in initial_gees[1:]:
                                            if station.tmean() >= initial_xy[0]:
                                                station_initial_g = initial_xy[1]
                                # Easy case: everything relative to the initial observation.
                                else:
                                    other_initial_g = initial_g[other1.station_name]
                                    station_initial_g = initial_g[station.station_name]

                                vert_lines.append([(station.tmean(), station.tmean()),
                                                   ((other_initial_g -
                                                     other1.gmean() -
                                                     (sta2_dg * time_prorate)) * -1,
                                                    station.gmean() - station_initial_g)])
                                deltas.append(delta)
                                # roman_dg_model.insertRows(delta, 0)
                        other1 = other2

    return deltas, vert_lines

    # If there is more than one delta-g between a given station pair, average them
    # Setup dict to store averages '(sta1, sta2)':[g]
    avg_dg = dict()
    unique_pairs = set()
    for i in range(roman_dg_model.rowCount()):
        delta = roman_dg_model.data(roman_dg_model.index(i, 0), role=QtCore.Qt.UserRole)
        delta_key1 = (delta.station1.station_name, delta.station2[0].station_name)
        delta_key2 = (delta.station2[0].station_name, delta.station1.station_name)
        if delta_key1 not in unique_pairs and delta_key2 not in unique_pairs:
            unique_pairs.add(delta_key1)
            avg_dg[delta_key1] = [delta]
            for ii in range(i + 1, roman_dg_model.rowCount()):
                testdelta = roman_dg_model.data(roman_dg_model.index(ii, 0), role=QtCore.Qt.UserRole)
                testdelta_key1 = (testdelta.station1.station_name, testdelta.station2[0].station_name)
                testdelta_key2 = (testdelta.station2[0].station_name, testdelta.station1.station_name)
                if delta_key1 == testdelta_key1 or delta_key1 == testdelta_key2:
                    avg_dg[delta_key1].append(testdelta)

    roman_avg_dg_model = DeltaTableModel()
    for station_pair in avg_dg.items():
        # just send list of deltas, not key (station info is already in the deltas)
        avg_delta = Delta.from_list(station_pair[1])
        roman_avg_dg_model.insertRows(avg_delta, 0)

    return roman_dg_model, roman_avg_dg_model, vert_lines