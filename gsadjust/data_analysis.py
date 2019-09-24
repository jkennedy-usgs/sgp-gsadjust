import numpy as np
import logging
import datetime as dt

def compute_gravity_change(obstreemodel, full_table=False):
    """
    Shows PyQt table of gravity change.

    :param full_table: if True, entire tabular output file for data release is shown. Includes station
    coordinates, g, standard deviation, and gravity change.
    """

    # Check that all values are positive (it should work either way, but it avoids confusion)
    for i in range(obstreemodel.invisibleRootItem().rowCount()):
        survey = obstreemodel.invisibleRootItem().child(i)
        for ii in range(survey.results_model.rowCount()):
            adj_station = survey.results_model.data(survey.results_model.index(ii, 0),
                                                    role=256)  # 256=QtCore.Qt.UserRole
            if adj_station.g < 0:
                return False

    compare_station, initial_station, iteration_station, iteration_name = None, None, None, None
    logging.info('Calculating gravity change')
    first = True
    unique_station_names = set()
    unique_stations = list()
    for i in range(obstreemodel.invisibleRootItem().rowCount()):
        survey = obstreemodel.invisibleRootItem().child(i)
        for ii in range(survey.rowCount()):
            loop = survey.child(ii)
            for iii in range(loop.rowCount()):
                station = loop.child(iii)
                unique_station_names.add(station.station_name)
                unique_stations.append(station)
    unique_station_names = sorted(unique_station_names)
    out_table_iteration, out_table_cumulative = [], []
    header1, header2 = [], []
    lat, lon, elev, all_g = [], [], [], []
    if full_table:
        for station in unique_station_names:
            station_g = []
            g_header = []
            coords = obstreemodel.station_coords[station]
            lat.append(coords[0])
            lon.append(coords[1])
            elev.append(coords[2])
            for i in range(obstreemodel.invisibleRootItem().rowCount()):
                survey = obstreemodel.invisibleRootItem().child(i)
                station_found = False
                g_header.append(survey.name)
                g_header.append(survey.name + '_sd')
                for ii in range(survey.results_model.rowCount()):
                    adj_station = survey.results_model.data(survey.results_model.index(ii, 0),
                                                            role=256)  # 256=QtCore.Qt.UserRole
                    if adj_station.station[:6] == station[:6]:
                        station_found = True
                        break
                if station_found:
                    station_g.append('{:0.1f}'.format(adj_station.g))
                    station_g.append('{:0.1f}'.format(adj_station.sd))
                else:
                    station_g.append('-999')
                    station_g.append('-999')
            all_g.append(station_g)

    dates = []
    for i in range(obstreemodel.invisibleRootItem().rowCount()):
        survey = obstreemodel.invisibleRootItem().child(i)
        dates.append(dt.datetime.strptime(survey.name, '%Y-%m-%d'))
        diff_cumulative = []
        diff_iteration = []
        if full_table:
            diff_cumulative_sd, diff_iteration_sd = [], []
        if first:
            # Calculate both the between-survey change and the change from the initial survey
            initial_survey = survey.results_model
            iteration_reference = initial_survey
            reference_name = survey.name
            iteration_name = reference_name
            first = False
        else:
            if not full_table:
                header1.append(str(iteration_name) + '_to_' + str(survey.name))
                header2.append(str(reference_name) + '_to_' + str(survey.name))
            else:
                header1.append('dH2O_' + str(iteration_name) + '_to_' + str(survey.name))
                header1.append('dH2O_sd_' + str(iteration_name) + '_to_' + str(survey.name))
                header2.append('dH2O_' + str(reference_name) + '_to_' + str(survey.name))
                header2.append('dH2O_sd_' + str(reference_name) + '_to_' + str(survey.name))

            compare_survey = survey.results_model
            for station_name in unique_station_names:

                for ii in range(initial_survey.rowCount()):
                    initial_station = initial_survey.data(initial_survey.index(ii, 0),
                                                          role=256)  # 256=QtCore.Qt.UserRole
                    # Iterate through, look for matching station. 'if' statements deal with Gravnet, which truncates
                    # station names to 6 characters
                    if len(initial_station.station) > 6 and len(station_name) > 6:
                        if initial_station.station == station_name:
                            break
                    elif len(initial_station.station) > 6:
                        if initial_station.station[:6] == station_name:
                            break
                    elif initial_station.station == station_name:
                        break
                    else:
                        initial_station = None
                for ii in range(iteration_reference.rowCount()):
                    iteration_station = iteration_reference.data(iteration_reference.index(ii, 0),
                                                                 role=256)  # 256=QtCore.Qt.UserRole
                    # Iterate through, look for matching station
                    # if iteration_station.station == station_name[:6]:
                    #     break
                    # else:
                    #     iteration_station = None

                    if len(iteration_station.station) > 6 and len(station_name) > 6:
                        if iteration_station.station == station_name:
                            break
                    elif len(iteration_station.station) > 6:
                        if iteration_station.station[:6] == station_name:
                            break
                    elif iteration_station.station == station_name:
                        break
                    else:
                        iteration_station = None

                for ii in range(compare_survey.rowCount()):
                    compare_station = compare_survey.data(compare_survey.index(ii, 0),
                                                          role=256)  # 256=QtCore.Qt.UserRole
                    if len(compare_station.station) > 6 and len(station_name) > 6:
                        if compare_station.station == station_name:
                            break
                    elif len(compare_station.station) > 6:
                        if compare_station.station[:6] == station_name:
                            break
                    elif compare_station.station == station_name:
                        break
                    else:
                        compare_station = None

                if initial_station is not None and compare_station is not None:
                    if not full_table:
                        diff_cumulative.append("{:0.1f}".format(compare_station.g - initial_station.g))
                    else:
                        diff_cumulative.append("{:0.2f}".format((compare_station.g - initial_station.g) / 41.9))
                        var = np.sqrt(compare_station.sd ** 2 + initial_station.sd ** 2) / 41.9
                        if np.isnan(var):
                            diff_cumulative_sd.append('-999')
                        else:
                            diff_cumulative_sd.append("{:0.2f}".format(var))
                else:
                    diff_cumulative.append("-999")
                    if full_table:
                        diff_cumulative_sd.append("-999")  # for sd column
                if iteration_station is not None and compare_station is not None:
                    if not full_table:
                        diff_iteration.append("{:0.1f}".format(compare_station.g - iteration_station.g))
                    else:
                        diff_iteration.append("{:0.2f}".format((compare_station.g - iteration_station.g) / 41.9))
                        var = np.sqrt(compare_station.sd ** 2 + iteration_station.sd ** 2) / 41.9
                        if np.isnan(var):
                            diff_iteration_sd.append('-999')
                        else:
                            diff_iteration_sd.append("{:0.2f}".format(var))
                else:
                    diff_iteration.append("-999")
                    if full_table:
                        diff_iteration_sd.append("-999")  # for sd column
            out_table_iteration.append(diff_iteration)
            out_table_cumulative.append(diff_cumulative)
            if full_table:
                out_table_iteration.append(diff_iteration_sd)
                out_table_cumulative.append(diff_cumulative_sd)
            iteration_reference = compare_survey
            iteration_name = survey.name
    out_table = [list(unique_station_names)] + out_table_iteration + out_table_cumulative

    if not full_table:
        header = ['station'] + header1 + header2
        table = out_table
        return (header, table, dates)
    else:
        header = ['Station', 'Longitude', 'Latitude', 'Elevation'] + g_header + header1 + header2
        # transpose table
        g = [list(i) for i in zip(*all_g)]
        table = [unique_station_names, lat, lon, elev]
        table += g
        table += out_table_iteration
        table += out_table_cumulative
        # transpose back
        table = [list(i) for i in zip(*table)]
        # table = [header] + table
        return (header, table, dates)


def adjusted_vs_observed_datum_analysis(self):
    """
    Leave one out analysis. For each datum station, this repeats the network adjustment for each survey, but with
    the datum station excluded. Results are sent to plot_LOO_analysis, which generates a line plot of the measured
    datum time series and the corresponding adjusted time series.

    TODO: This isn't quite right, it should be comparing (at each datum station) the adjusted value with the station
    included in the adjustment, with the adjusted value with the station excluded.
    :return:
    """
    # Get list of datum stations for all surveys
    pbar = ProgressBar(total=self.obsTreeModel.invisibleRootItem().rowCount(), textmess='Adjusting surveys')
    pbar.show()

    datums = self.obsTreeModel.datums()
    # Loop over surveys
    x_all, adj_g_all, obs_g_all = [], [], []
    ctr = 0
    for station in datums:
        pbar.progressbar.setValue(ctr)
        ctr += 1
        QtWidgets.QApplication.processEvents()
        station_adj_g = []
        station_obs_g = []
        x_data = []

        # Iterate through surveys
        for i in range(self.obsTreeModel.invisibleRootItem().rowCount()):

            # If datum is in survey, uncheck it and do inversion
            obstreesurvey = self.obsTreeModel.invisibleRootItem().child(i)
            done = False
            adj_g = None
            for ii in range(obstreesurvey.datum_model.rowCount()):
                if not done:
                    idx = obstreesurvey.datum_model.index(ii, 0)
                    datum = obstreesurvey.datum_model.data(idx, role=QtCore.Qt.UserRole)
                    if datum.station == station:
                        checkstate = obstreesurvey.datum_model.data(idx, role=QtCore.Qt.CheckStateRole)
                        obstreesurvey.datum_model.setData(idx, 0, QtCore.Qt.CheckStateRole)
                        if self.menus.mnAdjPyLSQ.isChecked():
                            adj_type = 'PyLSQ'
                        else:
                            adj_type = 'Gravnet'
                        obstreesurvey.run_inversion(adj_type)
                        # Restore check state
                        obstreesurvey.datum_model.setData(idx, checkstate, QtCore.Qt.CheckStateRole)
                        for iii in range(obstreesurvey.results_model.rowCount()):
                            idx = obstreesurvey.results_model.index(iii, 0)
                            adj_station = obstreesurvey.results_model.data(idx, role=QtCore.Qt.UserRole)

                            if adj_station.station == station:
                                adj_g = adj_station.g + (datum.gradient * datum.meas_height)
                                done = True
                                break
            if adj_g:
                station_adj_g.append(adj_g)
                station_obs_g.append(datum.g)
                x_data.append(dt.datetime.strptime(obstreesurvey.name, '%Y-%m-%d'))

        # Store results
        x_all.append(x_data)
        adj_g_all.append(station_adj_g)
        obs_g_all.append(station_obs_g)
    self.plot_LOO_analysis(x_all, adj_g_all, obs_g_all, datums)