"""
Many of these functions involve looping over the survey > loop > station data.

It's tempting to replace the for i in range... statements with __iter__ functions in obstreemodel, obstreesurvey, etc.

I tried this and it kind of worked, but it interfered with appendRows from PyQt = its apparently expecting the
default iterator which is presumable inherited from QtGui.QStandardItem.
"""

def init_cal_coeff_dict(obstreemodel):
    try:
        meter_list = {}
        for i in range(obstreemodel.invisibleRootItem().rowCount()):
            survey = obstreemodel.invisibleRootItem().child(i)
            for ii in range(survey.rowCount()):
                loop = survey.child(ii)
                if loop.meter not in meter_list:
                    meter_list[loop.meter] = 1.000
        return meter_list
    except Exception:
        return None

def return_delta_given_key(key, deltas):
    """
    Return a delta based on the station name and station count of the comprising stations.

    DEPRECATED, ONLY NEEDED TO OPEN OLD .P FILES

    :param keys: a tuple or list, depending on the type of delta
    :param deltas: List of deltas, as returned from assemble_all_deltas()
    :return:
    """
    for delta in deltas:
        # Sometimes they differ by a little bit (1e-8)
        if delta.key[0:len(key) - 5] == key[0:len(key) - 5]:
            return delta
    return None


def assemble_all_deltas(obstreemodel):
    """
    Get all deltas from loop delta models.

    DEPRECATED, ONLY NEEDED TO OPEN OLD .P FILES


    :return: One long list of all deltas in a campaign.
    """
    deltas = []
    for i in range(obstreemodel.invisibleRootItem().rowCount()):
        survey = obstreemodel.invisibleRootItem().child(i)
        for ii in range(survey.rowCount()):
            loop = survey.child(ii)
            for iii in range(loop.delta_model.rowCount()):
                delta = loop.delta_model.data(loop.delta_model.index(iii, 0), role=32)  # 32=QtCore.Qt.UserRole
                deltas.append(delta)
    return deltas

def init_station_coords_dict(obstreemodel):
    """
    Stores a single set of coordinates for each station with the obsTreeModel object. The coordinates of the last
    Station in the Survey > Loop > Station hierarchy will be used.
    """
    station_coords = dict()
    for i in range(obstreemodel.invisibleRootItem().rowCount()):
        survey = obstreemodel.invisibleRootItem().child(i)
        for ii in range(survey.rowCount()):
            loop = survey.child(ii)
            for iii in range(loop.rowCount()):
                station = loop.child(iii)
                station_coords[station.station_name] = (station.long[0], station.lat[0], station.elev[0])
    return station_coords

