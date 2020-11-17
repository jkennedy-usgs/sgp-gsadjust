from gsadjust.obstree import ObsTreeSurvey, ObsTreeLoop, ObsTreeStation, ObsTreeModel

from PyQt5 import QtGui
from utils import *


def test_ObsTreeSurvey(obstreesurvey):
    assert obstreesurvey.rowCount() == 1
    assert obstreesurvey.name == '2020-05-11'
    assert type(obstreesurvey.child(0)) == ObsTreeLoop
    assert obstreesurvey.child(0).rowCount() == 46
    assert obstreesurvey.populate_delta_model() is True


def test_create_ObsTreeStation(channellist):
    a = ObsTreeStation(channellist, 'teststa', '1')
    assert a.station_name == 'teststa'


def test_obstreemodel(qtmodeltester, channellist):
    model = ObsTreeModel()
    surv1 = ObsTreeSurvey('surv1')
    loop1 = ObsTreeLoop('loop1')
    sta1 = ObsTreeStation(channellist, 'sta1', '1')
    sta2 = ObsTreeStation(channellist, 'sta2', '2')
    sta3 = ObsTreeStation(channellist, 'sta3', '3')
    loop1.appendRow([sta1, QtGui.QStandardItem('a'), QtGui.QStandardItem('a')])
    loop1.appendRow([sta2, QtGui.QStandardItem('a'), QtGui.QStandardItem('a')])
    loop1.appendRow([sta3, QtGui.QStandardItem('a'), QtGui.QStandardItem('a')])
    surv1.appendRow([loop1, QtGui.QStandardItem('a'), QtGui.QStandardItem('a')])
    model.appendRow([surv1, QtGui.QStandardItem('a'), QtGui.QStandardItem('a')])

    qtmodeltester.check(model)

#
# def test_deltatablemodel(qtmodeltester, list_of_deltas):
#     model = pyqt_models.DeltaTableModel()
#     for delta in list_of_deltas:
#         model.insertRows(delta, 0)
#
#     qtmodeltester.check(model)


def test_initcalcoeffdict(obstreemodel):
    cal_dict = init_cal_coeff_dict(obstreemodel)
    assert type(cal_dict) is dict
    assert cal_dict['B44'] == 1.0


def test_init_coords(obstreemodel):
    coords = init_station_coords_dict(obstreemodel)
    assert len(coords) == 16
    assert coords['rg37'][0] == -106.5
    assert coords['rg37'][2] == 1600.0
