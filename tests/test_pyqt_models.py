import pytest
import pytestqt
import pyqt_models
from test_fixture_pyqt import channellist
from PyQt5 import QtGui

def test_create_ObsTreeSurvey():
    a = pyqt_models.ObsTreeSurvey('test')
    assert type(a) is pyqt_models.ObsTreeSurvey
    assert a.name == 'test'

def test_create_ObsTreeStation(channellist):
    a = pyqt_models.ObsTreeStation(channellist, 'teststa', '1')
    assert a.name == 'teststa'

def test_obstreemodel(qtmodeltester, channellist):
    model = pyqt_models.ObsTreeModel()
    surv1 = pyqt_models.ObsTreeSurvey('surv1')
    loop1 = pyqt_models.ObsTreeLoop('loop1')
    sta1 = pyqt_models.ObsTreeStation(channellist, 'sta1', '1')
    sta2 = pyqt_models.ObsTreeStation(channellist, 'sta2', '2')
    sta3 = pyqt_models.ObsTreeStation(channellist, 'sta3', '3')
    loop1.appendRow([sta1, QtGui.QStandardItem('a'), QtGui.QStandardItem('a')])
    loop1.appendRow([sta2, QtGui.QStandardItem('a'), QtGui.QStandardItem('a')])
    loop1.appendRow([sta3, QtGui.QStandardItem('a'), QtGui.QStandardItem('a')])
    surv1.appendRow([loop1, QtGui.QStandardItem('a'), QtGui.QStandardItem('a')])
    model.appendRow([surv1, QtGui.QStandardItem('a'), QtGui.QStandardItem('a')])

    qtmodeltester.check(model)