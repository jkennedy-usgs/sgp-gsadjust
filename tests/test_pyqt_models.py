import pytest
import pytestqt
import pyqt_models
import GSadjust
from test_fixture_pyqt import channellist, obstreesurvey
from PyQt5 import QtGui

def test_ObsTreeSurvey(obstreesurvey):
    assert obstreesurvey.rowCount() == 1
    assert obstreesurvey.name == 'test'
    assert type(obstreesurvey.child(0)) == pyqt_models.ObsTreeLoop
    assert obstreesurvey.child(0).rowCount() == 46
    assert obstreesurvey.populate_delta_model() == True

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