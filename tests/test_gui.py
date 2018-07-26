import pytest
import pytestqt
import GSadjust
import os
import gui_objects
from time import sleep
from PyQt5 import QtCore, Qt

def test_gui(qtbot, monkeypatch):
    # Not sure why, but need to store and restore the path after this test
    pwd = os.getcwd()
    window = GSadjust.MainProg()
    # window.show()
    qtbot.addWidget(window)
    window.show()
    qtbot.wait(1000)

    # Open data
    window.open_raw_data(r"E:\Shared\current\python\sgp-gsadjust\tests\test_BurrisData2.txt", 'Burris')
    assert window.obsTreeModel.rowCount() == 1

    # Divide into loops
    window.divide_survey(8/24)
    qtbot.wait(2000)
    survey = window.obsTreeModel.invisibleRootItem().child(0)
    loop = survey.child(0)
    assert loop.rowCount() == 12
    assert survey.rowCount() == 3

    # Change to Drift tab
    window.tab_widget.setCurrentIndex(1)

    # Step through loops
    window.activate_survey_or_loop(survey.child(1).index())
    qtbot.wait(1000)
    window.activate_survey_or_loop(survey.child(2).index())
    qtbot.wait(1000)
    window.tab_widget.setCurrentIndex(2)

    # Populate delta and datum tables
    qtbot.keyClick(window, 'a', modifier=QtCore.Qt.ControlModifier)
    qtbot.wait(1000)
    monkeypatch.setattr(gui_objects.AddDatumFromList, 'add_datum', classmethod(lambda *args: 'CDOT'))
    qtbot.keyClick(window, 'd', modifier=QtCore.Qt.ControlModifier)
    window.adjust_network()

    # Verify gravnet input
    assert survey.results_model.rowCount() == 30
    assert len(survey.adjustment.adjustmentresults.text) == 14  # number of lines in Numpy output
    assert survey.adjustment.adjustmentresults.n_deltas == 83
    assert survey.adjustment.adjustmentresults.n_datums == 1

    test_workspace = 'test1.p'
    success = window.obsTreeModel.save_workspace(test_workspace)
    assert success == True

    window.workspace_clear()
    assert window.obsTreeModel.rowCount() == 0

    window.workspace_open(test_workspace)
    survey = window.obsTreeModel.invisibleRootItem().child(0)
    loop = survey.child(0)
    assert loop.rowCount() == 12
    assert survey.rowCount() == 3

    window.menus.mnAdjPyLSQ.setChecked(True)
    window.adjust_network()

    window.adjust_network()
    for line in survey.adjustment.adjustmentresults.text:
        elems = line.split(' ')
        if elems[0] == 'SD':
            sd0 = float(elems[-1])

    # Disable some observations, save workspace, clear and reload, verify that we get the same adjustment results
    rows = [1,3,5]
    for row in rows:
        survey.delta_model.setData(survey.delta_model.index(row, 0), QtCore.Qt.Unchecked, role=QtCore.Qt.CheckStateRole)

    window.adjust_network()
    for line in survey.adjustment.adjustmentresults.text:
        elems = line.split(' ')
        if elems[0] == 'SD':
            sd1 = float(elems[-1])

    # Adjustment results should be different with some observations disabled
    assert abs(sd0 - sd1) > 0.01

    success = window.workspace_save()
    assert success == True

    window.workspace_clear()
    assert window.obsTreeModel.rowCount() == 0

    window.workspace_open(test_workspace)
    window.adjust_network()
    for line in survey.adjustment.adjustmentresults.text:
        elems = line.split(' ')
        if elems[0] == 'SD':
            sd2 = float(elems[-1])

    assert abs(sd1 - sd2) < 0.000001
    os.remove(test_workspace)
    window.close()
    os.chdir(pwd)