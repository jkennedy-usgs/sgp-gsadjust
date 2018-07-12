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
    qtbot.keyClick(window, '1', modifier=QtCore.Qt.ControlModifier)
    qtbot.wait(1000)
    survey.msg.close()

    # Verify gravnet input
    assert survey.results_model.rowCount() == 30
    assert len(survey.adjustment.adjustmentresults.text) == 22
    assert survey.adjustment.adjustmentresults.n_deltas == 83
    assert survey.adjustment.adjustmentresults.n_datums == 1
    window.close()
    os.chdir(pwd)