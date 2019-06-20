import pytest
import pytestqt
from PyQt5 import QtCore
from test_fixture_pyqt import mainprog

from gui_objects import CoordinatesTable, AdjustOptions

def test_coordinates_dialog(qtbot):
    coords = dict()
    coords['1'] = (-101.0, 30.0, 1600.0)
    coords['2'] = (-102.1, 31.77, 1600.0)
    coords['3'] = (-103.0, 32.0, 1600.0)
    coords['4'] = (-104.0, 33.25, 1600.0)
    coords['5'] = (-105.4, 34.63, 1600.0)
    coords['6'] = (-106.0, 35.0, 1600.0)

    ct = CoordinatesTable(coords)
    assert type(ct) == CoordinatesTable
    assert ct.coords() == coords
    ct.show()
    ct.table.selectAll()
    qtbot.keyPress(ct, QtCore.Qt.Key_C, QtCore.Qt.ControlModifier)
    qtbot.wait(100)
    assert ct.sys_clip.text()[0:10] == '1\t30.0\t-10'
    new_text = ct.sys_clip.text()
    new_text = new_text.replace('30', '40')
    ct.sys_clip.setText(new_text)
    qtbot.keyPress(ct, QtCore.Qt.Key_V, QtCore.Qt.ControlModifier)
    qtbot.wait(100)
    assert ct.sys_clip.text()[0:10] == '1\t40.0\t-10'
    ct.close()

def test_adjustmentoptions_dialog(qtbot, mainprog):
    obstreesurvey = mainprog.obsTreeModel.invisibleRootItem().child(0)
    ao = AdjustOptions('test_survey', obstreesurvey.adjustment.adjustmentoptions, parent=mainprog)
    ao.show()
    assert ao.sigma_add_chk.checkState() == 0
    ao.calc_cal_coeff_checked_or_unchecked(2)
    assert ao.cal_coeff_table.isEnabled() == False
    ao.calc_cal_coeff_checked_or_unchecked(0)
    ao.specify_cal_coeff_checked_or_unchecked(2)
    ao.apply_current()