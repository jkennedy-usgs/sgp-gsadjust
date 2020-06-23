import pytest
import pytestqt
import os
from PyQt5 import QtCore
import numpy as np
from test_fixture_pyqt import mainprog

from data_analysis import compute_gravity_change
from gui_objects import (
    CoordinatesTable, AdjustOptions, SelectAbsg, DialogApplyTimeCorrection, DialogOverwrite,
    DialogMeterType, DialogLoopProperties, GravityChangeTable
)
from gsa_plots import PlotGravityChange

# def test_GravityChangeMap(qtbot):
#
#     gravity_change_map = GravityChangeMap(table, header, coords, full_table, full_header)


def test_PlotGravityChange(qtbot, mainprog):
    mainprog.workspace_clear()
    mainprog.workspace_open_json(os.path.join('test_data', 'complete_example2.gsa'))
    mainprog.adjust_network()
    data = compute_gravity_change(mainprog.obsTreeModel)
    change_table = GravityChangeTable(mainprog, data=data)
    change_plot = PlotGravityChange(change_table.dates, change_table.table)
    jeff = 1
    assert change_plot.figure is not None
    assert len(change_plot.figure.axes[0].get_children()[0].get_xdata()) == 2
    assert np.abs(change_plot.figure.axes[0].get_children()[0].get_ydata()[1] + 6.8) < 0.01

def test_ApplyTimeCorrection(qtbot):
    time_correction_dialog = DialogApplyTimeCorrection()

    def on_timeout():
        qtbot.keyClick(time_correction_dialog.msg, QtCore.Qt.Key_Enter)

    QtCore.QTimer.singleShot(1000, on_timeout)
    time_correction_dialog.msg.exec_()
    correction_type = time_correction_dialog.time_correction_type
    assert correction_type == 'station'

    time_correction_dialog = DialogApplyTimeCorrection()

    def on_timeout():
        qtbot.keyClick(time_correction_dialog.msg, QtCore.Qt.Key_Tab, delay=50)
        qtbot.keyClick(time_correction_dialog.msg, QtCore.Qt.Key_Tab, delay=50)
        qtbot.keyClick(time_correction_dialog.msg, QtCore.Qt.Key_Enter, delay=50)

    QtCore.QTimer.singleShot(1000, on_timeout)
    time_correction_dialog.msg.exec_()
    correction_type = time_correction_dialog.time_correction_type
    assert correction_type == 'survey'


def test_Overwrite(qtbot):
    overwrite_dialog = DialogOverwrite()

    def on_timeout():
        qtbot.keyClick(overwrite_dialog, QtCore.Qt.Key_Enter)

    QtCore.QTimer.singleShot(1000, on_timeout)
    result = overwrite_dialog.exec_()
    assert result == 0
    overwrite_dialog = DialogOverwrite()

    def on_timeout():
        qtbot.keyClick(overwrite_dialog, QtCore.Qt.Key_Tab, delay=50)
        qtbot.keyClick(overwrite_dialog, QtCore.Qt.Key_Enter, delay=50)

    QtCore.QTimer.singleShot(1000, on_timeout)
    result = overwrite_dialog.exec_()
    assert result == 1


def test_meter_type(qtbot, mainprog):
    meter_dialog = DialogMeterType()

    def on_timeout():
        qtbot.keyClick(meter_dialog, QtCore.Qt.Key_Enter)

    QtCore.QTimer.singleShot(1000, on_timeout)
    meter_dialog.exec_()
    assert meter_dialog.meter_type == 'CG5'
    meter_dialog = DialogMeterType()

    def on_timeout():
        qtbot.keyClick(meter_dialog, QtCore.Qt.Key_Tab, delay=50)
        qtbot.keyClick(meter_dialog, QtCore.Qt.Key_Tab, delay=50)
        qtbot.keyClick(meter_dialog, QtCore.Qt.Key_Enter, delay=50)

    QtCore.QTimer.singleShot(1000, on_timeout)
    meter_dialog.exec_()
    assert meter_dialog.meter_type == 'CG6Tsoft'


def test_loop_options(qtbot, mainprog):
    loops = [mainprog.obsTreeModel.invisibleRootItem().child(0).child(0)]
    options_dialog = DialogLoopProperties(loops)

    def on_timeout():
        qtbot.keyClick(options_dialog, QtCore.Qt.Key_Tab, delay=100)
        qtbot.keyClick(options_dialog, QtCore.Qt.Key_Tab, delay=100)
        qtbot.keyClicks(options_dialog.comment_edit, "USGS", delay=100)
        qtbot.mouseClick(options_dialog.ok_button, QtCore.Qt.LeftButton, delay=1000)

    QtCore.QTimer.singleShot(1000, on_timeout)
    options_dialog.exec_()

    for loop in loops:
        loop.oper = options_dialog.operator_edit.text()
        loop.meter = options_dialog.meter_edit.text()
        loop.comment = options_dialog.comment_edit.toPlainText()

    assert loop.oper == 'USGS'
    assert loop.comment == 'USGS'
    assert loop.meter == 'B44'


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


def test_selectabsg():
    sa = SelectAbsg('./test_data/field/Absolute/')
    sa.show()
    sa.load_button.click()
    assert sa.ProxyModel.data(sa.ProxyModel.index(0, 0), role=QtCore.Qt.CheckStateRole) == 0
    sa.ProxyModel.data(sa.ProxyModel.index(0, 0), role=QtCore.Qt.UserRole).checked = 2
    assert sa.ProxyModel.data(sa.ProxyModel.index(0, 0), role=QtCore.Qt.CheckStateRole) == 2
    sa.export_and_close()
    assert len(sa.new_datums) == 1
