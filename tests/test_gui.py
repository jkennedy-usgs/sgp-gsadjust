import pytest
import pytestqt
import GSadjust
from time import sleep
from PyQt5 import QtCore, Qt

def test_gui(qtbot):
    window = GSadjust.MainProg()
    # window.show()
    qtbot.addWidget(window)
    window.show()
    qtbot.wait(1000)
    window.open_raw_data(r"E:\Shared\current\python\sgp-gsadjust\tests\test_BurrisData.txt", 'Burris')
    assert window.obsTreeModel.rowCount() == 1
