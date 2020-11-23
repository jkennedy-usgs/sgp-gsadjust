"""
gui/messages.py
==============

GUI message boxes
-----------------

This software is preliminary, provisional, and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.
"""
import logging

from PyQt5.QtWidgets import QMessageBox


class MessageBox:
    @staticmethod
    def critical(title, message, **kwargs):
        logging.exception("%s: %s", title, message, exc_info=True)
        return QMessageBox.critical(None, title, message, **kwargs)

    @staticmethod
    def information(title, message, **kwargs):
        logging.info("%s: %s", title, message)
        return QMessageBox.information(None, title, message, **kwargs)

    @staticmethod
    def question(title, message, **kwargs):
        logging.info("%s: %s", title, message)
        result = QMessageBox.question(None, title, message, **kwargs)
        logging.info("[User answered] - %s", result)
        return result

    @staticmethod
    def warning(title, message, **kwargs):
        logging.warning("%s: %s", title, message, exc_info=True)
        return QMessageBox.warning(None, title, message, **kwargs)
