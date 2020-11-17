#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
gui/tabs/network.py
===============

PyQt graphical elements on the network adjustment tab of GSadjust.
--------------------------------------------------------------------------------

This software is preliminary, provisional, and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.
"""
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtCore import Qt


###########################################################################
# GSadjust adjustment tab
###########################################################################
# noinspection PyUnresolvedReferences
class TabAdjust(QtWidgets.QWidget):
    def __init__(self, parent):
        super(TabAdjust, self).__init__()
        # Delta table (top left)
        self.parent = parent
        self.delta_view = QtWidgets.QTableView()
        # self.delta_proxy_model = QtCore.QSortFilterProxyModel(self)
        # self.delta_view.setModel(self.delta_proxy_model)
        self.delta_view.setSortingEnabled(True)

        self.delta_view.setStyleSheet(
            """
            QTableView::indicator:checked{
               outline: 1px solid #1e5180;
               }"""
        )

        # Datum (Abs. g) table (bottom left)
        self.datum_view = QtWidgets.QTableView()
        self.datum_view.setContextMenuPolicy(Qt.CustomContextMenu)
        self.datum_view.customContextMenuRequested.connect(self.datum_context_menu)
        self.datum_proxy_model = QtCore.QSortFilterProxyModel(self)
        self.datum_view.setModel(self.datum_proxy_model)
        self.datum_view.setSortingEnabled(True)

        self.datum_popup_menu = QtWidgets.QMenu("Datum Popup Menu", self)
        self.mnDeleteDatum = QtWidgets.QAction('Delete datum', self)
        self.mnDeleteDatum.triggered.connect(self.parent.delete_datum)

        self.results_popup_menu = QtWidgets.QMenu("Results Popup Menu", self)
        self.mnCopyResults = QtWidgets.QAction('Copy to clipboard', self)
        self.mnCopyResults.triggered.connect(self.copy_results)

        # Results table (top right)
        self.results_view = QtWidgets.QTableView()
        self.results_view.setContextMenuPolicy(Qt.CustomContextMenu)
        self.results_view.customContextMenuRequested.connect(self.results_context_menu)
        self.results_proxy_model = QtCore.QSortFilterProxyModel(self)
        self.results_view.setModel(self.results_proxy_model)
        self.results_view.setSortingEnabled(True)

        # Main window
        layout_final = QtWidgets.QHBoxLayout()
        main_layout = QtWidgets.QSplitter(Qt.Horizontal, self)

        # Left subwindow
        layout_left = QtWidgets.QSplitter(Qt.Vertical, self)
        lbl = QtWidgets.QLabel("Relative-gravity differences (delta-g's)", self)
        lbl.setFixedHeight(30)
        layout_left.addWidget(lbl)
        layout_left.addWidget(self.delta_view)
        lbl = QtWidgets.QLabel("Datum observations", self)
        lbl.setFixedHeight(30)
        layout_left.addWidget(lbl)
        layout_left.addWidget(self.datum_view)
        main_layout.addWidget(layout_left)

        # Right subwindow
        layout_right = QtWidgets.QSplitter(Qt.Vertical)
        lbl = QtWidgets.QLabel("Adjusted station values")
        lbl.setFixedHeight(30)
        layout_right.addWidget(lbl)
        layout_right.addWidget(self.results_view)
        lbl = QtWidgets.QLabel("Least-squares statistics")
        lbl.setFixedHeight(30)
        layout_right.addWidget(lbl)
        self.stats_view = QtWidgets.QTableView()
        layout_right.addWidget(self.stats_view)
        main_layout.addWidget(layout_right)

        self.delta_view.horizontalHeader().sectionResized.connect(
            self.save_delta_column_widths_to_settings
        )
        self.datum_view.horizontalHeader().sectionResized.connect(
            self.save_datum_column_widths_to_settings
        )
        self.results_view.horizontalHeader().sectionResized.connect(
            self.save_results_column_widths_to_settings
        )

        layout_final.addWidget(main_layout)
        self.setLayout(layout_final)

    def invalidate_sort(self):
        self.datum_proxy_model.invalidate()

    def save_delta_column_widths_to_settings(self):
        col_widths = []
        for i in range(self.delta_view.model().columnCount()):
            col_widths.append(int(self.delta_view.columnWidth(i)))
        self.parent.settings.setValue('delta_table_column_widths', col_widths)

    def save_datum_column_widths_to_settings(self):
        col_widths = []
        for i in range(self.datum_view.model().columnCount()):
            col_widths.append(int(self.datum_view.columnWidth(i)))
        self.parent.settings.setValue('datum_table_column_widths', col_widths)

    def save_results_column_widths_to_settings(self):
        col_widths = []
        for i in range(self.results_view.model().columnCount()):
            col_widths.append(int(self.results_view.columnWidth(i)))
        self.parent.settings.setValue('results_table_column_widths', col_widths)

    def update_col_widths(self):
        cw = self.parent.settings.value('delta_table_column_widths')
        if cw:
            for i, w in enumerate(cw):
                self.delta_view.setColumnWidth(i, int(w))
        cw = self.parent.settings.value('datum_table_column_widths')
        if cw:
            for i, w in enumerate(cw):
                self.datum_view.setColumnWidth(i, int(w))
        cw = self.parent.settings.value('results_table_column_widths')
        if cw:
            for i, w in enumerate(cw):
                self.results_view.setColumnWidth(i, int(w))

    def datum_context_menu(self, point):
        selected = self.datum_view.selectedIndexes()
        if selected:
            self.datum_popup_menu.addAction(self.mnDeleteDatum)
            self.datum_popup_menu.exec_(self.datum_view.mapToGlobal(point))

    def results_context_menu(self, point):
        """
        Right-click context menu on results table
        :param point: PyQt reference to click point, determines where to show popup.
        """
        self.results_popup_menu.addAction(self.mnCopyResults)
        self.results_popup_menu.exec_(self.results_view.mapToGlobal(point))

    def copy_results(self):
        """
        Copies results table to clipboard
        """
        self.results_proxy_model.sourceModel().copyToClipboard()
