"""
gui_objects.py
==============

Miscellaneous GUI objects for GSadjust
--------------------------------------------------------------------------------

Major GUI objects (tabs, table views) are in the tab_... files. This module has
primarily pop-up dialogs used to set network adjustment settings, show gravity
change over time, etc. Major dialogs are written as classes and instantiated in
GSadjust.py. Minor dialogs are written as functions and called directly.

This software is preliminary, provisional, and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.
"""
from . import dialogs, menus, messages, tabs, utils, widgets
