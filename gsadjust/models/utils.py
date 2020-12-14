"""
models/utils.py
===============

PyQt model utility functions.
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


def format_numeric_column(column):
    """
    Format fn for simple numeric columns.
    Returns column (zero-indexed) +1 to give 1-indexed label.
    """
    return column + 1
