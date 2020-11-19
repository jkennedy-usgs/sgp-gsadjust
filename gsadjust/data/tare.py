"""
data/tare.py
===============

Represents an offset applied to the data

This software is preliminary, provisional, and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition tha
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.
"""


class Tare:
    """
    Tares are offsets in the meter reading.

    Tares apply to all data collected after the time of the Tare. Typically the time
    is assigned so that a Tare falls between two station occupations.

    Attributes
    ----------
    checked: int
        For PyQt tables. Checked = 2, Unchecked = 0
    datetime : datetime
        Python datetime, the time the Tare occurred
    tare : Float
        Tare value, in µGal

    """

    def __init__(self, datetime, tare):
        self.datetime = datetime
        self.tare = float(tare)
        self.checked = 2

    def __str__(self):
        if self.checked == 2:
            in_use = "x"
        else:
            in_use = "o"
        return f"{in_use} {self.datetime} {self.tare}"

    def to_json(self):
        """Method for serializing Tare.

        Python datetime objects can't be jsonified, therefore we convert to string.

        """
        return {
            "checked": self.checked,
            "datetime": str(self.datetime),
            "tare": self.tare,
        }
