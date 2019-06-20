import sys, os

code_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, code_path + '/../gsadjust')
sys.path.insert(0, code_path)

import pytest
import GSadjust
from data_import import import_abs_g_complete, import_abs_g_simple


def test_readfile_exceptions():
    bad_filename = './tests/test_doesnotexist.txt'
    bad_indexdata_filename = './tests/test_badindex_BurrisData.txt'
    bad_valuedata_filename = './tests/test_badvalue_BurrisData.txt'
    meter_type = 'Burris'
    with pytest.raises(IOError) as excinfo:
        data = GSadjust.MainProg.read_raw_data_file(bad_filename, meter_type)
    with pytest.raises(IndexError) as excinfo:
        data = GSadjust.MainProg.read_raw_data_file(bad_indexdata_filename, meter_type)
    with pytest.raises(ValueError) as excinfo:
        data = GSadjust.MainProg.read_raw_data_file(bad_valuedata_filename, meter_type)

def test_read_Burris():
    filename = './tests/test_BurrisData.txt'
    meter_type = 'Burris'
    data = GSadjust.MainProg.read_raw_data_file(filename, meter_type)
    assert len(data.dial) == 497
    assert len(data.raw_grav) == 497
    first_station = 'rg37'
    for idx, el in enumerate(data.station):
        if el != first_station:
            first_station_data = data.raw_grav[:idx]
            break
    mean = sum(first_station_data) / len(first_station_data)
    # This value calculated independently in Excel
    assert abs(mean - 2775777) < 0.1

def test_read_ScintrexCG6():
    filename = './tests/test_ScintrexCG5Data.txt'
    meter_type = 'Scintrex'
    data = GSadjust.MainProg.read_raw_data_file(filename, meter_type)
    first_station = '1'
    assert len(data.dial) == 0
    assert len(data.keepdata) == 2096
    assert len(data.raw_grav) == 2096
    for idx, el in enumerate(data.station):
        if el != first_station:
            first_station_data = data.raw_grav[:idx]
            break
    mean = sum(first_station_data) / len(first_station_data)
    # This value calculated independently in Excel
    assert abs(mean - 2639256.214) < 0.1

def testreadabsgcomplete():
    filename = './tests/test_Absg_complete.txt'
    datums = import_abs_g_complete(filename)
    assert type(datums) == list
    assert len(datums) == 5
    assert datums[0].sd == 13.44
    assert type(datums[1].date) == str
    assert 1==1

def test_importabsgsimple():
    filename = './tests/test_Absg_simple.txt'
    datums = import_abs_g_simple(filename)
    assert type(datums) == list
    assert len(datums) == 5
    assert datums[0].sd == 13.44
    assert type(datums[1].date) == str
