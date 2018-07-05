import pytest
import GSadjust


def test_readfile_exceptions():
    bad_filename = './test_doesnotexist.txt'
    bad_indexdata_filename = './test_badindex_BurrisData.txt'
    bad_valuedata_filename = './test_badvalue_BurrisData.txt'
    meter_type = 'Burris'
    with pytest.raises(IOError) as excinfo:
        data = GSadjust.MainProg.read_raw_data_file(bad_filename, meter_type)
    with pytest.raises(IndexError) as excinfo:
        data = GSadjust.MainProg.read_raw_data_file(bad_indexdata_filename, meter_type)
    with pytest.raises(ValueError) as excinfo:
        data = GSadjust.MainProg.read_raw_data_file(bad_valuedata_filename, meter_type)

def test_read_Burris():
    filename = './test_BurrisData.txt'
    meter_type = 'Burris'
    # try:
    data = GSadjust.MainProg.read_raw_data_file(filename, meter_type)
    # except Exception as e:
    #     print('Exception')

    assert len(data.dial) == 497
    assert len(data.raw_grav) == 497
    assert len(data.corr_g) == 0