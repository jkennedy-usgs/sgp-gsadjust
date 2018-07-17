from test_fixture_pyqt import test_channellist_obj
from data_objects import Delta
import pytest


@pytest.mark.parametrize("test_channellist_obj, mean_g, station_length, station_subset",[
                         ('channellist_burris.p', 3092327.454, 11, 10),
                         ('channellist_scintrex.p', 2639321.583, 2096, 11)
                         ], indirect=['test_channellist_obj'])
def test_channellist(test_channellist_obj, mean_g, station_length, station_subset):
    # channelist_testobj = test_channellist_obj(test_channellist_obj)
    assert len(test_channellist_obj.raw_grav) == station_length
    test_data1 = test_channellist_obj.extract_subset_idx(0, station_subset)
    assert len(test_data1.raw_grav) - 1 == station_subset
    mean = sum(test_data1.raw_grav) / len(test_data1.raw_grav)
    assert abs(mean - mean_g) < 0.001

@pytest.mark
def test_delta():
    sta1 = Station1
    d = Delta(sta1, sta2, driftcorr=0., ls_drift=None, delta_type='normal', checked=2,
    adj_sd=999):

# def test_channellist_iter(channellist):
#     for field in
