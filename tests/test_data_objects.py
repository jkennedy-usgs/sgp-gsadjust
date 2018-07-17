from test_fixture_pyqt import test_channellist_fixture, test_twostations_fixture
from data_objects import Delta
import pytest


@pytest.mark.parametrize('test_channellist_fixture, mean_g, station_length, station_subset',[
                         ('channellist_burris.p', 3092327.454, 11, 10),
                         ('channellist_scintrex.p', 2639321.583, 2096, 11)
                         ], indirect=['test_channellist_fixture'])
def test_channellist(test_channellist_fixture, mean_g, station_length, station_subset):
    assert len(test_channellist_fixture.raw_grav) == station_length
    test_data1 = test_channellist_fixture.extract_subset_idx(0, station_subset)
    assert len(test_data1.raw_grav) - 1 == station_subset
    mean = sum(test_data1.raw_grav) / len(test_data1.raw_grav)
    assert abs(mean - mean_g) < 0.001

@pytest.mark.parametrize('test_twostations_fixture, true_delta_g, true_sd',[
                         ('channellist_scintrex.p', 2159.339, 0.3387),
                         ], indirect=['test_twostations_fixture'])
def test_delta_normal(test_twostations_fixture, true_delta_g, true_sd):
    sta1 = test_twostations_fixture[0]
    sta2 = test_twostations_fixture[1]
    # Test that delta with default values is the same as one with defaults specified
    d1 = Delta(sta1, sta2, driftcorr=0., ls_drift=None, delta_type='normal', checked=2,adj_sd=999)
    d2 = Delta(sta1, sta2)
    assert d1.__dict__ == d2.__dict__
    assert abs(d1.dg - true_delta_g) < 0.001
    assert abs(d1.sd - true_sd) < 0.001

# def test_channellist_iter(channellist):
#     for field in
