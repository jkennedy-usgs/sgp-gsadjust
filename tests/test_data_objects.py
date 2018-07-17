from test_fixture_pyqt import channellist


def test_channellist(channellist):
    test_data1 = channellist.extract_subset_idx(0, 28)
    assert len(test_data1.raw_grav) == 11
    mean = sum(test_data1.raw_grav) / len(test_data1.raw_grav)
    assert abs(mean - 3092327.454) < 0.001

# def test_channellist_iter(channellist):
#     for field in
