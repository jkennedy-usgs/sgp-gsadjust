import pytest
import os
import gsadjust

from gsadjust.file.read import read_csv, read_cg6, read_cg6tsoft


def test_import_csv():
    data_file = os.path.join(
        os.getcwd(), 'test_data', 'field', 'CG-6', 'CG-6_TestData.dat'
    )
    with pytest.raises(IndexError):
        with open(data_file, 'r') as fh:
            data = read_csv(fh)
    data_file = os.path.join(os.getcwd(), 'tests', 'test_csv_data.csv')
    with open(data_file, 'r') as fh:
        data = read_csv(fh)
    assert len(data.raw_grav) == 95


def test_import_CG6():
    data_file = os.path.join(
        os.getcwd(), 'test_data', 'field', 'CG-6', 'CG-6_TestData.dat'
    )
    with open(data_file, 'r') as fh:
        data = read_cg6(fh)
    assert len(data.raw_grav) == 43

    data_file = os.path.join(os.getcwd(), 'tests', 'test_csv_data.csv')
    with pytest.raises(IndexError):
        with open(data_file, 'r') as fh:
            data = read_cg6(fh)


def test_import_CG6tsoft():
    data_file = os.path.join(
        os.getcwd(), 'test_data', 'field', 'CG-6', 'CG-6_TsoftFormat.DAT'
    )
    with open(data_file, 'r') as fh:
        data = read_cg6tsoft(fh)
    assert len(data.raw_grav) == 5515

    data_file = os.path.join(os.getcwd(), 'tests', 'test_csv_data.csv')
    with pytest.raises(ValueError):
        with open(data_file, 'r') as fh:
            data = read_cg6tsoft(fh)
