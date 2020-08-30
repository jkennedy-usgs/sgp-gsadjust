import os
import gsadjust

import pytest

from test_fixture_pyqt import (
    obstreesurvey,
    obstreesurvey_with_results,
    obstreemodel_with_results,
    obstreemodel,
)

from gsadjust.file.write import export_metadata, export_data, export_summary


def test_export_workspace(obstreemodel_with_results, obstreemodel):
    fn = export_metadata(obstreemodel_with_results, os.getcwd())
    assert os.path.exists(fn)
    assert os.path.getsize(fn) > 200
    os.remove(fn)


def test_export_data(obstreemodel_with_results):
    fn = export_summary(obstreemodel_with_results, os.getcwd())
    assert os.path.exists(fn)
    assert os.path.getsize(fn) > 500
    os.remove(fn)


def test_export_summary(obstreemodel_with_results):
    fn = export_data(obstreemodel_with_results, os.getcwd())
    assert os.path.exists(fn)
    assert os.path.getsize(fn) > 500
    os.remove(fn)
