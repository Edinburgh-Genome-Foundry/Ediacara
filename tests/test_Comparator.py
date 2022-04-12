import os
import pytest

import ediacara as edi

data_dir = os.path.join("tests", "data")


def test_comparator():
    edi.Comparator


def test_comparatorgroup():
    edi.ComparatorGroup


def write_sequencinggroup_report():
    edi.write_sequencinggroup_report
