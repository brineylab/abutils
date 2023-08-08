#!usr/env/python
# filename: test_utilities.py

from ..utils import utilities


def test_nested_dict_lookup():
    d = {"1a": {"2a": "2a"}, "1b": {"2b": {"3b": "3b"}}}
    a = utilities.nested_dict_lookup(d, ["1a", "2a"])
    b = utilities.nested_dict_lookup(d, ["1b", "2b", "3b"])
    assert a == "2a" and b == "3b"
