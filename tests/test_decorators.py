"""Tests for abutils.utils.decorators."""

import pytest

from abutils.utils.decorators import lazy_property


class _Counter:
    """Helper class that tracks how many times a lazy property is computed."""

    def __init__(self, value):
        self._value = value
        self.compute_count = 0

    @lazy_property
    def expensive(self):
        self.compute_count += 1
        return self._value * 2


class TestLazyProperty:
    def test_computes_on_first_access(self):
        obj = _Counter(5)
        assert obj.expensive == 10
        assert obj.compute_count == 1

    def test_caches_after_first_access(self):
        obj = _Counter(5)
        _ = obj.expensive
        _ = obj.expensive
        _ = obj.expensive
        assert obj.compute_count == 1

    def test_setter_overrides_value(self):
        obj = _Counter(5)
        assert obj.expensive == 10
        obj.expensive = 42
        assert obj.expensive == 42
        assert obj.compute_count == 1  # not recomputed

    def test_deleter_forces_recomputation(self):
        obj = _Counter(5)
        assert obj.expensive == 10
        del obj.expensive
        obj._value = 100
        assert obj.expensive == 200
        assert obj.compute_count == 2

    def test_independent_across_instances(self):
        a = _Counter(3)
        b = _Counter(7)
        assert a.expensive == 6
        assert b.expensive == 14
        assert a.compute_count == 1
        assert b.compute_count == 1

    def test_delete_before_access_is_safe(self):
        obj = _Counter(5)
        del obj.expensive  # should not raise
        assert obj.expensive == 10

    def test_works_with_none_return(self):
        class _NoneReturn:
            @lazy_property
            def value(self):
                return None

        obj = _NoneReturn()
        assert obj.value is None

    def test_preserves_across_multiple_properties(self):
        class _Multi:
            def __init__(self):
                self.a_count = 0
                self.b_count = 0

            @lazy_property
            def prop_a(self):
                self.a_count += 1
                return "a"

            @lazy_property
            def prop_b(self):
                self.b_count += 1
                return "b"

        obj = _Multi()
        assert obj.prop_a == "a"
        assert obj.prop_b == "b"
        # Accessing one doesn't affect the other
        _ = obj.prop_a
        assert obj.a_count == 1
        assert obj.b_count == 1
