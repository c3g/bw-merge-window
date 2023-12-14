import logging
import pytest

from bw_merge_window import merge as m

logger = logging.getLogger(__name__)


def test_int_cast():
    with pytest.raises(ValueError) as e:
        m.int_cast_or_raise_with_msg("x", "x")
    assert "x must be an integer" in str(e)


def test_parse_window():
    assert m.get_window_contig_start_end("chr1:0-10000", logger) == ("chr1", 0, 10000)
    assert m.get_window_contig_start_end("chr1:-5-10000", logger) == ("chr1", 0, 10000)

    with pytest.raises(ValueError):
        m.get_window_contig_start_end("chr1:x-y", logger)


def test_parse_range():
    assert m.get_range_values("0-1000") == (0, 1000)
    assert m.get_range_values("-5-1000") == (-5, 1000)
    assert m.get_range_values("500-1000") == (500, 1000)
    assert m.get_range_values("-10--5") == (-10, -5)

    with pytest.raises(ValueError):
        m.get_range_values("x-y")

    with pytest.raises(ValueError):
        m.get_range_values("5-y")

    with pytest.raises(ValueError):
        m.get_range_values("10-5")
