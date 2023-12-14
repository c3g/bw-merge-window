import logging
import pyBigWig
import pytest

from bw_merge_window import merge as m

from . import shared_data as sd

logger = logging.getLogger(__name__)


def test_int_cast():
    with pytest.raises(ValueError) as e:
        m.int_cast_or_raise_with_msg("x", "x")
    assert "x must be an integer" in str(e)


def test_parse_window():
    assert m.get_window_contig_start_end("chr1:0-10000", logger) == ("chr1", 0, 10000)
    assert m.get_window_contig_start_end("chr1:-5-10000", logger) == ("chr1", 0, 10000)
    assert m.get_window_contig_start_end("chr1:0-", logger) == ("chr1", 0, None)
    assert m.get_window_contig_start_end("chr1", logger) == ("chr1", 0, None)

    with pytest.raises(ValueError):
        m.get_window_contig_start_end("chr1:x-y", logger)

    with pytest.raises(ValueError):
        m.get_window_contig_start_end("", logger)

    with pytest.raises(ValueError):
        m.get_window_contig_start_end("chr1:-", logger)

    with pytest.raises(ValueError):
        m.get_window_contig_start_end("chr1:-10000", logger)


def test_parse_range():
    assert m.get_range_values(None) is None
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


def test_open_handles():
    assert m.get_bigwig_file_handles(()) == ()

    fhs = m.get_bigwig_file_handles(sd.INPUT_FILES)  # TODO: fixture
    try:
        assert len(fhs) == len(sd.INPUT_FILES)
    finally:
        for fh in fhs:
            fh.close()


def test_contig_presence():
    h = pyBigWig.open(str(sd.INPUT_FILES[0]))

    try:
        assert m.verify_contig_presence_and_get_length(h, 0, "chr1") == 249250621  # GRCh37 chr1 length

        with pytest.raises(ValueError):
            m.verify_contig_presence_and_get_length(h, 0, "chrDNE")

    finally:
        h.close()


def test_contig_consensus_length():
    with pytest.raises(ValueError) as e:
        m.get_consensus_contig_length_or_raise((), "chr1")
    assert "At least one bigWig file must be specified" in str(e)

    fhs = m.get_bigwig_file_handles(sd.INPUT_FILES)
    try:
        assert m.get_consensus_contig_length_or_raise(fhs, "chr1") == 249250621  # GRCh37 chr1 length
    finally:
        for fh in fhs:
            fh.close()
