import functools
import logging
import math
import numpy as np
import pyBigWig
import re

from numpy.typing import NDArray
from pathlib import Path

__all__ = ["merge_bigwigs"]


WINDOW_FORMAT = re.compile(r"^(\w+)(?::(-?\d+)-(\d+)?)?$")
RANGE_FORMAT = re.compile(r"^(-?\d+)-(-?\d+)$")


def int_cast_or_raise_with_msg(val: str, field: str) -> int:
    try:
        return int(val)
    except ValueError:
        raise ValueError(f"{field} must be an integer")


def get_window_contig_start_end(window: str, logger: logging.Logger) -> tuple[str, int, int | None]:
    if window_match := WINDOW_FORMAT.match(window):
        contig = window_match.group(1)

        start_group = window_match.group(2)
        start: int = int_cast_or_raise_with_msg(window_match.group(2), "start") if start_group is not None else 0
        if start < 0:
            logger.warning("start should not be below 0; clamping to 0")
            start = 0

        end_group = window_match.group(3)
        end = int_cast_or_raise_with_msg(end_group, "end") if end_group is not None else None

        return contig, start, end

    raise ValueError("Malformatted window: should be contig, contig:start-, or contig:start-end")


def get_range_values(output_range: str | None) -> tuple[int, int] | None:
    if range_match := (RANGE_FORMAT.match(output_range) if output_range else None):
        range_min = int_cast_or_raise_with_msg(range_match.group(1), "range min")
        range_max = int_cast_or_raise_with_msg(range_match.group(2), "range max")

        if range_min >= range_max:
            raise ValueError("Malformatted output range: min must be less than max")

        return range_min, range_max

    else:
        if output_range:  # could not match range from input value
            raise ValueError("Malformatted output range: should be min-max")
        return None


def get_bigwig_file_handles(files: tuple[Path, ...]) -> tuple[pyBigWig.pyBigWig, ...]:
    return tuple(pyBigWig.open(str(f)) for f in files)


def verify_contig_presence_and_get_length(h: pyBigWig.pyBigWig, h_idx: int, contig: str) -> int:
    # noinspection PyArgumentList
    if (cl := h.chroms(contig)) is not None:
        return int(cl)
    else:
        raise ValueError(f"Contig {contig} not found in file at index {h_idx}")


def get_consensus_contig_length_or_raise(handles: tuple[pyBigWig.pyBigWig, ...], contig: str):
    if len(handles) == 0:
        raise ValueError("At least one bigWig file must be specified")

    contig_length: int = verify_contig_presence_and_get_length(handles[0], 0, contig)
    for i, h in enumerate(handles[1:], 1):
        if (cl := verify_contig_presence_and_get_length(h, i, contig)) != contig_length:
            raise ValueError(
                f"Encountered two reference lengths for contig {contig}: {contig_length} and {cl}. Are "
                f"these bigWigs from different assemblies?"
            )

    return contig_length


def get_values_for_file(
    fh: pyBigWig.pyBigWig,
    contig: str,
    start: int,
    end: int,
    range_vals: tuple[int, int] | None,
) -> NDArray:
    bw_header = fh.header()
    vals = fh.values(contig, start, end, numpy=True)
    if range_vals:
        # If the user has specified an output range, re-value the averaged array to be in this range
        vals -= bw_header["minVal"]
        vals /= bw_header["maxVal"] - bw_header["minVal"]
        # now, vals is between 0 and 1 - get it back to be between range min/max
        vals *= range_vals[1] - range_vals[0]
        vals += range_vals[0]
    return vals


async def merge_bigwigs(
    window: str,
    files: tuple[Path, ...],
    output_path: Path,
    output_range: str | None,
    treat_missing_as_zero: bool,
    logger: logging.Logger,
):
    # Process input values ----------------------------------------------------------------
    contig, start, end = get_window_contig_start_end(window, logger)
    range_vals: tuple[int, int] | None = get_range_values(output_range)
    bw_handles: tuple[pyBigWig.pyBigWig, ...] = get_bigwig_file_handles(files)

    # noinspection PyArgumentList
    merged_bw_h: pyBigWig.pyBigWig = pyBigWig.open(str(output_path), "w")
    # -------------------------------------------------------------------------------------

    try:
        # First, get the contig length and ensure that the contig length is consistent across every bigWig file
        ref_contig_length: int = get_consensus_contig_length_or_raise(bw_handles, contig)

        if end is None:
            end = ref_contig_length
        elif end > ref_contig_length:
            logger.warning(f"end should not be above contig length, clamping to {ref_contig_length}")
            end = ref_contig_length

        files_values_matrix = np.array(
            tuple(
                # For each bigWig, compute and collect the (optionally range-ified) base-level values
                get_values_for_file(h, contig, start, end, range_vals)
                for h in bw_handles
            )
        )

        if treat_missing_as_zero:
            np.nan_to_num(files_values_matrix, copy=False, nan=0.0)

        # Calculate averages for each position. If a NaN is present at any point, it makes the whole average NaN:
        files_values_matrix_avg = np.mean(files_values_matrix, axis=0).tolist()

        # Write the output bigWig

        merged_bw_h.addHeader([(contig, ref_contig_length)])

        entry_starts: list[int] = []
        entry_ends: list[int] = []
        entry_values: list[float] = []

        current_start: int = start
        current_value: float = files_values_matrix_avg[0]

        @functools.cache
        def _is_a_gap(f: float) -> bool:
            return (treat_missing_as_zero and f == 0.0) or math.isnan(f)

        # start the enumeration at 1, since we already handled the first entry.
        # add a NaN on the end of the values list to force a final write.
        for i, v in enumerate((*files_values_matrix_avg[1:], math.nan), start + 1):
            if (_is_a_gap(v) or (not _is_a_gap(v) and v != current_value)) and not _is_a_gap(current_value):
                # transition from non-NaN block to NaN block
                #  - add entry values
                entry_starts.append(current_start)
                entry_ends.append(i)
                entry_values.append(current_value)
                #  - switch over current values
                current_start = i
                current_value = v
            elif _is_a_gap(current_value) and not _is_a_gap(v):
                # transition from NaN block to non-NaN block - don't add any entries, since our values are undefined for
                # the region we just passed through.
                #  - switch over current value
                current_start = i
                current_value = v

        merged_bw_h.addEntries(
            [contig] * len(entry_starts),
            entry_starts,
            ends=entry_ends,
            values=entry_values,
        )

    finally:
        for h in bw_handles:
            h.close()
        merged_bw_h.close()
