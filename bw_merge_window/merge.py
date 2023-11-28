import math
import numpy as np
import pyBigWig
import re
from pathlib import Path

__all__ = ["merge_bigwigs"]


WINDOW_FORMAT = re.compile(r"^(\w+):(\d+)-(\d+)$")
RANGE_FORMAT = re.compile(r"^(\d+)-(\d+)$")


async def merge_bigwigs(window: str, files: tuple[Path, ...], output_path: Path, output_range: str | None):
    window_match = WINDOW_FORMAT.match(window)

    if not window_match:
        raise ValueError("Malformatted window: should be contig:start-end")

    range_match = RANGE_FORMAT.match(output_range) if output_range else None
    range_vals: tuple[int, int] | None = None
    if output_range and not range_match:
        raise ValueError("Malformatted output range: should be min-max")

    try:
        if range_match:
            range_vals = (int(range_match.group(1)), int(range_match.group(2)))
    except ValueError:
        raise ValueError("Output range min/max must be integers")

    contig = window_match.group(1)

    try:
        start = int(window_match.group(2))
    except ValueError:
        raise ValueError("start must be an integer")

    try:
        end = int(window_match.group(3))
    except ValueError:
        raise ValueError("end must be an integer")

    bw_handles: tuple[pyBigWig.pyBigWig, ...] = tuple(pyBigWig.open(str(f)) for f in files)

    # noinspection PyArgumentList
    merged_bw_h: pyBigWig.pyBigWig = pyBigWig.open(str(output_path), "w")

    try:
        files_values = []
        ref_contig_length: int | None = None

        # For each bigWig, collect the base-level values
        for i, h in enumerate(bw_handles):
            # noinspection PyArgumentList
            contig_length = h.chroms(contig)

            if not contig_length:
                raise ValueError(f"Contig {contig} not found in file {files[i]}")

            int_cl = int(contig_length)

            if ref_contig_length is not None and int_cl != ref_contig_length:
                raise ValueError(
                    f"Encountered two reference lengths for contig {contig}: {ref_contig_length} and {int_cl}. Are "
                    f"these bigWigs from different assemblies?")
            elif ref_contig_length is None:
                ref_contig_length = int_cl

            bw_header = h.header()

            vals = h.values(contig, start, end, numpy=True)
            if range_vals:
                # If the user has specified an output range, re-value the averaged array to be in this range
                vals -= bw_header["minVal"]
                vals /= (bw_header["maxVal"] - bw_header["minVal"])
                # now, vals is between 0 and 1 - get it back to be between range min/max
                vals *= (range_vals[1] - range_vals[0])
                vals += range_vals[0]

            files_values.append(vals)

        files_values_matrix = np.array(files_values)

        # Calculate averages for each position. If a NaN is present at any point, it makes the whole average NaN:
        files_values_matrix_avg = np.mean(files_values_matrix, axis=0).tolist()

        # Write the output bigWig

        merged_bw_h.addHeader([(contig, ref_contig_length)])

        entry_starts: list[int] = []
        entry_ends: list[int] = []
        entry_values: list[float] = []

        current_start: int = start
        current_value: float = files_values_matrix_avg[0]

        for i, v in enumerate((*files_values_matrix_avg[1:], math.nan)):  # add a NaN on the end to force final write
            if (math.isnan(v) or (not math.isnan(v) and v != current_value)) and not math.isnan(current_value):
                # transition from non-NaN block to NaN block
                #  - add entry values
                entry_starts.append(current_start)
                entry_ends.append(start + i)
                entry_values.append(current_value)
                #  - switch over current values
                current_start = start + i
                current_value = v

        merged_bw_h.addEntries([contig] * len(entry_starts), entry_starts, ends=entry_ends, values=entry_values)

    finally:
        for h in bw_handles:
            h.close()
        merged_bw_h.close()
