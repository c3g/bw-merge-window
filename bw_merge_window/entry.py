import argparse
import asyncio
from pathlib import Path

from . import __version__
from .logger import get_logger
from .merge import merge_bigwigs


async def main(args: tuple[str, ...] | None):
    parser = argparse.ArgumentParser(
        prog="bw-merge-window",
        description="A command-line utility for merging a window of a set of bigWigs.",
    )

    parser.add_argument('--version', action="version", version=__version__)

    parser.add_argument("window", type=str, help="Window to merge: contig:start-end")
    parser.add_argument("bigWig", type=Path, nargs="+", help="One or more bigWig files to merge.")
    parser.add_argument("--output", "-o", type=Path, required=True, help="Path to write merged bigWig file to.")
    parser.add_argument(
        "--range", "-r",
        type=str,
        help="Range for output values, applied to each bigWig individually. If not set, raw values will be used.")
    parser.add_argument("--treat-missing-as-zero", "-t", action="store_true")

    logger = get_logger()

    args = parser.parse_args(args)
    await merge_bigwigs(args.window, tuple(args.bigWig), args.output, args.range, args.treat_missing_as_zero, logger)


def entry(args: tuple[str, ...] | None = None):
    asyncio.run(main(args))


if __name__ == "__main__":  # pragma: no cover
    entry()  # run without args specified (will use sys.argv)
