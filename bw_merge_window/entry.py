import argparse
import asyncio
from pathlib import Path

from . import __version__
from .merge import merge_bigwigs


async def main():
    parser = argparse.ArgumentParser(
        prog="bw-merge-window",
        description="A command-line utility for merging a window of a set of bigWigs.",
    )

    parser.add_argument('--version', action="version", version=__version__)

    parser.add_argument("window", type=str, help="Window to merge: contig:start-end")
    parser.add_argument("bigWig", type=Path, nargs="+", help="One or more bigWig files to merge.")
    parser.add_argument("--output", "-o", type=Path, required=True, help="Path to write merged bigWig file to.")
    parser.add_argument(
        "--range", "-r", type=str, help="Range for output values; if not set, raw values will be used.")

    args = parser.parse_args()
    await merge_bigwigs(args.window, tuple(args.bigWig), args.output, args.range)


def entry():
    asyncio.run(main())


if __name__ == "__main__":  # pragma: no cover
    entry()
