# bw-merge-window

[![PyPI version](https://badge.fury.io/py/bw-merge-window.svg)](https://badge.fury.io/py/bw-merge-window)
[![codecov](https://codecov.io/gh/c3g/bw-merge-window/graph/badge.svg?token=HU9UTAASQG)](https://codecov.io/gh/c3g/bw-merge-window)

A command-line utility for merging and averaging a window of a set of bigWigs, inspired by/tested against the previous 
C3G `bigWigMergePlus` tool developed by Romain Grégoire.

This is part of the core set of tools that powers the [EpiVar Browser](https://github.com/c3g/epivar-browser).

**Requires:** Python 3.10+


## License

Everything in this repository **EXCEPT FOR EVERYTHING** under the `tests/bin` and `tests/data` directories is licensed 
under the terms of the [GNU General Public License, v3.0](./LICENSE). 

&copy; McGill University 2023

Everything under the `tests/bin` directory is &copy; UCSC. Everything under the `tests/data` directory is from the
IHEC Epigenomes portal and the copyright is retained by the relevant holders.


## Installation

```bash
pip install bw-merge-window
```


## Usage

```bash
bw-merge-window chr1:100000-200000 file1.bw file2.bw --output average.bw [--range 0-1000] [--treat-missing-as-zero]
```

If the optional flag `--treat-missing-as-zero` is passed, `N/A` values in bigWigs become 0, and any 0-values in the
merged file will not be written. Otherwise, if a value is missing in one or more of the files, it'll be missing in the
merged file as well. This flag is **required** if you want full backwards-compatibility with `bigWigMergePlus`.

If you want to include a negative value for `--range`, use "equals" syntax instead, e.g., `--range='-500-1000'`.
