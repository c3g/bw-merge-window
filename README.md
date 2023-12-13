# bw-merge-window

A command-line utility for merging and averaging a window of a set of bigWigs, inspired by/tested against the previous 
C3G `bigWigMergePlus` tool developed by Romain Grégoire.

This is part of the core set of tools that powers the [EpiVar Browser](https://github.com/c3g/epivar-browser).


## Installation

```bash
pip install bw-merge-window
```


## Usage

```bash
bw-merge-window chr1:100000-200000 file1.bw file2.bw --output average.bw [--treat-missing-as-zero]
```

If the optional flag `--treat-missing-as-zero` is passed, `N/A` values in bigWigs become 0, and any 0-values in the
merged file will not be written. Otherwise, if a value is missing in one or more of the files, it'll be missing in the
merged file as well. This flag is **required** if you want full backwards-compatibility with `bigWigMergePlus`.
