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
bw-merge-window chr1:100000-200000 file1.bw file2.bw --output average.bw
```
