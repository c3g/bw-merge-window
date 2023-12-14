# Test against the software this replaces: bigWigMergePlus
import subprocess

from bw_merge_window.entry import entry

from . import shared_data as sd


def parse_bedgraph_line(line: str) -> tuple[str, int, int, float]:
    p = line.strip().split("\t")
    return p[0], int(p[1]), int(p[2]), float(p[3])


def test_against_bwmp():
    out_new = str(sd.TEST_OUT_DIR / "merged-new.bw")
    out_new_bedgraph = str(sd.TEST_OUT_DIR / "merged-new.bedGraph")
    entry(
        (
            sd.TEST_POS,
            *sd.INPUT_FILES_STR,
            "--treat-missing-as-zero",
            "--range",
            "0-1000",
            "--output",
            out_new,
        )
    )
    subprocess.run((str(sd.BIGWIG_TO_BEDGRAPH_BIN), out_new, out_new_bedgraph))

    out_old = str(sd.TEST_OUT_DIR / "merged-old.bw")
    out_old_bedgraph = str(sd.TEST_OUT_DIR / "merged-old.bedGraph")
    subprocess.run(
        (
            str(sd.BIGWIG_MERGE_PLUS_BIN),
            "-range=0-1000",
            "-compress",
            f"-position={sd.TEST_POS}",
            *sd.INPUT_FILES_STR,
            str(sd.TEST_OUT_DIR / "merged-old.bw"),
        )
    )
    subprocess.run((str(sd.BIGWIG_TO_BEDGRAPH_BIN), out_old, out_old_bedgraph))

    with open(out_old_bedgraph, "r") as fho:
        fho_lines = fho.readlines()
    with open(out_new_bedgraph, "r") as fhn:
        fhn_lines = fhn.readlines()

    for fho_line, fhn_line in zip(map(parse_bedgraph_line, fho_lines), map(parse_bedgraph_line, fhn_lines)):
        assert fho_line[:3] == fhn_line[:3]
        assert abs(fhn_line[3] - fho_line[3]) < 0.002  # difference is smaller than some epsilon
