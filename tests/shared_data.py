import pathlib

TEST_BIN_DIR = pathlib.Path(__file__).parent / "bin"
TEST_DATA_DIR = pathlib.Path(__file__).parent / "data"
TEST_OUT_DIR = pathlib.Path(__file__).parent / "out"

BIGWIG_TO_BEDGRAPH_BIN = TEST_BIN_DIR / "bigWigToBedGraph.linux.x86_64"
BIGWIG_MERGE_PLUS_BIN = TEST_BIN_DIR / "bigWigMergePlus"

INPUT_FILES = (
    TEST_DATA_DIR / "25574.Blueprint.ERS487305.RNA-Seq.signal_reverse.bigWig",
    TEST_DATA_DIR / "25575.Blueprint.ERS487305.RNA-Seq.signal_forward.bigWig",
    TEST_DATA_DIR / "25582.Blueprint.ERS487306.RNA-Seq.signal_reverse.bigWig",
)
INPUT_FILES_STR = [*map(str, INPUT_FILES)]

TEST_POS = "chr1:0-500000"
