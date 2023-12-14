import pytest

from bw_merge_window.entry import entry

from . import shared_data as sd

ALL_TEST_ARGS = (
    (sd.TEST_POS, *sd.INPUT_FILES_STR),
    (sd.TEST_POS, *sd.INPUT_FILES_STR, "--range", "0-1000"),
    (sd.TEST_POS, *sd.INPUT_FILES_STR, "--treat-missing-as-zero"),
    (sd.TEST_POS, *sd.INPUT_FILES_STR, "--treat-missing-as-zero", "--range", "0-1000"),
    (sd.TEST_POS, *sd.INPUT_FILES_STR, "--treat-missing-as-zero", "--range=-500-1000"),
    (sd.TEST_POS, *sd.INPUT_FILES_STR, "--treat-missing-as-zero", "--range", "10-50"),
    # TODO: these should be tested, but it's too slow...
    # ("chr1", *sd.INPUT_FILES_STR, "--treat-missing-as-zero"),
    # ("chr1:0-", *sd.INPUT_FILES_STR, "--treat-missing-as-zero"),
    # ("chr1:5000-", *sd.INPUT_FILES_STR, "--treat-missing-as-zero"),
)


@pytest.mark.parametrize("args", ALL_TEST_ARGS)
def test_merge_conditions(args):
    out_test = str(sd.TEST_OUT_DIR / "merged-test.bw")
    entry((*args, "--output", out_test))
