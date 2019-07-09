import gzip
import os
import subprocess


def test_regession():

    script_dir = os.path.realpath(os.path.dirname(__file__))
    base_dir = os.path.realpath(os.path.join(script_dir, '..', '..'))
    bin_dir = os.path.join(base_dir, 'ASMC_SRC', 'BIN')
    old_file = os.path.join(script_dir, 'data', 'regression_test_original.gz')
    assert os.path.isfile(old_file)

    # Find the executable in the BIN dir
    exe_names = ['ASMC_exe', 'ASMC_exe.exe']
    asmc_exe = None

    for exe_name in exe_names:
        exe = os.path.join(bin_dir, exe_name)
        if os.path.isfile(exe):
            asmc_exe = exe
            break

    assert asmc_exe is not None,\
        "ASMC_exe not found in {}. Ensure ASMC_exe target is built before running this test".format(bin_dir)

    # Old file contents are before OxfordRSE involvement in ASMC
    with gzip.open(old_file, 'rt') as gz_f:
        old_lines = gz_f.readlines()

    # New file contents are the result of running the example with the current ASMC source
    decoding_file = os.path.join(base_dir, 'FILES', 'DECODING_QUANTITIES', '30-100-2000.decodingQuantities.gz')
    haps_file = os.path.join(base_dir, 'FILES', 'EXAMPLE', 'exampleFile.n300.array')

    subprocess.call([
        asmc_exe,
        '--decodingQuantFile', decoding_file,
        '--hapsFileRoot', haps_file,
        '--posteriorSums',
    ])

    new_file = os.path.join(base_dir, 'FILES', 'EXAMPLE', 'exampleFile.n300.array.1-1.sumOverPairs.gz')
    assert os.path.isfile(new_file), \
        "No output file found at {}. Did the executable run as expected?".format(new_file)

    with gzip.open(new_file, 'rt') as gz_f:
        new_lines = gz_f.readlines()

    assert len(old_lines) == len(new_lines),\
        "The outputs have different numbers of lines ({} and {})".format(len(old_lines), len(new_lines))

    for i, (old, new) in enumerate(zip(old_lines, new_lines)):
        assert old == new, "The outputs first differ at line {}".format(i)
