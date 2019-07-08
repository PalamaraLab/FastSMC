import os


def pytest_addoption(parser):

    exe_locations = [
        os.path.join('build', 'ASMC_exe'),  # azure pipelines default
        os.path.join('Release', 'ASMC_exe'),
        os.path.join('Debug', 'ASMC_exe'),
    ]

    exe = None

    for location in exe_locations:
        if os.path.isfile(location):
            exe = location
            break

    parser.addoption("--asmc_exe", action="store", default=exe)


def pytest_generate_tests(metafunc):
    # This is called for every test. Only get/set command line arguments
    # if the argument is specified in the list of test "fixturenames".
    option_value = metafunc.config.option.asmc_exe
    if 'asmc_exe' in metafunc.fixturenames and option_value is not None:

        metafunc.parametrize("asmc_exe", [option_value])
