# Before any changes (Release):
# times: [51.93998432 52.30349731 51.96640849 52.12452006 51.67673945]
# median: 51.966408491134644
#
#
#
#
#

import pathlib
import subprocess
import sys
import time

import numpy as np

test_exe = pathlib.Path(__file__).parent / 'Release' / 'ASMC_regression'

if not test_exe.is_file():
    print(f'Expected {test_exe} to exist but it does not...')
    sys.exit(1)

num_repeats = 1
times = np.zeros(num_repeats)

for i, _ in enumerate(times):
    t_start = time.time()
    subprocess.call([test_exe, '[HMM_regression]'])
    times[i] = time.time() - t_start

print(f'times: {times}')
print(f'median: {np.median(times)}')
