import numpy as np

import gzip
import os
import sys

"""
Compare two files to check that they are the "same".

For lines containing only numbers, we compare with a delta using numpy.allclose.
For other lines, we compare the strings directly.

The files may be plaintext or gzipped, and this script handles both.
"""


def validate_argv():
    """
    Expects to be called with two arguments: the (relative) paths of the two files the compare
    """
    if len(sys.argv) != 3:
        print('Usage: python test_file_contents.py <file1> <file2>')
        sys.exit(1)

    if not os.path.isfile(sys.argv[1]):
        print(f'Expected a file, but got {sys.argv[1]}')
        sys.exit(1)

    if not os.path.isfile(sys.argv[2]):
        print(f'Expected a file, but got {sys.argv[2]}')
        sys.exit(1)


def is_float(value: str):
    """
    Check if a string can be converted to a float
    :param value: the string to try converting
    :return: whether the string could be converted to a float
    """
    try:
        float(value)
        return True
    except ValueError:
        return False


def line_is_numeric(line: str):
    """
    Check if a line consists entirely of numbers that can be converted to float
    :param line: a string to check
    :return: whether the line consists entirely of numbers
    """
    delim = ' ' if ' ' in line else '\t'

    for x in line.split(delim):
        if not is_float(x):
            return False
    return True


def two_lines_identical(line1, line2):
    """
    Check if two lines are identical.  If they contain only numbers, check with numpy.allclose, else simply compare
    the two strings.
    :param line1: the first line
    :param line2: the second line
    :return: whether the two lines are the same
    """
    if line_is_numeric(line1):
        l1_array = np.fromstring(line1, dtype='float', sep=' ')
        l2_array = np.fromstring(line2, dtype='float', sep=' ')
        return np.allclose(l1_array, l2_array)
    else:
        return line1 == line2


def extract_lines_from_file(file):
    """
    Get the lines from a file that may be gzipped.  Each line should become a stripped string.
    :param file: the file to open
    :return: a list of strings representing the lines in the file
    """
    lines_in_file = []
    if file.strip().endswith('.gz'):
        with gzip.open(file, 'rb') as f:
            for line in f.readlines():
                lines_in_file.append(line.decode('utf-8').strip())
    else:
        with open(file, 'r') as f:
            for line in f.readlines():
                lines_in_file.append(line.strip())

    return lines_in_file


def main():
    """
    Check whether the files contain the same contents (up to small floating point errors).
    """

    # Check files are valid
    validate_argv()
    file_1 = sys.argv[1]
    file_2 = sys.argv[2]

    lines_1 = extract_lines_from_file(file_1)
    lines_2 = extract_lines_from_file(file_2)

    if len(lines_1) != len(lines_2):
        print(f'Files {file_1} and {file_2} do not have the same number of lines'
              f' ({len(lines_1)} and {len(lines_2)} respectively)')
        sys.exit(1)

    for idx, [l1, l2] in enumerate(zip(lines_1, lines_2)):
        if not two_lines_identical(l1, l2):
            print(f'Line {1 + idx} of {file_1} and {file_2} do not match:\n{l1}{l2}')
            sys.exit(1)

    print(f'Files {file_1} and {file_2} match')
    sys.exit(0)


if __name__ == '__main__':
    main()
