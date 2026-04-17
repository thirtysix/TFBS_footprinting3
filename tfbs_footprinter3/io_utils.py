"""Pure utility helpers with no tfbs_footprinter3 internal dependencies.

Extracted from the monolith (previously tfbs_footprinter3.py:223-274, 758-766,
2108-2120). These are imported back into tfbs_footprinter3.py for
backward compatibility with callers that do
``from tfbs_footprinter3.tfbs_footprinter3 import load_json``.
"""
from __future__ import annotations

import json
import logging
import os
import socket
import sys

import msgpack


def signal_handler(signal, frame):
    print('You have manually stopped tfbs_footprinter3 with Ctrl+C')
    sys.exit(0)


def load_json(filename):
    if os.path.exists(filename):
        with open(filename) as open_infile:
            return json.load(open_infile)
    else:
        return None


def dump_json(filename, json_data):
    with open(filename, 'w') as open_outfile:
        json_data = json.dump(json_data, open_outfile)


def load_msgpack(object_filename):
    """unpack a msgpack file to object."""

    if os.path.exists(object_filename):
        with open(object_filename, 'rb') as object_file:
            return msgpack.unpack(object_file, max_array_len=200000, use_list=False, strict_map_key=False)
    else:
        return None


def directory_creator(directory_name):
    """
    Create directory if it does not already exist.
    """

    if not os.path.isdir(directory_name):
        os.mkdir(directory_name)


def is_online():
    """
    Test if the system is online.
    This breaks when TFBS_footprinter3 outlasts Google.
    """

    REMOTE_SERVER = "www.google.com"
    try:
        host = socket.gethostbyname(REMOTE_SERVER)
        s = socket.create_connection((host, 80), 2)
        return True

    except:
        logging.info(" ".join(["System does not appear to be connected to the internet."]))
        return False


def overlap_range(x, y):
    """
    Identify an overlap between two lists of two numbers.
    """

    x.sort()
    y.sort()

    return range(max(x[0], y[0]), min(x[-1], y[-1]) + 1)


def distance_solve(r1, r2):
    # sort the two ranges such that the range with smaller first element
    # is assigned to x and the bigger one is assigned to y
    r1.sort()
    r2.sort()
    x, y = sorted((r1, r2))

    # now if x[1] lies between x[0] and y[0] (x[1] != y[0] but can be equal to x[0])
    # then the ranges are not overlapping and return the difference of y[0] and x[1]
    # otherwise return 0
    if x[1] < y[0]:
        return y[0] - x[1]
    return 0
