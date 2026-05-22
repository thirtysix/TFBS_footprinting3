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
    Test if the hosts TFBS_footprinter3 actually depends on are reachable.

    Probes the S3 experimental-data bucket and the (optionally overridden)
    Ensembl REST endpoint on port 443. Returns True if *any* target is
    reachable, so users behind firewalls that block one host but not the
    other are not false-negatived. The previous implementation probed
    www.google.com, which made the tool unusable for users behind
    Google-blocking firewalls (see GitHub issue: "System does not appear
    to be connected to the internet" on networks where the S3 bucket is
    in fact reachable).
    """

    ensembl_host = os.environ.get(
        "TFBS_FOOTPRINTER3_ENSEMBL_REST",
        "https://oct2024.rest.ensembl.org",
    )
    # Strip scheme/path -> bare hostname for socket.gethostbyname.
    ensembl_hostname = ensembl_host.split("://", 1)[-1].split("/", 1)[0]

    targets = [
        "s3.us-east-2.amazonaws.com",
        ensembl_hostname,
    ]

    for hostname in targets:
        try:
            host = socket.gethostbyname(hostname)
            socket.create_connection((host, 443), 2).close()
            return True
        except Exception:
            continue

    logging.info(
        " ".join([
            "System does not appear to be connected to the internet.",
            "Tried:", ", ".join(targets) + ".",
            "If you are behind a firewall that blocks one of these,",
            "set TFBS_FOOTPRINTER3_ENSEMBL_REST to a reachable mirror.",
        ])
    )
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
