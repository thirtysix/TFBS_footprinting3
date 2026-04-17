"""Input-file parsing and JASPAR TF ID helpers.

Extracted from tfbs_footprinter3.py (previously lines 351-404, 1464-1499).
Pure parsing and name-cleaning utilities; no Ensembl/network dependencies.
"""
from __future__ import annotations

import csv
import logging
import os


def parse_transcript_ids(transcript_ids_filename):
    """
    If user has provided a file with Ensembl transcript ids, parse these to a list.
    """

    with open(transcript_ids_filename) as transcript_ids_file:
        transcript_ids_list = transcript_ids_file.read().splitlines()
        transcript_ids_list = [x for x in transcript_ids_list if len(x) > 0]

    return transcript_ids_list


def parse_tf_ids(target_tfs_filename):
    """
    If user has provided a file with Ensembl transcript ids, parse these to a list.
    """

    with open(target_tfs_filename) as target_tfs_file:
        target_tfs_list = target_tfs_file.read().splitlines()
        target_tfs_list = [x.upper() for x in target_tfs_list if len(x) > 0]

    return target_tfs_list


def file_to_datalist(data_filename, delimiter):
    """
    Starting with a filename, import and convert data to a list.
    """

    if os.path.exists(data_filename):
        with open(data_filename) as data_file:
            csv_reader = csv.reader(data_file, delimiter=delimiter)
            all_data = list(csv_reader)
    else:
        all_data = []
        print(data_filename, "does not exist.")

    return all_data


def compare_tfs_list_jaspar(target_tfs_list, TFBS_matrix_dict):
    """
    If user has provided a file containing Jaspar TF ids,
    compare candidate entries to those in the loaded dictionary of Jaspar PWMs.
    """

    jaspar_dict_keys = TFBS_matrix_dict.keys()
    erroneous = list(set(target_tfs_list) - set(jaspar_dict_keys))

    target_tfs_list = list(set(jaspar_dict_keys).intersection(target_tfs_list))
    if len(erroneous) > 0:
        logging.warning(" ".join(["the following tf ids are not in the Jaspar database:", ", ".join(erroneous)]))

    return target_tfs_list


def clean_jaspar_names(uncleaned_jaspar_ids):
    """
    Clean names of jaspar transcription factor names.
    MSX3 <- lost in humans.
    RHOX11 <- only present in 3 species.
    DUX <- mouse only gene.
    EWSR1 <- didn't end up in the Ensembl BioMart export.
    MIX-A <- jaspar says present in xenopus laevis, but not even showing
    in Ensembl as a gene for any species.
    """

    special_dict = {"EWSR1-FLI1": ["EWSR1", "FLI1"]}
    names_list = []

    # split the combined names
    for uncleaned_jaspar_id in uncleaned_jaspar_ids:
        uncleaned_jaspar_id = uncleaned_jaspar_id.upper()
        split_names = uncleaned_jaspar_id.split("::")
        for name in split_names:
            names_list.append(name)

    # replace variants
    for i, name in enumerate(names_list):
        names_list[i] = name.replace("(VAR.2)", "").replace("(VAR.3)", "")

    tmp_list = []
    for name in names_list:
        if name in special_dict:
            tmp_list += special_dict[name]
        else:
            tmp_list.append(name)

    names_list = list(set(tmp_list))
    names_list.sort()

    return names_list
