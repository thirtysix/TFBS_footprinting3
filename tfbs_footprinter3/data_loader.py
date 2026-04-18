"""Species-specific experimental data download and loading.

Responsible for two concerns:

1. Base-data and per-species download from the public S3 bucket
   (`tfbssexperimentaldata`) at first use or whenever `-update` is
   passed. See `experimentalDataUpdater` and `experimentaldata`.

2. Loading the set of JSON / msgpack / TSV files that `find_clusters`
   consumes during scoring — PWM thresholds, GERP conservation, CAGE
   peaks and correlations, GTEx variants, GTRD metaclusters, ATAC-Seq
   peaks, CpG weights, and the CAS p-values lookup. See
   `species_specific_data`.

Extracted from tfbs_footprinter3.py. The data directory lives at
`tfbs_footprinter3/data/` (sibling of this module); `script_dir` below
computes it via `os.path.dirname(__file__)` which resolves to the
package directory whether the caller is this module or the monolith.
"""
from __future__ import annotations

import logging
import os
import tarfile
import time

import numpy as np
import pandas as pd
import wget

from tfbs_footprinter3.io_utils import (
    directory_creator,
    dump_json,
    load_json,
    load_msgpack,
)

script_dir = os.path.dirname(__file__)
AWS_SERVER = "https://s3.us-east-2.amazonaws.com"


def experimentalDataUpdater(exp_data_update):
    """
    Update the experimental data by downloading it from the Amazon repository.
    Only activates if the user specifically calls for an update, or the data directory does not exist.
    """

    experimental_data_dir = os.path.join(script_dir, 'data')
    experimental_data_present = False

    # test if data dir exists
    if not os.path.exists(experimental_data_dir):
        directory_creator(experimental_data_dir)
        exp_data_update = True
        print("Data dir doesn't exist")

    # if the data dir exists, check to see that all of the required file patterns are present
    else:
        required_data_file_patterns = ["JASPAR_2026_pwms.json"]
        experimental_data_filenames = [x for x in os.listdir(experimental_data_dir) if os.path.isfile(os.path.join(experimental_data_dir, x))]
        all_patterns_matched = all([any([required_data_file_pattern in experimental_data_filename
                                          for experimental_data_filename in experimental_data_filenames])
                                         for required_data_file_pattern in required_data_file_patterns])

        if all_patterns_matched:
            experimental_data_present = True
        else:
            exp_data_update = True

    # perform an update of the base data dir
    if exp_data_update:
        experimental_data_url = "/".join([AWS_SERVER, "tfbssexperimentaldata", "data.tar.gz"])
        experimental_data_down_loc = os.path.join(script_dir, 'data.tar.gz')
        print("Downloading the most current experimental data")
        logging.info(" ".join(["Downloading most current experimental data."]))

        try:
            wget.download(experimental_data_url, out=experimental_data_down_loc)
            tar = tarfile.open(experimental_data_down_loc)
            tar.extractall(experimental_data_dir)
            experimental_data_present = True

        except:
            logging.warning(" ".join(["Error in downloading experimental data.  Check your internet connection."]))
            experimental_data_present = False

    return experimental_data_present


def experimentaldata(target_species):
    """
    Retrieve the experimental data for non-human species by downloading it from the Amazon repository.
    Activates if the current installation Data directory does not include data for the target species.
    Example location: https://s3.us-east-2.amazonaws.com/tfbssexperimentaldata/acanthochromis_polyacanthus.tar.gz
    """

    experimental_data_dir = os.path.join(script_dir, 'data')
    experimental_data_species_dir = os.path.join(experimental_data_dir, target_species)
    experimental_data_species_dir_tar = os.path.join(experimental_data_dir, ".".join([target_species, "tar.gz"]))

    if not os.path.exists(experimental_data_species_dir):
        experimental_data_species_url = "/".join([AWS_SERVER, "tfbssexperimentaldata", ".".join([target_species, "tar.gz"])])
        print(experimental_data_species_url)
        logging.info("Downloading most current experimental data for %s." % target_species)

        try:
            wget.download(experimental_data_species_url, out=experimental_data_species_dir_tar)
            tar = tarfile.open(experimental_data_species_dir_tar)
            tar.extractall(experimental_data_dir)

        except:
            logging.warning(" ".join(["Error in downloading experimental data.  Check your internet connection."]))


def experimentalDataUpdater_beta():
    """
    Update the experimental data by downloading it from the Amazon repository.
    Using a file which contains an dictionary of the most up to date exp. data filenames,
    activates if data directory does not exist, if the data directory does not contain the most recent files,
    or if it has been >= 60 days since last update.
    This version of the updater is perhaps too error prone for this stage of development.
    """

    ## download experimental data if not already present or if it is outdated
    current_version_url = "https://s3.us-east-2.amazonaws.com/tfbssexperimentaldata/experimental_data.current_versions.json"
    experimental_data_url = "https://s3.us-east-2.amazonaws.com/tfbssexperimentaldata/data.tar.gz"
    experimental_data_down_loc = os.path.join(script_dir, 'data.tar.gz')
    experimental_data_dir = os.path.join(script_dir, 'data')
    current_versions_file = os.path.join(experimental_data_dir, "experimental_data.current_versions.json")
    update_required = False

    if not os.path.exists(experimental_data_dir):
        directory_creator(experimental_data_dir)
        update_required = True
    else:
        if os.path.exists(current_versions_file):
            # check if all current versions are in the data dir
            current_versions = load_json(current_versions_file)
            current_versions_filenames = current_versions.values()
            owned_versions_filenames = os.listdir(experimental_data_dir)
            missing_files = [x for x in current_versions_filenames if x not in owned_versions_filenames]
            if len(missing_files) > 0:
                update_required = True

            # check if 60 days has passed since last check of versions
            if 'last_checked' in current_versions:
                current_versions_last_checked = current_versions['last_checked']
                if (time.time() - current_versions_last_checked) / (3600 * 24) >= 60:
                    update_required = True
        else:
            update_required = True

    # download the most current experimental data
    if update_required:
        print("Downloading the most current experimental data")
        logging.info(" ".join(["Downloading most current experimental data."]))

        try:
            wget.download(current_version_url, out=experimental_data_dir)
            wget.download(experimental_data_url, out=experimental_data_down_loc)
            tar = tarfile.open(experimental_data_down_loc)
            tar.extractall(experimental_data_dir)

        except:
            logging.warning(" ".join(["Error in downloading experimental data.  Check your internet connection."]))

    # update the current versions file last checked time with current time
    if os.path.exists(current_versions_file):
        current_versions = load_json(current_versions_file)
        current_versions['last_checked'] = time.time()
        dump_json(current_versions_file, current_versions)


def species_specific_data(target_species, chromosome, species_specific_data_dir):
    """
    Many datasets are species-specific.  If the current target species has species-specific datasets, load them.
    Altered to allow for multiple versions of data.  The matching files are sorted and those occuring later in the list, are presumably later editions.
    """

    logging.info(" ".join(["Species-specific data: loading"]))

    # load - PWM thresholds
    species_pwm_score_threshold_df = None
    species_pwm_score_threshold_dict_filenames = [os.path.join(species_specific_data_dir, x) for x in os.listdir(species_specific_data_dir) if ".tfs_thresholds." in x and target_species in x]
    if len(species_pwm_score_threshold_dict_filenames) > 0:
        species_pwm_score_threshold_dict_filenames.sort()
        species_pwm_score_threshold_dict_filename = species_pwm_score_threshold_dict_filenames[-1]
        species_pwm_score_threshold_df = pd.read_csv(species_pwm_score_threshold_dict_filename, sep="\t")

    # load GERP locations
    gerp_conservation_locations_dict = {}
    species_group = ""
    gerp_data_dir = os.path.join(species_specific_data_dir, "gerp_data")
    gerp_conservation_locations_dict_filenames = [os.path.join(gerp_data_dir, x) for x in os.listdir(gerp_data_dir) if ".gerp_loc_weight." in x and target_species in x]
    if len(gerp_conservation_locations_dict_filenames) > 0:
        gerp_conservation_locations_dict_filenames.sort()
        gerp_conservation_locations_dict_filename = gerp_conservation_locations_dict_filenames[-1]
        gerp_conservation_locations_dict = load_msgpack(gerp_conservation_locations_dict_filename)

    # load - CAGE data; human or non-human
    cage_dict = {}
    cage_data_dir = os.path.join(species_specific_data_dir, "cage_data")
    if os.path.exists(cage_data_dir):
        cage_dict_filenames = [os.path.join(cage_data_dir, x) for x in os.listdir(cage_data_dir) if ".CAGE.peak_dict." in x and target_species in x]
        if len(cage_dict_filenames) > 0:
            cage_dict_filename = sorted(cage_dict_filenames)[-1]
        if os.path.exists(cage_dict_filename):
            cage_dict = load_json(cage_dict_filename)

    # load CAGE locs occuring near promoters of TFs
    TF_cage_dict = {}
    cage_data_dir = os.path.join(species_specific_data_dir, "cage_data")
    if os.path.exists(cage_data_dir):
        TF_cage_dict_filename = os.path.join(cage_data_dir, ".".join(["homo_sapiens", "CAGE", "jasparTFs", "dict", "json"]))
        if os.path.exists(TF_cage_dict_filename):
            TF_cage_dict = load_json(TF_cage_dict_filename)

    # load CAGE dist weights
    cage_dist_weights_dict = {}
    if os.path.exists(cage_data_dir):
        cage_dist_weights_dict_filenames = [os.path.join(cage_data_dir, x) for x in os.listdir(cage_data_dir) if "cage_dist_weights" in x and target_species in x]
        if len(cage_dist_weights_dict_filenames) > 0:
            cage_dist_weights_dict_filenames.sort()
            cage_dist_weights_dict_filename = cage_dist_weights_dict_filenames[-1]
            cage_dist_weights_dict = load_json(cage_dist_weights_dict_filename)

    # load CAGE correlations
    cage_correlations_dict = {}
    cage_corr_data_dir = os.path.join(species_specific_data_dir, "cage_corr_data")
    if os.path.exists(cage_corr_data_dir):
        cage_correlations_dict_filename = os.path.join(cage_corr_data_dir, ".".join([target_species, "CAGE_corr", "Chr" + chromosome.upper(), "hg38", "msg"]))
        if os.path.exists(cage_correlations_dict_filename):
            cage_correlations_dict = load_msgpack(cage_correlations_dict_filename)
        else:
            print("cage_correlations_dict not loaded")

    # load CAGE correlation weights
    cage_corr_weights_dict = {}
    if os.path.exists(cage_corr_data_dir):
        cage_corr_weights_dict_filename = os.path.join(cage_corr_data_dir, ".".join([target_species, "CAGE_corr", "weight_dict", "hg38", "json"]))
        if os.path.exists(cage_corr_weights_dict_filename):
            cage_corr_weights_dict = load_json(cage_corr_weights_dict_filename)
            cage_corr_weights_dict = {float(k): v for k, v in cage_corr_weights_dict.items()}
        else:
            print("cage_corr_weights_dict not loaded")

    # load CpG score weights
    cpg_data_dir = os.path.join(species_specific_data_dir, "cpg_data")
    cpg_obsexp_weights_dict_filenames = [os.path.join(cpg_data_dir, x) for x in os.listdir(cpg_data_dir) if ".cpg_obsexp_weights." in x and target_species in x]
    if len(cpg_obsexp_weights_dict_filenames) > 0:
        cpg_obsexp_weights_dict_filenames.sort()
        cpg_obsexp_weights_dict_filename = cpg_obsexp_weights_dict_filenames[-1]
        cpg_obsexp_weights_dict = load_json(cpg_obsexp_weights_dict_filename)
        cpg_obsexp_weights_dict = {float(k): float(v) for k, v in cpg_obsexp_weights_dict.items()}
        cpg_obsexp_weights_dict_keys = list(cpg_obsexp_weights_dict.keys())
        cpg_obsexp_weights_dict_keys.sort()
    else:
        cpg_obsexp_weights_dict = {}
        cpg_obsexp_weights_dict_keys = []

    # load GTEx variants
    gtex_variants = {}
    gtex_data_dir = os.path.join(species_specific_data_dir, "gtex_data")
    if os.path.exists(gtex_data_dir):
        gtex_chrom_dict_filename = os.path.join(gtex_data_dir, ".".join([target_species, "gtex_v7", "Chr" + chromosome.upper(), "min_unique", "eqtls", "grch38", "msg"]))
        if os.path.exists(gtex_chrom_dict_filename):
            gtex_variants = load_msgpack(gtex_chrom_dict_filename)

    # load GTEx weights
    gtex_weights_dict = {}
    if os.path.exists(gtex_data_dir):
        gtex_weights_dict_filenames = [os.path.join(gtex_data_dir, x) for x in os.listdir(gtex_data_dir) if "gtex_weights" in x and target_species in x]
        if len(gtex_weights_dict_filenames) > 0:
            gtex_weights_dict_filenames.sort()
            gtex_weights_dict_filename = gtex_weights_dict_filenames[-1]
            gtex_weights_dict = load_json(gtex_weights_dict_filename)
            gtex_weights_dict = {float(k): float(v) for k, v in gtex_weights_dict.items()}

    # load meta clusters from GTRD project
    gtrd_metaclusters_dict = {}
    gtrd_data_dir = os.path.join(species_specific_data_dir, "gtrd_data")
    if os.path.exists(gtrd_data_dir):
        gtrd_metaclusters_chrom_dict_filename = os.path.join(gtrd_data_dir, ".".join([target_species, "metaclusters", "interval", "Chr" + chromosome.upper(), "clipped", "ordered", "tupled", "msg"]))
        if os.path.exists(gtrd_metaclusters_chrom_dict_filename):
            gtrd_metaclusters_dict = load_msgpack(gtrd_metaclusters_chrom_dict_filename)

    # load metacluster overlap weights
    metacluster_overlap_weights_dict = {}
    if os.path.exists(gtrd_data_dir):
        metacluster_overlap_weights_dict_filenames = [os.path.join(gtrd_data_dir, x) for x in os.listdir(gtrd_data_dir) if "metaclusters_overlap_weights_dict" in x and target_species in x]
        if len(metacluster_overlap_weights_dict_filenames) > 0:
            metacluster_overlap_weights_dict_filenames.sort()
            metacluster_overlap_weights_dict_filename = metacluster_overlap_weights_dict_filenames[-1]
            metacluster_overlap_weights_dict = load_json(metacluster_overlap_weights_dict_filename)

    # load ATAC-Seq from ChIP-Atlas
    atac_seq_dict = {}
    atac_seq_data_dir = os.path.join(species_specific_data_dir, "atac_data")
    if os.path.exists(atac_seq_data_dir):
        atac_seq_chrom_dict_filename = os.path.join(atac_seq_data_dir, ".".join([target_species, "ATAC_loc_weight", "Chr" + chromosome.upper(), "msg"]))

        if os.path.exists(atac_seq_chrom_dict_filename):
            atac_seq_dict = load_msgpack(atac_seq_chrom_dict_filename)

    # load pre-calculated CAS threshold table: {species}.CAS_thresholds.<jaspar>.tsv.gz
    # Produced by hpc/puhti/build_cas_distributions.py on the campaign output.
    # Format: header (tf_name, p_value, score) + one row per TF per target p-value
    # grid, rows sorted p_value descending / score ascending within each TF.
    # We ingest into {tf_name: {"scores": ndarray (asc), "pvalues": ndarray (desc)}}
    # for searchsorted lookup in scoring.calcCombinedAffinityPvalue.
    cas_pvalues_dict = {}
    cas_thresholds_candidates = sorted(
        os.path.join(species_specific_data_dir, x)
        for x in os.listdir(species_specific_data_dir)
        if ".CAS_thresholds." in x and x.endswith(".tsv.gz")
    )
    if cas_thresholds_candidates:
        # Pick the newest-jaspar version (lexicographic sort puts jaspar_2026
        # after jaspar_2018 / jaspar_2024).
        cas_thresholds_file = cas_thresholds_candidates[-1]
        df = pd.read_csv(cas_thresholds_file, sep="\t")
        for tf_name, grp in df.groupby("tf_name", sort=False):
            sorted_grp = grp.sort_values("score")
            cas_pvalues_dict[tf_name] = {
                "scores": sorted_grp["score"].to_numpy(dtype=np.float64),
                "pvalues": sorted_grp["p_value"].to_numpy(dtype=np.float64),
            }

    return species_pwm_score_threshold_df, gerp_conservation_locations_dict, species_group, cage_dict, TF_cage_dict, cage_dist_weights_dict, cage_correlations_dict, cage_corr_weights_dict, metacluster_overlap_weights_dict, cpg_obsexp_weights_dict, cpg_obsexp_weights_dict_keys, gtex_variants, gtex_weights_dict, gtrd_metaclusters_dict, atac_seq_dict, cas_pvalues_dict
