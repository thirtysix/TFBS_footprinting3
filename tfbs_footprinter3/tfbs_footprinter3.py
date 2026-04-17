#!/usr/bin/python
# Python vers. 3.8.0 ###########################################################
import argparse
import logging
import os
import signal

# Libraries ####################################################################
import tarfile
import textwrap
import time

import matplotlib
import pandas as pd
import wget
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from tfbs_footprinter3 import __version__

# Public re-exports from extracted modules. These keep
# `from tfbs_footprinter3.tfbs_footprinter3 import <func>` working for
# existing callers (tests, hpc/, downstream users) after the monolith
# split. When we eventually move main() into cli.py, this shim stays.
from tfbs_footprinter3.alignment import (  # noqa: F401
    fasta_writer,
    load_genome_aligned,
    remove_duplicate_species,
    remove_gap_only,
    remove_non_ACGT,
)
from tfbs_footprinter3.data_parsing import (  # noqa: F401
    clean_jaspar_names,
    compare_tfs_list_jaspar,
    file_to_datalist,
    parse_tf_ids,
    parse_transcript_ids,
)
from tfbs_footprinter3.ensembl import ensemblrest, ensemblrest_rate  # noqa: F401
from tfbs_footprinter3.io_utils import (  # noqa: F401
    directory_creator,
    distance_solve,
    dump_json,
    is_online,
    load_json,
    load_msgpack,
    overlap_range,
    signal_handler,
)
from tfbs_footprinter3.output import (  # noqa: F401
    sort_target_species_hits,
    target_species_hits_table_writer,
    top_greatest_hits,
    top_x_greatest_hits,
)
from tfbs_footprinter3.plotting import plot_promoter, plot_promoter_all  # noqa: F401
from tfbs_footprinter3.pwm import PWM_scorer, pwm_maker  # noqa: F401
from tfbs_footprinter3.scoring import (  # noqa: F401
    atac_weights_summing,
    cage_correlations_summing,
    cage_correlations_summing_preparation,
    cage_weights_summing,
    calcCombinedAffinityPvalue,
    cpg_weights_summing,
    eqtl_overlap_likelihood,
    eqtls_weights_summing,
    find_clusters,
    gerp_weights_summing,
    metacluster_weights_summing,
)
from tfbs_footprinter3.translators import (  # noqa: F401
    CpG,
    atac_pos_translate,
    cage_position_translate,
    gerp_positions_translate,
    gtex_position_translate,
    gtrd_positions_translate,
    reg_position_translate,
    start_end_found_motif,
    unaligned2aligned_indexes,
)

matplotlib.use('Agg')
from bisect import bisect_left
from operator import itemgetter

################################################################################
# Description ##################################################################
################################################################################
"""
Thresholds are allowed to be negative, as based on whole genome scoring.

Improvements to be made:
1) Change FANTOM CAGE correlation analysis to only include top CAGE peaks,
instead of all CAGES.  Should reduce size and focus on most relevant CAGEs.
2) Account for GTRD metacluster peak count.  Switch to GTRD peaks by TF if the
number of TFs available becomes larger than those offered by JASPAR.  Similarly,
if all JASPAR TFs become present in GTRD database, then switch from metaclusters
to peaks by TF.
3) Incorporate best combination of data sources by individual TF.
"""


################################################################################
# Functions ####################################################################
################################################################################
script_dir = os.path.dirname(__file__)
curdir = os.getcwd()

################################################################################
# Arguments ####################################################################
################################################################################

def get_args():
    """
    Retrieve arguments provided by the user.
    """

    # Instantiate the parser
    parser = argparse.ArgumentParser(
        prog="tfbs_footprinter3",
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            TFBS Footprinting - Identification of conserved vertebrate transcription factor binding sites (TFBSs).
            See https://github.com/thirtysix/TFBS_footprinting for additional usage instructions.

            ------------------------------------------------------------------------------------------------------
            Example Usage:

                simplest:
                tfbs_footprinter3 PATH_TO/sample_ensembl_ids.txt

                all arguments:
                tfbs_footprinter3 -t PATH_TO/sample_ensembl_ids.txt -tfs PATH_TO/sample_jaspar_tf_ids.txt -pb 900 -pa 100 -tx 10 -p 0.01 -update

                run the sample analysis:
                Option #1: tfbs_footprinter3 -t PATH_TO/sample_analysis/sample_analysis_list.csv
                Option #2: tfbs_footprinter3 -t PATH_TO/sample_analysis/sample_ensembl_ids.txt

                update the experimental data files (not needed often):
                tfbs_footprinter3 -update

                Results will be output to the current directory in a created directory named "tfbs_results"
            ------------------------------------------------------------------------------------------------------
            """))


    # Arguments
    parser.add_argument('--t_ids_file', '-t', metavar='', type=str,
                        help='Required for running an analysis.  Location of a file containing Ensembl target_species transcript ids.  Input options are either a text file of Ensembl transcript ids or a .csv file with individual values set for each parameter.')

    parser.add_argument('--tf_ids_file', '-tfs', metavar='', type=str, default = None, help='Optional: Location of a file containing a limited list of Jaspar TFs to use in scoring alignment \
                                                                                                (see sample file tf_ids.txt at https://github.com/thirtysix/TFBS_footprinting) [default: all Jaspar TFs]')

    parser.add_argument('--promoter_before_tss', '-pb', metavar='', choices = range(-100000, 100001), type=int, default=900,
                        help='(0-100,000) [default: 900] - Number (integer) of nucleotides upstream of TSS to include in analysis.  If this number is negative the start point will be downstream of the TSS, the end point will then need to be further downstream.')

    parser.add_argument('--promoter_after_tss', '-pa', metavar='', choices = range(-100000, 100001), type=int, default=100,
                        help='(0-100,000) [default: 100] - Number (integer) of nucleotides downstream of TSS to include in analysis.  If this number is negative the end point will be upstream of the TSS.  The start point will then need to be further upstream.')

    parser.add_argument('--top_x_tfs', '-tx', metavar='', choices = range(1, 21), type=int, default=10,
                        help='(1-20) [default: 10] - Number (integer) of unique TFs to include in output .svg figure.')

##    parser.add_argument('--output_dir', '-o', metavar='', type=str, default=os.path.join(curdir, "tfbs_results"),
##                        help=" ".join(['[default:', os.path.join(curdir, "tfbs_results"), '] - Full path of directory where result directories will be output.  Make sure that the root directory already exists.']))

    # for now pvalue refers to the PWM score, in the future it will need to relate to the combined affinity score
    parser.add_argument('--pval', '-p', type=float, default=1, help='P-value (float) for PWM score cutoff (range: 1 (all results) to 0.0000001; in divisions of 10 (i.e. 1, 0.1, 0.01, 0.001 etc.) [default: 0.01]')

    parser.add_argument('--pvalc', '-pc', type=float, default=1, help='P-value (float) for PWM score cutoff (range: 1 (all results) to 0.0000001; in divisions of 10 (i.e. 1, 0.1, 0.01, 0.001 etc.) [default: 0.01]')

    parser.add_argument('--exp_data_update', '-update', action="store_true", help='Download the latest experimental data files for use in analysis.  Will run automatically if the "data" directory does not already exist (e.g. first usage).')

    parser.add_argument('--nofig', '-no', action="store_true", help="Don't output a figure.")

    # Functionality to add later
    ##parser.add_argument('--noclean', '-nc', action = 'store_true', help='Optional: Don't clean retrieved alignment. Off by default.')

    # pre-processing the arguments
    args = parser.parse_args()
    args_lists = []
    transcript_ids_filename = args.t_ids_file
    exp_data_update = args.exp_data_update
    nofigure = args.nofig

    if transcript_ids_filename:
        filename, file_extension = os.path.splitext(transcript_ids_filename)

        if file_extension == ".csv" or file_extension == ".tsv":

            if file_extension == ".csv":
                parsed_arg_lines = file_to_datalist(transcript_ids_filename, delimiter = ",")[1:]

            elif file_extension == ".tsv":
                parsed_arg_lines = file_to_datalist(transcript_ids_filename, delimiter = "\t")[1:]

            # If the user has provided a .csv file with the required parameters defined for each Ensembl transcript id
            # this can be parsed to run unique analyses for each.

            for i, parsed_arg_line in enumerate(parsed_arg_lines):
                if len(parsed_arg_line) < 6:
                    print("Incomplete arguments in input file on line", i)
                    print(parsed_arg_line)

                else:
                    transcript_id, target_tfs_filename, promoter_before_tss, promoter_after_tss, top_x_tfs_count, pval, pvalc = parsed_arg_line

                    # promoter_before_tss/promoter_after_tss
                    try:
                        promoter_before_tss = int(promoter_before_tss)
                    except:
                        print("Entered promoter before TSS", promoter_before_tss, "in line", i, "is not an integer.  Defaulting to 900.")
                        promoter_before_tss = 900

                    try:
                        promoter_after_tss = int(promoter_after_tss)
                    except:
                        print("Entered promoter after TSS", promoter_after_tss, "in line", i, "is not an integer.  Defaulting to 100.")
                        promoter_after_tss = 100

                    # top_x_tfs_count
                    try:
                        top_x_tfs_count = int(top_x_tfs_count)
                    except:
                        print("Entered top x tfs count", top_x_tfs_count, "in line", i, "is not an integer.  Defaulting to 10.")
                        top_x_tfs_count = 10

                    # p-value PWM
                    try:
                        pval = float(pval)
                    except:
                        print("Entered p-value threshold", pval, "in line", i, "is not float.  Defaulting to 0.01.")
                        pval = 0.01

                    # p-value combined affinity
                    try:
                        pvalc = float(pvalc)
                    except:
                        print("Entered p-value threshold", pvalc, "in line", i, "is not float.  Defaulting to 0.01.")
                        pvalc = 0.01

                    # update exp data
                    exp_data_update = False

                    parsed_cleaned_arg_line = [transcript_id, target_tfs_filename, promoter_before_tss, promoter_after_tss, top_x_tfs_count, pval, pvalc]
                    args_lists.append([args, transcript_ids_filename] + parsed_cleaned_arg_line)

        else:
            # If the analysis does not require setting the parameters individually for each Ensembl transcript id then build
            # build a list which has all of the parameters set as the same, in this way there can be a similar input format
            # as a .tsv, and standardized handling in the rest of the analysis.
            target_tfs_filename = args.tf_ids_file
            promoter_before_tss = args.promoter_before_tss
            promoter_after_tss = args.promoter_after_tss
            top_x_tfs_count = args.top_x_tfs
            pval = args.pval
            pvalc = args.pvalc
            exp_data_update = args.exp_data_update
            nofigure = args.nofig

            transcript_ids_list = parse_transcript_ids(transcript_ids_filename)
            for transcript_id in transcript_ids_list:
                args_list = [args, transcript_ids_filename, transcript_id, target_tfs_filename, promoter_before_tss, promoter_after_tss, top_x_tfs_count, pval, pvalc]
                args_lists.append(args_list)

    return args_lists, exp_data_update, nofigure


################################################################################
# House-keeping functions ######################################################
################################################################################
# signal_handler, load_json, dump_json, load_msgpack, directory_creator,
# is_online moved to tfbs_footprinter3/io_utils.py (re-exported above).


# ensemblrest and ensemblrest_rate moved to tfbs_footprinter3/ensembl.py
# (re-exported above). The HPC cache patches tff.ensemblrest here in the
# monolith's namespace — see hpc/ensembl_cache.py::patch_tfbs_footprinter3().


# parse_transcript_ids, parse_tf_ids, file_to_datalist, compare_tfs_list_jaspar
# moved to tfbs_footprinter3/data_parsing.py (re-exported above).


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
        #required_data_file_patterns = ["pwms.json", "all_tfs_thresholds", "all_pwms_loglikelihood_dict"]
        required_data_file_patterns = ["pwms.json"]
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
        aws_server = "https://s3.us-east-2.amazonaws.com"
        experimental_data_url = "/".join([aws_server, "tfbssexperimentaldata", "data.tar.gz"])
        experimental_data_down_loc = os.path.join(script_dir,'data.tar.gz')
##        current_versions_file = os.path.join(experimental_data_dir, "experimental_data.current_versions.json")
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
        aws_server = "https://s3.us-east-2.amazonaws.com"
        experimental_data_species_url = "/".join([aws_server, "tfbssexperimentaldata", ".".join([target_species, "tar.gz"])])
        print(experimental_data_species_url)
        logging.info("Downloading most current experimental data for %s." %target_species)

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
    experimental_data_down_loc = os.path.join(script_dir,'data.tar.gz')
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
                if (time.time() - current_versions_last_checked)/(3600*24) >= 60:
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

##        os.remove(experimental_data_down_loc)

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

##    # load GERP locations
##    gerp_conservation_locations_dict = {}
##    species_group = ""
##    gerp_data_dir = os.path.join(species_specific_data_dir, "gerp_data")
##    gerp_conservation_locations_dict_filenames = [os.path.join(gerp_data_dir, x) for x in os.listdir(gerp_data_dir) if "gerp_conservation.locations_dict" in x and target_species in x]
##    if len(gerp_conservation_locations_dict_filenames) > 0:
##        gerp_conservation_locations_dict_filenames.sort()
##        gerp_conservation_locations_dict_filename = gerp_conservation_locations_dict_filenames[-1]
##        gerp_conservation_locations_dict = load_msgpack(gerp_conservation_locations_dict_filename)
##        species_group = gerp_conservation_locations_dict_filename.split(".")[3]
##
##    # load GERP conservation weights
##    gerp_conservation_weight_dict = {}
##    gerp_conservation_weight_dict_filenames = [os.path.join(gerp_data_dir, x) for x in os.listdir(gerp_data_dir) if "gerp_conservation.weight_dict" in x and target_species in x]
##    if len(gerp_conservation_weight_dict_filenames) > 0:
##        gerp_conservation_weight_dict_filenames.sort()
##        gerp_conservation_weight_dict_filename = gerp_conservation_weight_dict_filenames[-1]
##        gerp_conservation_weight_dict = load_msgpack(gerp_conservation_weight_dict_filename)

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
        #cage_dict_filename = os.path.join(cage_data_dir, ".".join([target_species, "CAGE", "peak_dict", "gene", "hg38", "json"]))
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
        cage_correlations_dict_filename = os.path.join(cage_corr_data_dir, ".".join([target_species, "CAGE_corr", "Chr"+chromosome.upper(), "hg38", "msg"]))
        if os.path.exists(cage_correlations_dict_filename):
            cage_correlations_dict = load_msgpack(cage_correlations_dict_filename)
        else:
            print("cage_correlations_dict not loaded")
##        cage_correlations_dict_filenames = [os.path.join(cage_corr_data_dir, x) for x in os.listdir(cage_corr_data_dir) if "rekeyed_combined_cage_corr_dict" in x and target_species in x]
##        if len(cage_correlations_dict_filenames) > 0:
##            cage_correlations_dict_filenames.sort()
##            cage_correlations_dict_filename = cage_correlations_dict_filenames[-1]
##            cage_correlations_dict = load_msgpack(cage_correlations_dict_filename)

    # load CAGE correlation weights
    cage_corr_weights_dict = {}
    if os.path.exists(cage_corr_data_dir):
        cage_corr_weights_dict_filename = os.path.join(cage_corr_data_dir, ".".join([target_species, "CAGE_corr", "weight_dict", "hg38", "json"]))
        if os.path.exists(cage_corr_weights_dict_filename):
            cage_corr_weights_dict = load_json(cage_corr_weights_dict_filename)
            cage_corr_weights_dict = {float(k):v for k,v in cage_corr_weights_dict.items()}
        else:
            print("cage_corr_weights_dict not loaded")

##        cage_corr_weights_dict_filenames = [os.path.join(cage_data_dir, x) for x in os.listdir(cage_data_dir) if "cage_corr_weights" in x and target_species in x]
##        if len(cage_corr_weights_dict_filenames) > 0:
##            cage_corr_weights_dict_filenames.sort()
##            cage_corr_weights_dict_filename = cage_corr_weights_dict_filenames[-1]
##            cage_corr_weights_dict = load_json(cage_corr_weights_dict_filename)
##            cage_corr_weights_dict = {float(k):v for k,v in cage_corr_weights_dict.iteritems()}

##    # load CAGE keys
##    cage_keys_dict = {}
##    if os.path.exists(cage_data_dir):
##        cage_keys_dict_filenames = [os.path.join(cage_data_dir, x) for x in os.listdir(cage_data_dir) if "cage_ids_key_dict" in x and target_species in x]
##        if len(cage_keys_dict_filenames) > 0:
##            cage_keys_dict_filenames.sort()
##            cage_keys_dict_filename = cage_keys_dict_filenames[-1]
##            cage_keys_dict = load_json(cage_keys_dict_filename)


##    # load JASPAR tfs to Ensembl transcript ids
##    jasparTFs_transcripts_dict_filenames = [os.path.join(species_specific_data_dir, x) for x in os.listdir(species_specific_data_dir) if "jasparTFs.transcripts.single_protein" in x and target_species in x]
##    if len(jasparTFs_transcripts_dict_filenames) > 0:
##        jasparTFs_transcripts_dict_filenames.sort()
##        jasparTFs_transcripts_dict_filename = jasparTFs_transcripts_dict_filenames[-1]
##        jasparTFs_transcripts_dict = load_json(jasparTFs_transcripts_dict_filename)
##    else:
##        jasparTFs_transcripts_dict = {}

    # load CpG score weights
    cpg_data_dir = os.path.join(species_specific_data_dir, "cpg_data")
    cpg_obsexp_weights_dict_filenames = [os.path.join(cpg_data_dir, x) for x in os.listdir(cpg_data_dir) if ".cpg_obsexp_weights." in x and target_species in x]
    if len(cpg_obsexp_weights_dict_filenames) > 0:
        cpg_obsexp_weights_dict_filenames.sort()
        cpg_obsexp_weights_dict_filename = cpg_obsexp_weights_dict_filenames[-1]
        cpg_obsexp_weights_dict = load_json(cpg_obsexp_weights_dict_filename)
        cpg_obsexp_weights_dict = {float(k):float(v) for k,v in cpg_obsexp_weights_dict.items()}
        cpg_obsexp_weights_dict_keys = list(cpg_obsexp_weights_dict.keys())
        cpg_obsexp_weights_dict_keys.sort()
    else:
        cpg_obsexp_weights_dict = {}
        cpg_obsexp_weights_dict_keys = []

    # load GTEx variants
    gtex_variants = {}
    gtex_data_dir = os.path.join(species_specific_data_dir, "gtex_data")
    if os.path.exists(gtex_data_dir):
        gtex_chrom_dict_filename = os.path.join(gtex_data_dir, ".".join([target_species, "gtex_v7", "Chr"+chromosome.upper(), "min_unique", "eqtls", "grch38","msg"]))
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
            gtex_weights_dict = {float(k):float(v) for k,v in gtex_weights_dict.items()}

    # load meta clusters from GTRD project
    gtrd_metaclusters_dict = {}
    gtrd_data_dir = os.path.join(species_specific_data_dir, "gtrd_data")
    if os.path.exists(gtrd_data_dir):
        gtrd_metaclusters_chrom_dict_filename = os.path.join(gtrd_data_dir, ".".join([target_species, "metaclusters", "interval", "Chr"+chromosome.upper(), "clipped", "ordered", "tupled", "msg"]))
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
            #metacluster_overlap_weights_dict = {float(k):float(v) for k,v in metacluster_overlap_weights_dict.items()}

##    # load ATAC-Seq from Encode project
##    atac_seq_dict = {}
##    atac_seq_data_dir = os.path.join(species_specific_data_dir, "atac_data")
##    if os.path.exists(atac_seq_data_dir):
##        atac_seq_chrom_dict_filename = os.path.join(atac_seq_data_dir, ".".join([target_species, "atac-seq", "Chr"+chromosome.upper(), "msg"]))
##
##        if os.path.exists(atac_seq_chrom_dict_filename):
##            atac_seq_dict = load_msgpack(atac_seq_chrom_dict_filename)
##    # load ATAC-Seq dist weights
##    atac_dist_weights_dict = {}
##    if os.path.exists(atac_seq_data_dir):
##        atac_dist_weights_dict_filenames = [os.path.join(atac_seq_data_dir, x) for x in os.listdir(atac_seq_data_dir) if "atac_dist_weights" in x and target_species in x]
##        if len(atac_dist_weights_dict_filenames) > 0:
##            atac_dist_weights_dict_filenames.sort()
##            atac_dist_weights_dict_filename = atac_dist_weights_dict_filenames[-1]
##            atac_dist_weights_dict = load_json(atac_dist_weights_dict_filename)

    # load ATAC-Seq from ChIP-Atlas
    atac_seq_dict = {}
    atac_seq_data_dir = os.path.join(species_specific_data_dir, "atac_data")
    if os.path.exists(atac_seq_data_dir):
        atac_seq_chrom_dict_filename = os.path.join(atac_seq_data_dir, ".".join([target_species, "ATAC_loc_weight", "Chr"+chromosome.upper(), "msg"]))

        if os.path.exists(atac_seq_chrom_dict_filename):
            atac_seq_dict = load_msgpack(atac_seq_chrom_dict_filename)

    # load pre-calculated combined affinity score, by tf, p-values
    cas_pvalues_dict = {}
    cas_pvalues_dict_filename = os.path.join(species_specific_data_dir, ".".join(["CAS_pvalues", "0.1", "tf_ls", "json"]))
    if os.path.exists(cas_pvalues_dict_filename):
        cas_pvalues_dict = load_json(cas_pvalues_dict_filename)


    #return species_pwm_score_threshold_df, gerp_conservation_locations_dict, gerp_conservation_weight_dict, species_group, cage_dict, TF_cage_dict, cage_dist_weights_dict, cage_correlations_dict, cage_corr_weights_dict, atac_dist_weights_dict, metacluster_overlap_weights_dict, cpg_obsexp_weights_dict, cpg_obsexp_weights_dict_keys, gtex_variants, gtex_weights_dict, gtrd_metaclusters_dict, atac_seq_dict, cas_pvalues_dict
    #return species_pwm_score_threshold_df, gerp_conservation_locations_dict, species_group, cage_dict, TF_cage_dict, cage_dist_weights_dict, cage_correlations_dict, cage_corr_weights_dict, atac_dist_weights_dict, metacluster_overlap_weights_dict, cpg_obsexp_weights_dict, cpg_obsexp_weights_dict_keys, gtex_variants, gtex_weights_dict, gtrd_metaclusters_dict, atac_seq_dict, cas_pvalues_dict
    return species_pwm_score_threshold_df, gerp_conservation_locations_dict, species_group, cage_dict, TF_cage_dict, cage_dist_weights_dict, cage_correlations_dict, cage_corr_weights_dict, metacluster_overlap_weights_dict, cpg_obsexp_weights_dict, cpg_obsexp_weights_dict_keys, gtex_variants, gtex_weights_dict, gtrd_metaclusters_dict, atac_seq_dict, cas_pvalues_dict


# overlap_range moved to io_utils; pwm_maker and PWM_scorer moved to pwm.py
# (all re-exported above).

################################################################################
# PWM analysis #################################################################
################################################################################


#def tfbs_finder(transcript_name, alignment, target_tfs_list, TFBS_matrix_dict, target_dir, pwm_score_threshold_dict, species_nt_freq_d, all_pwms_loglikelihood_dict, unaligned2aligned_index_dict, promoter_after_tss, pval, pvalc):
def tfbs_finder(transcript_name, alignment, target_tfs_list, TFBS_matrix_dict, target_dir, pwm_score_threshold_dict, species_nt_freq_d, unaligned2aligned_index_dict, promoter_after_tss, pval, pvalc):
    """
    1.Convert PFM to PWM for each TF in the Jaspar dictionary.
    2.Score all positions in the cleaned sequence
    3.If score is greater than or equal to precomputed threshold, then keep, otherwise set to zero
    4.Return dictionary of pwm scores for [species][tf_name][strand]
    """
    start_time = time.time()

    tfbss_found_dict_outfilename = os.path.join(target_dir, "TFBSs_found.all.json")

    # Determine if the analysis has been done already, load results if so
    if os.path.isfile(tfbss_found_dict_outfilename):
        logging.info(" ".join(["tfbss_found_dict already exists: loading"]))
        tfbss_found_dict = load_json(tfbss_found_dict_outfilename)

    # If results don't already exist, time to party
    else:
        tfbss_found_dict = {}
        align_chars = '-N .'
        mononuc_pwm_dict = {"A":0,"C":1,"G":2,"T":3}

        entry = alignment[0]
        species = entry['species']

        # Remove alignment/ambiguous characters from the sequences
        cleaned_seq = entry['seq']
        for char in align_chars:
            cleaned_seq = cleaned_seq.replace(char,"")
##        entry_seqrecord = SeqRecord(Seq(cleaned_seq, alphabet=IUPAC.unambiguous_dna), id=species)
        entry_seqrecord = SeqRecord(Seq(cleaned_seq), id=species)
        forward_seq = str(entry_seqrecord.seq)
        reverse_seq = str(entry_seqrecord.seq.reverse_complement())
        seq_dict = {"+1": forward_seq, "-1":reverse_seq}

        # generate background frequencies of each mono-nucleotide for forward and reverse strands
        bg_nuc_freq_dict = {}
        #neg_bg_nuc_freq_dict = {}

        # bool - check if species is human; use previously calculated nt freqs
        if species == "homo_sapiens":
            # https://arxiv.org/pdf/q-bio/0611041.pdf
            # empirical data from complete genome
            bg_nuc_freq_dict = {'A':0.292, 'C':0.207, 'G':0.207, 'T':0.292}
            #neg_bg_nuc_freq_dict = {'A':0.292, 'C':0.207, 'G':0.207, 'T':0.292}

        # bool - check if species is in our calculated species_nt_freq_d; use previously calculated nt freqs
        elif species in species_nt_freq_d:
            bg_nuc_freq_dict = species_nt_freq_d[species]

        # iterate through each tf_name and its motif
        for tf_name in target_tfs_list:
            if tf_name in TFBS_matrix_dict:
                tf_motif = TFBS_matrix_dict[tf_name]
                motif_length = len(tf_motif[0])

                if motif_length > 0:

                    # retrieve precomputed threshold and other information required for calculating the pvalue of the score
                    tf_pwm_score_threshold_dict = pwm_score_threshold_dict[tf_name]
                    pvals_scores_list = [[k,v] for k,v in tf_pwm_score_threshold_dict.items()]
                    pvals_scores_list_sorted = sorted(pvals_scores_list, key=itemgetter(1))
                    scores_list_sorted = [x[1] for x in pvals_scores_list_sorted]
                    pvals_list = [x[0] for x in pvals_scores_list_sorted]
                    pvals_list.sort()

                    if pval in tf_pwm_score_threshold_dict:
                        tf_pwm_score_threshold = tf_pwm_score_threshold_dict[pval]
                    else:
                        if pval == 1:
                            tf_pwm_score_threshold = -100000
                        else:
                            pval = pvals_list[bisect_left(pvals_list, pval)]
                            tf_pwm_score_threshold = tf_pwm_score_threshold_dict[pval]

                    tfbss_found_dict[tf_name] = []

                    # iterate through the forward and reverse strand sequences
                    for strand, seq in seq_dict.items():
                        pwm = pwm_maker(strand, motif_length, tf_motif, bg_nuc_freq_dict)
                        #pwm = pwm_maker(strand, motif_length, tf_motif, bg_nuc_freq_dict, neg_bg_nuc_freq_dict)

                        seq_length = len(seq)
                        # iterate through the nt sequence, extract a current frame based on the motif size and score
                        for i in range(0, seq_length - motif_length):
                            current_frame = seq[i:i+motif_length]
                            current_frame_score = PWM_scorer(current_frame, pwm, mononuc_pwm_dict, 'mono')

                            # keep results that are above the precomputed threshold
                            if current_frame_score >= tf_pwm_score_threshold:
                                current_frame_score = round(current_frame_score, 2)

                                pval_index = bisect_left(scores_list_sorted, current_frame_score)
                                if pval_index >= len(pvals_scores_list_sorted):
                                    pval_index = -1

                                else:
                                    pval_index -= 1

                                # account for pvalues which are larger than the current largest, so that they can be distinguished appropriately in the result table.
                                if pval_index < 0:
                                    current_frame_score_pvalue = ">" + str(pvals_scores_list_sorted[0][0])
                                else:
                                    current_frame_score_pvalue = str(pvals_scores_list_sorted[pval_index][0])

                                hit_loc_start, hit_loc_end, hit_loc_before_TSS_start, hit_loc_before_TSS_end = start_end_found_motif(i, strand, seq_length, promoter_after_tss, motif_length)

                                # identify position in alignment from start of found motif in unaligned sequence
                                aligned_position = unaligned2aligned_index_dict[species][hit_loc_start]

                                # add to results dictionary by tf_name
                                tfbss_found_dict[tf_name].append([current_frame, strand, hit_loc_start, hit_loc_end, hit_loc_before_TSS_start, hit_loc_before_TSS_end, current_frame_score, current_frame_score_pvalue])

    end_time = time.time()
    logging.info(" ".join(["total time for tfbs_finder() for this transcript:", str(end_time - start_time), "seconds"]))

    return tfbss_found_dict


# Scoring (CAS engine) moved to tfbs_footprinter3/scoring.py:
#   calcCombinedAffinityPvalue, find_clusters, eqtl_overlap_likelihood,
#   eqtls_weights_summing, cage_correlations_summing_preparation,
#   cage_correlations_summing, cage_weights_summing, atac_weights_summing,
#   metacluster_weights_summing, gerp_weights_summing, cpg_weights_summing.


# clean_jaspar_names moved to tfbs_footprinter3/data_parsing.py (re-exported above).


################################################################################
# Alignment Manipulation #######################################################
################################################################################
def retrieve_genome_aligned(target_species, chromosome, strand, promoter_start, promoter_end):
    """
    Takes as input target_species CCDS start position and size of promoter to be extracted.  Retrieves genome aligned,
    corresponding regions in all orthologs.
    ***Alignment no longer necessary as the use of pre-computed GERP scores means that no conservation calculation needs to be made.
    Additionally, the newest implementation of Ensembl REST alignment will not retrieve target species sequence,
    if there is no alignment with other species at that location.
    For example, a request for homo_sapiens alignment from chr1:1-10,000 will only return the locations where an alignment
    exists with the target species group (e.g. mammals_low).  This may only exist at chr1:2000-9500.
    """

    query_type = "/sequence/region/"
    pre_options = target_species + "/" + chromosome + ":" + str(promoter_start) + "-" + str(promoter_end) + ":" + str(strand)
    options = pre_options + "?content-type=application/json"
    target_only_decoded = ensemblrest(query_type, options, 'json', "", log=True)
    if 'seq' in target_only_decoded:
        target_only_decoded['species'] = target_species
        alignment = [target_only_decoded]

    else:
        alignment = []

    return alignment


# fasta_writer, remove_non_ACGT, remove_gap_only, remove_duplicate_species,
# load_genome_aligned moved to tfbs_footprinter3/alignment.py (re-exported above).


def alignment_tools(ensembl_aligned_filename, cleaned_aligned_filename, target_species, chromosome, strand, promoter_start, promoter_end):
    """
    Return cleaned alignment for further analysis.
    """

    # if cleaned alignment file doesn't exist, or the size is zero.
    if not os.path.isfile(cleaned_aligned_filename) or (os.path.isfile(cleaned_aligned_filename) and os.path.getsize(cleaned_aligned_filename) == 0):

        # If uncleaned Ensembl alignment file doesn't exist, or the size is zero: retrieve from Ensembl, write to file.
        if not os.path.isfile(ensembl_aligned_filename) or (os.path.isfile(ensembl_aligned_filename) and os.path.getsize(ensembl_aligned_filename) == 0):
            alignment = retrieve_genome_aligned(target_species, chromosome, strand, promoter_start, promoter_end)
            fasta_writer(alignment, ensembl_aligned_filename)

        # If uncleaned Ensembl file exists and size is not zero: clean, write to cleaned filename.
        if os.path.isfile(ensembl_aligned_filename) and (os.path.isfile(ensembl_aligned_filename) and os.path.getsize(ensembl_aligned_filename) > 0):
            alignment = load_genome_aligned(ensembl_aligned_filename)
            alignment = remove_non_ACGT(alignment)
            alignment = remove_duplicate_species(alignment, target_species)
            alignment = remove_gap_only(alignment)
            fasta_writer(alignment, cleaned_aligned_filename)

        # Uncleaned alignment file still doesn't exist (or size is zero): note in logfile.
        else:
            logging.warning(" ".join(["No ensembl alignment, or size is zero"]))
            alignment = []

    # Cleaned alignment file exists and size is not zero: load cleaned alignment.
    else:
        alignment = load_genome_aligned(cleaned_aligned_filename)

    return alignment


def test_transcript_id(decoded_json_description, transcript_id):
    """
    Test if the dictionary for the target transcript id:
    Indicates it is a transcript.
    Doesn't contain errors/Is complete.
    """

    transcript_id_pass = False

    if 'error' not in decoded_json_description:
        if 'object_type' in decoded_json_description:
            if decoded_json_description['object_type'].lower() == 'transcript':
                transcript_id_pass = True
            else:
                logging.warning(" ".join([transcript_id, "This input does not appear to be a valid Ensembl transcript ID.  Ensembl REST defines it as:", decoded_json_description['object_type']]))
        else:
            logging.warning(" ".join([transcript_id, "This input does not appear to be a valid Ensembl transcript ID.  Please check it for errors."]))
    else:
        logging.warning(" ".join(["Ensembl REST responds with error:", decoded_json_description['error']]))

    return transcript_id_pass


def transfabulator(transcript, transcript_dict_filename):
    """
    Given a transcript ID, retrieve Ensembl descriptive data for that transcript.
    """

    retrieve_transcript_data = True
    # load transcript position data from json file if it already exists
    if os.path.isfile(transcript_dict_filename):
        if os.path.getsize(transcript_dict_filename) != 0:
            logging.info(" ".join(["transcript_dict already exists: loading"]))
            decoded_json_description = load_json(transcript_dict_filename)
            retrieve_transcript_data = False

    # retrieve transcript position data from json file if it does not exist
    # or contains no data.
    if retrieve_transcript_data:
        # Set parameters for retrieving Ensembl data via REST
        query_type = '/lookup/id/'
        options = '?feature=transcript;content-type=application/json'

        # populate 'transcript_dict' dictionary with sub-dictionaries.
        # key[transcript_id] = {chromosome, strand, start, end} for each ensembl transcript id
        decoded_json_description = ensemblrest(query_type, options, 'json', transcript, log=True)
        decoded_json_description = {k.lower():v for k,v in decoded_json_description.items()}

    return decoded_json_description


def transcript_data_retrieve(decoded_json_description, transcript_dict_filename, promoter_before_tss, promoter_after_tss):
    """
    The retrieved transcript data is assumed complete, extract important data.
    Write json data to file.
    Based on user-defined values for target region (referenced to TSS),
    calculate genomic coordinates of target region.
    """

    # Extract position data
    chromosome = decoded_json_description['seq_region_name']
    chr_start = decoded_json_description['start']
    chr_end = decoded_json_description['end']
    strand = decoded_json_description['strand']
    ens_gene_id = decoded_json_description['parent']
    decoded_id = decoded_json_description['id']
    target_species = decoded_json_description['species']

    if 'display_name' in decoded_json_description:
        transcript_name = decoded_json_description['display_name']
    else:
        transcript_name = "-".join([ens_gene_id, decoded_id])
    dump_json(transcript_dict_filename, decoded_json_description)

    if strand == 1:
        tss = chr_start
        #[promoter_start][promoter_end][TSS=chr_start][>GENE--->][chr_end]
        promoter_start = tss - promoter_before_tss
        promoter_end = tss - 1 + promoter_after_tss

    if strand == -1:
        tss = chr_end
        #[chr_start][<---GENE<][TSS=chr_end][promoter_start][promoter_end]
        promoter_start = tss + 1 - promoter_after_tss
        promoter_end = tss + promoter_before_tss

    return target_species, transcript_name, ens_gene_id, chromosome, tss, strand, promoter_start, promoter_end, chr_start, chr_end


def gene_data_retrieve(gene_dict_filename, ens_gene_id):
    """
    Retrieve gene data for the parent gene of the target transcript.
    """

    # determine likelihood of overlapping an eQTL at all.
    # Set parameters for retrieving Ensembl data via REST

    decoded_json_description = load_json(gene_dict_filename)

    if not decoded_json_description or len(decoded_json_description)==0:
        query_type = '/lookup/id/'
        options = '?feature=transcript;content-type=application/json'

        # retrieve_gene_len
        decoded_json_description = ensemblrest(query_type, options, 'json', ens_gene_id, log=True)
        decoded_json_description = {k.lower():v for k,v in decoded_json_description.items()}

        dump_json(gene_dict_filename, decoded_json_description)

    gene_start = decoded_json_description['start']
    gene_end = decoded_json_description['end']

    if 'display_name' in decoded_json_description:
        gene_name = decoded_json_description['display_name'].upper()
    else:
        gene_name = decoded_json_description['id']

    gene_len = gene_end - gene_start

    return gene_name, gene_len


################################################################################
# Regulatory & Conservation Features ###########################################
################################################################################
def retrieve_regulatory(chromosome, strand, promoter_start, promoter_end, regulatory_decoded_filename, target_species):
    """
    Retrieve ensembl JSON data for regulatory features within the coordinates provided.
    """

    # determine if the regulatory data has already been retrieved, if so load, if not retrieve.
    if os.path.isfile(regulatory_decoded_filename):
        logging.info(" ".join(["regulatory_decoded already exists: loading"]))
        regulatory_decoded = load_json(regulatory_decoded_filename)

    else:
        query_type = "/overlap/region/"
        pre_options = target_species + "/" + chromosome + ":" + str(promoter_start) + "-" + str(promoter_end) + ":" + str(strand)

        options = pre_options + "?feature=regulatory;content-type=application/json"
        regulatory_decoded = ensemblrest(query_type, options, 'json', "", log=True)

        # rename the Ensembl regulatory elements so that they don't overtake the space available for names.
        for reg in regulatory_decoded:
            if "description" in reg:
                if reg["description"] == "Transcription factor binding site":
                    reg["description"] = "Pred. TFBS"
                if reg["description"] == "Open chromatin region":
                    reg["description"] = "Open chromatin"
                if reg["description"] == "Predicted promoter":
                    reg["description"] = "Pred. promoter"

        dump_json(regulatory_decoded_filename, regulatory_decoded)

    return regulatory_decoded


# distance_solve moved to io_utils (re-exported above).


##    # variable x-ticks
##    dist = promoter_before_tss + promoter_after_tss
##    rough_interval = dist/10
##    power = int(np.log10(rough_interval))
##    xtick_jump = (rough_interval/(10**power)) * 10**power
##    ax3.set_xticks(range(-1 * promoter_before_tss, promoter_after_tss + 1, xtick_jump))


################################################################################
# Initiating Variables #########################################################
################################################################################

################################################################################
# Execution ####################################################################
################################################################################
signal.signal(signal.SIGINT, signal_handler)

def main():
    """
    All the things.
    """

    total_time_start = time.time()
    print("Executing tfbs_footprinter3 version %s." % __version__)

    # Create directory for results
    output_dir = os.path.join(curdir, "tfbs_results")
    directory_creator(output_dir)

    # begin timing and logging
    logging.basicConfig(filename=os.path.join(os.path.dirname(output_dir), 'TFBS_footprinter3.log'), level=logging.INFO, format='%(asctime)s:    [%(levelname)s]    %(message)s')
    logging.info(" ".join(["***NEW SET OF ANALYSES HAS BEGUN***"]))

    if is_online():
        args_lists, exp_data_update, nofigure = get_args()

        # if experimental data dir does not exist or user has requested an exp data update, then update.
        experimental_data_present = experimentalDataUpdater(exp_data_update)

        if experimental_data_present:
            if len(args_lists) > 0:
                # analysis variables
                # non-species-specific
                # dictionary of thresholds for each TF
                ##    # updated version, which requires the presence of a current versions file
                ##    pwm_score_threshold_dict_filename = os.path.join(experimental_data_dir, current_versions["jaspar_thresholds"])
##                pwm_score_threshold_dict_filename = os.path.join(script_dir, 'data/all_tfs_thresholds.jaspar_2018.1.json')
##                pwm_score_threshold_dicta = load_json(pwm_score_threshold_dict_filename)
##                pwm_score_threshold_dict = {}
##                for k,v in pwm_score_threshold_dicta.items():
##                    pwm_score_threshold_dict[k] = {float(kk):float(vv) for kk,vv in v.items()}

                # load mono-nuc PFMs
                TFBS_matrix_filename = os.path.join(script_dir, 'data/pwms.json')
                TFBS_matrix_dict = load_json(TFBS_matrix_filename)
                TFBS_matrix_dict = {k.upper():v for k,v in TFBS_matrix_dict.items()}

##                # load JASPAR PWM score weights
##                all_pwms_loglikelihood_dict_filename = os.path.join(script_dir, 'data/all_pwms_loglikelihood_dict.reduced.msg')
##                all_pwms_loglikelihood_dict = load_msgpack(all_pwms_loglikelihood_dict_filename)

                # load species nt frequencies
                species_nt_freq_fn = os.path.join(script_dir, 'data/species_nt_freq.json')
                species_nt_freq_d = load_json(species_nt_freq_fn)

                last_target_species = None
                last_chromosome = None

            for i, args_list in enumerate(args_lists):
                args, transcript_ids_filename, transcript_id, target_tfs_filename, promoter_before_tss, promoter_after_tss, top_x_tfs_count, pval, pvalc = args_list
                print("Ensembl transcript id:", transcript_id)
                logging.info("\n"+" ".join(["***ANALYSIS OF A NEW TRANSCRIPT HAS BEGUN:", transcript_id]))
                logging.info(" ".join(["Arguments used in this run:", str(args_list)]))

                # target dir naming
                start_end = "("+"_".join([str(promoter_before_tss), str(promoter_after_tss)])+")"
                target_dir_name = "_".join([transcript_id+start_end, str(pval)])
                target_dir = os.path.join(output_dir, target_dir_name)

                # declare all possible results filenames.
                tfbss_found_dict_outfilename = os.path.join(target_dir, "TFBSs_found.all.json")
                cluster_dict_filename = os.path.join(target_dir, "cluster_dict.json")
                ensembl_aligned_filename = os.path.join(target_dir, "alignment_uncleaned.fasta")
                cleaned_aligned_filename = os.path.join(target_dir, "alignment_cleaned.fasta")
                transcript_dict_filename = os.path.join(target_dir, "transcript_dict.json")
                gene_dict_filename = os.path.join(target_dir, "gene_dict.json")
                regulatory_decoded_filename = os.path.join(target_dir, "regulatory_decoded.json")
                sortedclusters_table_filename = os.path.join(target_dir, ".".join(["TFBSs_found", "sortedclusters", "csv"]))

                # check if results have been created for this query.
                required_results_filenames = [cluster_dict_filename, ensembl_aligned_filename, cleaned_aligned_filename, transcript_dict_filename, gene_dict_filename, regulatory_decoded_filename, sortedclusters_table_filename]
                results_files_exist = all([os.path.exists(x) for x in required_results_filenames])

                if not results_files_exist:
                    # identify if the target transcript id exists in Ensembl
                    decoded_json_description = transfabulator(transcript_id, transcript_dict_filename)
                    transcript_id_pass = test_transcript_id(decoded_json_description, transcript_id)

                    # parse target transcript id data from successful retrieval, and continue
                    if transcript_id_pass:
                        # create target output dir
                        directory_creator(target_dir)
                        logging.info(" ".join(["Results will be output to:", target_dir]))
                        target_species, transcript_name, ens_gene_id, chromosome, tss, strand, promoter_start, promoter_end, chr_start, chr_end = transcript_data_retrieve(decoded_json_description, transcript_dict_filename, promoter_before_tss, promoter_after_tss)
                        gene_name, gene_len = gene_data_retrieve(gene_dict_filename, ens_gene_id)

                        # species-specific
                        species_specific_data_dir = os.path.join(script_dir, 'data', target_species)
                        experimentaldata(target_species)
                        if target_species != last_target_species or chromosome != last_chromosome:
                            if os.path.exists(species_specific_data_dir):
                                #species_pwm_score_threshold_df, gerp_conservation_locations_dict, gerp_conservation_weight_dict, species_group, cage_dict, TF_cage_dict, cage_dist_weights_dict, cage_correlations_dict, cage_corr_weights_dict, atac_dist_weights_dict, metacluster_overlap_weights_dict, cpg_obsexp_weights_dict, cpg_obsexp_weights_dict_keys, gtex_variants, gtex_weights_dict, gtrd_metaclusters_dict, atac_seq_dict, cas_pvalues_dict = species_specific_data(target_species, chromosome, species_specific_data_dir)
                                #species_pwm_score_threshold_df, gerp_conservation_locations_dict, species_group, cage_dict, TF_cage_dict, cage_dist_weights_dict, cage_correlations_dict, cage_corr_weights_dict, atac_dist_weights_dict, metacluster_overlap_weights_dict, cpg_obsexp_weights_dict, cpg_obsexp_weights_dict_keys, gtex_variants, gtex_weights_dict, gtrd_metaclusters_dict, atac_seq_dict, cas_pvalues_dict = species_specific_data(target_species, chromosome, species_specific_data_dir)
                                species_pwm_score_threshold_df, gerp_conservation_locations_dict, species_group, cage_dict, TF_cage_dict, cage_dist_weights_dict, cage_correlations_dict, cage_corr_weights_dict, metacluster_overlap_weights_dict, cpg_obsexp_weights_dict, cpg_obsexp_weights_dict_keys, gtex_variants, gtex_weights_dict, gtrd_metaclusters_dict, atac_seq_dict, cas_pvalues_dict = species_specific_data(target_species, chromosome, species_specific_data_dir)
                            last_target_species = target_species
                            last_chromosome = chromosome

                        # format - target species tf:pvalue:score
                        pwm_score_threshold_dict = {}
                        if len(species_pwm_score_threshold_df) >0:
                            pwm_score_threshold_dict = species_pwm_score_threshold_df[["tf_name", "p_value", "score"]].groupby('tf_name')[["p_value", "score"]].apply(lambda x: dict(x.to_numpy())).to_dict()

                        # load target tfs
                        if target_tfs_filename == "" or target_tfs_filename is None:
                            target_tfs_filename = None
                            target_tfs_list = TFBS_matrix_dict.keys()

                        if target_tfs_filename is not None:
                            target_tfs_list = parse_tf_ids(target_tfs_filename)
                            target_tfs_list = compare_tfs_list_jaspar(target_tfs_list, TFBS_matrix_dict)

                        # filenames for alignment and ensembl regulatory data
                        alignment = alignment_tools(ensembl_aligned_filename, cleaned_aligned_filename, target_species, chromosome, strand, promoter_start, promoter_end)

                        # continue if there is an alignment from Ensembl, and after cleaning
                        if len(alignment) > 0:

                            target_species_row = alignment[0]
                            alignment_len = len(target_species_row['seq'].replace('-',''))

                            # retrieve regulatory
                            regulatory_decoded = retrieve_regulatory(chromosome, strand, promoter_start, promoter_end, regulatory_decoded_filename, target_species)
                            converted_reg_dict = reg_position_translate(tss,regulatory_decoded,promoter_start,promoter_end,strand,promoter_before_tss,promoter_after_tss)

                            # conservation
                            converted_gerps_in_promoter = gerp_positions_translate(target_dir, gerp_conservation_locations_dict, chromosome, strand, promoter_start, promoter_end, tss)

                            # identify information content of each column of the alignment
                            cpg_list = CpG(cleaned_aligned_filename)

                            # identify CAGEs in proximity to Ensembl TSS, convert for plotting
                            converted_cages = cage_position_translate(target_species, gene_name, ens_gene_id, transcript_id,tss,cage_dict,promoter_start,promoter_end,strand,promoter_before_tss,promoter_after_tss)

                            # identify eQTLs in proximity to Ensembl TSS, convert for plotting
                            converted_eqtls = gtex_position_translate(ens_gene_id,gtex_variants,tss,promoter_start,promoter_end,strand,promoter_before_tss,promoter_after_tss)

                            # GTRD metaclusters
                            #converted_metaclusters_in_promoter = gtrd_positions_translate(target_dir, gtrd_metaclusters_dict, chromosome, strand, promoter_start, promoter_end, tss)
                            converted_metaclusters_in_promoter, metacluster_in_promoter_counts = gtrd_positions_translate(target_dir, gtrd_metaclusters_dict, chromosome, strand, promoter_start, promoter_end, tss)

                            # ATAC-seq data
                            converted_atac_seqs_in_promoter = atac_pos_translate(atac_seq_dict, chromosome, strand, promoter_start, promoter_end, tss)

                            # create index of aligned to unaligned positions
                            unaligned2aligned_index_dict = unaligned2aligned_indexes(cleaned_aligned_filename)

                            if not (os.path.exists(cluster_dict_filename) and os.path.exists(sortedclusters_table_filename)):
                                # score alignment for tfbss
                                #tfbss_found_dict = tfbs_finder(transcript_name, alignment, target_tfs_list, TFBS_matrix_dict, target_dir, pwm_score_threshold_dict, species_nt_freq_d, all_pwms_loglikelihood_dict, unaligned2aligned_index_dict, promoter_after_tss, pval, pvalc)
                                tfbss_found_dict = tfbs_finder(transcript_name, alignment, target_tfs_list, TFBS_matrix_dict, target_dir, pwm_score_threshold_dict, species_nt_freq_d, unaligned2aligned_index_dict, promoter_after_tss, pval, pvalc)

                                # sort through scores, identify hits in target_species supported in other species
                                #cluster_dict = find_clusters(gene_name, ens_gene_id, chr_start, chr_end, alignment, target_species, chromosome, tfbss_found_dict, cleaned_aligned_filename, converted_gerps_in_promoter, gerp_conservation_weight_dict,  converted_cages, converted_metaclusters_in_promoter, converted_atac_seqs_in_promoter, converted_eqtls, gtex_weights_dict, transcript_id, cage_dict, TF_cage_dict, cage_dist_weights_dict, atac_dist_weights_dict, metacluster_overlap_weights_dict, cpg_list, cpg_obsexp_weights_dict, cpg_obsexp_weights_dict_keys, cage_correlations_dict, cage_corr_weights_dict, gtex_variants, gene_len, cas_pvalues_dict, pvalc)
                                #cluster_dict = find_clusters(gene_name, ens_gene_id, chr_start, chr_end, alignment, target_species, chromosome, tfbss_found_dict, cleaned_aligned_filename, converted_gerps_in_promoter, converted_cages, converted_metaclusters_in_promoter, converted_atac_seqs_in_promoter, converted_eqtls, gtex_weights_dict, transcript_id, cage_dict, TF_cage_dict, cage_dist_weights_dict, atac_dist_weights_dict, metacluster_overlap_weights_dict, cpg_list, cpg_obsexp_weights_dict, cpg_obsexp_weights_dict_keys, cage_correlations_dict, cage_corr_weights_dict, gtex_variants, gene_len, cas_pvalues_dict, pvalc)
                                cluster_dict = find_clusters(gene_name, ens_gene_id, chr_start, chr_end, alignment, target_species, chromosome, tfbss_found_dict, cleaned_aligned_filename, converted_gerps_in_promoter,  converted_cages, converted_metaclusters_in_promoter, metacluster_in_promoter_counts, converted_atac_seqs_in_promoter, converted_eqtls, gtex_weights_dict, transcript_id, cage_dict, TF_cage_dict, cage_dist_weights_dict, metacluster_overlap_weights_dict, cpg_list, cpg_obsexp_weights_dict, cpg_obsexp_weights_dict_keys, cage_correlations_dict, cage_corr_weights_dict, gtex_variants, gene_len, cas_pvalues_dict, pvalc)
                                tfbss_found_dict.clear()
                                dump_json(cluster_dict_filename, cluster_dict)

                            else:
                                cluster_dict = load_json(cluster_dict_filename)

                            # sort the target_species hits supported by other species
                            sorted_clusters_target_species_hits_list = sort_target_species_hits(cluster_dict)

                            target_species_hits_table_writer(sorted_clusters_target_species_hits_list, sortedclusters_table_filename)

                            if len(sorted_clusters_target_species_hits_list) > 0:

                                # extract the top x target_species hits supported by other species
                                top_x_greatest_hits_dict = top_x_greatest_hits(sorted_clusters_target_species_hits_list, top_x_tfs_count)
                                top_greatest_hits_dict = top_greatest_hits(sorted_clusters_target_species_hits_list, top_x_tfs_count)

                                # plot the top x target_species hits
                                if not nofigure:
                                    if len(top_x_greatest_hits_dict) > 0:
                                        plot_promoter(target_species, transcript_id, species_group, alignment, alignment_len, promoter_before_tss, promoter_after_tss, transcript_name, top_x_greatest_hits_dict, target_dir, converted_reg_dict, converted_gerps_in_promoter, cpg_list, converted_cages, converted_metaclusters_in_promoter, converted_atac_seqs_in_promoter, converted_eqtls)
    ##                                    plot_promoter_all(target_species, transcript_id, species_group, alignment, alignment_len, promoter_before_tss, promoter_after_tss, transcript_name, top_greatest_hits_dict, target_dir, converted_reg_dict, converted_gerps_in_promoter, cpg_list, converted_cages, converted_metaclusters_in_promoter, converted_atac_seqs_in_promoter, converted_eqtls)

            total_time_end = time.time()
            logging.info(" ".join(["Total time for", str(len(args_lists)), "transcripts:", str(total_time_end - total_time_start), "seconds"]) + "\n\n")

    else:
        print("System does not appear to be connected to the internet.  Exiting TFBS_footprinter3.")
