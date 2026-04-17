#!/usr/bin/python
# Python vers. 3.8.0 ###########################################################
import argparse
import logging
import os
import signal

# Libraries ####################################################################
import textwrap
import time

import matplotlib

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
from tfbs_footprinter3.data_loader import (  # noqa: F401
    experimentaldata,
    experimentalDataUpdater,
    experimentalDataUpdater_beta,
    species_specific_data,
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
from tfbs_footprinter3.pipeline import (  # noqa: F401
    alignment_tools,
    gene_data_retrieve,
    retrieve_genome_aligned,
    retrieve_regulatory,
    test_transcript_id,
    tfbs_finder,
    transcript_data_retrieve,
    transfabulator,
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


# overlap_range moved to io_utils; pwm_maker and PWM_scorer moved to pwm.py
# (all re-exported above).

################################################################################
# PWM analysis #################################################################
################################################################################


#def tfbs_finder(transcript_name, alignment, target_tfs_list, TFBS_matrix_dict, target_dir, pwm_score_threshold_dict, species_nt_freq_d, all_pwms_loglikelihood_dict, unaligned2aligned_index_dict, promoter_after_tss, pval, pvalc):


# Scoring (CAS engine) moved to tfbs_footprinter3/scoring.py:
#   calcCombinedAffinityPvalue, find_clusters, eqtl_overlap_likelihood,
#   eqtls_weights_summing, cage_correlations_summing_preparation,
#   cage_correlations_summing, cage_weights_summing, atac_weights_summing,
#   metacluster_weights_summing, gerp_weights_summing, cpg_weights_summing.


# clean_jaspar_names moved to tfbs_footprinter3/data_parsing.py (re-exported above).


# fasta_writer, remove_non_ACGT, remove_gap_only, remove_duplicate_species,
# load_genome_aligned moved to tfbs_footprinter3/alignment.py (re-exported above).


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
