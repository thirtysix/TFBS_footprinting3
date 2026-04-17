"""Command-line entry point for tfbs_footprinter3.

Wires together the pipeline modules: argument parsing, experimental-data
update, per-transcript orchestration (transcript metadata -> alignment
-> PWM scan -> CAS scoring -> output + optional figure).

`main()` is the console-script entry point declared in pyproject.toml:

    [project.scripts]
    tfbs_footprinter3 = "tfbs_footprinter3.tfbs_footprinter3:main"

The monolith shim `tfbs_footprinter3.tfbs_footprinter3` re-exports
`main` from here for backward compatibility with that entry point path.
"""
from __future__ import annotations

import argparse
import logging
import os
import signal
import textwrap
import time

from tfbs_footprinter3 import __version__
from tfbs_footprinter3.data_loader import (
    experimentaldata,
    experimentalDataUpdater,
    species_specific_data,
)
from tfbs_footprinter3.data_parsing import (
    compare_tfs_list_jaspar,
    file_to_datalist,
    parse_tf_ids,
    parse_transcript_ids,
)
from tfbs_footprinter3.io_utils import (
    directory_creator,
    dump_json,
    is_online,
    load_json,
    load_msgpack,  # noqa: F401  (re-exported for backward-compat)
    signal_handler,
)
from tfbs_footprinter3.output import (
    sort_target_species_hits,
    target_species_hits_table_writer,
    target_species_hits_table_writer_parquet,
    top_greatest_hits,
    top_x_greatest_hits,
)
from tfbs_footprinter3.pipeline import (
    alignment_tools,
    gene_data_retrieve,
    retrieve_regulatory,
    test_transcript_id,
    tfbs_finder,
    transcript_data_retrieve,
    transfabulator,
)
from tfbs_footprinter3.plotting import plot_promoter
from tfbs_footprinter3.scoring import find_clusters
from tfbs_footprinter3.translators import (
    CpG,
    atac_pos_translate,
    cage_position_translate,
    gerp_positions_translate,
    gtex_position_translate,
    gtrd_positions_translate,
    reg_position_translate,
    unaligned2aligned_indexes,
)

script_dir = os.path.dirname(__file__)
curdir = os.getcwd()


def _prefetch_and_sort_transcripts(args_lists, output_dir):
    """Group transcripts by (species, chromosome) so the per-transcript
    processing loop gets maximum cache hits on species_specific_data.

    Does a single up-front pass over args_lists that:
      * Calls transfabulator for each transcript (hits Ensembl REST once,
        or reads the transcript_dict.json cache if a prior run wrote it).
      * Saves transcript_dict.json under each target_dir ONLY when the
        response validates as a real transcript, so failed IDs don't
        litter the results directory with empty caches.
      * Returns args_lists reordered by (species, seq_region_name). Within
        each (species, chromosome) bucket the original relative order is
        preserved (Python's sort is stable).

    For a user who organized their input by chromosome this is a no-op.
    For anyone else it cuts species_specific_data calls from one-per-
    transcript down to one-per-unique-(species, chromosome) — and that
    load is the heaviest one-time cost in the pipeline.
    """
    if len(args_lists) <= 1:
        return args_lists

    annotated = []
    unknown_key = ("~~", "~~")  # sorts after real species/chromosomes
    for al in args_lists:
        # args_list layout matches get_args():
        #   [args_ns, transcript_ids_filename, transcript_id, target_tfs_filename,
        #    promoter_before_tss, promoter_after_tss, top_x_tfs_count, pval, pvalc]
        transcript_id = al[2]
        promoter_before_tss = al[4]
        promoter_after_tss = al[5]
        pval = al[7]

        start_end = "(" + "_".join([str(promoter_before_tss), str(promoter_after_tss)]) + ")"
        target_dir_name = "_".join([transcript_id + start_end, str(pval)])
        target_dir = os.path.join(output_dir, target_dir_name)
        transcript_dict_filename = os.path.join(target_dir, "transcript_dict.json")

        if os.path.isfile(transcript_dict_filename) and os.path.getsize(transcript_dict_filename) > 0:
            decoded = load_json(transcript_dict_filename)
        else:
            decoded = transfabulator(transcript_id, transcript_dict_filename)
            if test_transcript_id(decoded, transcript_id):
                directory_creator(target_dir)
                dump_json(transcript_dict_filename, decoded)

        if decoded and "species" in decoded and "seq_region_name" in decoded:
            key = (decoded["species"], str(decoded["seq_region_name"]))
        else:
            key = unknown_key
        annotated.append((key, al))

    annotated.sort(key=lambda pair: pair[0])
    return [al for _, al in annotated]


def get_args():
    """
    Retrieve arguments provided by the user.
    """

    # Instantiate the parser
    parser = argparse.ArgumentParser(
        prog="tfbs_footprinter3",
        formatter_class=argparse.RawDescriptionHelpFormatter,
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

    parser.add_argument('--tf_ids_file', '-tfs', metavar='', type=str, default=None, help='Optional: Location of a file containing a limited list of Jaspar TFs to use in scoring alignment \
                                                                                                (see sample file tf_ids.txt at https://github.com/thirtysix/TFBS_footprinting) [default: all Jaspar TFs]')

    parser.add_argument('--promoter_before_tss', '-pb', metavar='', choices=range(-100000, 100001), type=int, default=900,
                        help='(0-100,000) [default: 900] - Number (integer) of nucleotides upstream of TSS to include in analysis.  If this number is negative the start point will be downstream of the TSS, the end point will then need to be further downstream.')

    parser.add_argument('--promoter_after_tss', '-pa', metavar='', choices=range(-100000, 100001), type=int, default=100,
                        help='(0-100,000) [default: 100] - Number (integer) of nucleotides downstream of TSS to include in analysis.  If this number is negative the end point will be upstream of the TSS.  The start point will then need to be further upstream.')

    parser.add_argument('--top_x_tfs', '-tx', metavar='', choices=range(1, 21), type=int, default=10,
                        help='(1-20) [default: 10] - Number (integer) of unique TFs to include in output .svg figure.')

    # for now pvalue refers to the PWM score, in the future it will need to relate to the combined affinity score
    parser.add_argument('--pval', '-p', type=float, default=1, help='P-value (float) for PWM score cutoff (range: 1 (all results) to 0.0000001; in divisions of 10 (i.e. 1, 0.1, 0.01, 0.001 etc.) [default: 0.01]')

    parser.add_argument('--pvalc', '-pc', type=float, default=1, help='P-value (float) for PWM score cutoff (range: 1 (all results) to 0.0000001; in divisions of 10 (i.e. 1, 0.1, 0.01, 0.001 etc.) [default: 0.01]')

    parser.add_argument('--exp_data_update', '-update', action="store_true", help='Download the latest experimental data files for use in analysis.  Will run automatically if the "data" directory does not already exist (e.g. first usage).')

    parser.add_argument('--nofig', '-no', action="store_true", help="Don't output a figure.")

    parser.add_argument('--cache_intermediates', '-ci', action="store_true",
                        help="Write cluster_dict.json (a large internal cache of per-TF hits with CAS + weights) alongside the user-facing TFBSs_found.sortedclusters.csv. Disabled by default because this file can be 100s of MB per transcript at pvalc=1 and is only needed for re-running without recomputing. The CSV output is unaffected.")

    parser.add_argument('--output_format', '-of', choices=['csv', 'parquet', 'both'], default='csv',
                        help="Format for the sorted-clusters output table. 'csv' (default) is human-readable. 'parquet' is ~10x faster to write and ~10x smaller on disk but requires the optional pyarrow dep (install with: pip install 'tfbs_footprinting3[parquet]'); recommended for HPC/batch runs. 'both' writes both files.")

    # pre-processing the arguments
    args = parser.parse_args()
    args_lists = []
    transcript_ids_filename = args.t_ids_file
    exp_data_update = args.exp_data_update
    nofigure = args.nofig
    cache_intermediates = args.cache_intermediates
    output_format = args.output_format

    if transcript_ids_filename:
        filename, file_extension = os.path.splitext(transcript_ids_filename)

        if file_extension == ".csv" or file_extension == ".tsv":

            if file_extension == ".csv":
                parsed_arg_lines = file_to_datalist(transcript_ids_filename, delimiter=",")[1:]
            elif file_extension == ".tsv":
                parsed_arg_lines = file_to_datalist(transcript_ids_filename, delimiter="\t")[1:]

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
            # a list which has all of the parameters set as the same, in this way there can be a similar input format
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

    return args_lists, exp_data_update, nofigure, cache_intermediates, output_format


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
        args_lists, exp_data_update, nofigure, cache_intermediates, output_format = get_args()

        # if experimental data dir does not exist or user has requested an exp data update, then update.
        experimental_data_present = experimentalDataUpdater(exp_data_update)

        if experimental_data_present:
            if len(args_lists) > 0:
                # load mono-nuc PFMs
                TFBS_matrix_filename = os.path.join(script_dir, 'data/JASPAR_2026_pwms.json')
                TFBS_matrix_dict = load_json(TFBS_matrix_filename)
                TFBS_matrix_dict = {k.upper(): v for k, v in TFBS_matrix_dict.items()}

                # load species nt frequencies
                species_nt_freq_fn = os.path.join(script_dir, 'data/species_nt_freq.json')
                species_nt_freq_d = load_json(species_nt_freq_fn)

                # Group transcripts by (species, chromosome) so the
                # species_specific_data cache below fires maximally.
                # Amortizes one Ensembl /lookup/id call per transcript up
                # front in exchange for potentially dozens of skipped
                # species_specific_data reloads.
                args_lists = _prefetch_and_sort_transcripts(args_lists, output_dir)

                last_target_species = None
                last_chromosome = None

            for i, args_list in enumerate(args_lists):
                args, transcript_ids_filename, transcript_id, target_tfs_filename, promoter_before_tss, promoter_after_tss, top_x_tfs_count, pval, pvalc = args_list
                print("Ensembl transcript id:", transcript_id)
                logging.info("\n" + " ".join(["***ANALYSIS OF A NEW TRANSCRIPT HAS BEGUN:", transcript_id]))
                logging.info(" ".join(["Arguments used in this run:", str(args_list)]))

                # target dir naming
                start_end = "(" + "_".join([str(promoter_before_tss), str(promoter_after_tss)]) + ")"
                target_dir_name = "_".join([transcript_id + start_end, str(pval)])
                target_dir = os.path.join(output_dir, target_dir_name)

                # declare all possible results filenames.
                cluster_dict_filename = os.path.join(target_dir, "cluster_dict.json")
                ensembl_aligned_filename = os.path.join(target_dir, "alignment_uncleaned.fasta")
                cleaned_aligned_filename = os.path.join(target_dir, "alignment_cleaned.fasta")
                transcript_dict_filename = os.path.join(target_dir, "transcript_dict.json")
                gene_dict_filename = os.path.join(target_dir, "gene_dict.json")
                regulatory_decoded_filename = os.path.join(target_dir, "regulatory_decoded.json")
                sortedclusters_csv_filename = os.path.join(target_dir, "TFBSs_found.sortedclusters.csv")
                sortedclusters_parquet_filename = os.path.join(target_dir, "TFBSs_found.sortedclusters.parquet")
                # The resume check keys off whichever format(s) were requested.
                if output_format == "parquet":
                    sortedclusters_table_filename = sortedclusters_parquet_filename
                else:
                    sortedclusters_table_filename = sortedclusters_csv_filename

                # check if results have been created for this query.
                # cluster_dict.json is only required when -ci/--cache_intermediates is on;
                # by default the resume check keys off the user-facing CSV + the small
                # Ensembl/alignment caches (transcript_dict / gene_dict / regulatory_decoded / alignment.fasta).
                required_results_filenames = [ensembl_aligned_filename, cleaned_aligned_filename, transcript_dict_filename, gene_dict_filename, regulatory_decoded_filename, sortedclusters_table_filename]
                if cache_intermediates:
                    required_results_filenames.append(cluster_dict_filename)
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
                                species_pwm_score_threshold_df, gerp_conservation_locations_dict, species_group, cage_dict, TF_cage_dict, cage_dist_weights_dict, cage_correlations_dict, cage_corr_weights_dict, metacluster_overlap_weights_dict, cpg_obsexp_weights_dict, cpg_obsexp_weights_dict_keys, gtex_variants, gtex_weights_dict, gtrd_metaclusters_dict, atac_seq_dict, cas_pvalues_dict = species_specific_data(target_species, chromosome, species_specific_data_dir)
                            last_target_species = target_species
                            last_chromosome = chromosome

                        # format - target species tf:pvalue:score
                        pwm_score_threshold_dict = {}
                        if len(species_pwm_score_threshold_df) > 0:
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
                            alignment_len = len(target_species_row['seq'].replace('-', ''))

                            # retrieve regulatory
                            regulatory_decoded = retrieve_regulatory(chromosome, strand, promoter_start, promoter_end, regulatory_decoded_filename, target_species)
                            converted_reg_dict = reg_position_translate(tss, regulatory_decoded, promoter_start, promoter_end, strand, promoter_before_tss, promoter_after_tss)

                            # conservation
                            converted_gerps_in_promoter = gerp_positions_translate(target_dir, gerp_conservation_locations_dict, chromosome, strand, promoter_start, promoter_end, tss)

                            # identify information content of each column of the alignment
                            cpg_list = CpG(cleaned_aligned_filename)

                            # identify CAGEs in proximity to Ensembl TSS, convert for plotting
                            converted_cages = cage_position_translate(target_species, gene_name, ens_gene_id, transcript_id, tss, cage_dict, promoter_start, promoter_end, strand, promoter_before_tss, promoter_after_tss)

                            # identify eQTLs in proximity to Ensembl TSS, convert for plotting
                            converted_eqtls = gtex_position_translate(ens_gene_id, gtex_variants, tss, promoter_start, promoter_end, strand, promoter_before_tss, promoter_after_tss)

                            # GTRD metaclusters
                            converted_metaclusters_in_promoter, metacluster_in_promoter_counts = gtrd_positions_translate(target_dir, gtrd_metaclusters_dict, chromosome, strand, promoter_start, promoter_end, tss)

                            # ATAC-seq data
                            converted_atac_seqs_in_promoter = atac_pos_translate(atac_seq_dict, chromosome, strand, promoter_start, promoter_end, tss)

                            # create index of aligned to unaligned positions
                            unaligned2aligned_index_dict = unaligned2aligned_indexes(cleaned_aligned_filename)

                            # Decide whether to (re)compute or load an existing cluster_dict.
                            # With caching on: load from cluster_dict.json when it exists.
                            # Without caching: always recompute (the file won't be there).
                            can_load_cluster = cache_intermediates and os.path.exists(cluster_dict_filename) and os.path.exists(sortedclusters_table_filename)
                            if not can_load_cluster:
                                # score alignment for tfbss
                                tfbss_found_dict = tfbs_finder(transcript_name, alignment, target_tfs_list, TFBS_matrix_dict, target_dir, pwm_score_threshold_dict, species_nt_freq_d, unaligned2aligned_index_dict, promoter_after_tss, pval, pvalc)

                                # sort through scores, identify hits in target_species supported in other species
                                cluster_dict = find_clusters(gene_name, ens_gene_id, chr_start, chr_end, alignment, target_species, chromosome, tfbss_found_dict, cleaned_aligned_filename, converted_gerps_in_promoter, converted_cages, converted_metaclusters_in_promoter, metacluster_in_promoter_counts, converted_atac_seqs_in_promoter, converted_eqtls, gtex_weights_dict, transcript_id, cage_dict, TF_cage_dict, cage_dist_weights_dict, metacluster_overlap_weights_dict, cpg_list, cpg_obsexp_weights_dict, cpg_obsexp_weights_dict_keys, cage_correlations_dict, cage_corr_weights_dict, gtex_variants, gene_len, cas_pvalues_dict, pvalc)
                                tfbss_found_dict.clear()
                                if cache_intermediates:
                                    dump_json(cluster_dict_filename, cluster_dict)

                            else:
                                cluster_dict = load_json(cluster_dict_filename)

                            # sort the target_species hits supported by other species
                            sorted_clusters_target_species_hits_list = sort_target_species_hits(cluster_dict)

                            if output_format in ("csv", "both"):
                                target_species_hits_table_writer(sorted_clusters_target_species_hits_list, sortedclusters_csv_filename)
                            if output_format in ("parquet", "both"):
                                # Parquet path also applies sci-pval formatting in-place; if we
                                # already ran the CSV writer on this same list, the p-value columns
                                # are already formatted (sci_cache is a no-op the 2nd time). Safe.
                                target_species_hits_table_writer_parquet(sorted_clusters_target_species_hits_list, sortedclusters_parquet_filename)

                            if len(sorted_clusters_target_species_hits_list) > 0:
                                # extract the top x target_species hits supported by other species
                                top_x_greatest_hits_dict = top_x_greatest_hits(sorted_clusters_target_species_hits_list, top_x_tfs_count)
                                top_greatest_hits_dict = top_greatest_hits(sorted_clusters_target_species_hits_list, top_x_tfs_count)

                                # plot the top x target_species hits
                                if not nofigure:
                                    if len(top_x_greatest_hits_dict) > 0:
                                        plot_promoter(target_species, transcript_id, species_group, alignment, alignment_len, promoter_before_tss, promoter_after_tss, transcript_name, top_x_greatest_hits_dict, target_dir, converted_reg_dict, converted_gerps_in_promoter, cpg_list, converted_cages, converted_metaclusters_in_promoter, converted_atac_seqs_in_promoter, converted_eqtls)

            total_time_end = time.time()
            logging.info(" ".join(["Total time for", str(len(args_lists)), "transcripts:", str(total_time_end - total_time_start), "seconds"]) + "\n\n")

    else:
        print("System does not appear to be connected to the internet.  Exiting TFBS_footprinter3.")
