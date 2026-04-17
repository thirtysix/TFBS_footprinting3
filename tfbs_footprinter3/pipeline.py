"""Pipeline stages: transcript metadata retrieval + alignment + TFBS scanning.

Each transcript in a tfbs_footprinter3 run flows through these functions
in order, orchestrated by `main()` in the monolith:

    1. transfabulator         - fetch transcript metadata from Ensembl
    2. test_transcript_id     - validate the response
    3. transcript_data_retrieve - parse position + compute promoter window
    4. gene_data_retrieve     - fetch parent-gene metadata
    5. retrieve_regulatory    - fetch regulatory annotations
    6. alignment_tools        - fetch + clean the genome alignment
    7. tfbs_finder            - PWM-scan the cleaned sequence for TFBS hits

Callers reach Ensembl through `ensembl.ensemblrest(...)` via ATTRIBUTE
ACCESS at call time rather than a bare name imported at module load.
That matters for the HPC cache layer: `hpc/ensembl_cache.py::patch_tfbs_footprinter3()`
rebinds `tfbs_footprinter3.ensembl.ensemblrest`, and attribute lookup
picks up the patched version transparently. A bare
`from tfbs_footprinter3.ensembl import ensemblrest` at the top of this
module would bind the original function here at import time and miss
subsequent patches.
"""
from __future__ import annotations

import logging
import os
import time
from bisect import bisect_left
from operator import itemgetter

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from tfbs_footprinter3 import ensembl
from tfbs_footprinter3.alignment import (
    fasta_writer,
    load_genome_aligned,
    remove_duplicate_species,
    remove_gap_only,
    remove_non_ACGT,
)
from tfbs_footprinter3.io_utils import dump_json, load_json
from tfbs_footprinter3.pwm import PWM_scorer, pwm_maker
from tfbs_footprinter3.translators import start_end_found_motif


def tfbs_finder(transcript_name, alignment, target_tfs_list, TFBS_matrix_dict, target_dir, pwm_score_threshold_dict, species_nt_freq_d, unaligned2aligned_index_dict, promoter_after_tss, pval, pvalc):
    """
    1. Convert PFM to PWM for each TF in the Jaspar dictionary.
    2. Score all positions in the cleaned sequence.
    3. If score is greater than or equal to precomputed threshold, then keep, otherwise set to zero.
    4. Return dictionary of pwm scores for [species][tf_name][strand].
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
        mononuc_pwm_dict = {"A": 0, "C": 1, "G": 2, "T": 3}

        entry = alignment[0]
        species = entry['species']

        # Remove alignment/ambiguous characters from the sequences
        cleaned_seq = entry['seq']
        for char in align_chars:
            cleaned_seq = cleaned_seq.replace(char, "")
        entry_seqrecord = SeqRecord(Seq(cleaned_seq), id=species)
        forward_seq = str(entry_seqrecord.seq)
        reverse_seq = str(entry_seqrecord.seq.reverse_complement())
        seq_dict = {"+1": forward_seq, "-1": reverse_seq}

        # generate background frequencies of each mono-nucleotide for forward and reverse strands
        bg_nuc_freq_dict = {}

        # bool - check if species is human; use previously calculated nt freqs
        if species == "homo_sapiens":
            # https://arxiv.org/pdf/q-bio/0611041.pdf — empirical data from complete genome
            bg_nuc_freq_dict = {'A': 0.292, 'C': 0.207, 'G': 0.207, 'T': 0.292}

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
                    pvals_scores_list = [[k, v] for k, v in tf_pwm_score_threshold_dict.items()]
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

                        seq_length = len(seq)
                        # iterate through the nt sequence, extract a current frame based on the motif size and score
                        for i in range(0, seq_length - motif_length):
                            current_frame = seq[i:i + motif_length]
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
    target_only_decoded = ensembl.ensemblrest(query_type, options, 'json', "", log=True)
    if 'seq' in target_only_decoded:
        target_only_decoded['species'] = target_species
        alignment = [target_only_decoded]
    else:
        alignment = []

    return alignment


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

    # retrieve transcript position data from json file if it does not exist or contains no data.
    if retrieve_transcript_data:
        # Set parameters for retrieving Ensembl data via REST
        query_type = '/lookup/id/'
        options = '?feature=transcript;content-type=application/json'

        # populate 'transcript_dict' dictionary with sub-dictionaries.
        # key[transcript_id] = {chromosome, strand, start, end} for each ensembl transcript id
        decoded_json_description = ensembl.ensemblrest(query_type, options, 'json', transcript, log=True)
        decoded_json_description = {k.lower(): v for k, v in decoded_json_description.items()}

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

    if not decoded_json_description or len(decoded_json_description) == 0:
        query_type = '/lookup/id/'
        options = '?feature=transcript;content-type=application/json'

        # retrieve_gene_len
        decoded_json_description = ensembl.ensemblrest(query_type, options, 'json', ens_gene_id, log=True)
        decoded_json_description = {k.lower(): v for k, v in decoded_json_description.items()}

        dump_json(gene_dict_filename, decoded_json_description)

    gene_start = decoded_json_description['start']
    gene_end = decoded_json_description['end']

    if 'display_name' in decoded_json_description:
        gene_name = decoded_json_description['display_name'].upper()
    else:
        gene_name = decoded_json_description['id']

    gene_len = gene_end - gene_start

    return gene_name, gene_len


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
        regulatory_decoded = ensembl.ensemblrest(query_type, options, 'json', "", log=True)

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
