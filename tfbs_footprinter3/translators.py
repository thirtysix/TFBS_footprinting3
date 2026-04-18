"""Coordinate / position translators between genome space and plot space.

These convert genomic coordinates of regulatory features (CAGE, GTEX eQTLs,
GERP conservation, GTRD metaclusters, ATAC-Seq peaks, Ensembl regulatory
regions, CpG windows) into TSS-relative coordinates that the plotter and
the scoring code consume.

Extracted from tfbs_footprinter3.py. All functions here are pure
numerical transformations; no network or filesystem I/O beyond
`unaligned2aligned_indexes` (reads a cleaned FASTA) and `CpG` (reads
an alignment via Biopython).
"""
from __future__ import annotations

from operator import itemgetter

from Bio import AlignIO, SeqIO

from tfbs_footprinter3.io_utils import overlap_range


def start_end_found_motif(i, strand, seq_length, promoter_after_tss, motif_length):
    """
    Determine the start/end positions of the found motif.
    """
    if strand == "+1":
        hit_loc_start = i
        hit_loc_before_TSS_start = i - seq_length + promoter_after_tss
        hit_loc_end = i + motif_length
        hit_loc_before_TSS_end = i - seq_length + motif_length + promoter_after_tss
    if strand == "-1":
        hit_loc_start = seq_length - i - motif_length
        hit_loc_before_TSS_start = (seq_length - i - motif_length) - seq_length + promoter_after_tss
        hit_loc_end = seq_length - i
        hit_loc_before_TSS_end = (seq_length - i) - seq_length + promoter_after_tss

    return hit_loc_start, hit_loc_end, hit_loc_before_TSS_start, hit_loc_before_TSS_end


def unaligned2aligned_indexes(cleaned_aligned_filename):
    """
    Create a dictionary for mapping aligned positions to unaligned positions.
    """

    with open(cleaned_aligned_filename) as cleaned_aligned_file:
        aligned_entries_dict = SeqIO.to_dict(SeqIO.parse(cleaned_aligned_file, 'fasta'))

    unaligned2aligned_index_dict = {}

    for species, seqRecord in aligned_entries_dict.items():
        unaligned2aligned_index_dict[species] = {}
        seq_str = str(seqRecord.seq)
        for aligned_index in range(len(seq_str)):
            if seq_str[aligned_index] != "-":
                unaligned_index = aligned_index - seq_str[:aligned_index].count("-")
                unaligned2aligned_index_dict[species][unaligned_index] = aligned_index

    return unaligned2aligned_index_dict


def reg_position_translate(tss, regulatory_decoded, promoter_start, promoter_end, strand, promoter_before_tss, promoter_after_tss):
    """
    Convert positions of regulatory elements (json) to coordinates usable by the plot_promoter function.
    Requires supplying location relative to TSS (i.e. negative).
    For creating a left to right plot of the promoter, regardless of strand:
    Converted_reg_start is the leftmost regulatory position.
    Converted_reg_end is the rightmost regulatory position.
    """

    converted_reg_dict = {}
    for reg in regulatory_decoded:
        if "error" not in reg:
            reg_id = reg['id']
            reg_start = reg['start']
            reg_end = reg['end']
            description = reg['description']

            if strand == 1:
                #[promoter_start][reg_start][reg_end][promoter_end][chr_start][TSS>GENE--->][chr_end]
                converted_reg_start = (tss - reg_start) * -1
                converted_reg_end = (tss - reg_end) * -1
                if reg_start <= promoter_start:
                    converted_reg_start = (-1 * promoter_before_tss)
                if reg_end >= promoter_end:
                    converted_reg_end = promoter_after_tss - 0.001

            if strand == -1:
                #[chr_start][<---GENE<TSS][chr_end][promoter_start][reg_start][reg_end][promoter_end]
                converted_reg_start = (tss - reg_start)
                converted_reg_end = (tss - reg_end)

                if reg_start <= promoter_start:
                    converted_reg_start = promoter_after_tss - 0.001
                if reg_end >= promoter_end:
                    converted_reg_end = (-1 * promoter_before_tss + promoter_after_tss) + 0.001

            converted_reg_dict[reg_id] = {'converted_start': converted_reg_start, 'converted_end': converted_reg_end, 'description': description}

    return converted_reg_dict


def CpG(aligned_filename):
    """
    Score the CpG content of the target_species sequence over a 200 nt window.
    """

    alignment = AlignIO.read(aligned_filename, "fasta")
    target_species_row = alignment[0]
    cpg_list = []

    # [1 C, 1 if G, 1 if CPG, CorG, num_cpg, obs2exp]
    for i in range(0, len(target_species_row)):
        current_pos = target_species_row[i]
        if current_pos != '-':
            if i < len(target_species_row) - 1:
                next_pos = target_species_row[i + 1]
            else:
                next_pos = False

            if current_pos == 'C' and next_pos == 'G':
                cpg_list.append([1, 0, 1])

            elif current_pos == 'C' and next_pos != 'G':
                cpg_list.append([1, 0, 0])

            elif current_pos == 'G':
                cpg_list.append([0, 1, 0])

            else:
                cpg_list.append([0, 0, 0])

    for i in range(0, len(cpg_list)):
        if i < 100:
            rolling_island = cpg_list[:i] + cpg_list[i:i + 100]
        elif i > len(cpg_list) - 100:
            rolling_island = cpg_list[i - 100:i] + cpg_list[i:]
        else:
            rolling_island = cpg_list[i - 100:i + 100]

        Cs = sum(x[0] for x in rolling_island)
        Gs = sum(x[1] for x in rolling_island)
        CorG_ratio = (Cs + Gs) / len(rolling_island)
        num_cpg = sum(x[2] for x in rolling_island)

        obs = num_cpg / len(rolling_island)
        exp = (CorG_ratio / 2) ** 2
        if exp == 0:
            exp = 0.0000000001
        obs2exp = obs / exp
        cpg_list[i] = cpg_list[i] + [CorG_ratio, num_cpg, obs2exp]

    return cpg_list


def cage_position_translate(target_species, gene_name, ens_gene_id, transcript_id, tss, cage_dict, promoter_start, promoter_end, strand, promoter_before_tss, promoter_after_tss):
    """
    Convert the CAGE data genome positions to those which can be mapped into the final figure.
    """

    converted_cages = []

    if target_species == "homo_sapiens":
        cage_key = gene_name
    else:
        cage_key = ens_gene_id

    if cage_key in cage_dict:
        cages = cage_dict[cage_key]
        cages_peak_count_sum = float(sum(int(x[5]) for x in cages))

        for cage in cages:
            #ref-point
            cage_desc = cage[1]
            cage_start = int(cage[3])
            cage_end = cage_start + int(cage[4])
            cage_peak_count = int(cage[5])
            cage_peak_count_ratio = cage_peak_count / cages_peak_count_sum
            cage_strand = cage[6]

            if cage_strand == "+":
                #[promoter_start][cage_start][cage_end][promoter_end][chr_start][TSS>GENE--->][chr_end]
                converted_cage_start = (tss - cage_start) * -1
                converted_cage_end = (tss - cage_end) * -1

            if cage_strand == "-":
                #[chr_start][<---GENE<TSS][chr_end][promoter_start][cage_start][cage_end][promoter_end]
                converted_cage_start = (tss - cage_start)
                converted_cage_end = (tss - cage_end)

            converted_cage = [converted_cage_start, converted_cage_end, cage_desc, cage_peak_count_ratio]
            converted_cages.append(converted_cage)

        converted_cages = sorted(converted_cages, key=itemgetter(2))

    return converted_cages


def gtex_position_translate(ens_gene_id, gtex_variants, tss, promoter_start, promoter_end, strand, promoter_before_tss, promoter_after_tss):
    """
    Convert the GTEx data genome positions to those which can be mapped into the final figure.
    Reduce to those that are within range of the promoter before/after tss.
    """

    converted_eqtls = []
    if ens_gene_id in gtex_variants:
        eqtls = gtex_variants[ens_gene_id]

        for eqtl in eqtls:
            if len(eqtl) == 2:
                loc = eqtl[0]
                eqtl_length = 1
                eqtl_effect = eqtl[1]
            else:
                loc = eqtl[0]
                eqtl_length = eqtl[1]
                eqtl_effect = eqtl[2]

            if promoter_start <= loc <= promoter_end or promoter_start <= loc + eqtl_length <= promoter_end:
                #[promoter_start][eqtl_start][eqtl_end][promoter_end][chr_start][TSS>GENE--->][chr_end]
                if strand == 1:
                    converted_eqtl_start = (tss - loc) * -1
                    converted_eqtl_end = (tss - loc + eqtl_length) * -1

                #[chr_start][<---GENE<TSS][chr_end][promoter_start][eqtl_start][eqtl_end][promoter_end]
                if strand == -1:
                    converted_eqtl_start = (tss - loc)
                    converted_eqtl_end = (tss - loc + eqtl_length)

                # save to final list
                converted_eqtl = [converted_eqtl_start, converted_eqtl_end, eqtl_effect]
                converted_eqtls.append(converted_eqtl)

    return converted_eqtls


def gerp_positions_translate(target_dir, gerp_conservation_locations_dict, chromosome, strand, promoter_start, promoter_end, tss):
    """
    Identify GERP constrained conservation locations which occur within the defined promoter region.
    Convert positions of GERP elements (json) to coordinates usable by the plot_promoter function.
    Requires supplying location relative to TSS (i.e. negative).
    For creating a left to right plot of the promoter, regardless of strand:
    Converted_reg_start is the leftmost regulatory position.
    Converted_reg_end is the rightmost regulatory position.
    """

    gerps_in_promoter = []

    # because a prediction can occur at the start/end of a defined promoter
    extended_range = 1000

    if chromosome in gerp_conservation_locations_dict:
        for potential_gerp_in_promoter in gerp_conservation_locations_dict[chromosome]:
            if promoter_start - extended_range <= potential_gerp_in_promoter[0] <= promoter_end + extended_range or promoter_start - extended_range <= potential_gerp_in_promoter[0] + potential_gerp_in_promoter[1] <= promoter_end + extended_range:
                gerps_in_promoter.append(potential_gerp_in_promoter)

    # convert the positions of the in-promoter metaclusters to tss-relative
    converted_gerps_in_promoter = []
    for gerp_in_promoter in gerps_in_promoter:
        gerp_start = gerp_in_promoter[0]
        gerp_end = gerp_start + gerp_in_promoter[1]
        gerp_weight = gerp_in_promoter[2]

        if strand == 1:
            converted_gerp_start = (tss - gerp_start) * -1
            converted_gerp_end = (tss - gerp_end) * -1
        if strand == -1:
            # gerp_end = gerp_start + length (genomic), so gerp_end > gerp_start.
            # For reverse-strand genes we need to swap so the converted range
            # still has start <= end; otherwise downstream overlap checks
            # (which test max(point, start) <= min(point+1, end)) always fail.
            converted_gerp_start = (tss - gerp_end)
            converted_gerp_end = (tss - gerp_start)

        converted_gerp = [converted_gerp_start, converted_gerp_end, gerp_weight]
        converted_gerps_in_promoter.append(converted_gerp)

    return converted_gerps_in_promoter


def gtrd_positions_translate(target_dir, gtrd_metaclusters_dict, chromosome, strand, promoter_start, promoter_end, tss):
    """
    Identify GTRD metaclusters which occur within the defined promoter region.
    Convert positions of metaclusters (json) to coordinates usable by the plot_promoter function.
    Requires supplying location relative to TSS (i.e. negative).
    For creating a left to right plot of the promoter, regardless of strand:
    Converted_reg_start is the leftmost regulatory position.
    Converted_reg_end is the rightmost regulatory position.
    """

    potential_metaclusters_in_promoter = []
    promoter_start_millions = int(promoter_start / 1000000)
    promoter_end_millions = int(promoter_end / 1000000)

    # retrieve the metacluster peaks on which the chrom that the transcript is found
    # if the millions place is the same for each then the metaclusters come from a single
    # subdict entry
    if promoter_start_millions == promoter_end_millions:
        if promoter_start_millions in gtrd_metaclusters_dict:
            potential_metaclusters_in_promoter += gtrd_metaclusters_dict[promoter_start_millions]

    # have to account for the possibility that this location spans a millions place
    # e.g. from 999,000 - 1,001,000
    else:
        if promoter_start_millions in gtrd_metaclusters_dict:
            potential_metaclusters_in_promoter += gtrd_metaclusters_dict[promoter_start_millions]

        if promoter_end_millions in gtrd_metaclusters_dict:
            potential_metaclusters_in_promoter += gtrd_metaclusters_dict[promoter_end_millions]

    # identify if the metacluster occurs within user-defined promoter region
    metaclusters_in_promoter = []
    for potential_metacluster in potential_metaclusters_in_promoter:
        metacluster_start = potential_metacluster[0]
        metacluster_end = metacluster_start + potential_metacluster[1]

        overlap = overlap_range([promoter_start, promoter_end], [metacluster_start, metacluster_end])

        if len(overlap) > 0:
            metaclusters_in_promoter.append(potential_metacluster)

    # convert the positions of the in-promoter metaclusters to tss-relative
    converted_metaclusters_in_promoter = []

    for metacluster_in_promoter in metaclusters_in_promoter:
        metacluster_start = metacluster_in_promoter[0]
        metacluster_end = metacluster_start + metacluster_in_promoter[1]
        metacluster_peak_count = metacluster_in_promoter[2]

        if strand == 1:
            converted_metacluster_start = (tss - metacluster_start) * -1
            converted_metacluster_end = (tss - metacluster_end) * -1
        if strand == -1:
            converted_metacluster_start = (tss - metacluster_start)
            converted_metacluster_end = (tss - metacluster_end)

        converted_metacluster = [converted_metacluster_start, converted_metacluster_end, metacluster_peak_count]
        converted_metaclusters_in_promoter.append(converted_metacluster)

    # make - count list
    metacluster_in_promoter_counts = [0] * (promoter_end - promoter_start)

    # iterate - intervals
    for metacluster_in_promoter in metaclusters_in_promoter:
        metacluster_relative_start = metacluster_in_promoter[0] - promoter_start
        metacluster_relative_end = metacluster_in_promoter[0] + metacluster_in_promoter[1] - promoter_start

        # iterate - all positions in current interval
        for pos in range(metacluster_relative_start, metacluster_relative_end - 1):
            if pos >= 0 and pos < len(metacluster_in_promoter_counts):
                metacluster_in_promoter_counts[pos] += 1

    if strand == -1:
        metacluster_in_promoter_counts = metacluster_in_promoter_counts[::-1]

    return converted_metaclusters_in_promoter, metacluster_in_promoter_counts


def atac_pos_translate(atac_seq_dict, chromosome, strand, promoter_start, promoter_end, tss):
    """
    Identify merged ATAC-Seq peaks which occur within the defined promoter region.
    Convert positions of ATAC-Seq peaks (json) to coordinates usable by the plot_promoter function.
    Requires supplying location relative to TSS (i.e. negative).
    For creating a left to right plot of the promoter, regardless of strand:
    Converted_reg_start is the leftmost regulatory position.
    Converted_reg_end is the rightmost regulatory position.
    """

    potential_atac_seqs_in_promoter = []

    promoter_start_millions = int(promoter_start / 1000000)
    promoter_end_millions = int(promoter_end / 1000000)

    # retrieve the ATAC-Seq peaks on which the chrom that the transcript is found
    # if the millions place is the same for each then the atac-seqs come from a single subdict entry
    if promoter_start_millions == promoter_end_millions:
        if promoter_start_millions in atac_seq_dict:
            potential_atac_seqs_in_promoter += atac_seq_dict[promoter_start_millions]

    # have to account for the possibility that this location spans a millions place
    # e.g. from 999,000 - 1,001,000
    else:
        if promoter_start_millions in atac_seq_dict:
            potential_atac_seqs_in_promoter += atac_seq_dict[promoter_start_millions]

        if promoter_end_millions in atac_seq_dict:
            potential_atac_seqs_in_promoter += atac_seq_dict[promoter_end_millions]

    # identify if the ATAC-Seq peak occurs within user-defined promoter region
    atac_seqs_in_promoter = []
    for potential_atac_seq in potential_atac_seqs_in_promoter:
        atac_seq_start = potential_atac_seq[0]
        atac_seq_end = potential_atac_seq[0] + potential_atac_seq[1]
        atac_seq_weight = potential_atac_seq[2]

        overlap = overlap_range([promoter_start, promoter_end], [atac_seq_start, atac_seq_end])

        if len(overlap) > 0:
            atac_seqs_in_promoter.append([atac_seq_start, atac_seq_end, atac_seq_weight])

    # convert the positions of the in-promoter atac_seqs to tss-relative
    converted_atac_seqs_in_promoter = []
    for atac_seq_in_promoter in atac_seqs_in_promoter:
        atac_seq_start = atac_seq_in_promoter[0]
        atac_seq_end = atac_seq_in_promoter[1]
        atac_seq_weight = atac_seq_in_promoter[2]

        if strand == 1:
            converted_atac_seq_start = (tss - atac_seq_start) * -1
            converted_atac_seq_end = (tss - atac_seq_end) * -1
        if strand == -1:
            # Swap so start <= end after conversion; see the matching comment
            # in gerp_positions_translate above.
            converted_atac_seq_start = (tss - atac_seq_end)
            converted_atac_seq_end = (tss - atac_seq_start)

        converted_atac_seq = [converted_atac_seq_start, converted_atac_seq_end, atac_seq_weight]
        converted_atac_seqs_in_promoter.append(converted_atac_seq)

    return converted_atac_seqs_in_promoter
