"""Sequence alignment processing helpers.

Pure alignment transformations extracted from tfbs_footprinter3.py
(previously lines 1642-1737). These do not hit the network; the
Ensembl-facing `retrieve_genome_aligned` and the orchestrator
`alignment_tools` remain in the monolith until the ensembl.py
extraction lands, at which point they can move here too.
"""
from __future__ import annotations

import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def fasta_writer(alignment, outfile):
    """
    Write ensembl JSON alignment to fasta file.
    """

    if not os.path.isfile(outfile) or (os.path.isfile(outfile) and os.path.getsize(outfile) == 0):
        with open(outfile, "w") as aligned_file:
            for entry in alignment:
                record = SeqRecord(Seq(entry['seq']), id=entry['species'], description="")
                SeqIO.write(record, aligned_file, 'fasta')


def remove_non_ACGT(alignment):
    """
    Remove non alignment characters and ambiguous nucleotides.  should consider changing to replacing any non ACGT char to '-'.
    """

    # account for sequences which are non-standard code
    non_standard_dict = {'R': ['A', 'G'],
                         'Y': ['C', 'T'],
                         'S': ['G', 'C'],
                         'W': ['A', 'T'],
                         'K': ['G', 'T'],
                         'M': ['A', 'C'],
                         'B': ['C', 'G', 'T'],
                         'D': ['A', 'G', 'T'],
                         'H': ['A', 'C', 'T'],
                         'V': ['A', 'C', 'G']}

    non_alignment_chars = " .N"
    for entry in alignment:
        for non_alignment_char in non_alignment_chars:
            entry['seq'] = entry['seq'].replace(non_alignment_char, '-')

        for multi_char, replacement_list in non_standard_dict.items():
            entry['seq'] = entry['seq'].replace(non_alignment_char, replacement_list[0])

    return alignment


def remove_gap_only(alignment):
    """
    Find columns in the alignment where the entire column is '-',
        replace the '-' with 'P', then remove the '*'.
    """

    if len(alignment) > 0:
        for entry in alignment:
            entry['seq'] = list(entry['seq'])

        for i in range(0, len(alignment[0]['seq'])):
            col = [x['seq'][i] for x in alignment]
            if col.count('-') == len(col):
                for entry in alignment:
                    entry['seq'][i] = 'Z'
        for entry in alignment:
            entry['seq'] = "".join(entry['seq']).replace('Z', "")

    return alignment


def remove_duplicate_species(alignment, target_species):
    """
    If there are multiple entries for a single species in an alignment retrieved from Ensembl,
    keep the one which has more ACGT characters.
    """

    entry_ids = [x['species'] for x in alignment]
    duplicate_ids = list({x for x in entry_ids if entry_ids.count(x) > 1})
    non_duplicate_alignment = [x for x in alignment if x['species'] not in duplicate_ids]
    for duplicate_id in duplicate_ids:
        duplicate_seqs = [x for x in alignment if x['species'] == duplicate_id]
        duplicate_seqs_lens = [x['seq'].count('-') for x in duplicate_seqs]
        sorted_duplicate_seqs_lens = duplicate_seqs_lens[:]
        sorted_duplicate_seqs_lens.sort()
        longest_seq = sorted_duplicate_seqs_lens[0]
        longest_seq_index = duplicate_seqs_lens.index(longest_seq)
        kept_seq = duplicate_seqs[longest_seq_index]
        if duplicate_id == target_species:
            non_duplicate_alignment = [kept_seq] + non_duplicate_alignment
        else:
            non_duplicate_alignment.append(kept_seq)

    return non_duplicate_alignment


def load_genome_aligned(aligned_filename):
    """
    Load previously retrieved alignment fasta file into dictionary.
    """

    with open(aligned_filename) as alignment_handle:
        alignment_list = list(SeqIO.parse(alignment_handle, 'fasta'))
    alignment = [{'seq': str(entry.seq), 'species': entry.id} for entry in alignment_list if '[' not in entry.id]

    return alignment
