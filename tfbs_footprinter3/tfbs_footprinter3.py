"""Backward-compatibility shim for tfbs_footprinter3.

The 3,300-line monolith that once lived here has been split across
dedicated modules: cli, data_loader, data_parsing, ensembl, io_utils,
output, pipeline, plotting, pwm, scoring, translators, alignment.

This shim re-exports every public name from those modules so that
existing callers who do

    from tfbs_footprinter3.tfbs_footprinter3 import <anything>

keep working unchanged. The console-script entry point in
pyproject.toml also points here for the same reason:

    [project.scripts]
    tfbs_footprinter3 = "tfbs_footprinter3.tfbs_footprinter3:main"

New code should import directly from the specific module (e.g.
`from tfbs_footprinter3.scoring import find_clusters`) — this shim
adds no functionality, only compatibility.

The HPC cache layer (`hpc/ensembl_cache.py::patch_tfbs_footprinter3`)
rebinds `ensemblrest` on both this shim's namespace AND on
`tfbs_footprinter3.ensembl`, so pipeline callers that use attribute
access transparently pick up the cached wrapper.
"""
from __future__ import annotations

import os

from tfbs_footprinter3 import __version__  # noqa: F401

# Re-exports. Sorted by module name for readability.
from tfbs_footprinter3.alignment import (  # noqa: F401
    fasta_writer,
    load_genome_aligned,
    remove_duplicate_species,
    remove_gap_only,
    remove_non_ACGT,
)
from tfbs_footprinter3.cli import get_args, main  # noqa: F401
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

# Legacy module-level globals. A few callers (historically) reached
# into tff.script_dir / tff.curdir; preserved here for compatibility.
script_dir = os.path.dirname(__file__)
curdir = os.getcwd()
