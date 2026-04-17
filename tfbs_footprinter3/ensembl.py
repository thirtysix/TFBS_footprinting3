"""Ensembl REST client.

Thin wrapper around httplib2 with per-call retry and x-ratelimit-aware
sleeping. Used by transcript metadata retrieval, promoter sequence
retrieval, and regulatory overlap queries.

Monkey-patch target for the HPC cache layer:
`hpc/ensembl_cache.py::patch_tfbs_footprinter3()` assigns
`tfbs_footprinter3.tfbs_footprinter3.ensemblrest = <cached_version>`.
That assignment works because callers inside the monolith resolve
`ensemblrest` via the module's globals, and the re-export binding
exists there. If a future PR moves a caller OUT of the monolith into
another submodule, it must either import `ensemblrest` at call time
(via `from tfbs_footprinter3 import tfbs_footprinter3 as _tff;
_tff.ensemblrest(...)`) OR the HPC patch must be extended to also
rebind `tfbs_footprinter3.ensembl.ensemblrest` — pick one and document.

As of v0.0.7 the 'fasta' branch here is unreached; all four call sites
use output_type='json'. Retained verbatim during extraction so the
behavior is preserved — a follow-up can classify transient vs fatal
errors and drop the dead branch.
"""
from __future__ import annotations

import json
import logging
import time

import httplib2


def ensemblrest(query_type, options, output_type, ensembl_id=None, log=False):
    """
    Retrieve REST data from Ensembl using provided ID, query type, and options.
    """

    http = httplib2.Http()
    server = "http://rest.ensembl.org"
    full_query = server + query_type + ensembl_id + options

    if log:
        logging.info(" ".join(["Ensembl REST query made:", full_query]))

    success = False
    try_count = 0
    max_tries = 10
    fail_sleep_time = 120
    decoded_json = {}
    fasta_content = []

    if output_type == 'json':
        while not success and try_count < max_tries:
            try:
                resp, json_data = http.request(full_query, method="GET")
                decoded_json = json.loads(json_data)
                ensemblrest_rate(resp)
                try_count += 1
                success = True
                return decoded_json

            except:
                logging.info(" ".join(["Ensembl REST query unsuccessful, attempt:", "/".join([str(try_count), str(max_tries)]), "Sleeping:", str(fail_sleep_time), "seconds.", full_query]))
                try_count += 1
                print(" ".join(["Ensembl REST query unsuccessful, attempt:", "/".join([str(try_count), str(max_tries)]), "Sleeping:", str(fail_sleep_time), "seconds.", "See logfile for query."]))
                time.sleep(fail_sleep_time)

        # return empty decoded_json if max tries has elapsed
        return decoded_json

    if output_type == 'fasta':
        while not success and try_count < max_tries:
            try:
                resp, fasta_content = http.request(server, method="GET", headers={"Content-Type": "text/x-fasta"})
                ensemblrest_rate(resp)
                try_count += 1
                success = True
                return fasta_content

            except:
                logging.info(" ".join(["Ensembl REST query unsuccessful, attempt:", "/".join([str(try_count), str(max_tries)]), "Sleeping:", str(fail_sleep_time), "seconds.", full_query]))
                try_count += 1
                print(" ".join(["Ensembl REST query unsuccessful, attempt:", "/".join([str(try_count), str(max_tries)]), "Sleeping:", str(fail_sleep_time), "seconds.", "See logfile for query."]))
                time.sleep(fail_sleep_time)

        # return empty fasta_content if max tries has elapsed
        return fasta_content


def ensemblrest_rate(resp):
    """
    Read ensembl REST headers and determine if rate-limit has been exceeded, sleep appropriately if necessary.
    """

    if int(resp['x-ratelimit-remaining']) == 0:
        if 'Retry-After' in resp:
            sleep_time = int(resp['Retry-After'])
            logging.warning(" ".join(["Ensembl REST (Retry-After) requests sleeping for:", str(sleep_time)]))
            time.sleep(sleep_time)

        else:
            sleep_time = 60
            logging.warning(" ".join(["Ensembl REST requests sleeping for:", str(sleep_time)]))
            time.sleep(sleep_time)
