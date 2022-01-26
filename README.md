
Full documentation now available at
[tfbs-footprinting.readthedocs.io](tfbs-footprinting.readthedocs.io)

# 1. Introduction

## 1.1 Basics

**Primary Usage:** Identification of cis-regulatory elements initially identified by matrix scoring and then additionally scored on 7 other relevant contextual datapoints. Based on analysis of protein-coding transcripts in the Ensembl database.

![logo](https://raw.githubusercontent.com/thirtysix/TFBS_footprinting/master/tfbs_logo.png)* This work is a derivative of ["Transcription factors"](https://commons.wikimedia.org/wiki/File:Transcription_Factors.svg) by [kelvin13](https://commons.wikimedia.org/wiki/User:Kelvin13), used under [CC BY 3.0](https://creativecommons.org/licenses/by/3.0/)
* * *


## 1.2 Experimental data

The TFBS footprinting method computationally predicts transcription factor binding sites (TFBSs) in a target species (e.g., homo sapiens) using 575 position weight matrices (PWMs) based on binding data from the JASPAR database. Additional experimental data from a variety of sources is used to support or detract from these predictions:

-   DNA sequence conservation in homologous mammal species sequences
-   proximity to CAGE-supported transcription start sites (TSSs)
-   correlation of expression between target gene and predicted
    transcription factor (TF) across 1800+ samples
-   proximity to ChIP-Seq determined TFBSs (GTRD project)
-   proximity to qualitative trait loci (eQTLs) affecting expression of
    the target gene (GTEX project)
-   proximity to CpGs
-   proximity to ATAC-Seq peaks (ENCODE project)

## 1.3 Ensembl transcript-ids

In order to incorporate several of the experimental datapoints, and for ease of use, the Ensembl transcript ID was chosen as the primary identifier around which TFBS_footprinting analyses would revolve. This allows users to predict TFBSs in the promoters any of 1-80,000 human protein coding transcripts in the Ensembl database, and incorporate expression data from the FANTOM project, and eQTLS from the GTEx project, which has been assigned to these IDs. TFBS predictions can also be made for 87 unique non-human species (including model organisms such as mouse and zebrafish), present in the following groups:

-   70 Eutherian mammals
-   24 Primates
-   11 Fish
-   7 Sauropsids

# 2. Installation

## 2.1 Pypi installation

The TFBS_footprinting package can be installed directly to your linux system using [PIP](https://pip.pypa.io/en/stable/installing/) install.

`$ pip install tfbs_footprinting3`

## 2.2 Docker installation

TFBS_footprinting is now available in a [Docker image](https://hub.docker.com/r/thirtysix/tfbs_footprinting/) which already contains all of the software requirements. This can be used on both Linux and Windows systems.

`$ docker pull thirtysix/tfbs_footprinting`

The Docker installation will have a default RAM allocation that is too low (\~2GB) to run TFBS_footprinting. This setting should be changed to \>6144MB. In Windows this can be adjusted by navigating through: Docker system tray icon\>Settings\>Advanced\>Memory. After changing this value, Docker will restart, which could take several minutes.

## 2.3 Requirements

These requirements are automatically installed by the PIP installer, and are already present in the Docker image:
-   biopython
-   numpy
-   matplotlib
-   httplib2
-   msgpack
-   wget

# 3. Input

## 3.1 CSV of Ensembl transcript
### IDs and arguments 
In this first option a .csv table of Ensembl transcript IDs and optional arguments is used. Except for the first column indicating the Ensembl transcript id, any argument can be left blank which will result in the default value being used. An example of this file is included in the sample_analysis folder and can be downloaded from [Github](https://github.com/thirtysix/TFBS_footprinting/blob/master/sample_analysis/sample_analysis_list.csv).

![input_csv](https://raw.githubusercontent.com/thirtysix/TFBS_footprinting/master/sample_csv.png)

## 3.2 Text file of transcript IDs

In this second option a simple text-file of Ensembl transcript IDs is used. Whatever arguments are provided to the command-line will be applied to analysis of all transcripts indicated in this file. An example of this file is included in the sample_analysis folder and can be downloaded from [Github](https://github.com/thirtysix/TFBS_footprinting/blob/master/sample_analysis/sample_ensembl_ids.txt).

![input_txt](https://raw.githubusercontent.com/thirtysix/TFBS_footprinting/master/sample_ids.txt.png)

## 3.3 Text file TF IDs

A file of JASPAR TF IDs can be provided which will limit the analysis of TFs to just those contained within. If no file name is provided then the analysis will use all 575 JASPAR TFs in the analysis, in this case the results in the output table can be filtered to just those TFs you wish to focus on. An example of this file, which contains all JASPAR TF IDs, is included in the sample_analysis folder which can be downloaded from [Github](https://github.com/thirtysix/TFBS_footprinting/blob/master/sample_analysis/sample_jaspar_tf_ids.txt).

![input_sample_tfs](https://raw.githubusercontent.com/thirtysix/TFBS_footprinting/master/sample_tf_ids.txt.png)

# 4. Output

## 4.1 Promoter Figure

-   **Output name format**
 \[Ensembl_transcript_id\]\[nt_before_TSS\]\[nt_after_TSS\]\[species_group\]\[species_coverage\]\[pvalue\].Promoterhisto.svg.

![output_primary](https://raw.githubusercontent.com/thirtysix/TFBS_footprinting/master/sample_output_figure.png)

-   **Figure Title (left)** - Ensembl gene name and Transcript ID.
-   **Predicted TFBSs** - Depicts the region analyzed in this analysis.
    Bars displayed on the positive y-axis indicate TFBS predictions made on the positive (+) DNA strand. The legend of colors is located directly above, and correlates colors to the TF which is predicted to bind. Combined affinity score is the summation of scores for the seven characteristics scored by the TFBS Footprinting tool. Using the x-axis values at the bottom of the figure, the x-axis indicates nucleotide position relative to the TSS (e.g., -1000 is 1000 nt upstream of the TSS).
 -  **TF Expression Correlation** - For the TFs with the highest combined affinity scores, the correlation of expression with the target gene has been determined and is presented in descending order of log-likelihood score. The color legend is the same as the Predicted TFBSs subplot, and is located at the top of the figure.
-   **Conservation** - A measure of the conservation at each position in the alignment.
-   **CpG Observed/Expected** - Indidents of CpG (Cytosine-Phosphate-Guanine) are indicated with a vertical black line. A line plot in red indicates the ratio of actual to expected CpGs (by chance alone) in a 200 nt window.
-   **eQTLs** - Locations of eQTLs identified in the [gTEX project](https://www.gtexportal.org/home/) which affect the target gene. Green bars on the positive y-axis indicate eQTLs assocaiate with enhanced expression of the target gene, while red bars on the negative y-axis indicate those associated with reduced expression. The magnitude of the bars indicates the projected effect on gene expression.
-   **TFBS Meta Clusters** - Locations of TFBS metaclusters representing clusters of experimentally identifed (ENCODE and SRA databases) TFBSs which have been compiled by the [GTRD project](http://gtrd.biouml.org/).
-   **ATAC-Seq** - Locations of ATAC-Seq peaks identified by ENCODE experiments representing open chromatin.
-   **CAGE peaks (TSSs)** - Locations of TSSs identified by CAGE in the [FANTOM project](http://fantom.gsc.riken.jp/).

## 4.2 Table of predictions

All predicted TFBSs for target transcript, filtered by p-value threshold and sorted by combined affinity score (TFBSs_found.sortedclusters.csv).

-   **binding prot** - Gene name of protein predicted to bind.
-   **species** - Name of target species.
-   **motif** - DNA sequence of target species which is predicted to be bound by the binding protein.
-   **strand** - Positive or negative strand of DNA.
-   **start** - Location in unaligned target species sequence of the start of the motif.
-   **end** - Location in unaligned target species sequence of the end of the motif.
-   **TSS-relative start** - Location relative to the TSS of the start of the motif.
-   **TSS-relative end** - Location relative to the TSS of the end of the motif.
-   **frame score** - Log-likelihood score of the prediction of the binding of the TF.
-   **p-value** - P-value corresponding to the frame score.
-   **pos in align.** - Location in the aligned target species sequence corresponding the start of the motif.
-   **support** - Number of species in alignment.
-   **combined affinity score** - Summation of all scores for this motif/prediction.
-   **species weights sum** - Information content of the alignment at the locations occupied by the TFBS.
-   **cage weights sum** - Summation of log-likehood scores for the CAGE peaks near the predicted TFBS.
-   **eqtls weights sum** - Summation of log-likehood scores for the eQTLs peaks near the predicted TFBS.
-   **atac weights sum** - Summation of log-likehood scores for the ATAC-Seq peaks near the predicted TFBS.
-   **metacluster weights sum** - Summation of log-likehood scores for the metacluster peaks near the predicted TFBS.
-   **cpg weight** - Log-likehood score for the CpG obs/exp ratio at the center of the predicted TFBS.
-   **corr weight sum** - Summation of log-likehood scores for correlations of expression between transcripts of the predicted TF and those of the target gene.

## 4.3 Additional outputs

-   Original alignment as retrieved from Ensembl (alignment_uncleaned.fasta).
-   Cleaned alignment (alignment_cleaned.fasta).
-   Regulatory information for the target transcripts user-defined promoter region (regulatory_decoded.json).
-   Transcript properties for target transcript (transcript_dict.json).
-   All predicted TFBSs for the target species which satisfy p-value threshold (TFBSs_found.all.json).

# 5. Examples

## 5.1 TFBS_footprinter3 Use Examples
**Running the sample analyses**
- Run the sample analysis using a .csv of arguments:

    > `$ tfbs_footprinter3 -t PATH_TO/sample_analysis/sample_analysis_list.csv`

-   Run the sample analysis using a .txt of Ensembl transcript ids, and minimal arguments:

    > `$ tfbs_footprinter3 -t PATH_TO/sample_analysis/sample_ensembl_ids.txt`

-   Run the sample analysis using a .txt of Ensembl transcript ids, and all arguments:

    > `$ tfbs_footprinter3 -t PATH_TO/sample_ensembl_ids.txt -tfs PATH_TO/sample_jaspar_tf_ids.txt -pb 900 -pa 100 -tx 10 -p 0.01 -update`

## 5.2 TFBS_footprinter3 Use Examples (within Docker)
**Running the sample analyses**
1. Within Docker we first need to mount a volume so that the results of the analyis can be viewed on our host computer. It is recommended that you create an empty directory on your host computer:

    > `$ docker run -v /ABSOLUTE_PATH_TO/EMPTY_DIR_ON_HOST:/home/sample_analysis/tfbs_results -it thirtysix/tfbs_footprinting bash`

2.  Then we move into the pre-existing sample analysis directory in the Docker container to perform the analysis there so that the results generated there will automatically appear in the designated location on our host computer:

    > `$ cd ./sample_analysis`

3.  Then we can run the sample analysis in Docker in the same way that we would normally use tfbs_footprinter3 (above), e.g., using a .csv of arguments:

    > `$ tfbs_footprinter3 -t ./sample_analysis_list.csv`

-   Or (again, as above) using a .txt of Ensembl transcript ids, and minimal arguments:

    > `$tfbs_footprinter3 -t ./sample_ensembl_ids.txt`

-   Or (again, as above) using a .txt of Ensembl transcript ids, and multiple arguments:

    > `$ tfbs_footprinter3 -t ./sample_ensembl_ids.txt -tfs ./sample_jaspar_tf_ids.txt -pb 900 -pa 100 -tx 10 -update`

### Example using user-defined files/arguments

1.  Within Docker we first need to mount a volume so that we can load your analysis files from your host computer into docker AND save the results of the analysis on our host computer:

    > `$ docker run -v /ABSOLUTE_PATH_TO/DIR_ON_HOST/CONTAINING_ANALYSIS_FILES:/home/analysis_dir -it thirtysix/tfbs_footprinting bash`

2.  Then we move into your analysis directory in the Docker container to perform the analysis there so that the results generated there will automatically appear in the designated location on our host computer:

    > `$ cd ./analysis_dir`

3.  Then we can run the sample analysis in Docker in the same way that we would normally use tfbs_footprinter3 (above), e.g., using a .csv of arguments:

    > `$ tfbs_footprinter3 -t ./USER_TABLE_OF_ENSEMBL_IDS_AND_ARGS.csv`

-   Or (again, as above) using a .txt of Ensembl transcript ids, and minimal arguments:

    > `$ tfbs_footprinter3 -t ./USER_LIST_OF_ENSEMBL_IDS.txt`

-   Or (again, as above) using a .txt of Ensembl transcript ids, and multiple arguments:

    > `$ tfbs_footprinter3 -t ./USER_LIST_OF_ENSEMBL_IDS.txt -tfs ./USER_LIST_OF_TF_NAMES.txt -s homo_sapiens -g mammals -e low -pb 900 -pa 100 -tx 10 -o PATH_TO/Results/ -update`

## 5.3 Update experimental data files

Check the documentation on [ReadTheDocs](https://tfbs-footprinting.readthedocs.io/en/latest/) for new releases.

> `$ tfbs_footprinter3 -update`

# 6. Parameters

**arguments:**
- **\-h, --help**
show this help message and exit.

-  **\-\-t_ids_file , -t**
Required for running an analysis. Location of a file containing Ensembl target_species transcript ids. Input options are either a text file of Ensembl transcript ids or a .csv file with individual values set for each parameter.

-  **\-\-tf_ids_file , -tfs**
Optional: Location of a file containing a limited list of Jaspar TFs to use in scoring alignment (see sample file tf_ids.txt at https://github.com/thirtysix/TFBS_footprinting) [default: all Jaspar TFs].

-  **\-\-promoter_before_tss , -pb** 
(0-100,000) [default: 900] - Number (integer) of nucleotides upstream of TSS to include in analysis. If this number is negative the start point will be downstream of the TSS, the end point will then need to be further downstream.

-  **\-\-promoter_after_tss , -pa** 
(0-100,000) [default: 100] - Number (integer) of nucleotides downstream of TSS to include in analysis. If this number is negative the end point will be upstream of the TSS. The start point will then need to be further upstream.

-  **\-\-top_x_tfs , -tx**     
(1-20) [default: 10] - Number (integer) of unique TFs to include in output .svg figure.

-  **\-\-pval, -p** 
PVAL  P-value (float) for PWM score cutoff (range: 1 (all results) to 0.0000001; in divisions of 10 (i.e. 1, 0.1, 0.01, 0.001 etc.) [default: 0.01].

-  **\-\-pvalc, -pc**
P-value (float) for PWM score cutoff (range: 1 (all results) to 0.0000001; in divisions of 10 (i.e. 1, 0.1, 0.01, 0.001 etc.) [default: 0.01].

-  **\-\-exp_data_update, -update**
Download the latest experimental data files for use in analysis. Will run automatically if the "data" directory does not already exist (e.g., first usage).

-  **\-\-nofig, -no**
Don't output a figure.


# 7.  Process

## 7.1 Steps
### Iterate through each user provided Ensembl transcript id:

**1.  Retrieve EPO aligned orthologous sequences from Ensembl** database for user-defined species group (mammals, primates, fish, sauropsids) for promoter of user-provided transcript id, between user-defined TSS-relative start/stop sites.

**2. Edit retrieved alignment:**

- Replace characters not corresponding to nucleotides (ACGT), with gaps characters "-".

- Remove gap-only columns from alignment.

**3. Generate position weight matrices (PWMs)** from Jaspar position frequency matrices (PFMs). 

**4. Score target species sequence** using either all or a user-defined  list of PWMs.

**5. Keep predictions with a log-likelihood score greater** than score threshold corresponding to p-value of 0.01, or user-defined  p-value. 

**6. Score each of the following for the target sequence region**:

- DNA sequence conservation in homologous mammal species sequences

- proximity to CAGE-supported transcription start sites (TSSs)

- correlation of expression between target gene and predicted transcription factor (TF) across 1800 samples

- proximity to ChIP-Seq determined TFBSs (GTRD project)

- proximity to qualitative trait loci (eQTLs) affecting expression of the target gene (GTEX project)

- proximity to CpGs

- proximity to ATAC-Seq peaks (ENCODE project) 

**7. Compute 'combined affinity score'** as a sum of scores for all  experimental data. 

**8. Output predictions table** 

**9. Output a vector graphics figure** showing the top 10 (or user-defined) unique TFs mapped onto the promoter of the target transcript, and additional output as described below.



## 7.2 Flowchart

![image](https://raw.githubusercontent.com/thirtysix/TFBS_footprinting/master/flowchart.png)

# 8. Ensembl available species
**Table of species which can be analyzed by tfbs_footprinter3 at this time**

|EPO_LOW mammals           |EPO_LOW fish          |EPO_LOW sauropsids |EPO mammals          |EPO primates       |EPO fish              |EPO sauropsids     |
|--------------------------|----------------------|-------------------|---------------------|-------------------|----------------------|-------------------|
|ailuropoda_melanoleuca    |astyanax_mexicanus    |anas_platyrhynchos |bos_taurus           |callithrix_jacchus |danio_rerio           |anolis_carolinensis|
|bos_taurus                |danio_rerio           |anolis_carolinensis|callithrix_jacchus   |chlorocebus_sabaeus|gasterosteus_aculeatus|gallus_gallus      |
|callithrix_jacchus        |gadus_morhua          |ficedula_albicollis|canis_familiaris     |gorilla_gorilla    |lepisosteus_oculatus  |meleagris_gallopavo|
|canis_familiaris          |gasterosteus_aculeatus|gallus_gallus      |chlorocebus_sabaeus  |homo_sapiens       |oryzias_latipes       |taeniopygia_guttata|
|cavia_porcellus           |lepisosteus_oculatus  |meleagris_gallopavo|equus_caballus       |macaca_mulatta     |tetraodon_nigroviridis|                   |
|chlorocebus_sabaeus       |oreochromis_niloticus |pelodiscus_sinensis|felis_catus          |pan_troglodytes    |                      |                   |
|choloepus_hoffmanni       |oryzias_latipes       |taeniopygia_guttata|gorilla_gorilla      |papio_anubis       |                      |                   |
|dasypus_novemcinctus      |poecilia_formosa      |                   |homo_sapiens         |pongo_abelii       |                      |                   |
|dipodomys_ordii           |takifugu_rubripes     |                   |macaca_mulatta       |                   |                      |                   |
|echinops_telfairi         |tetraodon_nigroviridis|                   |mus_musculus         |                   |                      |                   |
|equus_caballus            |xiphophorus_maculatus |                   |oryctolagus_cuniculus|                   |                      |                   |
|erinaceus_europaeus       |                      |                   |ovis_aries           |                   |                      |                   |
|felis_catus               |                      |                   |pan_troglodytes      |                   |                      |                   |
|gorilla_gorilla           |                      |                   |papio_anubis         |                   |                      |                   |
|homo_sapiens              |                      |                   |pongo_abelii         |                   |                      |                   |
|ictidomys_tridecemlineatus|                      |                   |rattus_norvegicus    |                   |                      |                   |
|loxodonta_africana        |                      |                   |sus_scrofa           |                   |                      |                   |
|macaca_mulatta            |                      |                   |                     |                   |                      |                   |
|microcebus_murinus        |                      |                   |                     |                   |                      |                   |
|mus_musculus              |                      |                   |                     |                   |                      |                   |
|mustela_putorius_furo     |                      |                   |                     |                   |                      |                   |
|myotis_lucifugus          |                      |                   |                     |                   |                      |                   |
|nomascus_leucogenys       |                      |                   |                     |                   |                      |                   |
|ochotona_princeps         |                      |                   |                     |                   |                      |                   |
|oryctolagus_cuniculus     |                      |                   |                     |                   |                      |                   |
|otolemur_garnettii        |                      |                   |                     |                   |                      |                   |
|ovis_aries                |                      |                   |                     |                   |                      |                   |
|pan_troglodytes           |                      |                   |                     |                   |                      |                   |
|papio_anubis              |                      |                   |                     |                   |                      |                   |
|pongo_abelii              |                      |                   |                     |                   |                      |                   |
|procavia_capensis         |                      |                   |                     |                   |                      |                   |
|pteropus_vampyrus         |                      |                   |                     |                   |                      |                   |
|rattus_norvegicus         |                      |                   |                     |                   |                      |                   |
|sorex_araneus             |                      |                   |                     |                   |                      |                   |
|sus_scrofa                |                      |                   |                     |                   |                      |                   |
|tarsius_syrichta          |                      |                   |                     |                   |                      |                   |
|tupaia_belangeri          |                      |                   |                     |                   |                      |                   |
|tursiops_truncatus        |                      |                   |                     |                   |                      |                   |
|vicugna_pacos             |                      |                   |                     |                   |                      |                   |


# 9. Troubleshooting

## 9.1 Log file

The first step in troubleshooting any possible errors that might arise is to check the log file. The log file is created in the directory where you initiated the analysis, and it is named \'TFBS_footprinter3.log\'.

**Many relevant events are logged there instead of output to terminal:**
-   start time of analysis
-   arguments/settings used in analysis
-   full queries made to the Ensembl REST system
-   warnings from Ensembl REST system regarding exceeding of rate limit
-   if given transcript ids are not in the Ensembl system (possibly misformed or deprecated)
-   if experimental data used in TFBS prediction has been downloaded
-   if there was an error in downloading experimental data
-   if results already exist for the current analysis, which are then loaded
-   total time of analysis
-   if there was an error in retrieving an alignment for the defined region
-   if the transcript information file already exists
-   if the transcript regulatory information file already exists

## 9.2 Docker

The Docker installation will have a default RAM allocation that is too low (\~2GB) to run TFBS_footprinting. This setting should be changed to \>6144MB.

In Windows this can be adjusted by navigating through:

**Docker system tray icon\>Settings\>Advanced\>Memory**.

After changing this value Docker will restart, which could take several minutes. If the allocated memory is too low, then the Docker image will terminate and you will see the `$ killed` message in the console.
