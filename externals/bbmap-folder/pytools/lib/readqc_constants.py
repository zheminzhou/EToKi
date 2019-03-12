#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Constants specific for RQC ReadQC pipeline (Illumina).
Global constants are defined in a file under the directory lib/rqc_constants.py .

Created: 7.29.2013

sulsj (ssul@lbl.gov)

"""

import os

# from rqc_constants import *
from common import get_run_path
# from os_utility import get_tool_path

# RQCReadQcConfig, RQCReadQc, ReadqcStats, RQCReadQcCommands

# global configuration variables
class RQCReadQcConfig:
    CFG = {}
    CFG["output_path"] = ""
    CFG["status_file"] = ""
    CFG["files_file"] = ""
    CFG["stats_file"] = ""
    CFG["log_file"] = ""
    CFG["no_cleanup"] = ""
    CFG["run_path"] = ""


# ReadQC constants
class RQCReadQc:
    ILLUMINA_SAMPLE_COUNT = 50000           ## Used in STEP13
    ILLUMINA_SAMPLE_PERCENTAGE = 0.01       ## Used in STEP1 and STEP13
    ILLUMINA_MER_SAMPLE_MER_SIZE = 25       ## Used in STEP2 (20 ==> 25 RQC-823 08102016)
    ILLUMINA_MER_SAMPLE_REPORT_FRQ = 25000  ## Used in STEP2 for mer sampling
    ILLUMINA_MIN_NUM_READS = 1000           ## Used in STEP1 to check min number of reads needed for processing


#Reference databases
class RQCReadQcReferenceDatabases:
    ## STEP10
    ARTIFACT_FILE = "/global/dna/shared/rqc/ref_databases/qaqc/databases/illumina.artifacts/Illumina.artifacts.2012.04.fa"

    ## STEP12 for bbduk
    SCICLONE_RNA2 = "/global/projectb/sandbox/gaag/bbtools/data/sciclone_rna.fa"
    SCICLONE_DNA2 = "/global/projectb/sandbox/gaag/bbtools/data/sciclone_dna.fa"

    ## STEP17
    END_OF_READ_ILLUMINA_ADAPTER_CHECK_DB = "/global/dna/shared/rqc/ref_databases/qaqc/databases/Artifacts.adapters_primers_only.fa"


class RQCContamDb:
    # ARTIFACT_FILE_NO_SPIKEIN  = '/global/dna/shared/rqc/ref_databases/qaqc/databases/illumina.artifacts/Illumina.artifacts.2012.10.no_DNA_RNA_spikeins.fa' #this db includes no DNA/RNA spike-in sequences
    ARTIFACT_FILE_NO_SPIKEIN  = '/global/dna/shared/rqc/ref_databases/qaqc/databases/illumina.artifacts/Illumina.artifacts.2013.12.no_DNA_RNA_spikeins.fa'
    ARTIFACT_FILE_DNA_SPIKEIN = '/global/dna/shared/rqc/ref_databases/qaqc/databases/illumina.artifacts/DNA_spikeins.artifacts.2012.10.fa.bak' #this db has just DNA spike-in sequences
    ARTIFACT_FILE_RNA_SPIKEIN = '/global/dna/shared/rqc/ref_databases/qaqc/databases/illumina.artifacts/RNA_spikeins.artifacts.2012.10.NoPolyA.fa' #this db has just RNA spike-in sequences
    CONTAMINANTS              = '/global/dna/shared/rqc/ref_databases/qaqc/databases/JGIContaminants.fa'  ## '/home/blast_db2_admin/qaqc_db/2010-11-19/JGIContaminants.fa'
    FOSMID_VECTOR             = '/global/dna/shared/rqc/ref_databases/qaqc/databases/pCC1Fos.ref.fa'
    PHIX                      = '/global/dna/shared/rqc/ref_databases/qaqc/databases/phix174_ill.ref.fa'
    GENERAL_RRNA_FILE         = '/global/dna/shared/rqc/ref_databases/qaqc/databases/rRNA.fa' ## including Chlamy

    CHLOROPLAST_NCBI_REFSEQ   = '/global/dna/shared/rqc/ref_databases/ncbi/CURRENT/refseq.plastid/refseq.plastid'
    MITOCHONDRION_NCBI_REFSEQ = '/global/dna/shared/rqc/ref_databases/ncbi/CURRENT/refseq.mitochondrion/refseq.mitochondrion'

    MICROBES = "/global/projectb/sandbox/gaag/bbtools/commonMicrobes/fusedERPBBmasked.fa.gz" ## non-synthetic contaminants
    SYNTHETIC = "/global/projectb/sandbox/gaag/bbtools/data/Illumina.artifacts.2013.12.no_DNA_RNA_spikeins.fa.gz" ## synthetic contaminants
    ADAPTERS = "/global/projectb/sandbox/gaag/bbtools/data/adapters.fa"

    ILLUMINA_READ_PERCENT_CONTAMINATION_ARTIFACT = "illumina read percent contamination artifact"
    ILLUMINA_READ_PERCENT_CONTAMINATION_ARTIFACT_50BP = "illumina read percent contamination artifact 50bp" ## 20131203 Added for 50bp contam
    ILLUMINA_READ_PERCENT_CONTAMINATION_ARTIFACT_20BP = "illumina read percent contamination artifact 20bp" ## 11092015 Added for 20bp contam for smRNA
    ILLUMINA_READ_PERCENT_CONTAMINATION_DNA_SPIKEIN = "illumina read percent contamination DNA spikein"
    ILLUMINA_READ_PERCENT_CONTAMINATION_RNA_SPIKEIN = "illumina read percent contamination RNA spikein"
    ILLUMINA_READ_PERCENT_CONTAMINATION_CONTAMINANTS = "illumina read percent contamination contaminants"
    # ILLUMINA_READ_PERCENT_CONTAMINATION_ECOLI_COMBINED = "illumina read percent contamination ecoli combined"
    ILLUMINA_READ_PERCENT_CONTAMINATION_FOSMID = "illumina read percent contamination fosmid"
    ILLUMINA_READ_PERCENT_CONTAMINATION_MITOCHONDRION = "illumina read percent contamination mitochondrion"
    ILLUMINA_READ_PERCENT_CONTAMINATION_PHIX = "illumina read percent contamination phix"
    ILLUMINA_READ_PERCENT_CONTAMINATION_PLASTID = "illumina read percent contamination plastid"
    ILLUMINA_READ_PERCENT_CONTAMINATION_RRNA = "illumina read percent contamination rrna"
    ILLUMINA_READ_PERCENT_CONTAMINATION_MICROBES = "illumina read percent contamination microbes" ## non-synthetic
    ILLUMINA_READ_PERCENT_CONTAMINATION_SYNTH = "illumina read percent contamination adapters"
    ILLUMINA_READ_PERCENT_CONTAMINATION_ADAPTERS = "illumina read percent contamination adapters"

    CONTAM_DBS = {}
    CONTAM_DBS['artifact'] = ARTIFACT_FILE_NO_SPIKEIN
    CONTAM_DBS['artifact_50bp'] = ARTIFACT_FILE_NO_SPIKEIN ## 20131203 Added for 50bp contam
    CONTAM_DBS['DNA_spikein'] = ARTIFACT_FILE_DNA_SPIKEIN
    CONTAM_DBS['RNA_spikein'] = ARTIFACT_FILE_RNA_SPIKEIN
    CONTAM_DBS['contaminants'] = CONTAMINANTS
    # CONTAM_DBS['ecoli_combined'] = ECOLI_COMBINED
    CONTAM_DBS['fosmid'] = FOSMID_VECTOR
    CONTAM_DBS['mitochondrion'] = MITOCHONDRION_NCBI_REFSEQ
    CONTAM_DBS['phix'] = PHIX
    CONTAM_DBS['plastid'] = CHLOROPLAST_NCBI_REFSEQ
    CONTAM_DBS['rrna'] = GENERAL_RRNA_FILE
    CONTAM_DBS['microbes'] = MICROBES ## non-synthetic
    CONTAM_DBS['synthetic'] = SYNTHETIC
    CONTAM_DBS['adapters'] = ADAPTERS

    CONTAM_KEYS = {}
    CONTAM_KEYS['artifact'] = ILLUMINA_READ_PERCENT_CONTAMINATION_ARTIFACT
    CONTAM_KEYS['artifact_50bp'] = ILLUMINA_READ_PERCENT_CONTAMINATION_ARTIFACT_50BP ## 12032013 Added for 50bp contam
    CONTAM_KEYS['artifact_20bp'] = ILLUMINA_READ_PERCENT_CONTAMINATION_ARTIFACT_20BP ## 11092015 Added for 20bp contam for smRNA
    CONTAM_KEYS['DNA_spikein'] = ILLUMINA_READ_PERCENT_CONTAMINATION_DNA_SPIKEIN
    CONTAM_KEYS['RNA_spikein'] = ILLUMINA_READ_PERCENT_CONTAMINATION_RNA_SPIKEIN
    CONTAM_KEYS['contaminants'] = ILLUMINA_READ_PERCENT_CONTAMINATION_CONTAMINANTS
    # CONTAM_KEYS['ecoli_combined'] = ILLUMINA_READ_PERCENT_CONTAMINATION_ECOLI_COMBINED
    CONTAM_KEYS['fosmid'] = ILLUMINA_READ_PERCENT_CONTAMINATION_FOSMID
    CONTAM_KEYS['mitochondrion'] = ILLUMINA_READ_PERCENT_CONTAMINATION_MITOCHONDRION
    CONTAM_KEYS['phix'] = ILLUMINA_READ_PERCENT_CONTAMINATION_PHIX
    CONTAM_KEYS['plastid'] = ILLUMINA_READ_PERCENT_CONTAMINATION_PLASTID
    CONTAM_KEYS['rrna'] = ILLUMINA_READ_PERCENT_CONTAMINATION_RRNA
    CONTAM_KEYS['microbes'] = ILLUMINA_READ_PERCENT_CONTAMINATION_MICROBES ## non-synthetic
    CONTAM_KEYS['synthetic'] = ILLUMINA_READ_PERCENT_CONTAMINATION_SYNTH
    CONTAM_KEYS['adapters'] = ILLUMINA_READ_PERCENT_CONTAMINATION_ADAPTERS

## Constants for reporting
class ReadqcStats:
    ILLUMINA_READ_Q20_READ1 = "read q20 read1"
    ILLUMINA_READ_Q20_READ2 = "read q20 read2"
    ILLUMINA_READ_QHIST_TEXT = "read qhist text"
    ILLUMINA_READ_QHIST_PLOT = "read qhist plot"
    ILLUMINA_READ_QHIST_D3_HTML_PLOT = "ILLUMINA_READ_QHIST_D3_HTML_PLOT"
    ILLUMINA_READ_QUAL_POS_PLOT_1 = "read qual pos plot 1"
    ILLUMINA_READ_QUAL_POS_PLOT_2 = "read qual pos plot 2"
    ILLUMINA_READ_QUAL_POS_PLOT_MERGED = "read qual pos plot merged"
    ILLUMINA_READ_QUAL_POS_PLOT_MERGED_D3_HTML_PLOT = "ILLUMINA_READ_QUAL_POS_PLOT_MERGED_D3_HTML_PLOT"
    ILLUMINA_READ_QUAL_POS_QRPT_1 = "read qual pos qrpt 1"
    ILLUMINA_READ_QUAL_POS_QRPT_2 = "read qual pos qrpt 2"
    ILLUMINA_READ_QUAL_POS_QRPT_MERGED = "ILLUMINA_READ_QUAL_POS_QRPT_MERGED"
    ILLUMINA_READ_QUAL_POS_QRPT_BOXPLOT_1 = "ILLUMINA_READ_QUAL_POS_QRPT_BOXPLOT_1"
    ILLUMINA_READ_QUAL_POS_QRPT_BOXPLOT_2 = "ILLUMINA_READ_QUAL_POS_QRPT_BOXPLOT_2"
    ILLUMINA_READ_QUAL_POS_QRPT_D3_HTML_BOXPLOT_1 = "ILLUMINA_READ_QUAL_POS_QRPT_D3_HTML_BOXPLOT_1"
    ILLUMINA_READ_QUAL_POS_QRPT_D3_HTML_BOXPLOT_2 = "ILLUMINA_READ_QUAL_POS_QRPT_D3_HTML_BOXPLOT_2"
    ILLUMINA_READ_QUAL_POS_QRPT_BOXPLOT_TEXT = "ILLUMINA_READ_QUAL_POS_QRPT_BOXPLOT_TEXT"

    ILLUMINA_READ_BASE_COUNT_TEXT_1 = "read base count text 1"
    ILLUMINA_READ_BASE_COUNT_TEXT_2 = "read base count text 2"
    ILLUMINA_READ_BASE_COUNT_PLOT_1 = "read base count plot 1"
    ILLUMINA_READ_BASE_COUNT_PLOT_2 = "read base count plot 2"
    ILLUMINA_READ_BASE_COUNT_D3_HTML_PLOT_1 = "ILLUMINA_READ_BASE_COUNT_D3_HTML_PLOT_1" # cyclone nucleotide comp plot
    ILLUMINA_READ_BASE_COUNT_D3_HTML_PLOT_2 = "ILLUMINA_READ_BASE_COUNT_D3_HTML_PLOT_2" # cyclone nucleotide comp plot

    ILLUMINA_READ_BASE_PERCENTAGE_TEXT_1 = "read base percentage text 1"
    ILLUMINA_READ_BASE_PERCENTAGE_TEXT_2 = "read base percentage text 2"
    ILLUMINA_READ_BASE_PERCENTAGE_PLOT_1 = "read base percentage plot 1"
    ILLUMINA_READ_BASE_PERCENTAGE_PLOT_2 = "read base percentage plot 2"
    ILLUMINA_READ_BASE_PERCENTAGE_D3_HTML_PLOT_1 = "ILLUMINA_READ_BASE_PERCENTAGE_D3_HTML_PLOT_1"
    ILLUMINA_READ_BASE_PERCENTAGE_D3_HTML_PLOT_2 = "ILLUMINA_READ_BASE_PERCENTAGE_D3_HTML_PLOT_2"

    ILLUMINA_READ_20MER_SAMPLE_SIZE = "read 20mer sample size"
    ILLUMINA_READ_20MER_PERCENTAGE_STARTING_MERS = "read 20mer percentage starting mers"
    ILLUMINA_READ_20MER_PERCENTAGE_RANDOM_MERS = "read 20mer percentage random mers"
    ILLUMINA_READ_20MER_UNIQUENESS_TEXT = "read 20mer uniqueness text"
    ILLUMINA_READ_20MER_UNIQUENESS_PLOT = "read 20mer uniqueness plot"
    ILLUMINA_READ_BWA_ALIGNED = "read bwa aligned"
    ILLUMINA_READ_BWA_ALIGNED_DUPLICATE = "read bwa aligned duplicate"
    ILLUMINA_READ_BWA_ALIGNED_DUPLICATE_PERCENT = "read bwa aligned duplicate percent"

    ILLUMINA_READ_N_FREQUENCE = "read N frequence"
    ILLUMINA_READ_N_PATTERN = "read N pattern"
    ILLUMINA_READ_GC_MEAN = "read GC mean"
    ILLUMINA_READ_GC_STD = "read GC std"
    ILLUMINA_READ_GC_MED = "read GC median"
    ILLUMINA_READ_GC_MODE = "read GC mode"

    ILLUMINA_READ_GC_PLOT = "read GC plot"
    ILLUMINA_READ_GC_D3_HTML_PLOT = "ILLUMINA_READ_GC_D3_HTML_PLOT"

    ILLUMINA_READ_GC_TEXT = "read GC text hist"
    ILLUMINA_READ_LENGTH_1 = "read length 1"
    ILLUMINA_READ_LENGTH_2 = "read length 2"
    ILLUMINA_READ_BASE_COUNT = "read base count"
    ILLUMINA_READ_COUNT = "read count"

    ILLUMINA_READ_TOPHIT_FILE = "read tophit file of"
    ILLUMINA_READ_TOP100HIT_FILE = "read top100hit file of"
    ILLUMINA_READ_TAXLIST_FILE = "read taxlist file of"

    ILLUMINA_READ_TAX_SPECIES = "read tax species of"
    ILLUMINA_READ_TOP_HITS = "read top hits of"
    ILLUMINA_READ_TOP_100HITS = "read top 100 hits of"
    ILLUMINA_READ_PARSED_FILE = "read parsed file of"
    ILLUMINA_READ_MATCHING_HITS = "read matching hits of"
    ILLUMINA_READS_NUMBER = "reads number"

    ILLUMINA_READ_DEMULTIPLEX_STATS = "demultiplex stats"
    ILLUMINA_READ_DEMULTIPLEX_STATS_PLOT = "demultiplex stats plot"
    ILLUMINA_READ_DEMULTIPLEX_STATS_D3_HTML_PLOT = "ILLUMINA_READ_DEMULTIPLEX_STATS_D3_HTML_PLOT"

    ILLUMINA_READ_BASE_QUALITY_STATS = "read base quality stats"
    ILLUMINA_READ_BASE_QUALITY_STATS_PLOT = "read base quality stats plot"

    ILLUMINA_READ_PERCENT_CONTAMINATION_ARTIFACT = "illumina read percent contamination artifact"
    ILLUMINA_READ_PERCENT_CONTAMINATION_ARTIFACT_50BP = "illumina read percent contamination artifact 50bp" ## 20131203 Added for 50bp contam
    ILLUMINA_READ_PERCENT_CONTAMINATION_ARTIFACT_20BP = "illumina read percent contamination artifact 20bp" ## 11092015 Added for 20bp contam for smRNA
    ILLUMINA_READ_PERCENT_CONTAMINATION_DNA_SPIKEIN = "illumina read percent contamination DNA spikein"
    ILLUMINA_READ_PERCENT_CONTAMINATION_RNA_SPIKEIN = "illumina read percent contamination RNA spikein"
    ILLUMINA_READ_PERCENT_CONTAMINATION_CONTAMINANTS = "illumina read percent contamination contaminants"
    ILLUMINA_READ_PERCENT_CONTAMINATION_ECOLI_B = "illumina read percent contamination ecoli b"
    ILLUMINA_READ_PERCENT_CONTAMINATION_ECOLI_K12 = "illumina read percent contamination ecoli k12"
    # ILLUMINA_READ_PERCENT_CONTAMINATION_ECOLI_COMBINED = "illumina read percent contamination ecoli combined"
    ILLUMINA_READ_PERCENT_CONTAMINATION_FOSMID = "illumina read percent contamination fosmid"
    ILLUMINA_READ_PERCENT_CONTAMINATION_MITOCHONDRION = "illumina read percent contamination mitochondrion"
    ILLUMINA_READ_PERCENT_CONTAMINATION_PLASTID = "illumina read percent contamination plastid"
    ILLUMINA_READ_PERCENT_CONTAMINATION_PHIX = "illumina read percent contamination phix"
    ILLUMINA_READ_PERCENT_CONTAMINATION_RRNA = "illumina read percent contamination rrna"
    ILLUMINA_READ_SCICLONE_DNA_COUNT_FILE = "illumina read sciclone DNA count file"
    ILLUMINA_READ_SCICLONE_RNA_COUNT_FILE = "illumina read sciclone RNA count file"

    ILLUMINA_READ_SCICLONE_RNA_COUNT_TOTAL = "ILLUMINA_READ_SCICLONE_RNA_COUNT_TOTAL"
    ILLUMINA_READ_SCICLONE_RNA_COUNT_MATCHED = "ILLUMINA_READ_SCICLONE_RNA_COUNT_MATCHED"
    ILLUMINA_READ_SCICLONE_RNA_COUNT_MATCHED_PERC = "ILLUMINA_READ_SCICLONE_RNA_COUNT_MATCHED_PERC"
    ILLUMINA_READ_SCICLONE_DNA_COUNT_TOTAL = "ILLUMINA_READ_SCICLONE_DNA_COUNT_TOTAL"
    ILLUMINA_READ_SCICLONE_DNA_COUNT_MATCHED = "ILLUMINA_READ_SCICLONE_DNA_COUNT_MATCHED"
    ILLUMINA_READ_SCICLONE_DNA_COUNT_MATCHED_PERC = "ILLUMINA_READ_SCICLONE_DNA_COUNT_MATCHED_PERC"

    ILLUMINA_BASE_Q30 = 'base Q30'
    ILLUMINA_BASE_Q25 = 'base Q25'
    ILLUMINA_BASE_Q20 = 'base Q20'
    ILLUMINA_BASE_Q15 = 'base Q15'
    ILLUMINA_BASE_Q10 = 'base Q10'
    ILLUMINA_BASE_Q5 = 'base Q5'

    ILLUMINA_BASE_C30 = 'base C30'
    ILLUMINA_BASE_C25 = 'base C25'
    ILLUMINA_BASE_C20 = 'base C20'
    ILLUMINA_BASE_C15 = 'base C15'
    ILLUMINA_BASE_C10 = 'base C10'
    ILLUMINA_BASE_C5 = 'base C5'

    ILLUMINA_READ_Q30 = 'read Q30'
    ILLUMINA_READ_Q25 = 'read Q25'
    ILLUMINA_READ_Q20 = 'read Q20'
    ILLUMINA_READ_Q15 = 'read Q15'
    ILLUMINA_READ_Q10 = 'read Q10'
    ILLUMINA_READ_Q5 = 'read Q5'

    ILLUMINA_BASE_Q30_SCORE_MEAN = 'Q30 bases Q score mean'
    ILLUMINA_BASE_Q30_SCORE_STD = 'Q30 bases Q score std'
    ILLUMINA_BASE_OVERALL_BASES_Q_SCORE_MEAN = 'overall bases Q score mean'
    ILLUMINA_BASE_OVERALL_BASES_Q_SCORE_STD = 'overall bases Q score std'

    ## read qc step 1
    ILLUMINA_TOO_SMALL_NUM_READS = "ILLUMINA_TOO_SMALL_NUM_READS" ## flag to notify too low number of reads

    ## read qc step 2
    ILLUMINA_READ_20MER_UNIQUENESS_D3_HTML_PLOT = "ILLUMINA_READ_20MER_UNIQUENESS_D3_HTML_PLOT"

    ## read qc step 13, 14, 15
    ILLUMINA_SKIP_BLAST = "ILLUMINA_SKIP_BLAST" ## to record the blast step was skipped

    ## read qc step 17
    ILLUMINA_READ_END_OF_READ_ADAPTER_CHECK_DATA = "ILLUMINA_READ_END_OF_READ_ADAPTER_CHECK_DATA"
    ILLUMINA_READ_END_OF_READ_ADAPTER_CHECK_PLOT = "ILLUMINA_READ_END_OF_READ_ADAPTER_CHECK_PLOT"
    ILLUMINA_READ_END_OF_READ_ADAPTER_CHECK_D3_HTML_PLOT = "ILLUMINA_READ_END_OF_READ_ADAPTER_CHECK_D3_HTML_PLOT"

    ## read qc step 18
    ILLUMINA_READ_INSERT_SIZE_HISTO_PLOT = "ILLUMINA_READ_INSERT_SIZE_HISTO_PLOT"
    ILLUMINA_READ_INSERT_SIZE_HISTO_DATA = "ILLUMINA_READ_INSERT_SIZE_HISTO_DATA"
    ILLUMINA_READ_INSERT_SIZE_HISTO_D3_HTML_PLOT = "ILLUMINA_READ_INSERT_SIZE_HISTO_D3_HTML_PLOT"

    ILLUMINA_READ_INSERT_SIZE_TOTAL_TIME = "ILLUMINA_READ_INSERT_SIZE_TOTAL_TIME"
    ILLUMINA_READ_INSERT_SIZE_NUM_READS = "ILLUMINA_READ_INSERT_SIZE_NUM_READS"
    ILLUMINA_READ_INSERT_SIZE_JOINED_NUM = "ILLUMINA_READ_INSERT_SIZE_JOINED_NUM"
    ILLUMINA_READ_INSERT_SIZE_JOINED_PERC = "ILLUMINA_READ_INSERT_SIZE_JOINED_PERC"

    ILLUMINA_READ_INSERT_SIZE_AMBIGUOUS_NUM = "ILLUMINA_READ_INSERT_SIZE_AMBIGUOUS_NUM"
    ILLUMINA_READ_INSERT_SIZE_AMBIGUOUS_PERC = "ILLUMINA_READ_INSERT_SIZE_AMBIGUOUS_PERC"
    ILLUMINA_READ_INSERT_SIZE_NO_SOLUTION_NUM = "ILLUMINA_READ_INSERT_SIZE_NO_SOLUTION_NUM"
    ILLUMINA_READ_INSERT_SIZE_NO_SOLUTION_PERC = "ILLUMINA_READ_INSERT_SIZE_NO_SOLUTION_PERC"
    ILLUMINA_READ_INSERT_SIZE_TOO_SHORT_NUM = "ILLUMINA_READ_INSERT_SIZE_TOO_SHORT_NUM"
    ILLUMINA_READ_INSERT_SIZE_TOO_SHORT_PERC = "ILLUMINA_READ_INSERT_SIZE_TOO_SHORT_PERC"

    ILLUMINA_READ_INSERT_SIZE_AVG_INSERT = "ILLUMINA_READ_INSERT_SIZE_AVG_INSERT"
    ILLUMINA_READ_INSERT_SIZE_STD_INSERT = "ILLUMINA_READ_INSERT_SIZE_STD_INSERT"
    ILLUMINA_READ_INSERT_SIZE_MODE_INSERT = "ILLUMINA_READ_INSERT_SIZE_MODE_INSERT"

    ILLUMINA_READ_INSERT_SIZE_INSERT_RANGE_START = "ILLUMINA_READ_INSERT_SIZE_INSERT_RANGE_START"
    ILLUMINA_READ_INSERT_SIZE_INSERT_RANGE_END = "ILLUMINA_READ_INSERT_SIZE_INSERT_RANGE_END"

    ILLUMINA_READ_INSERT_SIZE_90TH_PERC = "ILLUMINA_READ_INSERT_SIZE_90TH_PERC"
    ILLUMINA_READ_INSERT_SIZE_50TH_PERC = "ILLUMINA_READ_INSERT_SIZE_50TH_PERC"
    ILLUMINA_READ_INSERT_SIZE_10TH_PERC = "ILLUMINA_READ_INSERT_SIZE_10TH_PERC"
    ILLUMINA_READ_INSERT_SIZE_75TH_PERC = "ILLUMINA_READ_INSERT_SIZE_75TH_PERC"
    ILLUMINA_READ_INSERT_SIZE_25TH_PERC = "ILLUMINA_READ_INSERT_SIZE_25TH_PERC"

    ## step 19
    GC_DIVERGENCE_CSV_FILE = "GC_DIVERGENCE_CSV_FILE"
    GC_DIVERGENCE_PLOT_FILE = "GC_DIVERGENCE_PLOT_FILE"
    GC_DIVERGENCE_COEFFICIENTS_CSV_FILE = "GC_DIVERGENCE_COEFFICIENTS_CSV_FILE"

    #GC_DIVERGENCE_VAL = "GC_DIVERGENCE_VAL"
    GC_DIVERGENCE_COEFF_R1_AT = "GC_DIVERGENCE_COEFF_R1_AT"
    GC_DIVERGENCE_COEFF_R1_ATCG = "GC_DIVERGENCE_COEFF_R1_ATCG"
    GC_DIVERGENCE_COEFF_R1_CG = "GC_DIVERGENCE_COEFF_R1_CG"
    GC_DIVERGENCE_COEFF_R2_AT = "GC_DIVERGENCE_COEFF_R2_AT"
    GC_DIVERGENCE_COEFF_R2_ATCG = "GC_DIVERGENCE_COEFF_R2_ATCG"
    GC_DIVERGENCE_COEFF_R2_CG = "GC_DIVERGENCE_COEFF_R2_CG"


## EOF
