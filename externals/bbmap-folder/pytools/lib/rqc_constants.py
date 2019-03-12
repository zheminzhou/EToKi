#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Constants that are global for all RQC pipelines.

"""

from os import path as ospath
from common import get_run_path

class RQCConstants:
    UMASK = 0o022 # 755

# Exit codes
class RQCExitCodes:
    # used by run_blastplus_taxserver, run_diamond, readqc, pacbio, dapseq, assemblyqc
    JGI_SUCCESS = 0
    JGI_FAILURE = 1

# The locations of different directories are hard-coded
class RQCPaths:
    RQC_SW_BIN_DIR = ospath.join(get_run_path(), "tools")

# General commands
class RQCCommands:
    GNUPLOT_CMD = 'gnuplot'

    CAT_CMD = "cat"
    BZCAT_CMD = "bzcat"

    RM_CMD = "rm -f" # readqc utils
    SETE_CMD = "set -e" # readqc_utils

    # Read QC
    BBTOOLS_REFORMAT_CMD = "module load bbtools; reformat.sh"


# Blast specific constants
class RQCBlast:
    MEGAN_CMD = RQCPaths.RQC_SW_BIN_DIR + "/megan.pl" # used in assemblyqc.py and pacbioqc2_utils.py

# Reference databases
class RQCReferenceDatabases:
    ## NCBI databases
    ## NERSC's DB location: /scratch/blastdb/global/dna/shared/rqc/ref_databases/ncbi/CURRENT"

    NR                   = "/global/dna/shared/rqc/ref_databases/ncbi/CURRENT/nr/nr"
    REFSEQ_ARCHAEA       = "/global/dna/shared/rqc/ref_databases/ncbi/CURRENT/refseq.archaea/refseq.archaea"
    REFSEQ_BACTERIA      = "/global/dna/shared/rqc/ref_databases/ncbi/CURRENT/refseq.bacteria/refseq.bacteria"
    REFSEQ_FUNGI         = "/global/dna/shared/rqc/ref_databases/ncbi/CURRENT/refseq.fungi/refseq.fungi"
    REFSEQ_MITOCHONDRION = "/global/dna/shared/rqc/ref_databases/ncbi/CURRENT/refseq.mitochondrion/refseq.mitochondrion"
    REFSEQ_PLANT         = "/global/dna/shared/rqc/ref_databases/ncbi/CURRENT/refseq.plant/refseq.plant"
    REFSEQ_PLASMID       = "/global/dna/shared/rqc/ref_databases/ncbi/CURRENT/refseq.plasmid/refseq.plasmid"
    REFSEQ_PLASTID       = "/global/dna/shared/rqc/ref_databases/ncbi/CURRENT/refseq.plastid/refseq.plastid"
    REFSEQ_VIRAL         = "/global/dna/shared/rqc/ref_databases/ncbi/CURRENT/refseq.viral/refseq.viral"
    NT                   = "/global/dna/shared/rqc/ref_databases/ncbi/CURRENT/nt/nt"

    # should use this by default for NT
    NT_maskedYindexedN_BB = "/global/dna/shared/rqc/ref_databases/ncbi/CURRENT/nt/bbtools_dedupe_mask/nt_bbdedupe_bbmasked_formatted"
    SAGTAMINANTS = "/global/projectb/sandbox/rqc/qcdb/blast+/sagtaminants/sagtaminants.fa"
    GREEN_GENES  = "/global/dna/shared/rqc/ref_databases/misc/CURRENT/green_genes16s.insa_gg16S.fasta" # same file as above but has db index
    CONTAMINANTS = "/global/dna/shared/rqc/ref_databases/misc/CURRENT/JGIContaminants.fa" # april 2015
    COLLAB16S = "/global/projectb/sandbox/rqc/qcdb/collab16s/collab16s.fa" # recent
    LSU_REF = "/global/dna/shared/rqc/ref_databases/silva/CURRENT/LSURef_tax_silva"
    SSU_REF = "/global/dna/shared/rqc/ref_databases/silva/CURRENT/SSURef_tax_silva"
    LSSU_REF = "/global/dna/shared/rqc/ref_databases/silva/CURRENT/LSSURef_tax_silva"

    ## Alignment
    ARTIFACT_REF = "/global/dna/shared/rqc/ref_databases/qaqc/databases/Artifacts.adapters_primers_only.fa"

    ## PHiX Pipeline
    PHIX_REF        = "/global/dna/shared/rqc/ref_databases/qaqc/databases/phix174_ill.ref.fa"
    PHIX_CIRCLE_REF = "/global/dna/shared/rqc/ref_databases/qaqc/databases/phix174.ilmn_ref_concat.fa"
    FRAG_ADAPTERS   = "/global/projectb/sandbox/gaag/bbtools/data/adapters.fa"

## EOF
