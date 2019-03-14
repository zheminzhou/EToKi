COMMAND:
usage: EToKi.py prepare [-h] [--pe PE] [--se SE] [-p PREFIX] [-q READ_QUAL]
    [-b MAX_BASE] [-m MEMORY] [--noTrim] [--merge]
    [--noRename]

EToKi.py prepare
(1) Concatenates reads of the same library together.
(2) Trims sequences based on base-qualities.
(3) Removes potential adapters and barcodes.
(4) Limits total amount of reads to be used.
(5) Renames reads using their indexes.

optional arguments:
    -h, --help            show this help message and exit
    --pe PE               sets of PE reads from the same library, delimited by commas.
        e.g.
        If you have two sets of reads from a same library, it writes like:
        --pe a_R1.fq.gz,a_R2.fq.gz,b_R1.fq.gz,b_R2.fq.gz
        Specify this multiple times if there are multiple libraries.
    --se SE               sets of SE reads from the same library, delimited by commas.
        e.g.
        If you have two sets of reads from a same library, it writes like:
        --se a.fq.gz,b.fq.gz
        Specify this multiple times if there are multiple libraries.
    -p PREFIX, --prefix PREFIX
        prefix for the outputs. Default: EToKiPrepare
    -q READ_QUAL, --read_qual READ_QUAL
        Minimum quality for trimming. Default: 6
    -b MAX_BASE, --max_base MAX_BASE
        Total amount of bases (in BPs) to be kept. Set -1 to no restriction. Use ~100X coverage when you know the size of genome. Default: -1
    -m MEMORY, --memory MEMORY
        maximum amount of memory to be used in bbduk2. Default: 30g
    --noTrim              Do not do quality trim using bbduk2
    --merge               Try to merge PE reads by their overlaps
    --noRename            Do not rename reads

EXAMPLE 1:
python EToKi.py prepare --pe examples/A_R1.fastq.gz,examples/A_R2.fastq.gz -p examples/prep_out

OUTPUT 1:
2019-03-14 15:39:49.797692      Load in 2 read files from 1 libraries
2019-03-14 15:39:51.239857      Obtained 273596 bases in 3003 reads after Trimming in Lib 0
--se examples/prep_out_L1_SE.fastq.gz --pe examples/prep_out_L1_R1.fastq.gz,examples/prep_out_L1_R2.fastq.gz


EXAMPLE 2:
python EToKi.py prepare --pe examples/OAGR_ModernL7_10K_R1.fastq.gz,examples/OAGR_ModernL7_10K_R2.fastq.gz -p examples/QAGR_ModernL7_trim --noRename --merge

OUTPUT 2:
2019-03-11 21:37:22.361375      Load in 2 read files from 1 libraries
2019-03-11 21:37:29.218261      Obtained 1370354 bases in 19895 reads after Trimming in Lib 0
--se QAGR_Modern_L7_L1_MP.fastq.gz --se QAGR_Modern_L7_L1_SE.fastq.gz --pe QAGR_Modern_L7_L1_R1.fastq.gz,QAGR_Modern_L7_L1_R2.fastq.gz
