
# EToKi (Enterobase Tool Kit)
all methods related to Enterobase

# INSTALLATION:

EToKi was developed and tested in both python 2.7 and 3.5. EToKi depends on several python libraries: 
~~~~~~~~~~
ete3
numba
numpy
pandas
sklearn
~~~~~~~~~~

All libraries can be installed using pip: 

~~~~~~~~~~
pip install ete3 numba numpy pandas sklearn
~~~~~~~~~~
EToKi also calls many 3rd party programs for different pipelines. 

~~~~~~~~~~
raxml
fasttree
rapidnj
bbmap
mmseqs
ncbi-blast
usearch
spades
megahit
samtools
pilon
gatk
bwa
bowtie2
minimap2
kraken2 & minikraken2
pilercr
trf
~~~~~~~~~~

All 3rd party programs except for usearch can be automatically installed using *configure* command:
~~~~~~~~~~
python EToKi.py configure --install --download_krakenDB
~~~~~~~~~~

NOTE: This has only been tested in Ubutu 16.06 but is expected to be run in other 64-bit Linux systems. 
 
Usearch is a commercial program and allows free uses of 32-bit version for individuals. Please download it from [https://www.drive5.com/usearch/](https://www.drive5.com/usearch/)

After it is downloaded. pass its executable file to EToKi using **--usearch**

~~~~~~~~~~
python EToKi.py configure --usearch /path/to/usearch
~~~~~~~~~~

 You can also run both **--install** and **--usearch** at the same time:
~~~~~~~~~~
python EToKi.py configure --install --download_krakenDB --usearch /path/to/usearch
~~~~~~~~~~

Note that **--download_krakenDB** will download the minikraken2 database, which is about 8GB in size. Alternatively, you can use **--link_krakenDB** to pass a different database into EToKi
~~~~~~~~~~
python EToKi.py configure --install --link_krakenDB /path/to/krakenDB --usearch /path/to/usearch
~~~~~~~~~~

You can also use pre-installed 3rd party programs in EToKi, by passing their absolute paths into the program using **--path**. This argument can be specified multiple times:
~~~~~~~~~~
python EToKi.py configure --path fasttree=/path/to/fasttree --path raxml=/path/to/raxml
~~~~~~~~~~
  

# Quick Start (with examples)

### Trim genomic reads
~~~~~~~~~~~
python EToKi.py prepare --pe examples/A_R1.fastq.gz,examples/A_R2.fastq.gz -p examples/prep_out
~~~~~~~~~~~
### Merge and trim metagenomic reads
~~~~~~~~~~~
python EToKi.py prepare --pe examples/OAGR_ModernL7_10K_R1.fastq.gz,examples/OAGR_ModernL7_10K_R2.fastq.gz -p examples/QAGR_ModernL7_trim --noRename --merge
~~~~~~~~~~~
### Assemble genomic reads using SPAdes
~~~~~~~~~~~
python EToKi.py assemble --pe examples/prep_out_L1_R1.fastq.gz,examples/prep_out_L1_R2.fastq.gz --se examples/prep_out_L1_SE.fastq.gz -p examples/asm_out
~~~~~~~~~~~
### Assemble genomic reads using MEGAHIT
~~~~~~~~~~~
python EToKi.py assemble --pe examples/prep_out_L1_R1.fastq.gz,examples/prep_out_L1_R2.fastq.gz --se examples/prep_out_L1_SE.fastq.gz -p examples/asm_out2 --assembler megahit
~~~~~~~~~~~
### Find orthologous groups from public annotations
~~~~~~~~~~~
python EToKi.py ortho -P examples/GCF_000010485.combined.gff.gz -p examples/ST131 examples/*.gff.gz
~~~~~~~~~~~
### Prepare reference alleles and a local database for 7 Gene MLST scheme
~~~~~~~~~~~
python EToKi.py MLSTdb -i examples/Escherichia.Achtman.alleles.fasta -r examples/Escherichia.Achtman.references.fasta -d examples/Escherichia.Achtman.convert.tab
~~~~~~~~~~~
### Calculate 7 Gene MLST genotype for a queried genome
~~~~~~~~~~~
python EToKi.py MLSType -i examples/GCF_001566635.1_ASM156663v1_genomic.fna -r examples/Escherichia.Achtman.references.fasta -k G749 -o stdout -d examples/Escherichia.Achtman.convert.tab
~~~~~~~~~~~
### Construct HierCC (hierarchical clustering of cgMLST) for Yersinia cgMLST
~~~~~~~~~~~
python EToKi.py hierCC -p examples/Yersinia.cgMLST.profile.gz --o examples/Yersinia.cgMLST.hierCC
~~~~~~~~~~~
### Run EBEis (EnteroBase Escherichia in silico Serotyping)
~~~~~~~~~~~
python EToKi.py EBEis -t Escherichia -q examples/GCF_000010485.1_ASM1048v1_genomic.fna -p SE15
~~~~~~~~~~~
### Cluster sequences into similarity-based groups 
~~~~~~~~~~~
python EToKi.py clust -p examples/Escherichia.Achtman.alleles_clust -i examples/Escherichia.Achtman.alleles.fasta -d 0.95 -c 0.95
~~~~~~~~~~~
### Do a joint BLASTn-like search using BLASTn, uSearch (uBLASTp), Mimimap and mmseqs
~~~~~~~~~~~
python EToKi.py uberBlast -q examples/Escherichia.Achtman.alleles.fasta -r examples/GCF_001566635.1_ASM156663v1_genomic.fna -o examples/G749_7Gene.bsn --blastn --ublast --minimap --mmseq -s 2 -f
~~~~~~~~~~~
### Build ML tree using RAxML and place all SNPs onto branches in the tree
~~~~~~~~~~~
python EToKi.py phylo -t all -p phylo_out -m examples/phylo_rec.fasta
~~~~~~~~~~~
### Identify recombination sketches from the SNP matrix, and revise the branch lengths of the tree
~~~~~~~~~~~
python EToKi.py RecHMM -d phylo_out.mutations.gz -p examples/rec_out
~~~~~~~~~~~
### Strip out recombinant SNPs from a SNP matrix
~~~~~~~~~~~
python EToKi.py RecFilter -s phylo_out.matrix.gz -t phylo_out.labelled.nwk -r examples/rec_out.recombination.region -p examples/rec_out.nonrec
~~~~~~~~~~~

# USAGE:
The first argument passed into EToKi specifies the command to be called and the rests are the parameters for each command. To see all the commands available in EToKi, use
> python EToKi.py -h

And to see the parameters for each command, use:
> EToKi.py \<command\> -h

## configure - install and/or configure 3rd party programs
See the INSTALL section or the help page below.
~~~~~~~~~~~~~~
usage: EToKi.py configure [-h] [--install] [--usearch USEARCH]
                          [--download_krakenDB]
                          [--link_krakenDB KRAKEN_DATABASE] [--path PATH]

Install or modify 3rd party programs.

optional arguments:
  -h, --help            show this help message and exit
  --install             install 3rd party programs
  --usearch USEARCH     usearch is required for ortho and MLSType. A 32-bit
                        version of usearch can be downloaded from
                        https://www.drive5.com/usearch/.
  --download_krakenDB   When specified, miniKraken2 (8GB) will be downloaded
                        into the EToKi folder. You can also use
                        --link_krakenDB to use a pre-installed kraken2
                        database.
  --link_krakenDB KRAKEN_DATABASE
                        Kraken is optional in the assemble module. You can
                        specify your own database here
  --path PATH, -p PATH  Specify path to the 3rd party programs manually. format: 
                        <program>=<path>. This parameter can be specified
                        multiple times`
~~~~~~~~~~~~~~~~~

## prepare - trim, collapse, downsize and rename the short reads
~~~~~~~~~~~~~
usage: EToKi.py prepare [-h] [--pe PE] [--se SE] [-p PREFIX] [-q READ_QUAL]
                        [-b MAX_BASE] [-m MEMORY] [--noTrim] [--merge]
                        [--noRename]

EToKi.py prepare
(1) Concatenates reads of the same library together.
(2) Merge pair-end sequences for metagenomic reads (bbmap).
(3) Trims sequences based on base-qualities (bbduk).
(4) Removes potential adapters and barcodes (bbduk).
(5) Limits total amount of reads to be used.
(6) Renames reads using sequential numbers.

optional arguments:
  -h, --help            show this help message and exit
  --pe PE               comma delimited files of PE reads from the same library.
                        e.g. --pe a_R1.fq.gz,a_R2.fq.gz,b_R1.fq.gz,b_R2.fq.gz
                        This can be specified multiple times for different libraries.
  --se SE               comma delimited files of SE reads from the same library.
                        e.g. --se c_SE.fq.gz,d_SE.fq.gz
                        This can be specified multiple times for different libraries.
  -p PREFIX, --prefix PREFIX
                        prefix for the outputs. Default: EToKi_prepare
  -q READ_QUAL, --read_qual READ_QUAL
                        Minimum quality to be kept in bbduk. Default: 6
  -b MAX_BASE, --max_base MAX_BASE
                        Total amount of bases (in BPs) to be kept.
                        Default as -1 for no restriction.
                        Suggest to use ~100X coverage for de novo assembly.
  -m MEMORY, --memory MEMORY
                        maximum amount of memory to be used in bbduk. Default: 30g
  --noTrim              Do not do quality trim using bbduk
  --merge               Try to merge PE reads by their overlaps using bbmap
  --noRename            Do not rename reads
~~~~~~~~~~~~~~~~

## assemble - *de novo* or reference-guided assembly for genomic or metagenomic reads
**EToKi assemble** is a joint methods for both *de novo* assembly and reference-guided assembly. 
* *de novo* assembly approach calls either SPAdes (default) or MEGAHIT (default for metagenomic data) on the short reads that have been cleaned up using **EToKi prepare**, and uses Pilon to polish the assembled scaffolds and evaluate the reliability of consensus bases of the scaffolds. 

* Reference-guided assembly is also called as "reference mapping". The short reads are aligned onto a user-specified reference genome using minimap2. The nucleotide bases of the reference genome are updated using Pilon, according to the consensus base callings of the covered reads. Non-specific metagenomic reads of closely related species can sometimes also align onto the reference genome and confuse the consensus callings. Two arguments, **--outgroup** and **--ingroup**, are given to pre-filter these non-specific reads and obtain clean SNP callings. 
~~~~~~~~~~~~~~~~~
usage: EToKi.py assemble [-h] [--pe PE] [--se SE] [-p PREFIX] [-a ASSEMBLER]
                         [-k KMERS] [-m MAPPER] [-d MAX_DIFF] [-r REFERENCE]
                         [-i INGROUP] [-o OUTGROUP] [-S SNP] [-c CONT_DEPTH]
                         [--excluded EXCLUDED] [--metagenome] [--reassemble]
                         [--noPolish] [--onlySNP] [--noQuality] [--onlyEval]
                         [--kraken]

EToKi.py assemble
(1.1) Assembles short reads into assemblies, or
(1.2) Maps them onto a reference.
And
(2) Polishes consensus using polish,
(3) Removes low level contaminations.
(4) Estimates the base quality of the consensus.
(5) Predicts taxonomy using Kraken.

optional arguments:
  -h, --help            show this help message and exit
  --pe PE               comma delimited two files of PE reads.
  --se SE               a file of SE read.
  -p PREFIX, --prefix PREFIX
                        prefix for the outputs. Default: EToKi_assemble
  -a ASSEMBLER, --assembler ASSEMBLER
                        Assembler used for de novo assembly. 
                        Disabled if you specify a reference.
                        Default: spades for single colony isolates, megahit for metagenome
  -r REFERENCE, --reference REFERENCE
                        Reference for read mapping. Specify this for reference mapping module.
  -k KMERS, --kmers KMERS
                        relative lengths of kmers used in SPAdes. Default: 30,50,70,90
  -m MAPPER, --mapper MAPPER
                        aligner used for read mapping.
                        options are: miminap (default), bwa and bowtie2
  -d MAX_DIFF, --max_diff MAX_DIFF
                        Maximum proportion of variations allowed for a aligned reads. 
                        Default: 0.1 for single isolates, 0.05 for metagenome
  -i INGROUP, --ingroup INGROUP
                        Additional references presenting intra-population genetic diversities.
  -o OUTGROUP, --outgroup OUTGROUP
                        Additional references presenting genetic diversities outside of the studied population. 
                        Reads that are more similar to outgroups will be excluded from analysis.
  -S SNP, --SNP SNP     Exclusive set of SNPs. This will overwrite the polish process.
                        Required format:
                        <cont_name> <site> <base_type>
                        ...
  -c CONT_DEPTH, --cont_depth CONT_DEPTH
                        Allowed range of read depth variations relative to average value. 
                        Default: 0.2,2.5
                        Contigs with read depths outside of this range will be removed from the final assembly.
  --excluded EXCLUDED   A name of the file that contains reads to be excluded from the analysis.
  --metagenome          Reads are from metagenomic samples
                        this will set --cont_depth 0.0001,10000 --assembler megahit --max_diff 0.05
  --noPolish            Do not do PILON polish.
  --reassemble          Do local re-assembly in PILON
  --onlySNP             Only modify substitutions during the PILON polish.
  --noQuality           Do not estimate base qualities.
  --onlyEval            Do not run assembly/mapping. Only evaluate assembly status.
  --kraken              Run kmer based species prediction on contigs.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ortho  - pan-genome (and wgMLST scheme) prediction
**EToKi ortho** aligns reference genes onto genomic sequences, to generate similarity-based gene predictions on these genomes. Genes that are similar to the same reference gene are assigned as an orthologous group, which is further refined using a phylogeny-aware algorithm. The workflow is:
1. Retrieve **seed genes** from genome annotations of GFF format. These can be from published annotations or *ab initio* predictions using prodigal or prokka. 
2. Cluster seed genes into groups using MMseqs linclust. A collection of **exemplar genes** is assembled by selecting the centroid sequence of each group. 
3. exemplar genes are aligned onto queried genomes using both BLASTn and uSearch (optional), to identify homologous regions that are similar to each exemplar gene. 
4.  The sets of homologous regions with potential paralogs are identified if there are at duplicate matched within any single genome. These regions are iteratively sub-clustered based on phylogenetic topology. 
    * Firstly, each set of sets of homologous regions are aligned together. 
    * The resulting alignment are used to generate a neighbor-joining tree (RapidNJ; default) or Maximum likelihood tree (FastTree). 
    * The ETE3 package are used to bipartition the tree to maximise the nucleotide diversity (at least 5%) between the subtrees. Each of the resulted subtrees was evaluated iteratively until no two regions came from the same genome in the same subtree, or the maximum inter-subtree diversity is less than 5%. 
    * Then we replace the original set of homolog regions with all of its sub-trees.
6. After the division process, all the homolog sets were scored and ranked according to the summarised alignment scores of their homolog regions. Homolog sets were discarded if they had regions which overlapped with the regions within other sets that had greater scores.

~~~~~~~~~~~~~
usage: EToKi.py ortho [-h] [-g GENES] [-P PRIORITY] [-p PREFIX] [-o ORTHOLOGY]
                      [-t N_THREAD] [--min_cds MIN_CDS]
                      [--incompleteCDS INCOMPLETECDS]
                      [--clust_identity CLUST_IDENTITY]
                      [--clust_match_prop CLUST_MATCH_PROP] [--fast]
                      [--match_identity MATCH_IDENTITY]
                      [--match_prop MATCH_PROP] [--match_len MATCH_LEN]
                      [--match_prop1 MATCH_PROP1] [--match_len1 MATCH_LEN1]
                      [--match_prop2 MATCH_PROP2] [--match_len2 MATCH_LEN2]
                      [--match_frag_prop MATCH_FRAG_PROP]
                      [--match_frag_len MATCH_FRAG_LEN]
                      [--synteny_gap SYNTENY_GAP]
                      [--synteny_diff SYNTENY_DIFF]
                      [--allowed_variation ALLOWED_VARIATION] [--metagenome]
                      [N [N ...]]

EToKi.py ortho
(1) Retieves genes and genomic sequences from GFF files and FASTA files.
(2) Groups genes into clusters using mmseq.
(3) Maps gene clusters back to genomes.
(4) Discard paralogous alignments.
(5) Discard orthologous clusters if they had regions which overlapped with the regions within other sets that had greater scores.
(6) Re-annotate genomes using the remained of orthologs.

positional arguments:
  N                     GFF files containing both annotations and sequences.

optional arguments:
  -h, --help            show this help message and exit
  -g GENES, --genes GENES
                        Comma delimited files for additional genes.
  -P PRIORITY, --priority PRIORITY
                        Comma delimited, ordered list of filenames that contain genes with reliable starts and ends.
                        Genes listed in these files are preferred in all stages.
  -p PREFIX, --prefix PREFIX
                        prefix for the outputs. Default: EToKi_ortho
  -o ORTHOLOGY, --orthology ORTHOLOGY
                        Method to define orthologous groups.
                        nj [default], ml (for small dataset) or rapid (extremely large datasets)
  -t N_THREAD, --n_thread N_THREAD
                        Number of threads. Default: 30
  --min_cds MIN_CDS     Minimum length of a reference CDS. Default: 150.
  --incompleteCDS INCOMPLETECDS
                        Allowed types of imperfection for reference genes. Default: ''.
                        's': allows unrecognized start codon.
                        'e': allows unrecognized stop codon.
                        'i': allows stop codons in the coding region.
                        'f': allows frameshift in the coding region.
                        Multiple keywords can be used together. e.g., use 'sife' to allow random sequences.
  --clust_identity CLUST_IDENTITY
                        minimum identities of mmseqs clusters. Default: 0.95
  --clust_match_prop CLUST_MATCH_PROP
                        minimum matches in mmseqs clusters. Default: 0.85
  --fast                disable uBLAST search. Fast but less sensitive when nucleotide identities < 0.9
  --match_identity MATCH_IDENTITY
                        minimum identities in BLAST search. Default: 0.5
  --match_prop MATCH_PROP
                        minimum match proportion for normal genes in BLAST search. Default: 0.65
  --match_len MATCH_LEN
                        minimum match length for normal genes in BLAST search. Default: 300
  --match_prop1 MATCH_PROP1
                        minimum match proportion for short genes in BLAST search. Default: 0.8
  --match_len1 MATCH_LEN1
                        minimum match length for short genes in BLAST search. Default: 0
  --match_prop2 MATCH_PROP2
                        minimum match proportion for long genes in BLAST search. Default: 0.5
  --match_len2 MATCH_LEN2
                        minimum match length for long genes in BLAST search. Default: 500
  --match_frag_prop MATCH_FRAG_PROP
                        Min proportion of each fragment for fragmented matches. Default: 0.3
  --match_frag_len MATCH_FRAG_LEN
                        Min length of each fragment for fragmented matches. Default: 60
  --synteny_gap SYNTENY_GAP
                        Consider two fragmented matches within N bases as a synteny block. Default: 300
  --synteny_diff SYNTENY_DIFF
                        Form a synteny block when the covered regions in the reference gene
                        and the queried genome differed by no more than this value. Default: 1.2
  --allowed_variation ALLOWED_VARIATION
                        Allowed relative variation level compare to global.
                        The larger, the more variations are kept. Default: 1.
  --metagenome          Set to metagenome mode. equals to
                        "--fast --incompleteCDS sife --clust_identity 0.99 --clust_match_prop 0.8 --match_identity 0.98 --orthology rapid"
~~~~~~~~~~~~~

## MLSTdb - Set up exemplar alleles and database for MLST schemes
**EToKi MLSTdb** converts existing allelic sequences into two files: (1) multi-fasta file of exemplar allelic sequences and (2) a lookup table for **EToKi MLSType** method. 
* The exemplar alleles are defined as: 
   1. Over 40% identities to the allelic sequences of a reference genome specified by **--refstrain**
   2. Less than 90% identities between different exemplar sequences of the same locus
   3. Identities to sequences of any different locus are at least 10% less than the similarity to sequences of the same locus
~~~~~~~~~~~
usage: EToKi.py MLSTdb [-h] -i ALLELEFASTA [-r REFSET] [-d DATABASE]
                       [-s REFSTRAIN] [-x MAX_IDEN] [-m MIN_IDEN] [-p PARALOG]
                       [-c COVERAGE] [-e]

MLSTdb. Create reference sets of alleles for nomenclature.

optional arguments:
  -h, --help            show this help message and exit
  -i ALLELEFASTA, --input ALLELEFASTA
                        [REQUIRED] A single file contains all known alleles in
                        a MLST scheme.
  -r REFSET, --refset REFSET
                        [DEFAULT: No ref allele] Output - Reference alleles
                        used for MLSType.
  -d DATABASE, --database DATABASE
                        [DEFAULT: No allele DB] Output - A lookup table of all
                        alleles.
  -s REFSTRAIN, --refstrain REFSTRAIN
                        [DEFAULT: None] A single file contains alleles from
                        the reference genome.
  -x MAX_IDEN, --max_iden MAX_IDEN
                        [DEFAULT: 0.9 ] Maximum identities between resulting
                        refAlleles.
  -m MIN_IDEN, --min_iden MIN_IDEN
                        [DEFAULT: 0.6 ] Minimum identities between refstrain
                        and resulting refAlleles.
  -p PARALOG, --paralog PARALOG
                        [DEFAULT: 0.1 ] Minimum differences between difference
                        loci.
  -c COVERAGE, --coverage COVERAGE
                        [DEFAULT: 0.7 ] Proportion of aligned regions between
                        alleles.
  -e, --relaxEnd        [DEFAULT: False ] Allow changed ends (for pubmlst).
~~~~~~~~~~~

## MLSType - MLST nomenclature using a local set of references
**EToKi MLSType** identities allelic sequences in a queried genome, by comparing it with the exemplar alleles generated by **MLSTdb**. Allelic designations are also assigned when the sequences are known.

~~~~~~~~~~~
usage: EToKi.py MLSType [-h] -i GENOME -r REFALLELE -k UNIQUE_KEY
                        [-d DATABASE] [-o OUTPUT] [-q] [-f] [-m MIN_IDEN]
                        [-p MIN_FRAG_PROP] [-l MIN_FRAG_LEN] [-x INTERGENIC]
                        [--merging_prop MERGING_PROP]
                        [--merging_error MERGING_ERROR] [--max_dist MAX_DIST]
                        [--diag_diff DIAG_DIFF] [--max_diff MAX_DIFF]

MLSType. Find and designate MLST alleles from a queried assembly.

optional arguments:
  -h, --help            show this help message and exit
  -i GENOME, --genome GENOME
                        [REQUIRED] Input - filename for genomic assembly.
  -r REFALLELE, --refAllele REFALLELE
                        [REQUIRED] Input - fasta file for reference alleles.
  -k UNIQUE_KEY, --unique_key UNIQUE_KEY
                        [REQUIRED] An unique identifier for the assembly.
  -d DATABASE, --database DATABASE
                        [OPTIONAL] Input - lookup table of existing alleles.
  -o OUTPUT, --output OUTPUT
                        [DEFAULT: No output] Output - filename for the
                        generated alleles. Specify to STDOUT for screen
                        output.
  -q, --query_only      [DEFAULT: False] Do not submit new allele, only query.
  -f, --force           [DEFAULT: False] Force to accept low quality alleles.
  -m MIN_IDEN, --min_iden MIN_IDEN
                        [DEFAULT: 0.65 ] Minimum identities between refAllele
                        and genome.
  -p MIN_FRAG_PROP, --min_frag_prop MIN_FRAG_PROP
                        [DEFAULT: 0.6 ] Minimum covereage of a fragment.
  -l MIN_FRAG_LEN, --min_frag_len MIN_FRAG_LEN
                        [DEFAULT: 50 ] Minimum length of a fragment.
  -x INTERGENIC, --intergenic INTERGENIC
                        [DEFAULT: 50,500 ] Call alleles in intergenic region
                        if to closely located loci differ by a distance
                        between two numbers.
  --merging_prop MERGING_PROP
                        [DEFAULT: 0.5 ] Two hits are conflicted if they cover
                        by this proportion.
  --merging_error MERGING_ERROR
                        [DEFAULT: 0.05 ] Remove a secondary hit if its
                        similarity is lower than another overlapped region by
                        this value.
  --max_dist MAX_DIST   [DEFAULT: 300 ] Synteny block: Ignore if two
                        alignments seperate by at least this value.
  --diag_diff DIAG_DIFF
                        [DEFAULT: 1.2 ] Synteny block: Ignore if the lengths
                        of the resulted block differ by X fold between qry and
                        ref.
  --max_diff MAX_DIFF   [DEFAULT: 200 ] Synteny block: Ignore if the lengths
                        of the resulted block differ by this value between qry
                        and ref.
~~~~~~~~~~~

## hierCC - generate hierarchical clusters from cgMLST profiles
Almost all STs of cgMLST schemes contain some missing genes because they are called from draft genomes consisting of multiple contigs. As a result, almost every genome results in a unique cgST number, many of which only differ from other cgSTs by missing data. It is not obviously from the profiles to evaluate the genetic relationships of genomes. **EToKi hierCC** offers a finest clustering of cgMLST profiles into 1,000's of different levels, using single-linkage cluster algorithm. Clusters present natural populations can be extracted as sub-sets of the hierCC scheme later on. 
~~~~~~~~~~~
usage: EToKi.py hierCC [-h] -p PROFILE -o OUTPUT [-i INCREMENTAL] [-d DELTA]
                       [--immutable]

hierCC takes allelic profile (as in https://pubmlst.org/data/) and
work out specialised single linkage clustering result of all the profiles in the list.

optional arguments:
  -h, --help            show this help message and exit
  -p PROFILE, --profile PROFILE
                        [INPUT; REQUIRED] name of the profile file. Can be GZIPed.
  -o OUTPUT, --output OUTPUT
                        [OUTPUT; REQUIRED] Prefix for the output files. These include a NUMPY and TEXT verions of the same clustering result
  -i INCREMENTAL, --incremental INCREMENTAL
                        [INPUT; optional] The NUMPY version of an old clustering result
  -d DELTA, --delta DELTA
                        [optional] comma delimited list of threshold (delta). All values are included by default.
  --immutable           [optional] Use a immutable clustering system. The designations of old profiles are immutable. Faster but leads to non-optimal assignment.
~~~~~~~~~~~

## align - align multiple queried genomes to a single reference
~~~~~~~~~~~
usage: EToKi.py align [-h] -r REFERENCE [-p PREFIX] [-a] [-m] [-c CORE]
                      [-n N_PROC]
                      queries [queries ...]

Align multiple genomes onto a single reference.

positional arguments:
  queries               queried genomes. Use <Tag>:<Filename> format to feed
                        in a tag for each genome. Otherwise filenames will be
                        used as tags for genomes.

optional arguments:
  -h, --help            show this help message and exit
  -r REFERENCE, --reference REFERENCE
                        [REQUIRED; INPUT] reference genomes to be aligned
                        against. Use <Tag>:<Filename> format to assign a tag
                        to the reference.
  -p PREFIX, --prefix PREFIX
                        [OUTPUT] prefix for all outputs.
  -a, --alignment       [OUTPUT] Generate core genomic alignments in FASTA
                        format
  -m, --matrix          [OUTPUT] Do not generate core SNP matrix
  -c CORE, --core CORE  [PARAM] percentage of presences for core genome.
                        [DEFAULT: 0.95]
  -n N_PROC, --n_proc N_PROC
                        [PARAM] number of processes to use. [DEFAULT: 5]
~~~~~~~~~~~

## phylo - infer phylogeny and ancestral states from genomic alignments 
~~~~~~~~~~~
usage: EToKi.py phylo [-h] [--tasks TASKS] --prefix PREFIX
                      [--alignment ALIGNMENT] [--snp SNP] [--tree TREE]
                      [--ancestral ANCESTRAL] [--core CORE] [--n_proc N_PROC]

EToKi phylo runs to:
(1) Generate SNP matrix from alignment (-t matrix)
(2) Calculate ML phylogeny from SNP matrix using RAxML (-t phylogeny)
(3) Workout the nucleotide sequences of internal nodes in the tree using ML estimation (-t ancestral or -t ancestral_proportion for ratio frequencies)
(4) Place mutations onto branches of the tree (-t mutation)

optional arguments:
  -h, --help            show this help message and exit
  --tasks TASKS, -t TASKS
                        Tasks to call. Allowed tasks are:
                        matrix: generate SNP matrix from alignment.
                        phylogeny: generate phylogeny from SNP matrix.
                        ancestral: generate AS (ancestral state) matrix from SNP matrix and phylogeny
                        ancestral_proportion: generate possibilities of AS for each site
                        mutation: assign SNPs into branches from AS matrix

                        You can run multiple tasks by sending a comma delimited task list.
                        There are also some pre-defined task combo:
                        all: matrix,phylogeny,ancestral,mutation
                        aln2phy: matrix,phylogeny [default]
                        snp2anc: phylogeny,ancestral
                        mat2mut: ancestral,mutation
  --prefix PREFIX, -p PREFIX
                        prefix for all outputs.
  --alignment ALIGNMENT, -m ALIGNMENT
                        aligned sequences in either fasta format or Xmfa format. Required for "matrix" task.
  --snp SNP, -s SNP     SNP matrix in specified format. Required for "phylogeny" and "ancestral" if alignment is not given
  --tree TREE, -z TREE  phylogenetic tree. Required for "ancestral" task
  --ancestral ANCESTRAL, -a ANCESTRAL
                        Inferred ancestral states in a specified format. Required for "mutation" task
  --core CORE, -c CORE  Core genome proportion. Default: 0.95
  --n_proc N_PROC, -n N_PROC
                        Number of processes. Default: 7.
~~~~~~~~~~~


## RecHMM -  identify recombination sketches from a SNP matrix
**RecHMM** was initially implemented in R and described in supplementary material of [Zhou, Z. et al., PNAS, 2014](https://www.pnas.org/content/111/33/12199). RecHMM implements a modified version of Hidden Markov Model (HMM) to separates sites imported by recombination from those vertically inherited. 

Here is a 2nd version of RecHMM implemented in EToKi package. The changes are:
1. The phylogeny and mutation placement module is separated into **EToKi phylo**. 
2. The codes are re-implemented in Python Numpy and much more efficient. 
3. The EM procedure is distributed into multiple processes. 
4. HMM model is redesigned to include two extra states (4 in total). The new implementation can identify not only recombination from external sources, but also intra-population transfers. 
~~~~~~~~~~~
usage: EToKi.py RecHMM [-h] --data DATA [--model MODEL] [--task TASK]
                       [--init INIT] [--prefix PREFIX] [--cool_down COOL_DOWN]
                       [--n_proc N_PROC] [--bootstrap BOOTSTRAP] [--report]
                       [--marginal MARGINAL] [--tree TREE] [--clean]
                       [--local_r LOCAL_R] [--local_nu LOCAL_NU]
                       [--local_delta LOCAL_DELTA]

Parameters for RecHMM.

optional arguments:
  -h, --help            show this help message and exit
  --data DATA, -d DATA  A list of mutations generated by EnPhyl
  --model MODEL, -m MODEL
                        Read a saved best model.
  --task TASK, -t TASK  task to run.
                        0: One rec category from external sources.
                        1: Three rec categories considering internal, external and mixed sources [default].
  --init INIT, -i INIT  Initiate models with guesses of recombinant proportions.
                        Default: 0.05,0.5,0.95
  --prefix PREFIX, -p PREFIX
                        Prefix for all the outputs
  --cool_down COOL_DOWN, -c COOL_DOWN
                        Delete the worst model every N iteration. Default:5
  --n_proc N_PROC, -n N_PROC
                        Number of processes. Default: 5.
  --bootstrap BOOTSTRAP, -b BOOTSTRAP
                        Number of Randomizations for confidence intervals.
                        Default: 1000.
  --report, -r          Only report the model and do not calculate external sketches.
  --marginal MARGINAL, -M MARGINAL
                        Find recombinant regions using marginal likelihood rather than [DEFAULT] maximum likelihood method.
                        [DEFAULT] 0 to use Viterbi algorithm to find most likely path.
                         Otherwise (0, 1) use forward-backward algorithm, and report regions with >= M posterior likelihoods as recombinant sketches.
  --tree TREE, -T TREE  [INPUT, OPTIONAL] A labelled tree. Only used to generate corresponding mutational tree.
  --clean, -v           Do not show intermediate results during the iterations.
  --local_r LOCAL_R, -lr LOCAL_R
                        Specify a comma-delimited list of branches that share a different R/theta (Frequency of rec) ratio than the global consensus. Can be specified multiple times.
                        Use "*" to assign different value for each branch.
  --local_nu LOCAL_NU, -ln LOCAL_NU
                        Specify a comma-delimited list of branches that share a different Nu (SNP density in rec) than the global consensus. Can be specified multiple times.
                        Use "*" to assign different value for each branch.
  --local_delta LOCAL_DELTA, -ld LOCAL_DELTA
                        Specify a comma-delimited list of branches that share a different Delta (Length of rec sketches) than the global consensus. Can be specified multiple times.
                        Use "*" to assign different value for each branch.
~~~~~~~~~~~


## RecFilter - Remove recombination sketches from a SNP matrix
**EToKi RecFilter** automatically remove SNPs from an input matrix that were brought in by homologous recombinations, according to the results of recombination detection tool. It currently supports outputs of RecHMM, [ClonalFrameML](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004041) and [Gubbins](https://academic.oup.com/nar/article/43/3/e15/2410982). It also accepts simulated results from [SimBac](https://www.ncbi.nlm.nih.gov/pubmed/27713837). 

Traditionally, SNPs affected by recombination are removed from the SNP matrix straightly. This is fine when recombination only covers low proprotion of the core genome. However, when large proportion of the genome is covered by recombination, removing without counting the region as "missing" significantly shortens the branch lengths of the tree. **RecFilter** counts the removed recombinant sketches as missing data, and increases the weights of the remaining mutational SNPs of the same branch, and thus gets better estimation of its length. 
~~~~~~~~~~~
usage: EToKi.py RecFilter [-h] --prefix PREFIX --snp MATRIX --tree TREE --rec
                          REC [--prob PROB] [--clonalframeml] [--simbac]
                          [--gubbins]

Generate a matrix of only vertically inherited SNPs.

optional arguments:
  -h, --help            show this help message and exit
  --prefix PREFIX, -p PREFIX
                        prefix for the output
  --snp MATRIX, -s MATRIX
                        SNP matrix
  --tree TREE, -t TREE  Labeled tree
  --rec REC, -r REC     Recombinant sketches
  --prob PROB, -b PROB  Minimum probability for rec sketches. Default: 0.5
  --clonalframeml       The recombinant sketches are in ClonalFrameML format.
  --simbac              The recombinant sketches are in SimBac break format.
  --gubbins             convert VCF format to SNP matrix format.
~~~~~~~~~~~


## EB*Eis* - *in silico* serotype prediction for *Escherichia coli* & *Shigella spp.*
**EB*Eis*** is a BLASTn based prediction tool for O and H antigens of *Escherichia coli* and *Shigella*. It uses essential genes (*wzx, wzy, wzt & wzm* for O; *fliC* for H) as markers. The database its using consists of two sources:
1. [SeroTypeFinder ](https://bitbucket.org/genomicepidemiology/serotypefinder_db/src)
2. O-antigen gene sequences reported in [DebRoy et al., PLoS ONE, 2016](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0147434#pone.0147434.ref011)
~~~~~~~~~~~
usage: EToKi.py EBEis [-h] -q QUERY [-t TAXON] [-p PREFIX]

EnteroBase Escherichia in silico serotyping

optional arguments:
  -h, --help            show this help message and exit
  -q QUERY, --query QUERY
                        file name for the queried assembly in multi-FASTA format.
  -t TAXON, --taxon TAXON
                        Taxon database to compare with. 
                        Only support Escherichia (default) for the moment.
  -p PREFIX, --prefix PREFIX
                        prefix for intermediate files. Default: EBEis
~~~~~~~~~~~

## isCRISPOL - *in silico* prediction of CRISPOL array for *Salmonella enterica* serovar Typhimurium
CRISPOL is an oligo based Typhimurium sub-typing method described in ([Fabre et al., PLoS ONE, 2012](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0036995)). We use the direct repeats (DRs) and spacers in the Typhimurium CPRISR array to predict CRISPOL types from genomic assemblies.
~~~~~~~~~~~
usage: EToKi.py isCRISPOL [-h] [N [N ...]]

in silico Typhimurium subtyping using CRISPOL scheme (Fabre et al., PLoS ONE, 2012)

positional arguments:
  N           FASTA files containing assemblies of S. enterica Typhimurium.

optional arguments:
  -h, --help  show this help message and exit
~~~~~~~~~~~

## uberBlast - Use BLASTn, uBLASTp, minimap2 and/or mmseqs to identify similar sequences
**EToKi uberBlast** is also internally called by **EToKi ortho** to align exemplar genes onto queried genomes, using both BLASTn and uSearch-uBLASTp. Amino acid alignments are converted back to its original nucleotide sequences, such that the coordinates are consistent across different methods. 

* minimap2 --- Fastest alignment in nucleotide level. High accuracy in identities >= 90%, but lose sensitivity quickly for lower identities. 
* blastn --- Fast alignment in nucleotide level.  Lose sensitivity for identities < 80%
* mmseqs --- Amino acid based alignment for identities >= 70% (open source)
* uBLASTp --- Amino acid based alignment for identities < 50% (commercial software)
~~~~~~~~~~~
usage: EToKi.py uberBlast [-h] -r REFERENCE -q QUERY [-o OUTPUT] [--blastn]
                          [--ublast] [--ublastSELF] [--minimap] [--minimapASM]
                          [--mmseq] [--min_id MIN_ID] [--min_cov MIN_COV]
                          [--min_ratio MIN_RATIO] [-s RE_SCORE] [-f]
                          [--filter_cov FILTER_COV]
                          [--filter_score FILTER_SCORE] [-m]
                          [--merge_gap MERGE_GAP] [--merge_diff MERGE_DIFF]
                          [-O] [--overlap_length OVERLAP_LENGTH]
                          [--overlap_proportion OVERLAP_PROPORTION]
                          [-e FIX_END] [-t N_THREAD] [-p]

Five different alignment methods.

optional arguments:
  -h, --help            show this help message and exit
  -r REFERENCE, --reference REFERENCE
                        [INPUT; REQUIRED] filename for the reference. This is
                        normally a genomic assembly.
  -q QUERY, --query QUERY
                        [INPUT; REQUIRED] filename for the query. This can be
                        short-reads or genes or genomic assemblies.
  -o OUTPUT, --output OUTPUT
                        [OUTPUT; Default: None] save result to a file or to
                        screen (stdout). Default do nothing.
  --blastn              Run BLASTn. Slowest. Good for identities between [70,
                        100]
  --ublast              Run uBLAST on tBLASTn mode. Fast. Good for identities
                        between [30-100]
  --ublastSELF          Run uBLAST on tBLASTn mode. Fast. Good for identities
                        between [30-100]
  --minimap             Run minimap. Fast. Good for identities between
                        [90-100]
  --minimapASM          Run minimap on assemblies. Fast. Good for identities
                        between [90-100]
  --mmseq               Run mmseq2 on tBLASTn mode. Fast. Good for identities
                        between [60-100]
  --min_id MIN_ID       [DEFAULT: 0.3] Minimum identity before reScore for an
                        alignment to be kept
  --min_cov MIN_COV     [DEFAULT: 40] Minimum length for an alignment to be
                        kept
  --min_ratio MIN_RATIO
                        [DEFAULT: 0.05] Minimum length for an alignment to be
                        kept, proportional to the length of the query
  -s RE_SCORE, --re_score RE_SCORE
                        [DEFAULT: 0] Re-interpret alignment scores and
                        identities. 0: No rescore; 1: Rescore with
                        nucleotides; 2: Rescore with amino acid; 3: Rescore
                        with codons
  -f, --filter          [DEFAULT: False] Remove secondary alignments if they
                        overlap with any other regions
  --filter_cov FILTER_COV
                        [DEFAULT: 0.9]
  --filter_score FILTER_SCORE
                        [DEFAULT: 0]
  -m, --linear_merge    [DEFAULT: False] Merge consective alignments
  --merge_gap MERGE_GAP
                        [DEFAULT: 300]
  --merge_diff MERGE_DIFF
                        [DEFAULT: 1.2]
  -O, --return_overlap  [DEFAULT: False] Report overlapped alignments
  --overlap_length OVERLAP_LENGTH
                        [DEFAULT: 300] Minimum overlap to report
  --overlap_proportion OVERLAP_PROPORTION
                        [DEFAULT: 0.6] Minimum overlap proportion to report
  -e FIX_END, --fix_end FIX_END
                        [FORMAT: L,R; DEFAULT: 0,0] Extend alignment to the
                        edges if the un-aligned regions are <= [L,R]
                        basepairs.
  -t N_THREAD, --n_thread N_THREAD
                        [DEFAULT: 8] Number of threads to use.
  -p, --process         [DEFAULT: False] Use processes instead of threads.
~~~~~~~~~~~

## clust - linear-time clustering of short sequences using mmseqs linclust
**EToKi clust** is also called by **EToKi ortho** internally to cluster seed genes into gene clusters. Given its linear-time complexity, it can cluster millions of gene sequences in minutes. 
~~~~~~~~~~~
usage: EToKi.py clust [-h] -i INPUT -p PREFIX [-d IDENTITY] [-c COVERAGE]
                      [-t N_THREAD]

Get clusters and exemplars of clusters from gene sequences using mmseqs linclust.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        [INPUT; REQUIRED] name of the file containing gene sequneces in FASTA format.
  -p PREFIX, --prefix PREFIX
                        [OUTPUT; REQUIRED] prefix of the outputs.
  -d IDENTITY, --identity IDENTITY
                        [PARAM; DEFAULT: 0.9] minimum intra-cluster identity.
  -c COVERAGE, --coverage COVERAGE
                        [PARAM; DEFAULT: 0.9] minimum intra-cluster coverage.
  -t N_THREAD, --n_thread N_THREAD
                        [PARAM; DEFAULT: 8] number of threads to use.
~~~~~~~~~~~
