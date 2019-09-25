
# EToKi (Enterobase Tool Kit)
All methods related to Enterobase data analysis pipelines.

# INSTALLATION:

EToKi was developed and tested in both Python 2.7 and Python 3.5. EToKi depends on several Python libraries: 
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
EToKi also calls the following 3rd party programs for different pipelines:

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

NOTE: This has only been tested in Ubutu 16.06 but is expected to run on other 64-bit Linux systems. 
 
Usearch is a commercial program and allows free use of the 32-bit version for individuals. Please download it from [https://www.drive5.com/usearch/](https://www.drive5.com/usearch/)

After it is downloaded, pass its executable file to EToKi using **--usearch**

~~~~~~~~~~
python EToKi.py configure --usearch /path/to/usearch
~~~~~~~~~~

 You can also run both **--install** and **--usearch** at the same time:
~~~~~~~~~~
python EToKi.py configure --install --download_krakenDB --usearch /path/to/usearch
~~~~~~~~~~

Note that **--download_krakenDB** will download the minikraken2 database, which is about 8GB in size. Alternatively, you can use **--link_krakenDB** to pass a different Kraken database to EToKi.
~~~~~~~~~~
python EToKi.py configure --install --link_krakenDB /path/to/krakenDB --usearch /path/to/usearch
~~~~~~~~~~

You can also use pre-installed 3rd party programs in EToKi, by passing their absolute paths into the program using **--path**. This argument can be specified multiple times in the same command:
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
### Run EBEis (EnteroBase Escherichia in silico serotyping)
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
### Identify recombination stretches from the SNP matrix, and revise the branch lengths of a tree
~~~~~~~~~~~
python EToKi.py RecHMM -d phylo_out.mutations.gz -p examples/rec_out
~~~~~~~~~~~
### Strip out recombinant SNPs from a SNP matrix
~~~~~~~~~~~
python EToKi.py RecFilter -s phylo_out.matrix.gz -t phylo_out.labelled.nwk -r examples/rec_out.recombination.region -p examples/rec_out.nonrec
~~~~~~~~~~~

# USAGE:
The first argument passed into EToKi specifies the command to be called and the rest are the parameters for that command. To see all the commands available in EToKi, use
> python EToKi.py -h

And to see the parameters for an individual command, use:
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
**EToKi assemble** is a joint method for both *de novo* assembly and reference-guided assembly. 
* *de novo* assembly approach calls either SPAdes (default) or MEGAHIT (default for metagenomic data) on short reads that have been cleaned up using **EToKi prepare**, and uses Pilon to polish the assembled scaffolds and evaluate the reliability of consensus bases of the scaffolds. 

* Reference-guided assembly is also called "reference mapping". Short reads are aligned to a user-specified reference genome using minimap2. Nucleotide bases of the reference genome are updated using Pilon, according to the consensus base calls of the covered reads. Non-specific metagenomic reads of closely related species can sometimes also align to the reference genome and confuse consensus calling. Two arguments, **--outgroup** and **--ingroup**, are given to pre-filter these non-specific reads and obtain clean SNP calls. 
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
**EToKi ortho** has now been migrated to a [separate repository](https://github.com/zheminzhou/PEPPA) and renamed as **PEPPA**. 

## MLSTdb - Set up exemplar alleles and database for MLST schemes
**EToKi MLSTdb** converts existing allelic sequences into two files: (1) a multi-fasta file of exemplar allelic sequences and (2) a lookup table for the **EToKi MLSType** method. 
* The exemplar alleles are defined as: 
   1. Over 40% identity to the allelic sequences of a reference genome specified by **--refstrain**
   2. Less than 90% identity between different exemplar sequences of the same locus
   3. Identity to sequences of any different locus that is at least 10% less than the similarity to sequences of the same locus.
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
                        [DEFAULT: 0.4 ] Minimum identities between refstrain
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
**EToKi MLSType** identities allelic sequences in a queried genome, by comparing it with the exemplar alleles generated by **MLSTdb**. 
 ~~~~~~~~~~
usage: EToKi.py MLSType [-h] -i GENOME -r REFALLELE -k UNIQUE_KEY
                        [-d DATABASE] [-o OUTPUT] [-q] [-f] [-m MIN_IDEN]
                        [-p MIN_FRAG_PROP] [-l MIN_FRAG_LEN] [-x INTERGENIC]
                        [--overlap_prop OVERLAP_PROP]
                        [--overlap_iden OVERLAP_IDEN] [--max_dist MAX_DIST]
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
                        [DEFAULT: -1,-1 ] Call alleles in intergenic region if
                        the distance between two closely located loci fall
                        within the range defined by the two numbers. Suggest
                        to use 50,500. This is diabled by default with minus
                        numbers.
  --overlap_prop OVERLAP_PROP
                        [DEFAULT: 0.5 ] Given two hits, if <overlap_prop> of
                        their regions overlap, and the sequence identities of
                        one hits is <overlap_iden> lower than the other. The
                        hit with lower identities will be removed.
  --overlap_iden OVERLAP_IDEN
                        [DEFAULT: 0.05 ] Given two hits, if <overlap_prop> of
                        their regions overlap, and the sequence identities of
                        one hits is <overlap_iden> lower than the other. The
                        hit with lower identities will be removed.
  --max_dist MAX_DIST   [DEFAULT: 300 ] Consider two closely located hits as a
                        synteny block if their coordinates in both queried
                        genomes and reference gene are seperated by no more
                        than <max_dist> bps.
  --diag_diff DIAG_DIFF
                        [DEFAULT: 1.2 ] Consider two closely located hits as a
                        synteny block if, after merged, its covered region in
                        the queried genome is no more than <diag_diff> folds
                        of the region in the reference gene.
  --max_diff MAX_DIFF   [DEFAULT: 200 ] Consider two closely located hits as a
                        synteny block if, after merged, the lengths of its
                        covered regions in the queried genome and the
                        reference gene are differed by no more than <max_diff>
                        bps.
 ~~~~~~~~~~

## hierCC - generate hierarchical clusters from cgMLST profiles
Almost all STs called from cgMLST schemes contain some missing genes because they are called from draft genomes consisting of multiple contigs. As a result, almost every new genome results in a unique cgST number, many of which only differ from other cgSTs by missing data. It is therefore not obvious how to evaluate the genetic relationships of genomes using cgMLST profiles. **EToKi hierCC** offers hierarchical clustering of cgMLST profiles into thousands of different levels, using a single-linkage cluster algorithm. Clusters representing natural divisions of bacterial population structure can be extracted as sub-sets of the hierCC scheme later on. 
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


## RecHMM -  identify recombination stretches from a SNP matrix
**RecHMM** was initially implemented in R and described in supplementary material of [Zhou, Z. et al., PNAS, 2014](https://www.pnas.org/content/111/33/12199). RecHMM implements a modified version of Hidden Markov Model (HMM) to separate sites imported by recombination from those which are vertically inherited. 

Here we have implemented RecHMM in EToKi. The changes are:
1. The phylogeny and mutation placement module is separated into **EToKi phylo**. 
2. The code is re-implemented in Python Numpy and are much more efficient. 
3. The EM procedure is distributed across multiple processes. 
4. The HMM model is redesigned to include two extra states (4 in total). This new implementation can identify not only recombination from external sources, but also intra-population transfers. 
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
  --report, -r          Only report the model and do not calculate external stretches.
  --marginal MARGINAL, -M MARGINAL
                        Find recombinant regions using marginal likelihood rather than [DEFAULT] maximum likelihood method.
                        [DEFAULT] 0 to use Viterbi algorithm to find most likely path.
                         Otherwise (0, 1) use forward-backward algorithm, and report regions with >= M posterior likelihoods as recombinant stretches.
  --tree TREE, -T TREE  [INPUT, OPTIONAL] A labelled tree. Only used to generate corresponding mutational tree.
  --clean, -v           Do not show intermediate results during the iterations.
  --local_r LOCAL_R, -lr LOCAL_R
                        Specify a comma-delimited list of branches that share a different R/theta (Frequency of rec) ratio than the global consensus. Can be specified multiple times.
                        Use "*" to assign different value for each branch.
  --local_nu LOCAL_NU, -ln LOCAL_NU
                        Specify a comma-delimited list of branches that share a different Nu (SNP density in rec) than the global consensus. Can be specified multiple times.
                        Use "*" to assign different value for each branch.
  --local_delta LOCAL_DELTA, -ld LOCAL_DELTA
                        Specify a comma-delimited list of branches that share a different Delta (Length of rec stretches) than the global consensus. Can be specified multiple times.
                        Use "*" to assign different value for each branch.
~~~~~~~~~~~


## RecFilter - Remove recombinant regions from a SNP matrix
**EToKi RecFilter** automatically remove tracts of SNPs from an input matrix that are identified as having been imported by homologous recombinations. It currently supports as input outputs from RecHMM, [ClonalFrameML](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004041) and [Gubbins](https://academic.oup.com/nar/article/43/3/e15/2410982). It also accepts simulated results from [SimBac](https://www.ncbi.nlm.nih.gov/pubmed/27713837). 

Traditionally, SNPs affected by recombination are removed from the SNP matrix straight away. This is fine when recombination only affects a low proportion of the core genome. However, when a large proportion of the genome is affected by recombination, removing sites without counting the region as "missing" significantly shortens the branch lengths of a downstream phylogeny. **RecFilter** counts removed recombinant regions as missing data and increases the weights of the remaining mutational SNPs on the same branch, thus improves branch-length estimation. 
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
  --rec REC, -r REC     Recombinant stretches
  --prob PROB, -b PROB  Minimum probability for rec stretches. Default: 0.5
  --clonalframeml       The recombinant stretches are in ClonalFrameML format.
  --simbac              The recombinant stretches are in SimBac break format.
  --gubbins             convert VCF format to SNP matrix format.
~~~~~~~~~~~


## EB*Eis* - *in silico* serotype prediction for *Escherichia coli* & *Shigella spp.*
**EB*Eis*** is a BLASTn based prediction tool for the O and H antigens of *Escherichia coli* and *Shigella*. It uses essential genes (*wzx, wzy, wzt & wzm* for O; *fliC* for H) as markers. **EB*Eis*** uses a database built from two sources:
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
**EToKi uberBlast** is also internally called by **EToKi ortho** to align exemplar genes to queried genomes, using both BLASTn and uSearch-uBLASTp. Amino acid alignments are converted back to nucleotide sequences, meaning that genome coordinates remain consistent across different methods. 

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
  --blastn              Run BLASTn. Slowest. Good for identities between [80,
                        100]
  --ublast              Run uBLAST in tBLASTn mode. Fast. Good for identities
                        between [30-100]
  --ublastSELF          Run uBLAST in tBLASTn mode. Fast. Good for identities
                        between [30-100]
  --minimap             Run minimap. Fast. Good for identities between
                        [90-100]
  --minimapASM          Run minimap on assemblies. Fast. Good for identities
                        between [90-100]
  --mmseq               Run mmseq2 in tBLASTn mode. Fast. Good for identities
                        between [70-100]
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
**EToKi clust** is called internally by **EToKi ortho** to cluster seed genes into gene clusters. Given its linear-time complexity, it can cluster millions of gene sequences in minutes. 
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
