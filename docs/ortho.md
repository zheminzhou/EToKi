usage: EToKi.py ortho [-h] [-g GENES] [-P PRIORITY] [-p PREFIX] [-o ORTHOLOGY]
                      [-t N_THREAD] [--min_cds MIN_CDS]
                      [--clust_identity CLUST_IDENTITY]
                      [--clust_match_prop CLUST_MATCH_PROP]
                      [--match_identity MATCH_IDENTITY]
                      [--match_prop MATCH_PROP] [--match_len MATCH_LEN]
                      [--match_prop2 MATCH_PROP2] [--match_len2 MATCH_LEN2]
                      [--match_frag_prop MATCH_FRAG_PROP]
                      [--match_frag_len MATCH_FRAG_LEN]
                      [--synteny_gap SYNTENY_GAP]
                      [--synteny_diff SYNTENY_DIFF]
                      [--synteny_ovl_prop SYNTENY_OVL_PROP]
                      [--synteny_ovl_len SYNTENY_OVL_LEN]
                      [--edge_rescue EDGE_RESCUE]
                      [--mutation_variation MUTATION_VARIATION]
                      [--incompleteCDS] [--metagenome]
                      [--old_prediction OLD_PREDICTION] [--clust CLUST]
                      [--map_bsn MAP_BSN] [--self_bsn SELF_BSN]
                      [--conflicts CONFLICTS] [--global GLOBAL]
                      [--prediction PREDICTION]
                      [N [N ...]]

EToKi.py ortho
(1) Retieves genes and genomic sequences from GFF files and FASTA files.
(2) Groups genes into clusters using mmseq.
(3) Maps gene clusters back to genomes.
(4) Filters paralogous cluster alignments.
(5) identify a set of most probable non-overlapping orthologs.
(6) Re-annotate genomes using the new set of orthologs.

positional arguments:
  N                     GFF files containing both annotations and sequences.

optional arguments:
  -h, --help            show this help message and exit
  -g GENES, --genes GENES
                        Comma delimited files for additional genes.
  -P PRIORITY, --priority PRIORITY
                        Comma delimited filenames that contain highly confident genes.
  -p PREFIX, --prefix PREFIX
                        prefix for the outputs. Default: EToKi
  -o ORTHOLOGY, --orthology ORTHOLOGY
                        Method to define orthologous groups. nj [default], ml or rapid (for extremely large datasets)
  -t N_THREAD, --n_thread N_THREAD
                        Number of threads. Default: 30
  --min_cds MIN_CDS     Minimum length of a reference CDS. Default: 150.
  --clust_identity CLUST_IDENTITY
                        minimum identities in mmseq clusters. Default: 0.9
  --clust_match_prop CLUST_MATCH_PROP
                        minimum matches in mmseq clusters. Default: 0.9
  --match_identity MATCH_IDENTITY
                        minimum identities in BLAST search. Default: 0.6
  --match_prop MATCH_PROP
                        minimum match proportion for short genes in BLAST search. Default: 0.7
  --match_len MATCH_LEN
                        minimum match proportion for short genes in BLAST search. Default: 300
  --match_prop2 MATCH_PROP2
                        minimum match proportion for long genes in BLAST search. Default: 0.5
  --match_len2 MATCH_LEN2
                        minimum match proportion for long genes in BLAST search. Default: 500
  --match_frag_prop MATCH_FRAG_PROP
                        Min proportion of each fragment for fragmented matches. Default: 0.4
  --match_frag_len MATCH_FRAG_LEN
                        Min length of each fragment for fragmented matches. Default: 90
  --synteny_gap SYNTENY_GAP
                        Consider two fragmented matches within N bases as a synteny block. Default: 200
  --synteny_diff SYNTENY_DIFF
                        . Default: 1.2
  --synteny_ovl_prop SYNTENY_OVL_PROP
                        Max proportion of overlaps between two fragments in a synteny block. Default: 0.7
  --synteny_ovl_len SYNTENY_OVL_LEN
                        Max length of overlaps between two fragments in a synteny block. Default: 300
  --edge_rescue EDGE_RESCUE
                        Consider fragments that are within N bases of contig edges as part of a synteny block. Default: 150
  --mutation_variation MUTATION_VARIATION
                        Relative variation level in an ortholog group. Default: 2.
  --incompleteCDS       Do not do CDS checking for the reference genes. Default: False.
  --metagenome          Set to metagenome mode. equals to "--incompleteCDS --clust_identity 0.99 --clust_match_prop 0.8 --match_identity 0.98 --prefilter reference"
  --old_prediction OLD_PREDICTION
                        development param
  --clust CLUST         development param
  --map_bsn MAP_BSN     development param
  --self_bsn SELF_BSN   development param
  --conflicts CONFLICTS
                        development param
  --global GLOBAL       development param
  --prediction PREDICTION
                        development param

EXAMPLES:

