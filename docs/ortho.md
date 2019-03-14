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


EXAMPLE:
python EToKi.py ortho -P examples/GCF_000010485.combined.gff.gz -p examples/ST131 examples/*.gff.gz


OUTPUT:
2019-03-12 22:56:34.820236      Run BLASTn starts
2019-03-12 22:56:43.781332      Run BLASTn finishes. Got 8852 alignments
2019-03-12 22:56:43.784246      Run uBLAST starts
2019-03-12 22:56:55.120453      Run uBLAST finishes. Got 9136 alignments
2019-03-12 22:56:56.786114      Update scores: 0 / 17988
2019-03-12 22:56:58.652004      Update scores: 10000 / 17988
2019-03-12 22:57:00.054938      Run filtering. Start with 17988 hits.
2019-03-12 22:57:00.438833      Done filtering. End with 9592 hits.
2019-03-12 22:57:01.234523      Run BLASTn starts
2019-03-12 22:57:01.247737      Run BLASTn starts
2019-03-12 22:57:01.264930      Run BLASTn starts
2019-03-12 22:57:01.292299      Run BLASTn starts
2019-03-12 22:57:17.119423      Run BLASTn finishes. Got 6417 alignments
2019-03-12 22:57:17.121188      Run uBLAST starts
2019-03-12 22:57:25.064429      Run BLASTn finishes. Got 7119 alignments
2019-03-12 22:57:25.067563      Run uBLAST starts
2019-03-12 22:57:26.635847      Run BLASTn finishes. Got 7491 alignments
2019-03-12 22:57:26.640165      Run uBLAST starts
2019-03-12 22:57:29.263940      Run BLASTn finishes. Got 7159 alignments
2019-03-12 22:57:29.268327      Run uBLAST starts
2019-03-12 22:57:46.429826      Run uBLAST finishes. Got 6673 alignments
2019-03-12 22:57:47.861054      Update scores: 0 / 13090
2019-03-12 22:57:49.606199      Update scores: 10000 / 13090
2019-03-12 22:57:50.156121      Run filtering. Start with 13090 hits.
2019-03-12 22:57:50.301241      Run uBLAST finishes. Got 7441 alignments
2019-03-12 22:57:50.363712      Done filtering. End with 7077 hits.
2019-03-12 22:57:50.366075      Start merging neighboring regions.
2019-03-12 22:57:51.719056      Finish merging neighboring regions.
2019-03-12 22:57:51.775346      Calculate overlaps.
2019-03-12 22:57:51.805216      Searching 0 / 7077 tabs
2019-03-12 22:57:52.394031      Identified 5019 overlaps.
2019-03-12 22:57:52.671129      Run uBLAST finishes. Got 7758 alignments
2019-03-12 22:57:53.429443      Update scores: 0 / 14560
2019-03-12 22:57:54.233304      Update scores: 0 / 15249
2019-03-12 22:57:55.005147      Run uBLAST finishes. Got 7457 alignments
2019-03-12 22:57:55.166076      Update scores: 10000 / 14560
2019-03-12 22:57:55.999315      Run filtering. Start with 14560 hits.
2019-03-12 22:57:56.148743      Update scores: 10000 / 15249
2019-03-12 22:57:56.228733      Done filtering. End with 7820 hits.
2019-03-12 22:57:56.231901      Start merging neighboring regions.
2019-03-12 22:57:57.097486      Update scores: 0 / 14616
2019-03-12 22:57:57.117133      Run filtering. Start with 15249 hits.
2019-03-12 22:57:57.354025      Done filtering. End with 8188 hits.
2019-03-12 22:57:57.357395      Start merging neighboring regions.
2019-03-12 22:57:57.816805      Finish merging neighboring regions.
2019-03-12 22:57:57.879342      Calculate overlaps.
2019-03-12 22:57:57.915496      Searching 0 / 7820 tabs
2019-03-12 22:57:58.520919      Identified 5646 overlaps.
2019-03-12 22:57:58.854594      Update scores: 10000 / 14616
2019-03-12 22:57:59.036556      Finish merging neighboring regions.
2019-03-12 22:57:59.112952      Calculate overlaps.
2019-03-12 22:57:59.151602      Searching 0 / 8188 tabs
2019-03-12 22:57:59.682447      Run filtering. Start with 14616 hits.
2019-03-12 22:57:59.864508      Identified 6017 overlaps.
2019-03-12 22:57:59.909194      Done filtering. End with 7838 hits.
2019-03-12 22:57:59.912038      Start merging neighboring regions.
2019-03-12 22:58:01.446068      Finish merging neighboring regions.
2019-03-12 22:58:01.505332      Calculate overlaps.
2019-03-12 22:58:01.541629      Searching 0 / 7838 tabs
2019-03-12 22:58:02.119294      Identified 5548 overlaps.
2019-03-12 22:58:34.417865      finding ANIs between genomes. 0/2940
2019-03-12 22:58:34.765815      finding ANIs between genomes. 100/2940
2019-03-12 22:58:35.088129      finding ANIs between genomes. 200/2940
2019-03-12 22:58:35.861494      finding ANIs between genomes. 300/2940
2019-03-12 22:58:36.208333      finding ANIs between genomes. 400/2940
2019-03-12 22:58:36.528049      finding ANIs between genomes. 500/2940
2019-03-12 22:58:36.850539      finding ANIs between genomes. 600/2940
2019-03-12 22:58:37.175355      finding ANIs between genomes. 700/2940
2019-03-12 22:58:37.523101      finding ANIs between genomes. 800/2940
2019-03-12 22:58:37.865397      finding ANIs between genomes. 900/2940
2019-03-12 22:58:38.184739      finding ANIs between genomes. 1000/2940
2019-03-12 22:58:38.528896      finding ANIs between genomes. 1100/2940
2019-03-12 22:58:38.892841      finding ANIs between genomes. 1200/2940
2019-03-12 22:58:39.261528      finding ANIs between genomes. 1300/2940
2019-03-12 22:58:39.658862      finding ANIs between genomes. 1400/2940
2019-03-12 22:58:40.047839      finding ANIs between genomes. 1500/2940
2019-03-12 22:58:40.492494      finding ANIs between genomes. 1600/2940
2019-03-12 22:58:40.854988      finding ANIs between genomes. 1700/2940
2019-03-12 22:58:41.198490      finding ANIs between genomes. 1800/2940
2019-03-12 22:58:41.550186      finding ANIs between genomes. 1900/2940
2019-03-12 22:58:41.871337      finding ANIs between genomes. 2000/2940
2019-03-12 22:58:42.221496      finding ANIs between genomes. 2100/2940
2019-03-12 22:58:42.543091      finding ANIs between genomes. 2200/2940
2019-03-12 22:58:42.862191      finding ANIs between genomes. 2300/2940
2019-03-12 22:58:43.221999      finding ANIs between genomes. 2400/2940
2019-03-12 22:58:43.590109      finding ANIs between genomes. 2500/2940
2019-03-12 22:58:43.918722      finding ANIs between genomes. 2600/2940
2019-03-12 22:58:44.227243      finding ANIs between genomes. 2700/2940
2019-03-12 22:58:44.506163      finding ANIs between genomes. 2800/2940
2019-03-12 22:58:44.750844      finding ANIs between genomes. 2900/2940
2019-03-12 22:58:56.193431      100 / 4053: pan gene "GCF_000010485:ECSF_RS00675" : "GCF_000010485:ECSF_RS00675" picked from rank 0 and score 10392.0
2019-03-12 22:58:56.686455      200 / 4053: pan gene "GCF_000010485:ECSF_RS22855" : "GCF_000010485:ECSF_RS22855" picked from rank 0 and score 8591.98324022
2019-03-12 22:58:57.076103      300 / 4053: pan gene "GCF_000010485:ECSF_RS00660" : "GCF_000010485:ECSF_RS00660" picked from rank 0 and score 7568.9952381
2019-03-12 22:58:57.455614      400 / 4053: pan gene "GCF_000010485:ECSF_RS10955" : "GCF_000010485:ECSF_RS10955" picked from rank 0 and score 6768.0
2019-03-12 22:58:57.902834      500 / 4053: pan gene "GCF_000010485:ECSF_RS15800" : "GCF_000010485:ECSF_RS15800" picked from rank 0 and score 6288.0
2019-03-12 22:58:58.249284      600 / 4053: pan gene "GCF_000010485:ECSF_RS00370" : "GCF_000010485:ECSF_RS00370" picked from rank 0 and score 5942.98181818
2019-03-12 22:58:58.675149      700 / 4053: pan gene "GCF_000010485:ECSF_RS10195" : "GCF_000010485:ECSF_RS10195" picked from rank 0 and score 5666.91147235
2019-03-12 22:58:59.046396      800 / 4053: pan gene "GCF_000010485:ECSF_RS01165" : "GCF_000010485:ECSF_RS01165" picked from rank 0 and score 5436.0
2019-03-12 22:58:59.502807      900 / 4053: pan gene "GCF_000010485:ECSF_RS15210" : "GCF_000010485:ECSF_RS15210" picked from rank 0 and score 5216.99308756
2019-03-12 22:58:59.884479      1000 / 4053: pan gene "GCF_000010485:ECSF_RS01295" : "GCF_000010485:ECSF_RS01295" picked from rank 0 and score 5016.0
2019-03-12 22:59:00.224274      1100 / 4053: pan gene "GCF_000010485:ECSF_RS10555" : "GCF_000010485:ECSF_RS10555" picked from rank 0 and score 4811.74961599
2019-03-12 22:59:00.520172      1200 / 4053: pan gene "GCF_000010485:ECSF_RS00920" : "GCF_000010485:ECSF_RS00920" picked from rank 0 and score 4632.0
2019-03-12 22:59:00.804831      1300 / 4053: pan gene "GCF_000010485:ECSF_RS01710" : "GCF_000010485:ECSF_RS01710" picked from rank 0 and score 4440.0
2019-03-12 22:59:01.077531      1400 / 4053: pan gene "GCF_000010485:ECSF_RS02900" : "GCF_000010485:ECSF_RS02900" picked from rank 0 and score 4236.0
2019-03-12 22:59:01.384687      1500 / 4053: pan gene "GCF_000010485:ECSF_RS07655" : "GCF_000010485:ECSF_RS07655" picked from rank 0 and score 4080.0
2019-03-12 22:59:01.722279      1600 / 4053: pan gene "GCF_000010485:ECSF_RS03220" : "GCF_000010485:ECSF_RS03220" picked from rank 0 and score 3948.0
2019-03-12 22:59:02.026207      1700 / 4053: pan gene "GCF_000010485:ECSF_RS03800" : "GCF_000010485:ECSF_RS03800" picked from rank 0 and score 3828.0
2019-03-12 22:59:02.362247      1800 / 4053: pan gene "GCF_000010485:ECSF_RS00720" : "GCF_000010485:ECSF_RS00720" picked from rank 0 and score 3708.0
2019-03-12 22:59:02.627106      1900 / 4053: pan gene "GCF_000010485:ECSF_RS00590" : "GCF_000010485:ECSF_RS00590" picked from rank 0 and score 3576.0
2019-03-12 22:59:02.887095      2000 / 4053: pan gene "GCF_000010485:ECSF_RS01565" : "GCF_000010485:ECSF_RS01565" picked from rank 0 and score 3456.0
2019-03-12 22:59:03.192145      2100 / 4053: pan gene "GCF_000010485:ECSF_RS00515" : "GCF_000010485:ECSF_RS00515" picked from rank 0 and score 3324.0
2019-03-12 22:59:03.444977      2200 / 4053: pan gene "GCF_000010485:ECSF_RS00865" : "GCF_000010485:ECSF_RS00865" picked from rank 0 and score 3192.0
2019-03-12 22:59:03.689575      2300 / 4053: pan gene "GCF_000010485:ECSF_RS01215" : "GCF_000010485:ECSF_RS01215" picked from rank 0 and score 3072.0
2019-03-12 22:59:03.933907      2400 / 4053: pan gene "GCF_000010485:ECSF_RS02795" : "GCF_000010485:ECSF_RS02795" picked from rank 0 and score 2988.0
2019-03-12 22:59:04.131175      2500 / 4053: pan gene "GCF_000010485:ECSF_RS03450" : "GCF_000010485:ECSF_RS03450" picked from rank 0 and score 2868.0
2019-03-12 22:59:04.324957      2600 / 4053: pan gene "GCF_000010485:ECSF_RS01875" : "GCF_000010485:ECSF_RS01875" picked from rank 0 and score 2760.0
2019-03-12 22:59:04.528909      2700 / 4053: pan gene "GCF_000010485:ECSF_RS03920" : "GCF_000010485:ECSF_RS03920" picked from rank 0 and score 2640.0
2019-03-12 22:59:04.771687      2800 / 4053: pan gene "GCF_000010485:ECSF_RS02375" : "GCF_000010485:ECSF_RS02375" picked from rank 0 and score 2508.0
2019-03-12 22:59:04.940249      2900 / 4053: pan gene "GCF_000010485:ECSF_RS01450" : "GCF_000010485:ECSF_RS01450" picked from rank 0 and score 2376.0
2019-03-12 22:59:05.100765      3000 / 4053: pan gene "GCF_000010485:ECSF_RS00055" : "GCF_000010485:ECSF_RS00055" picked from rank 0 and score 2268.0
2019-03-12 22:59:05.259951      3100 / 4053: pan gene "GCF_000010485:ECSF_RS00255" : "GCF_000010485:ECSF_RS00255" picked from rank 0 and score 2124.0
2019-03-12 22:59:05.412704      3200 / 4053: pan gene "GCF_000010485:ECSF_RS14875" : "GCF_000010485:ECSF_RS14875" picked from rank 0 and score 1796.89552239
2019-03-12 22:59:05.541675      3300 / 4053: pan gene "GCF_000010485:ECSF_RS22100" : "GCF_000010485:ECSF_RS22100" picked from rank 0 and score 564.0
2019-03-12 22:59:05.756434      3400 / 4053: pan gene "GCF_001566635:AVR74_RS28510" : "GCF_001566635:AVR74_RS28510" picked from rank 1 and score 3432.0
2019-03-12 22:59:05.953180      3500 / 4053: pan gene "GCF_000214765:ECNA114_RS26550" : "GCF_000214765:ECNA114_RS26550" picked from rank 1 and score 2325.0
2019-03-12 22:59:06.160713      3600 / 4053: pan gene "GCF_000214765:ECNA114_RS05940" : "GCF_000214765:ECNA114_RS05940" picked from rank 1 and score 1710.0
2019-03-12 22:59:06.302989      3700 / 4053: pan gene "GCF_001566635:AVR74_RS24900" : "GCF_001566635:AVR74_RS24900" picked from rank 1 and score 1242.0
2019-03-12 22:59:06.430804      3800 / 4053: pan gene "GCF_000214765:ECNA114_RS06680" : "GCF_000214765:ECNA114_RS06680" picked from rank 1 and score 957.0
2019-03-12 22:59:06.571282      3900 / 4053: pan gene "GCF_001566635:AVR74_RS15980" : "GCF_000214765:ECNA114_RS17090" picked from rank 1 and score 735.0
2019-03-12 22:59:06.751163      4000 / 4053: pan gene "GCF_001566635:AVR74_RS24970" : "GCF_001566635:AVR74_RS24970" picked from rank 1 and score 576.0
2019-03-12 22:59:16.947706      Pan genome annotations have been saved in examples/ST131.EToKi.gff
2019-03-12 22:59:16.947795      Gene allelic sequences have been saved in examples/ST131.allele.fna
XAMPLES:

