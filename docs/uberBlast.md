COMMAND: 
usage: EToKi.py uberBlast [-h] -r REFERENCE -q QUERY [-o OUTPUT] [--blastn]
                          [--ublast] [--ublastSELF] [--minimap] [--minimapASM]
                          [--mmseq] [--min_id MIN_ID] [--min_cov MIN_COV]
                          [-s RE_SCORE] [-f] [--filter_cov FILTER_COV]
                          [--filter_score FILTER_SCORE] [-m]
                          [--merge_gap MERGE_GAP] [--merge_diff MERGE_DIFF]
                          [--return_overlap] [--overlap_length OVERLAP_LENGTH]
                          [--overlap_proportion OVERLAP_PROPORTION]
                          [-e FIX_END] [-t N_THREAD]

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
                        [OUTPUT] filename for output file. Default output to
                        sceen.
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
  --return_overlap      [DEFAULT: False] Report overlapped alignments
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

EXAMPLE:
python EToKi.py uberBlast -q examples/Escherichia.Achtman.alleles.fasta -r examples/GCF_001566635.1_ASM156663v1_genomic.fna -o examples/G749_7Gene.bsn --blastn --ublast --minimap --mmseq -s 2 -f
2019-03-11 13:38:12.919783      Run BLASTn starts
2019-03-11 13:38:17.447691      Run BLASTn finishes. Got 6568 alignments
2019-03-11 13:38:17.449272      Run uBLAST starts
2019-03-11 13:38:24.553244      Run uBLAST finishes. Got 8025 alignments
2019-03-11 13:38:24.559933      Run Minimap starts
2019-03-11 13:38:25.454961      Run Minimap finishes. Got 5529 alignments
2019-03-11 13:38:25.455421      Run MMSeqs starts
2019-03-11 13:38:38.872538      Run MMSeqs finishes. Got 10061 alignments
2019-03-11 13:38:39.714206      Update scores: 0 / 30183
2019-03-11 13:38:41.272581      Update scores: 10000 / 30183
2019-03-11 13:38:42.805440      Update scores: 20000 / 30183
2019-03-11 13:38:44.043799      Update scores: 30000 / 30183
2019-03-11 13:38:44.066854      Run filtering. Start with 30183 hits.
2019-03-11 13:38:44.540898      Done filtering. End with 9249 hits.

head examples/G749_7Gene.bsn
adk_1   NZ_CP014488.1   0.9606741573033708      536     0       0       1       536     485572  486107  0.0     891.0   536     4897758 536M    14593
adk_10  NZ_CP014488.1   0.9438202247191011      536     0       0       1       536     485572  486107  0.0     860.0   536     4897758 536M    14602
adk_100 NZ_CP014488.1   0.9606741573033708      536     0       0       1       536     485572  486107  0.0     887.0   536     4897758 536M    14691
adk_101 NZ_CP014488.1   0.898876404494382       536     0       0       1       536     485572  486107  0.0     813.0   536     4897758 536M    14692
adk_103 NZ_CP014488.1   0.9550561797752809      536     0       0       1       536     485572  486107  0.0     883.0   536     4897758 536M    14693
adk_108 NZ_CP014488.1   0.949438202247191       536     0       0       1       536     485572  486107  0.0     875.0   536     4897758 536M    14694
adk_109 NZ_CP014488.1   0.949438202247191       536     0       0       1       536     485572  486107  0.0     870.0   536     4897758 536M    14695
adk_11  NZ_CP014488.1   0.9550561797752809      536     0       0       1       536     485572  486107  0.0     886.0   536     4897758 536M    14603
adk_111 NZ_CP014488.1   0.9438202247191011      536     0       0       1       536     485572  486107  0.0     860.0   536     4897758 536M    14696
adk_112 NZ_CP014488.1   0.9438202247191011      536     0       0       1       536     485572  486107  0.0     869.0   536     4897758 536M    14697
