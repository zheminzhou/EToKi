#!/bin/bash
set -e
set -u

### Prepare reference alleles and a local database for 7 Gene MLST scheme
python EToKi.py MLSTdb -i examples/Escherichia.Achtman.alleles.fasta -r examples/Escherichia.Achtman.references.fasta -d examples/Escherichia.Achtman.convert.tab

### Calculate 7 Gene MLST genotype for a queried genome
gzip -cd examples/GCF_001566635.1_ASM156663v1_genomic.fna.gz > examples/GCF_001566635.1_ASM156663v1_genomic.fna && \
python EToKi.py MLSType -i examples/GCF_001566635.1_ASM156663v1_genomic.fna -r examples/Escherichia.Achtman.references.fasta -k G749 -o stdout -d examples/Escherichia.Achtman.convert.tab

### Run EBEis (EnteroBase Escherichia in silico Serotyping)
python EToKi.py EBEis -t Escherichia -q examples/GCF_000010485.1_ASM1048v1_genomic.fna -p SE15

### Cluster sequences into similarity-based groups 
python EToKi.py clust -p examples/Escherichia.Achtman.alleles_clust -i examples/Escherichia.Achtman.alleles.fasta -d 0.95 -c 0.95

### Do a joint BLASTn-like search using BLASTn, uSearch (uBLASTp), Mimimap and mmseqs
python EToKi.py uberBlast -q examples/Escherichia.Achtman.alleles.fasta -r examples/GCF_001566635.1_ASM156663v1_genomic.fna -o examples/G749_7Gene.bsn --blastn --ublast --minimap --mmseq -s 2 -f

### align multiple genomes onto one reference
python EToKi.py align -r GCF_000010485:examples/GCF_000010485.1_ASM1048v1_genomic.fna.gz -p examples/phylo_out \
GCF_000005845:examples/GCF_000005845.2_ASM584v2_genomic.fna.gz \
GCF_000214765:examples/GCF_000214765.2_ASM21476v3_genomic.fna.gz \
GCF_001566635:examples/GCF_001566635.1_ASM156663v1_genomic.fna.gz

### Build ML tree using RAxML and place all SNPs onto branches in the tree
cd examples && python ../EToKi.py phylo -t snp2mut -p phylo_out -s phylo_out.matrix.gz && cd ..

