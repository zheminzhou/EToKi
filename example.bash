### Trim genomic reads
python EToKi.py prepare --pe examples/A_R1.fastq.gz,examples/A_R2.fastq.gz -p examples/prep_out

### Merge and trim metagenomic reads
python EToKi.py prepare --pe examples/OAGR_ModernL7_10K_R1.fastq.gz,examples/OAGR_ModernL7_10K_R2.fastq.gz -p examples/QAGR_ModernL7_trim --noRename --merge

### Assemble genomic reads using SPAdes
python EToKi.py assemble --pe examples/prep_out_L1_R1.fastq.gz,examples/prep_out_L1_R2.fastq.gz --se examples/prep_out_L1_SE.fastq.gz -p examples/asm_out

### Assemble genomic reads using MEGAHIT
python EToKi.py assemble --pe examples/prep_out_L1_R1.fastq.gz,examples/prep_out_L1_R2.fastq.gz --se examples/prep_out_L1_SE.fastq.gz -p examples/asm_out2 --assembler megahit

### Find orthologous groups from public annotations
python EToKi.py ortho -P examples/GCF_000010485.combined.gff.gz --min_cds 60 --fast --incompleteCDS s -p examples/ST131 examples/*.gff.gz

### Prepare reference alleles and a local database for 7 Gene MLST scheme
python EToKi.py MLSTdb -i examples/Escherichia.Achtman.alleles.fasta -r examples/Escherichia.Achtman.references.fasta -d examples/Escherichia.Achtman.convert.tab

### Calculate 7 Gene MLST genotype for a queried genome
python EToKi.py MLSType -i examples/GCF_001566635.1_ASM156663v1_genomic.fna -r examples/Escherichia.Achtman.references.fasta -k G749 -o stdout -d examples/Escherichia.Achtman.convert.tab

### Construct HierCC (hierarchical clustering of cgMLST) for Yersinia cgMLST
python EToKi.py hierCC -p examples/Yersinia.cgMLST.profile.gz --o examples/Yersinia.cgMLST.hierCC

### Run EBEis (EnteroBase Escherichia in silico Serotyping)
python EToKi.py EBEis -t Escherichia -q examples/GCF_000010485.1_ASM1048v1_genomic.fna -p SE15

### Cluster sequences into similarity-based groups 
python EToKi.py clust -p examples/Escherichia.Achtman.alleles_clust -i examples/Escherichia.Achtman.alleles.fasta -d 0.95 -c 0.95

### Do a joint BLASTn-like search using BLASTn, uSearch (uBLASTp), Mimimap and mmseqs
python EToKi.py uberBlast -q examples/Escherichia.Achtman.alleles.fasta -r examples/GCF_001566635.1_ASM156663v1_genomic.fna -o examples/G749_7Gene.bsn --blastn --ublast --minimap --mmseq -s 2 -f

### Build ML tree using RAxML and place all SNPs onto branches in the tree
python EToKi.py phylo -t all -p phylo_out -m examples/phylo_rec.fasta

### Identify recombination sketches from the SNP matrix, and revise the branch lengths of the tree
python EToKi.py RecHMM -d phylo_out.mutations.gz -p examples/rec_out

### Strip out recombinant SNPs from a SNP matrix
python EToKi.py RecFilter -s phylo_out.matrix.gz -t phylo_out.labelled.nwk -r examples/rec_out.recombination.region -p examples/rec_out.nonrec
