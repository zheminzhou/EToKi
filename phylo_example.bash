cd examples
# get phylogeny
python ../EToKi.py phylo -t all -p phylo_out -m phylo_rec.fasta
# run RecHMM
python ../EToKi.py RecHMM -d phylo_out.mutations.gz -p rec_out
