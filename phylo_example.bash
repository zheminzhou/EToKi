cd examples
# get phylogeny
python ../EToKi.py phylo -t all -p phylo_out -m phylo_rec.fasta
# run RecHMM
python ../EToKi.py RecHMM -d phylo_out.mutations.gz -p rec_out
# run RecFilter
python ../EToKi.py RecFilter -s phylo_out.matrix.gz -p rec_out -r rec_out.recombination.region -t phylo_out.labeled.nwk
# get phylogeny with filtered SNPs
python ../EToKi.py phylo -t phylogeny -p phyloFilter_out -s rec_out.filtered.gz
