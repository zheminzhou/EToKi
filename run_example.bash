cd examples
# run EnPhyl
python ../EToKi.py phylo -t all -p sample_out -m sample.fasta
# run RecHMM
python ../EToKi.py RecHMM -d sample_out.mutations.gz -p sample_out
