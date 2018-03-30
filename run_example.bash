cd examples
python ../EToKi.py EnPhyl -t all -p sample_out -m sample.fasta
python ../EToKi.py RecHMM -d sample_out.mutations.gz -p sample_out
