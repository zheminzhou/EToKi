## the sample reads were copied from samples in SPAdes 3.10

cd examples
# prepare reads
python ../EToKi.py prepare --pe A_R1.fastq.gz,A_R2.fastq.gz -p prep_out
# assemble reads
python ../EToKi.py assemble --pe prep_out_L1_R1.fastq.gz,prep_out_L1_R2.fastq.gz --se prep_out_L1_R3.fastq.gz -p asm_out
