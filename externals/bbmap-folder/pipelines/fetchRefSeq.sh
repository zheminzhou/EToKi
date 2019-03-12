#!/bin/bash
set -e

#Written by Brian Bushnell
#Last updated July 19, 2018

#Fetches and renames RefSeq.
#Be sure the taxonomy server is updated first, or run with local taxonomy data!
#To use this script outside of NERSC when using local taxonomy data,
#add "taxpath=/path/to/taxonomy_directory/" to each BBTools command


#Ensure necessary executables are in your path
#module load bbtools
module load pigz

#Fetch RefSeq
#time wget -nv ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/*genomic.fna.gz
wget -q -O - ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/*genomic.fna.gz | gi2taxid.sh -Xmx1g in=stdin.fa.gz out=renamed.fa.gz pigz=32 unpigz zl=9 tree=null table=null accession=null server ow

#Sort by taxonomy.
#This makes sketching by taxa use much less memory because sketches can be written to disk as soon as they are finished.
time sortbyname.sh in=renamed.fa.gz memmult=0.33 out=sorted.fa.gz zl=8 pigz=64 taxa tree=auto gi=ignore fastawrap=1023 minlen=60 readbufferlen=2 readbuffers=1
