#!/bin/bash
set -e

#Written by Brian Bushnell
#Last updated May 31, 2018

#Sketches RefSeq.
#This script should be run after fetchRefSeq.sh.
#Be sure taxonomy is updated first!
#To use this script outside of NERSC,
#add "taxpath=/path/to/taxonomy_directory/" to each BBTools command


#Ensure necessary executables are in your path
#module load bbtools
#module load pigz

#Make a blacklist of kmers occuring in at least 250 different species.
time sketchblacklist.sh -Xmx63g in=sorted.fa.gz prepasses=1 tree=auto taxa taxlevel=species ow out=blacklist_refseq_species_250.sketch mincount=250 k=31,24 sizemult=2

#Generate 31 sketch files, with one sketch per species.
time bbsketch.sh -Xmx63g in=sorted.fa.gz out=taxa#.sketch mode=taxa tree=auto accession=null gi=null files=31 ow unpigz minsize=400 prefilter autosize blacklist=blacklist_refseq_species_250.sketch k=31,24 depth sizemult=2

#A query such as contigs.fa can now be compared to the new reference sketches like this:
#comparesketch.sh in=contigs.fa k=31,24 tree=auto taxa*.sketch blacklist=blacklist_refseq_species_250.sketch

#On NERSC systems, you can then set the default path to nt by pointing /global/projectb/sandbox/gaag/bbtools/refseq/current at the path to the new sketches.
#Then you can use the default set of refseq sketches like this:
#comparesketch.sh in=contigs.fa refseq tree=auto
#That command automatically adds the default path to the sketches, the blacklist, and the correct values for K.
