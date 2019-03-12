#!/bin/bash
#Run this on gpweb34
#module unload oracle-jdk
#module load oracle-jdk/1.8_64bit

nohup /global/projectb/sandbox/gaag/bbtools/jgi-bbtools/taxserver.sh -da -Xmx28g prealloc=0.9 port=3072 verbose tree=auto sizemult=2 sketchonly index domain=https://refseq-sketch.jgi-psf.org killcode=xxxx oldcode=xxxx oldaddress=https://refseq-sketch.jgi-psf.org/kill/ RefSeq k=31,24 1>>refseqlog22.o 2>&1 &

#simple mode, for testing:
#nohup /global/projectb/sandbox/gaag/bbtools/jgi-bbtools/taxserver.sh -ea -Xmx28g port=3072 verbose tree=auto sizemult=2 sketchonly RefSeq k=31,24 index=f
