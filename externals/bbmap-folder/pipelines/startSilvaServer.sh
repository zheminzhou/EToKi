#!/bin/bash
#Run this on gpweb34
#module unload oracle-jdk
#module load oracle-jdk/1.8_64bit

nohup /global/projectb/sandbox/gaag/bbtools/jgi-bbtools/taxserver.sh -da -Xmx10g port=3073 verbose tree=auto sketchonly index whitelist domain=https://ribo-sketch.jgi-psf.org killcode=xxxx oldcode=xxxx oldaddress=https://ribo-sketch.jgi-psf.org/kill/ ref=/global/projectb/sandbox/gaag/bbtools/silva/latest/both_seq#.sketch dbname=Silva blacklist=silva k=31,0 1>>silvalog22.o 2>&1 &

#simple mode, for testing:
#/global/projectb/sandbox/gaag/bbtools/jgi-bbtools/taxserver.sh -ea -Xmx10g port=3073 verbose tree=auto sketchonly silva k=31 index=f
