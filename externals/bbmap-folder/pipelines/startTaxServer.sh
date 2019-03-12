#!/bin/bash
#Run this on gpweb25

#curl http://localhost:3068/kill/TimeForANewServer
module unload oracle-jdk
module load oracle-jdk/1.7_64bit
module load pigz
sleep 2

LOG=taxlog35.txt
nohup /global/projectb/sandbox/gaag/bbtools/jgi-bbtools/taxserver.sh -da -Xmx31g port=3068 verbose accession=auto tree=auto table=auto size=auto img=auto pattern=auto prealloc domain=https://taxonomy.jgi-psf.org killcode=xxxx oldcode=xxxx oldaddress=https://taxonomy.jgi-psf.org/kill/ 1>>$LOG 2>&1 &

#simple mode, for testing:
#/global/projectb/sandbox/gaag/bbtools/jgi-bbtools/taxserver.sh -ea -Xmx8g port=3068 verbose accession=null tree=auto table=null
