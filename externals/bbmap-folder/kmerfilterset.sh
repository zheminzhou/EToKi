#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 26, 2018

Description:  Generates a set of kmers such that every input sequence will
contain at least one kmer in the set.  This is a greedy algorithm which
retains the top X most common kmers each pass, and removes the sequences
matching those kmers, so subsequent passes are faster.

This will not generate an optimally small set but the output will be
quite small.  The file size may be further decreased with kcompress.sh.

Usage:  kmerfilterset.sh in=<input file> out=<output file> k=<integer>

File parameters:
in=<file>       Primary input.
out=<file>      Primary output.
temp=<file>     Temporary file pattern (optional).  Must contain # symbol.
initial=<file>  Initial kmer set (optional).  This can be used to accelerate
                the process.

Processing parameters:
k=31            Kmer length.
rcomp=t         Consider forward and reverse-complement kmers identical.
minkpp=1        (minkmersperpass) Retain at least this many kmers per pass.
                Higher is faster but results in a larger set.
maxkpp=2        (maxkmersperpass) Retain at most this many kmers per pass;
                0 means unlimited.
mincount=1      Ignore kmers seen fewer than this many times in this pass.
maxpasses=4000  Don't run more than this many passes.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                specify 200 megs. The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
-da             Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

#This block allows symlinked shellscripts to correctly set classpath.
pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx1000m"
z2="-Xms1000m"
EA="-ea"
EOOM=""
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	parseXmx "$@"
}
calcXmx "$@"

kmerlimit() {
	if [[ $SHIFTER_RUNTIME == 1 ]]; then
		#Ignore NERSC_HOST
		shifter=1
	elif [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_144_64bit
		module load pigz
	elif [[ $NERSC_HOST == denovo ]]; then
		module unload java
		module load java/1.8.0_144
		module load pigz
	elif [[ $NERSC_HOST == cori ]]; then
		module use /global/common/software/m342/nersc-builds/denovo/Modules/jgi
		module use /global/common/software/m342/nersc-builds/denovo/Modules/usg
		module unload java
		module load java/1.8.0_144
		module load pigz
	fi
	local CMD="java $EA $EOOM $z -cp $CP jgi.KmerFilterSetMaker $@"
	echo $CMD >&2
	eval $CMD
}

kmerlimit "$@"
