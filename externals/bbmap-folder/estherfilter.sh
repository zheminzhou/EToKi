#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified January 21, 2015

Description:  BLASTs queries against reference, and filters out hits with
              scores less than 'cutoff'.  The score is taken from column 12
              of the BLAST output.  The specific BLAST command is:
              blastall -p blastn -i QUERY -d REFERENCE -e 0.00001 -m 8

Usage:  estherfilter.sh <query> <reference> <cutoff>

For example:

estherfilter.sh reads.fasta genes.fasta 1000 > results.txt

'fasta' can be used as a fourth argument to get output in Fasta format.  Requires more memory.

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

z="-Xmx3200m"
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

estherfilter() {
	if [[ $SHIFTER_RUNTIME == 1 ]]; then
		#Ignore NERSC_HOST
		shifter=1
	elif [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_144_64bit
		module load blast
		module load pigz
	elif [[ $NERSC_HOST == denovo ]]; then
		module unload java
		module load java/1.8.0_144
		module load blast
		module load pigz
	elif [[ $NERSC_HOST == cori ]]; then
		module use /global/common/software/m342/nersc-builds/denovo/Modules/jgi
		module use /global/common/software/m342/nersc-builds/denovo/Modules/usg
		module unload java
		module load java/1.8.0_144
		module load blast
		module load pigz
	fi
	local CMD="java $EA $EOOM $z -cp $CP driver.EstherFilter $@"
	echo $CMD >&2
	eval $CMD
}

estherfilter "$@"
