#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified August 24, 2015

Description:  Filters a sam file to select only reads with substitution errors
for bases with quality scores in a certain interval.  Used for manually 
examining specific reads that may have incorrectly calibrated quality scores.

Usage:  filtersubs.sh in=<file> out=<file> minq=<number> maxq=<number>

Parameters:
in=<file>       Input sam or bam file.
out=<file>      Output file.
minq=0          Keep only reads with substitutions of at least this quality.
maxq=99         Keep only reads with substitutions of at most this quality.
countindels=t   Also keep reads with indels in the quality range. 
keepperfect=f   Also keep error-free reads.

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

z="-Xmx120m"
EA="-ea"
EOOM=""
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

filtersubs() {
	if [[ $SHIFTER_RUNTIME == 1 ]]; then
		#Ignore NERSC_HOST
		shifter=1
	elif [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_144_64bit
		module load samtools/1.4
		module load pigz
	elif [[ $NERSC_HOST == denovo ]]; then
		module unload java
		module load java/1.8.0_144
		module load PrgEnv-gnu/7.1
		module load samtools/1.4
		module load pigz
	elif [[ $NERSC_HOST == cori ]]; then
		module use /global/common/software/m342/nersc-builds/denovo/Modules/jgi
		module use /global/common/software/m342/nersc-builds/denovo/Modules/usg
		module unload java
		module load java/1.8.0_144
		module unload PrgEnv-intel
		module load PrgEnv-gnu/7.1
		module load samtools/1.4
		module load pigz
	fi
	local CMD="java $EA $EOOM $z -cp $CP jgi.FilterReadsWithSubs $@"
#	echo $CMD >&2
	eval $CMD
}

filtersubs "$@"
