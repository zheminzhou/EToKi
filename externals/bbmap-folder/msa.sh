#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 15, 2018

Description:  Aligns a query sequence to reference sequences.
Outputs the best matching position per reference sequence.
If there are multiple queries, only the best-matching query will be used.
MSA in this context stands for MultiStateAligner, not Multiple Sequence Alignment.

Usage:
msa.sh in=<file> out=<file> literal=<literal,literal,...>
or
msa.sh in=<file> out=<file> ref=<lfile>

Parameters:
in=<file>       File containing reads.
out=<file>      Sam output file.
literal=        A sequence of bases to match, or a comma-delimited list.
ref=<file>      A fasta file of bases to match.  Please set either ref
                or literal, not both.
rcomp=t         Also look for reverse-complements of the sequences.
replicate=t     Make copies of sequences with undefined bases for every
                possible combination.  For example, ATN would expand to
                ATA, ATC, ATG, and ATT.
cutoff=0        Ignore alignments with identity below this (range 0-1).

Java Parameters:
-Xmx            This will set Java's memory usage, overriding automatic
                memory detection. -Xmx20g will specify 
                20 gigs of RAM, and -Xmx200m will specify 200 megs.  
                The max is typically 85% of physical memory.
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

z="-Xmx1g"
z2="-Xms1g"
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
	if [[ $set == 1 ]]; then
		return
	fi
	freeRam 2000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

msa() {
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
	local CMD="java $EA $EOOM $z -cp $CP jgi.FindPrimers $@"
	echo $CMD >&2
	eval $CMD
}

msa "$@"
