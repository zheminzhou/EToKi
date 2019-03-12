#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified August 25, 2016

Description:  Estimates cardinality of unique kmers in sequence data.
See also kmercountmulti.sh.

Usage:  loglog.sh in=<file> k=<31>

Parameters:
in=<file>       (in1) Input file, or comma-delimited list of files.
in2=<file>      ptional second file for paired reads.
k=31            Use this kmer length for counting.
buckets=1999    Use this many buckets for counting; higher decreases variance, for large datasets.
bits=8          Hash this many bits per cycle.  More bits use more memory.
seed=-1         Use this seed for hash functions.  A negative number forces a random seed.
minprob=0       Set to a value between 0 and 1 to exclude kmers with a lower probability of being correct.

Shortcuts:
The # symbol will be substituted for 1 and 2.
For example:
loglog.sh in=read#.fq
...is equivalent to:
loglog.sh in1=read1.fq in2=read2.fq

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
-da             Disable assertions.

Supported input formats are fastq, fasta, scarf, sam, and bam.
Supported compression formats are gzip and bz2.
To read from stdin, set 'in=stdin'.  The format should be specified with an extension, like 'in=stdin.fq.gz'

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

z="-Xmx200m"
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

function loglog() {
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
	local CMD="java $EA $EOOM $z -cp $CP jgi.LogLog $@"
	echo $CMD >&2
	eval $CMD
}

loglog "$@"
