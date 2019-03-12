#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 17, 2015

Description:  Randomly adds adapters to a file, or grades a trimmed file.
The input is a set of reads, paired or unpaired.
The output is those same reads with adapter sequence replacing some of the bases in some reads.
For paired reads, adapters are located in the same position in read1 and read2.
This is designed for benchmarking adapter-trimming software, and evaluating methodology.
randomreads.sh is better for paired reads, though, as it actually adds adapters at the correct location,
so that overlap may be used for adapter detection.

Usage:  addadapters.sh in=<file> in2=<file2> out=<outfile> out2=<outfile2> adapters=<file>

in2 and out2 are for paired reads and are optional.
If input is paired and there is only one output file, it will be written interleaved.

Parameters:
ow=f                (overwrite) Overwrites files that already exist.
int=f               (interleaved) Determines whether INPUT file is considered interleaved.
qin=auto            ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto.
qout=auto           ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input).
add                 Add adapters to input files.  Default mode.
grade               Evaluate trimmed input files.
adapters=<file>     Fasta file of adapter sequences.
literal=<sequence>  Comma-delimited list of adapter sequences.
left                Adapters are on the left (3') end of the read.
right               Adapters are on the right (5') end of the read.  Default mode.
adderrors=t         Add errors to adapters based on the quality scores.
addpaired=t         Add adapters to the same location for read 1 and read 2.
arc=f               Add reverse-complemented adapters as well as forward.
rate=0.5            Add adapters to this fraction of reads.

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

function addadapters() {
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
	local CMD="java $EA $EOOM $z -cp $CP jgi.AddAdapters $@"
	echo $CMD >&2
	eval $CMD
}

addadapters "$@"
