#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified January 28, 2019

Description:  Merges multiple sketches into a single sketch.

Please read bbmap/docs/guides/BBSketchGuide.txt for more information.

Usage:           mergesketch.sh in=a.sketch,b.sketch out=c.sketch
With wildcards:  mergesketch.sh *.sketch out=c.sketch

Standard parameters:
in=<file>       Input sketches or fasta files; may be a comma-delimited
                list.  in= is optional so wildcards may be used.
out=<file>      Output sketch.
amino=f         Use amino acid mode.

Sketch-making parameters:
mode=single     Possible modes, for fasta input:
                   single: Generate one sketch per file.
                   sequence: Generate one sketch per sequence.
autosize=t      Produce an output sketch of whatever size the union 
                happens to be.
size=           Restrict output sketch to this upper bound of size.
k=31            Kmer length, 1-32.  To maximize sensitivity and 
                specificity, dual kmer lengths may be used:  k=31,24
                Dual kmers are fastest if the shorter is a multiple 
                of 4.  Query and reference k must match.
                Currently, JGI-hosted nt and RefSeq use k=31,24
                while JGI-hosted Silva uses k=31.
keyfraction=0.2 Only consider this upper fraction of keyspace.
minkeycount=1   Ignore kmers that occur fewer times than this.  Values
                over 1 can be used with raw reads to avoid error kmers.
depth=f         Retain kmer counts if available.

Metadata parameters: (if blank the values of the first sketch will be used)
taxid=-1        Set the NCBI taxid.
imgid=-1        Set the IMG id.
spid=-1         Set the JGI sequencing project id.
name=           Set the name (taxname).
name0=          Set name0 (normally the first sequence header).
fname=          Set fname (normally the file name).
meta_=          Set an arbitrary metadata field.
                For example, meta_Month=March.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
-da             Disable assertions.

For more detailed information, please read /bbmap/docs/guides/BBSketchGuide.txt.
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

z="-Xmx4g"
z2="-Xms4g"
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
	freeRam 3200m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

sendsketch() {
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
	local CMD="java $EA $EOOM $z -cp $CP sketch.MergeSketch $@"
#	echo $CMD >&2
	eval $CMD
}

sendsketch "$@"
