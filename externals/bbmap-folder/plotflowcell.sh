#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified August 9, 2018

Description:  Generates statistics about flowcell positions.

Usage:	plotflowcell.sh in=<input> out=<output>

Input parameters:
in=<file>           Primary input file.
in2=<file>          Second input file for paired reads in two files.
indump=<file>       Specify an already-made dump file to use instead of
                    analyzing the input reads.
reads=-1            Process this number of reads, then quit (-1 means all).
interleaved=auto    Set true/false to override autodetection of the
                    input file as paired interleaved.

Output parameters:
out=<file>          Output file for filtered reads.
dump=<file>         Write a summary of quality information by coordinates.

Tile parameters:
xsize=500           Initial width of micro-tiles.
ysize=500           Initial height of micro-tiles.
size=               Allows setting xsize and ysize tot he same value.
target=800          Iteratively increase the size of micro-tiles until they
                    contain an average of at least this number of reads.

Other parameters:
trimq=-1            If set to a positive number, trim reads to that quality
                    level instead of filtering them.
qtrim=r             If trimq is positive, to quality trimming on this end
                    of the reads.  Values are r, l, and rl for right,
                    left, and both ends.

Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will specify 20 GB of RAM; -Xmx200m will specify 
                    200 MB.  The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an
                    out-of-memory exception occurs.  Requires Java 8u92+.
-da                 Disable assertions.

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

z="-Xmx8g"
z2="-Xms8g"
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

plotflowcell() {
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
	local CMD="java $EA $EOOM $z $z2 -cp $CP hiseq.PlotFlowCell $@"
	echo $CMD >&2
	eval $CMD
}

plotflowcell "$@"
