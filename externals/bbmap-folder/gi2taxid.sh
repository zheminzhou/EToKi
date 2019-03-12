#!/bin/bash

usage(){
echo "
Written by Brian Bushnell.
Last modified April 19, 2018

Description:  Renames fasta sequences to indicate their NCBI taxIDs.
The headers must be in NCBI or Silva format with gi numbers,
accessions, or organism names.  Only supports fasta sequence files.

Usage:  gi2taxid.sh in=<file> out=<file> gi=<file>

Parameters:
in=<file>       Input sequences; required parameter.  Must be fasta.
                This can alternatively be a comma-delimited list,
                or just a bunch of space-delimited filenames, e.g.:
                gi2taxid.sh x.fa y.fa z.fa out=tid.fa tree=auto table=auto
out=<file>      Destination for renamed sequences.
invalid=<file>  Destination for headers with no taxid.
keepall=t       Keep sequences with no taxid in normal output.
prefix=t        Append the taxid as a prefix to the old header, but keep
                the old header.
title=tid       Set the title of the new number (e.g. ncbi, taxid, tid).
ziplevel=2      (zl) Compression level for gzip output.
pigz=f          Spawn a pigz (parallel gzip) process for faster 
                compression than Java.  Requires pigz to be installed.
silva=f         Parse headers in Silva format.
shrinknames=f   Replace multiple concatenated headers with the first.

Taxonomy file flags:
server=f        Use the taxonomy server instead of local files.
                Server mode only works for accessions (like RefSeq).
tree=           Specify a taxtree file.  On Genepool, use 'auto'.
gi=             Specify a gitable file.  On Genepool, use 'auto'.
accession=      Specify one or more comma-delimited NCBI accession to
                taxid files.  On Genepool, use 'auto'.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx800m will specify 800 megs.
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
NATIVELIBDIR="$DIR""jni/"

z="-Xmx7g"
z2="-Xms7g"
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
	freeRam 7000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"


gi2taxid() {
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
	local CMD="java $EA $EOOM $z $z2 -cp $CP tax.RenameGiToNcbi $@"
	echo $CMD >&2
	eval $CMD
}

gi2taxid "$@"
