#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 28, 2018

Description:  Cuts out features defined by a gff file, and writes them
to a new fasta.  Features are output in their sense strand.

Usage:  cutgff.sh in=<fna file> gff=<gff file> out=<fna file>

Peak-calling parameters:
in=<file>           Input FNA (fasta) file.
gff=<file>          Input GFF file (optional).
out=<file>          Output FNA file.
types=CDS           Types of features to cut.
invert=false        Invert selection.
attributes=         A comma-delimited list of strings.  If present, one of
                    these strings must be in the gff line attributes.
bannedattributes=   A comma-delimited list of banned strings.
banpartial=t        Ignore lines with 'partial=true' in attributes.
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

gff() {
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
	local CMD="java $EA $EOOM $z -cp $CP prok.CutGff $@"
#	echo $CMD >&2
	eval $CMD
}

gff "$@"
