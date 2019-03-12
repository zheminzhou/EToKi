#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified September 17, 2018
This script requires at least 10GB RAM.
It is designed for NERSC and uses hard-coded paths.

Description:  Removes all reads that map to selected common microbial contaminant genomes.
Removes approximately 98.5% of common contaminant reads, with zero false-positives to non-bacteria.
NOTE!  This program uses hard-coded paths and will only run on Nersc systems.

Usage:  removemicrobes.sh in=<input file> outu=<clean output file>

Input may be fasta or fastq, compressed or uncompressed.

Parameters:
in=<file>           Input reads.  Should already be adapter-trimmed.
outu=<file>         Destination for clean reads.
outm=<file>         Optional destination for contaminant reads.
threads=auto        (t) Set number of threads to use; default is number of logical processors.
overwrite=t         (ow) Set to false to force the program to abort rather than overwrite an existing file.
interleaved=auto    (int) If true, forces fastq input to be paired and interleaved.
trim=t              Trim read ends to remove bases with quality below minq.
                    Values: t (trim both ends), f (neither end), r (right end only), l (left end only).
untrim=t            Undo the trimming after mapping.
minq=4              Trim quality threshold.
ziplevel=6          (zl) Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster.

build=1             Choses which masking mode was used:  
                    1 is most stringent and should be used for bacteria.
                    2 uses fewer bacteria for masking (only RefSeq references).
                    3 is only masked for plastids and entropy, for use on anything except bacteria.
                    4 is unmasked.

***** All BBMap parameters can be used; run bbmap.sh for more details. *****

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

z="-Xmx6000m"
z2="-Xms6000m"
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

function removemicrobes() {
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
	local CMD="java $EA $z $z2 -cp $CP align2.BBMap strictmaxindel=4 bwr=0.16 bw=12 ef=0.001 minhits=2 path=/global/projectb/sandbox/gaag/bbtools/commonMicrobes pigz unpigz zl=6 qtrim=r trimq=10 untrim idtag printunmappedcount ztd=2 kfilter=25 maxsites=1 k=13 minid=0.95 idfilter=0.95 minhits=2 build=1 bloomfilter $@"
	echo $CMD >&2
	eval $CMD
}

removemicrobes "$@"
