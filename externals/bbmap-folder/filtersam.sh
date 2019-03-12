#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified March 21, 2018

Description:  Filters a sam file to remove reads with substitutions
unsupported by other reads (bad subs).  For particularly bad data,
it may be advisable to iteratively re-call variants and re-run FilterSam.
Calling variants may be performed like this:

callvariants.sh in=mapped.sam out=vars.vcf clearfilters minreads=2

Usage:  filtersam.sh in=<file> out=<file> vcf=<file>

Parameters:
in=<file>       Input sam or bam file.
ref=<file>      Optional fasta reference file.
out=<file>      Output file for good reads.
outb=<file>     Output file for bad reads.
vcf=<file>      VCF file of variants called from these reads.
vars=<file>     Alternatively, variants can be provided in CallVariants'
                native output format.
maxbadsubs=2    A read will be discarded if it has more than maxbadsubs.
mbsad=2         (maxbadsuballeledepth) A substitution is bad if the 
                allele depth is at most this much. 
mbsrd=2         (minbadsubreaddepth) Substitutions may only be considered
                bad if the total read depth spanning the variant is
                at least this much.

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

filtersam() {
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
	local CMD="java $EA $EOOM $z $z2 -cp $CP var2.FilterSam $@"
#	echo $CMD >&2
	eval $CMD
}

filtersam "$@"
