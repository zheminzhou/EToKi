#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 1, 2017

Description:  Renames reads to <prefix>_<number> where you specify the prefix and the numbers are ordered.

Usage:  rename.sh in=<file> in2=<file2> out=<outfile> out2=<outfile2> prefix=<>

in2 and out2 are for paired reads and are optional.
If input is paired and there is only one output file, it will be written interleaved.

Parameters:
prefix=             The string to prepend to existing read names.
ow=f                (overwrite) Overwrites files that already exist.
zl=4                (ziplevel) Set compression level, 1 (low) to 9 (max).
int=f               (interleaved) Determines whether INPUT file is considered interleaved.
fastawrap=70        Length of lines in fasta output.
minscaf=1           Ignore fasta reads shorter than this.
qin=auto            ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto.
qout=auto           ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input).
ignorebadquality=f  (ibq) Fix out-of-range quality values instead of crashing with a warning.

Renaming modes (if not default):
renamebyinsert=f    Rename the read to indicate its correct insert size.
renamebymapping=f   Rename the read to indicate its correct mapping coordinates.
renamebytrim=f      Rename the read to indicate its correct post-trimming length.
addprefix=f         Rename the read by prepending the prefix to the existing name.
prefixonly=f        Only use the prefix; don't add _<number>
addunderscore=t     Add an underscore after the prefix (if there is a prefix).

Sampling parameters:
reads=-1            Set to a positive number to only process this many INPUT reads (or pairs), then quit.

Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
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

z="-Xmx1g"
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

function rename() {
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
	local CMD="java $EA $EOOM $z -cp $CP jgi.RenameReads $@"
	echo $CMD >&2
	eval $CMD
}

rename "$@"
