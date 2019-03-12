#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 25, 2014

Description:  Generates an identity matrix via all-to-all alignment.

*** WARNING: This program may produce incorrect results in some cirumstances.
*** It is not advisable to use until fixed.

Usage:	idmatrix.sh in=<file> out=<file>

Parameters:
in=<file>       File containing reads. in=stdin.fa will pipe from stdin.
out=<file>      Matrix output. out=stdout will pipe to stdout.
threads=auto    (t) Set number of threads to use; default is number of 
                logical processors.
percent=f       Output identity as percent rather than a fraction.
edits=          Allow at most this much edit distance.  Default is the
                length of the longest input sequence. Lower is faster.
width=          Alignment bandwidth, lower is faster.  Default: 2*edits+1.
usejni=f        (jni) Do alignments faster, in C code.  Requires 
                compiling the C code; details are in /jni/README.txt.

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
NATIVELIBDIR="$DIR""jni/"

z="-Xmx2g"
z2="-Xms2g"
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

idmatrix() {
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
	local CMD="java -Djava.library.path=$NATIVELIBDIR $EA $z -cp $CP jgi.IdentityMatrix $@"
	echo $CMD >&2
	eval $CMD
}

idmatrix "$@"
