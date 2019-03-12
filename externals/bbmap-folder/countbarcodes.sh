#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 16, 2015

Description: Counts the number of reads with each barcode.

Usage:   countbarcodes.sh in=<file> counts=<file>

Input may be stdin or a fasta or fastq file, raw or gzipped.
If you pipe via stdin/stdout, please include the file type; e.g. for gzipped fasta input, set in=stdin.fa.gz

Input parameters:
in=<file>           Input reads, whose names end in a colon then barcode.
counts=<file>       Output of counts.
interleaved=auto    (int) If true, forces fastq input to be paired and interleaved.
qin=auto            ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto.
unpigz=t            Use pigz to decompress.
expected=           Comma-delimited list of expected bar codes.
valid=              Comma-delimited list of valid bar codes.
countundefined=t    Count barcodes that contain non-ACGT symbols.
printheader=t       Print a header.
maxrows=-1          Optionally limit the number of rows printed.

Output parameters:
out=<file>          Write bar codes and counts here.  'out=stdout' will pipe to standard out.

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

countbarcodes() {
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
	local CMD="java $EA $EOOM $z -cp $CP jgi.CountBarcodes $@"
	echo $CMD >&2
	eval $CMD
}

countbarcodes "$@"
