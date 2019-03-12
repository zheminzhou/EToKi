#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 18, 2016

Description:   Filters NCBI assembly summaries according to their taxonomy.
The specific files are available here:

ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt
or ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt
ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
or ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt

Usage:  filterassemblysummary.sh in=<input file> out=<output file> tree=<tree file> table=<table file> ids=<numbers> level=<name or number>

Standard parameters:
in=<file>       Primary input.
out=<file>      Primary output.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Processing parameters:
level=          Taxonomic level, such as phylum.  Filtering will operate on
                sequences within the same taxonomic level as specified ids.
reqlevel=       Require nodes to have ancestors at these levels.  For example,
                reqlevel=species,genus would ban nodes that are not defined
                at both the species and genus levels.
ids=            Comma-delimited list of NCBI numeric IDs.
names=          Alternately, a list of names (such as 'Homo sapiens').
                Note that spaces need special handling.
include=f       'f' will discard filtered sequences, 't' will keep them.
tree=           A taxonomic tree made by TaxTree, such as tree.taxtree.gz.
table=          A table translating gi numbers to NCBI taxIDs.
                Only needed if gi numbers will be used.
* Note *
Tree and table files are in /global/projectb/sandbox/gaag/bbtools/tax
For non-Genepool users, or to make new ones, use taxtree.sh and gitable.sh

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                specify 200 megs. The max is typically 85% of physical memory.
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
	freeRam 1000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

filterassemblysummary() {
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
	local CMD="java $EA $EOOM $z -cp $CP driver.FilterAssemblySummary $@"
	echo $CMD >&2
	eval $CMD
}

filterassemblysummary "$@"
