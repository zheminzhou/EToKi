#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified August 16, 2017

Description:  Demultiplexes sequences into multiple files based on their names,
substrings of their names, or prefixes or suffixes of their names.
Opposite of muxbyname.

Usage:
demuxbyname.sh in=<file> in2=<file2> out=<outfile> out2=<outfile2> names=<string,string,string...>

Alternately:
demuxbyname.sh in=<file> out=<outfile> delimiter=whitespace prefixmode=f
This will demultiplex by the substring after the last whitespace.

demuxbyname.sh in=<file> out=<outfile> length=8 prefixmode=t
This will demultiplex by the first 8 characters of read names.

in2 and out2 are for paired reads and are optional.
If input is paired and there is only one output file, it will be written interleaved.

Parameters:
in=<file>       Input file.
out=<file>      Output files for reads with matched headers (must contain % symbol).
outu=<file>     Output file for reads with unmatched headers.
prefixmode=t    (pm) Match prefix of read header.  If false, match suffix of read header.
column=-1       If positive, split the header on a delimiter and match that column (1-based).
                For example, using this header:
                NB501886:61:HL3GMAFXX:1:11101:10717:1140 1:N:0:ACTGAGC+ATTAGAC
                You could split by tile (11101) using 'delimiter=: column=5'
substring=f     Names can be substrings of read headers.
names=          List of strings (or files containing strings) to parse from read names.
                This is optional.
length=0        If positive, use a suffix or prefix of this length from read name instead of or in addition to the list of names.
                For example, you could create files based on the first 8 characters of read names.
delimiter=      For prefix or suffix mode, specifying a delimiter will allow exact matches even if the length is variable.
                This allows demultiplexing based on names that are found without specifying a list of names.
                In suffix mode, for example, everything after the last delimiter will be used.
                Normally the delimiter will be used as a literal string (a Java regular expression); for example, ':' or 'HISEQ'.
                But there are some special delimiters which will be replaced by the symbol they name, 
                because they are reserved in some operating systems or cause other problems.
                These are provided for convenience due to possible OS conflicts:
                   space, tab, whitespace, pound, greaterthan, lessthan, equals,
                   colon, semicolon, bang, and, quote, singlequote
                These are provided because they interfere with Java regular expression syntax:
                   backslash, hat, dollar, dot, pipe, questionmark, star,
                   plus, openparen, closeparen, opensquare, opencurly
                In other words, to match '.', you should set 'delimiter=dot'.

Common parameters:
ow=t            (overwrite) Overwrites files that already exist.
app=f           (append) Append to files that already exist.
zl=1            (ziplevel) Set compression level, 1 (low) to 9 (max).
int=auto        (interleaved) Determines whether INPUT file is considered interleaved.
qin=auto        ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto.
qout=auto       ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input).
                    

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
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

z="-Xmx1200m"
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

function demuxbyname() {
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
	local CMD="java $EA $EOOM $z -cp $CP jgi.DemuxByName $@"
	echo $CMD >&2
	eval $CMD
}

demuxbyname "$@"
