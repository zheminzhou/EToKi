#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified January 28, 2019

Description:  Creates one or more sketches from a fasta file,
optionally annotated with taxonomic information.

Please read bbmap/docs/guides/BBSketchGuide.txt for more information.

Usage:  sketch.sh in=<fasta file> out=<sketch file>

Standard parameters:
in=<file>           A fasta file containing one or more sequences.
out=<file>          Output filename.  If multiple files are desired it must
                    contain the # symbol.
blacklist=<file>    Ignore keys in this sketch file.  Additionaly, there are
                    built-in blacklists that can be specified:
                       nt:      Blacklist for nt
                       refseq:  Blacklist for Refseq
                       silva:   Blacklist for Silva
                       img:     Blacklist for IMG
files=1             Number of parallel files, for faster loading.  This is
                    independent of the number of sketches produced; sketches
                    will be randomly distributed between files.
size=10000          Desired size of sketches (if not using autosize).
maxfraction=0.01    (mgf) Max fraction of genomic kmers to use.
minsize=100         Do not generate sketches for genomes smaller than this.
autosize=t          Use flexible sizing instead of fixed-length.
sizemult=1          Multiply the default size of sketches by this factor.
                    Normally a bacterial-size genome will get a sketch size
                    of around 10000; if autosizefactor=2, it would be ~20000.
k=31                Kmer length, 1-32.  To maximize sensitivity and 
                    specificity, dual kmer lengths may be used:  k=31,24
                    Dual kmers are fastest if the shorter is a multiple 
                    of 4.  Query and reference k must match.
rcomp=t             Look at reverse-complement kmers also.
amino=f             Use amino acid mode.  Input should be amino acids.
translate=f         Call genes and translate to proteins.  Input should be
                    nucleotides.  Designed for prokaryotes.
mode=single         Possible modes:
                       single: Write one sketch.
                       sequence: Write one sketch per sequence.
                       taxa: Write one sketch per taxonomic unit.
                          Requires more memory, and taxonomic annotation.
                       img: Write one sketch per IMG id.
delta=t             Delta-compress sketches.
a48=t               Encode sketches as ASCII-48 rather than hex.
depth=f             Track the number of times kmers appear.
                    Required for the depth2 field in comparisons.
entropy=0.66        Ignore sequence with entropy below this value.

Metadata flags (optional; intended for single-sketch mode):
taxid=-1            Set the NCBI taxid.
imgid=-1            Set the IMG id.
spid=-1             Set the JGI sequencing project id.
name=               Set the name (taxname).
name0=              Set name0 (normally the first sequence header).
fname=              Set fname (normally the file name).
meta_=              Set an arbitrary metadata field.
                    For example, meta_Month=March.

Taxonomy-specific flags:
tree=               Specify a taxtree file.  On Genepool, use 'auto'.
gi=                 Specify a gitable file.  On Genepool, use 'auto'.
accession=          Specify one or more comma-delimited NCBI accession to
                    taxid files.  On Genepool, use 'auto'.
imgdump=            Specify an IMG dump file containing NCBI taxIDs,
                    for IMG mode.
taxlevel=subspecies Taxa hits below this rank will be promoted and merged
                    with others.
prefilter=f         For huge datasets full of junk like nt, this flag
                    will save memory by ignoring taxa smaller than minsize.
                    Requires taxonomic information (tree and gi).
tossjunk=t          For taxa mode, discard taxonomically uninformative
                    sequences.  This includes sequences with no taxid,
                    with a tax level NO_RANK, of parent taxid of LIFE.
silva=f             Parse headers using Silva or semicolon-delimited syntax.


Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an
                    out-of-memory exception occurs.  Requires Java 8u92+.
-da                 Disable assertions.

For more detailed information, please read /bbmap/docs/guides/BBSketchGuide.txt.
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
	freeRam 4000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

sketch() {
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
	local CMD="java $EA $EOOM $z $z2 -cp $CP sketch.SketchMaker $@"
	echo $CMD >&2
	eval $CMD
}

sketch "$@"
