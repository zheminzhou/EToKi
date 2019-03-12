COMMAND:
    usage: EToKi.py configure [-h] [--usearch USEARCH]
        [--kraken_database KRAKEN_DATABASE]
    
    Specify links to kraken database and usearch program.
    
    optional arguments:
        -h, --help            show this help message and exit
        --usearch USEARCH     usearch is required for ortho and MLSType. Download
            the 32-bit version from
            https://www.drive5.com/usearch/.
        --kraken_database KRAKEN_DATABASE
            Kraken is optional in the assemble module. You can
            specify your own database or use MiniKraken2: https://
            ccb.jhu.edu/software/kraken2/dl/minikraken2_v2_8GB.tgz
    
EXAMPLE:
    python EToKi.py configure --usearch /path/to/usearch --kraken_database /path/to/minikraken2_v2

OUTPUT:
    2019-03-09 16:38:25.333438      adapters ("/path/to/EToKi/externals/adapters.fa") is present.
    2019-03-09 16:38:25.348602      blastn ("/path/to/EToKi/externals/blastn") is present.
    2019-03-09 16:38:25.409206      bowtie2 ("/path/to/EToKi/externals/bowtie2") is present.
    2019-03-09 16:38:25.513261      bowtie2build ("/path/to/EToKi/externals/bowtie2-build") is present.
    2019-03-09 16:38:25.518220      bwa ("/path/to/EToKi/externals/bwa") is present.
    2019-03-09 16:38:26.414283      enbler_filter ("/path/to/EToKi/modules/_EnFlt.py") is present.
    2019-03-09 16:38:30.665654      gatk ("/path/to/EToKi/externals/gatk-package-4.1.0.0-local.jar") is present.
    2019-03-09 16:38:30.719354      kraken2 ("/path/to/EToKi/externals/kraken2") is present.
    2019-03-09 16:38:30.719536      kraken_database ("/path/to/minikraken2_v2") is present.
    2019-03-09 16:38:30.732788      makeblastdb ("/path/to/EToKi/externals/makeblastdb") is present.
    2019-03-09 16:38:30.892571      megahit ("/path/to/EToKi/externals/megahit") is present.
    2019-03-09 16:38:30.897455      minimap2 ("/path/to/EToKi/externals/minimap2") is present.
    2019-03-09 16:38:30.902475      mmseqs ("/path/to/EToKi/externals/mmseqs") is present.
    2019-03-09 16:38:31.686454      pilon ("/path/to/EToKi/externals/pilon-1.23.jar") is present.
    2019-03-09 16:38:31.703595      samtools ("/path/to/EToKi/externals/samtools") is present.
    2019-03-09 16:38:31.900430      spades ("/path/to/EToKi/externals/spades.py") is present.
    2019-03-09 16:38:31.906154      usearch ("/path/to/usearch") is present.
    2019-03-09 16:38:31.911018      Configuration complete.
