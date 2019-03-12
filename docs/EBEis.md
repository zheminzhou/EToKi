COMMAND:
    EnteroBase Escherichia in silico serotyping
    
    optional arguments:
        -h, --help            show this help message and exit
        -q QUERY, --query QUERY
            file name for the queried assembly in multi-FASTA format.
        -t TAXON, --taxon TAXON
            Taxon database to compare with. Only support Escherichia
        -p PREFIX, --prefix PREFIX
            prefix for intermediate files.

EXAMPLE:
python EToKi.py EBEis -t Escherichia -q examples/GCF_000010485.1_ASM1048v1_genomic.fna -p SE15
python EToKi.py EBEis -t Escherichia -q examples/GCF_000214765.2_ASM21476v3_genomic.fna -p NA114

OUTPUT:
{"H": "H5", "O": "O16", "prefix": "SE15", "query": "examples/GCF_000010485.1_ASM1048v1_genomic.fna"}
{"H": "H4", "O": "O25", "prefix": "NA114", "query": "examples/GCF_000214765.2_ASM21476v3_genomic.fna"}
