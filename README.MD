# Summary
Python program to do 5' and 3' untranslated region (UTR) extraction. 

## UTR_Extractor.py
<ins>Language</ins>: Python3

<ins>Summary</ins>: Extracts all possible UTRs (greater than 1bp by default) from transcripts defined in a gff file. Should work on all species (tested on mouse and human). Sequences are always extracted in the 5' to 3' direction regardless of gene orientation. Requires a gff in standard format, and an indexed fasta file (index using samtools faidx). Outputs to stdout the gene, transcript, location, and sequence of the associated 5' and 3' UTRs. 

#### Usage: UTR_Extractor.py [gff file] [fasta file]


## UTR_Extractor_Bed.py
<ins>Language</ins>: Python3

<ins>Summary</ins>: Derived from UTR_Extractor.py. This one only requires a gff file and outputs the genomic intervals containing 5' and 3' UTRs in a modified bed format. Prints to stdout the Chromosome/Contig; the start and end coords, the strand the gene is found in, the gene name, transcript accession, and UTR type.

#### Usage: UTR_Extractor_Bed.py [gff file] 
