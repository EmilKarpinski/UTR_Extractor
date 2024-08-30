# Summary
Python program to do 5' and 3' untranslated region (UTR) extraction. 

## UTR_Extractor.py
<ins>Language</ins>: Python3

<ins>Summary</ins>: Extracts all possible UTRS (greater than 1bp by default) from transcripts defined in a gff file. Should work on all species (tested on mouse and human). Sequences are always extracted in the 5' to 3' direction regardless of gene orientation. Requires a gff in standard format, and an indexed fasta file (index using samtools faidx). Outputs a file containing the gene, transcript, location, and sequence of the associated 5' and 3' UTRs. 

###Usage: UTR_Extractor.py [gff file] [fasta file]
