#!/usr/bin/env python3
#Written by: Emil Karpinski (2024-10-07)

#Curr version: 2B

#Extracts UTRs from a genome given a gff file. 
#Will also extract the sequences if the genome is available and indexed. 
#Note will only extract UTRs from protein coding genes as per definition (5' UTR = region from TSS to start codon; 3' URT = region from stop codon to end of transcript)
#V1 - This works and extracts the "UTR" from the gene start to gene end. However, it only extracts the minimal UTR (i.e. the shortest span between the TSS and a coding exon). More importantly it also extracts raw DNA sequence not spliced DNA sequence. So if the mRNA is spliced outside of the coding region it doesn't capture that. 
#V2 - This one returns 5' and 3' UTRs on a transcript level and take's mRNA splicing outside of the coding regions into account. 
#V2B - this is a modified version of UTR_Extractor to output a bed file of the 5' and 3' UTR regions for later variant filtering.


#Importing one package to do the sequence extraction and one to catch command line args
#To install: pip install pyfaidx  
#from pyfaidx import Fasta  # type: ignore
import sys

#Function to do printing
def Printing(Name, TransSc, Chr, GS, GE, Dir, Flag, Orientation):
    if Flag == 0:
        print("Chromosome", "Start", "Stop", "Strand",  "Gene", "Transcript", "UTR_Type")
        Flag = 1
    #Adding one here to report accurate values for printing.
    if Dir == "+":
        if Orientation == "Up":
            print(Chr, GS+1, GE+1, Dir, Name, TransSc, "5p-UTR", sep='\t')
        elif Orientation == "Down":
            print(Chr, GS+1, GE, Dir, Name, TransSc, "3p-UTR", sep='\t')
    elif Dir == "-":
        if Orientation == "Up":
            print(Chr, GS+1, GE+1, Dir, Name, TransSc, "3p-UTR", sep='\t')
        if Orientation == "Down":     
            print(Chr, GS+1, GE, Dir, Name, TransSc, "5p-UTR", sep='\t')
    #Returning the value of Flag/Header here otherwise it won't update outside of here and it'll print the header everytime it prints a new gene. 
    return(Flag)

##Variables
#Stores information from the GFF File.
#Notably, we only store information from the gene record since that's the most expansive one. 
#For the CDS records we only record each start and end. 
#Setting these all to 0/empty right now since we'll be updating them as we go.
GeneName = ""
TranscriptID = ""
GeneStart = 0
GeneEnd = 0
Contig = ""
CDS_Start = []
CDS_End = []
Strand = ""
Exons = []

#A flag only used once in the Printing function to print the header row.
Header = 0

#Getting the name of the GFF file to parse through and the name of the genome fasta file.
#The genome fasta file has to be indexed beforehand.
#Both of these are passed via the command line.
GFF_File = sys.argv[1]

#Opening the GFF File for reading. 
with open(GFF_File, "r+") as GFF:
    #Loops through and reads one line at a time. Using a seperate readline() here causes it to skip lines. 
    for line in GFF:

        #Splits the current line of input on tabs since GFFs are tab-delimited
        CurrLine=line.split('\t')

        #Passing over the headers which all start with "#"
        if CurrLine[0][0] == "#":
            pass

        #Checking to see if the line contains a new mRNA or gene as these are the two lines that denote we're going to something else where we need to store new info and maybe do some printing.
        elif CurrLine[2] == "mRNA" or CurrLine[2] == "gene":
            #Checks if there's something to print.
            if TranscriptID != "" and len(CDS_Start) != 0 and (min(CDS_Start)-GeneStart != 0 or GeneEnd - max(CDS_End) != 0):
                Point1 = min(CDS_Start)
                Point2 = max(CDS_End)
                for i in range(0,len(Exons)):
                    #The pluses and minuses here work with the ones in the printing function for specific cases. Don't touch these.
                    if Exons[i][0] < Point1 and Exons[i][1] < Point1:
                        Header = Printing(GeneName, TranscriptID, Contig, Exons[i][0], Exons[i][1], Strand, Header, "Up")
                    elif Exons[i][0] < Point1 and Exons[i][1] > Point1:
                        Header = Printing(GeneName, TranscriptID, Contig, Exons[i][0], Point1-1, Strand, Header, "Up")
                    #Had to update this to a seperate if statment (not an elif) to get around the edge case where mRNA encodes a single exon that contains a CDS sequence. Otherwise it would only print a 5' UTR. 
                    if Exons[i][0] > Point2 and Exons[i][1] > Point2:
                        Header = Printing(GeneName, TranscriptID, Contig, Exons[i][0], Exons[i][1]+1, Strand, Header, "Down")
                    elif Exons[i][0] < Point2 and Exons[i][1] > Point2:
                        Header = Printing(GeneName, TranscriptID, Contig, Point2+1, Exons[i][1]+1, Strand, Header, "Down")

            #Modified checks to see if the line is a gene or transcript as we need to collect different information at different points.
            if CurrLine[2] == "gene":
                GeneName = CurrLine[8][8:CurrLine[8].index(";")]
            elif CurrLine[2] == "mRNA":
                #Storing relevant variables from the current transcript
                #Note: GeneStart and GeneEnd are the start and end of the current transcript (not gene) as of V2
                GeneStart = int(CurrLine[3])-1
                GeneEnd = int(CurrLine[4])-1
                Contig = CurrLine[0]
                Strand = CurrLine[6]
                TranscriptID = CurrLine[8][7:CurrLine[8].index(";")]
                
            #Clearing these so they don't inherit the values for the previous genes.
            CDS_Start.clear()
            CDS_End.clear()
            Exons.clear()

        #Need to record the exon junctions here so we can splice together the mRNA sequence correctly.
        elif CurrLine[2] == "exon" and TranscriptID in CurrLine[8]:
            Exons.append([(int(CurrLine[3])-1), (int(CurrLine[4])-1)])
        
        #Else checking to see if the info in column 3 indicates that the target is a CDS and as a safety check contains the info for GeneName we stored
        #If so appenending the values to a list so we can select the appropriate ones later. 
        elif CurrLine[2] == "CDS" and TranscriptID in CurrLine[8]:
            CDS_Start.append(int(CurrLine[3])-1)
            CDS_End.append(int(CurrLine[4])-1)


#One final print statment to print the last thing stuck in memory in case it ends with a valid entry. 
if TranscriptID != "" and len(CDS_Start) != 0 and (min(CDS_Start)-GeneStart != 0 or GeneEnd - max(CDS_End) != 0):
    #Copy and pasted the code above for the sequence acquisition. This should probably be in a function for cleanliness though.
    Point1 = min(CDS_Start)
    Point2 = max(CDS_End)
    for i in range(0,len(Exons)):
        #The pluses and minuses here work with the ones in the printing function for specific cases. Don't touch these.
        if Exons[i][0] < Point1 and Exons[i][1] < Point1:
            Header = Printing(GeneName, TranscriptID, Contig, Exons[i][0], Exons[i][1], Strand, Header, "Up")
        elif Exons[i][0] < Point1 and Exons[i][1] > Point1:
            Header = Printing(GeneName, TranscriptID, Contig, Exons[i][0], Point1-1, Strand, Header, "Up")
        #Had to update this to a seperate if statment (not an elif) to get around the edge case where mRNA encodes a single exon that contains a CDS sequence. Otherwise it would only print a 5' UTR. 
        if Exons[i][0] > Point2 and Exons[i][1] > Point2:
            Header = Printing(GeneName, TranscriptID, Contig, Exons[i][0], Exons[i][1]+1, Strand, Header, "Down")
        elif Exons[i][0] < Point2 and Exons[i][1] > Point2:
            Header = Printing(GeneName, TranscriptID, Contig, Point2+1, Exons[i][1]+1, Strand, Header, "Down")












