#!/usr/bin/env python3
#Written by: Emil Karpinski (2024-08-21)

#V2

#Extracts UTRs from a genome given a gff file. 
#Will also extract the sequences if the genome is available and indexed. 
#Note will only extract UTRs from protein coding genes as per definition (5' UTR = region from TSS to start codon; 3' URT = region from stop codon to end of transcript)
#V1 - This works and extracts the "UTR" from the gene start to gene end. However, it only extracts the minimal UTR (i.e. the shortest span between the TSS and a coding exon). More importantly it also extracts raw DNA sequence not spliced DNA sequence. So if the mRNA is spliced outside of the coding region it doesn't capture that. 
#V2 - This one returns 5' and 3' UTRs on a transcript level and take's mRNA splicing outside of the coding regions into account. 

#Importing one package to do the sequence extraction and one to catch command line args
#To install: pip install pyfaidx  
from pyfaidx import Fasta  # type: ignore
import sys

#Gets the requested sequence in either forward or reverse complement direction (returns always 5' to 3')
def GetSeq(Chr, Start, End, Dir):
    Seq=""
    if Dir == "+":
            #Don't need to add the .seq here to the end of these since apparently assigning the output to a variable automatically returns only the sequence.
            Seq = str(genes[Chr][Start:End])
    elif Dir == "-":
            Seq = str(-genes[Chr][Start:End])
    return (Seq.upper())

#The below function (GetUTRs) is decrepit as of this version. Left here only for historical sake.
'''
#Function to get the UTR sequnces from the file. 
def GetUTRs(GS, GE, CS, CE, Chr, Dir):
    Seq1 = ""
    Seq2 = ""
    if Seq_File != "":
        if Dir == "+":
            #Don't need to add the .seq here to the end of these since apparently assigning the output to a variable automatically returns only the sequence.
            Seq1 = str(genes[Chr][GS:CS])
            Seq2 = str(genes[Chr][CE+1:GE+1])
        elif Dir == "-":
            Seq1 = str(-genes[Chr][CE+1:GE+1])
            Seq2 = str(-genes[Chr][GS:CS])
    return(Seq1.upper(), Seq2.upper())
'''

#Header = Printing(GeneName, TranscriptID, Contig, GeneStart, GeneEnd, min(CDS_Start), max(CDS_End), UTR_5_Seq, UTR_3_Seq, Strand, Header)
#Only true in the forward direction otherwise reversed. 

#Function to do printing
def Printing(Name, TransSc, Chr, GS, GE, CS, CE, U5S, U3S, Dir, Flag):
    if Flag == 0:
        print("Gene_Name", "Transcript", "UTR_Type", "Locus", "Start_Position", "End_Position", "Sequence")
        Flag = 1
    #Adding one here to report accurate values for printing. Otherwise they're all 
    #Note that in this version the 5' and 3' UTR are relative to the DNA and strand agnostic. So they need to be fliped if the gene is on the reverse strand.
    if Dir == "+":
        if len(U5S) >= 1:
            print(Name, TransSc, "5p-UTR", Chr, GS+1, (CS-1)+1, U5S, sep='\t')
        if len(U3S) >= 1:
            print(Name, TransSc, "3p-UTR", Chr, (CE+1)+1, GE+1, U3S, sep='\t')
    elif Dir == "-":
        if len(U3S) >= 1:
            print(Name, TransSc, "5p-UTR", Chr, GE+1, (CE+1)+1, U3S, sep='\t')
        if len(U5S) >= 1:     
            print(Name, TransSc, "3p-UTR", Chr, (CS-1)+1, GS+1, U5S, sep='\t')
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
UTR_5_Seq = ""
UTR_3_Seq = ""

#A flag only used once in the Printing function to print the header row.
Header = 0

#Getting the name of the GFF file to parse through and the name of the genome fasta file.
#The genome fasta file has to be indexed beforehand.
#Both of these are passed via the command line.
GFF_File = sys.argv[1]
Seq_File = sys.argv[2]

#This doesn't work. If a Seq_File is not given then you get an error at line 68 (above - during variable assignment.)
#If a seq file is given, extract the names. 
#Could probably also check that this is indexed, but that should really be on the end user.
if Seq_File != "":
    genes = Fasta(Seq_File)

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
            #print(CurrLine[2])
            if TranscriptID != "" and len(CDS_Start) != 0 and (min(CDS_Start)-GeneStart != 0 or GeneEnd - max(CDS_End) != 0):
                Point1 = min(CDS_Start)
                Point2 = max(CDS_End)
                for i in range(0,len(Exons)):
                    #Going to need to move these pluses and minuses into the function up top or they won't be correct for the different strands.
                    #Not actually sure this is true. I think it should work, but need to test with reverse stuff.
                    if Exons[i][0] < Point1 and Exons[i][1] < Point1:
                        UTR_5_Seq += GetSeq(Contig, Exons[i][0], Exons[i][1]+1, Strand)
                    elif Exons[i][0] < Point1 and Exons[i][1] > Point1:
                        UTR_5_Seq += GetSeq(Contig, Exons[i][0], Point1, Strand)
                    #Had to update this to a seperate if statment (not an elif) to get around the edge case where mRNA encodes a single exon that contains a CDS sequence. Otherwise it would only print a 5' UTR. 
                    if Exons[i][0] > Point2 and Exons[i][1] > Point2:
                        UTR_3_Seq += GetSeq(Contig, Exons[i][0], Exons[i][1]+1, Strand)
                    elif Exons[i][0] < Point2 and Exons[i][1] > Point2:
                        UTR_3_Seq += GetSeq(Contig, Point2+1, Exons[i][1]+1, Strand)

                #Sending everything off for printing and catching if the header has already been printed or not.
                Header = Printing(GeneName, TranscriptID, Contig, GeneStart, GeneEnd, Point1, Point2, UTR_5_Seq, UTR_3_Seq, Strand, Header)
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
            UTR_5_Seq = ""
            UTR_3_Seq = ""

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
    #Doing seq retrieval and printing one final time to get the last stored entry (else the last gene in the document would be cut off. )
    #Copy and pasted the code above for the sequence acquisition. This should probably be in a function for cleanliness though.
    Point1 = min(CDS_Start)
    Point2 = max(CDS_End)
    for i in range(0,len(Exons)):
        #Going to need to move these pluses and minuses into the function up top or they won't be correct for the different strands.
        #Not actually sure this is true. I think it should work, but need to test with reverse stuff.
        if Exons[i][0] < Point1 and Exons[i][1] < Point1:
            UTR_5_Seq += GetSeq(Contig, Exons[i][0], Exons[i][1]+1, Strand)
        elif Exons[i][0] < Point1 and Exons[i][1] > Point1:
            UTR_5_Seq += GetSeq(Contig, Exons[i][0], Point1, Strand)
        elif Exons[i][0] > Point2 and Exons[i][1] > Point2:
            UTR_3_Seq += GetSeq(Contig, Exons[i][0], Exons[i][1]+1, Strand)
        elif Exons[i][0] < Point2 and Exons[i][1] > Point2:
            UTR_3_Seq += GetSeq(Contig, Point2+1, Exons[i][1]+1, Strand)

    #Sending everything off for printing and catching if the header has already been printed or not.
    Header = Printing(GeneName, TranscriptID, Contig, GeneStart, GeneEnd, Point1, Point2, UTR_5_Seq, UTR_3_Seq, Strand, Header)












