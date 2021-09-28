from Bio import SeqIO
import sys

LISTAveragecontigs=[]
LISTtotalcontigs=[]
LISTshortestcontigs=[]
LISTlongestcontigs=[]
LISTN50 =[]
LISTBP300N50=[]
LISTAverageforwardreads=[]
LISTAveragereversereads=[]

for inputpath in snakemake.input:

    #first input has to be config.fasta file
    file1 = snakemake.input[0]
    untrimmed=list(SeqIO.parse(file1, "fasta"))
    file2= snakemake.input[1]
    trimmed=list(SeqIO.parse(file2, "fasta"))

    archives = [untrimmed, trimmed]

    Averagecontigs= 0
# total length of contigs, average contig length, shortest contigs and longest contigs
    for number in archives:
        Numbercontigs= len(number)
        Totallengthcontigs=0
        Totallength300contigs = 0
        contigslength = []
        shortestcontigs= len(number[0])
        longestcontigs = len(number[0])
        contigs_above300 = []
        for record in number:
            contigslength.append(len(record))
            Totallengthcontigs += len(record)
            if len(record) < shortestcontigs:
                shortestcontigs = len(record)
            if len(record) > longestcontigs:
                longestcontigs = len(record)
            if len(record) > 300:
                contigs_above300.append(len(record))
                Totallength300contigs += len(record)

        Averagecontigs = Totallengthcontigs/Numbercontigs


        # N50 of all contigs
        contigslength.sort()
        numb = contigslength[0]
        num = 0
        while numb < (Totallengthcontigs/2):
            num += 1
            numb += contigslength[num]
        N50contigs = contigslength[num]

        #N50 of all the contigs longer than 300 bp
        contigs_above300.sort()
        numb300 = contigs_above300[0]
        num300 = 0
        while numb300 < (Totallength300contigs/2):
            num300 += 1
            numb300 += contigs_above300[num300]
        contigs_above300N50 = contigs_above300[num300]

        LISTAveragecontigs.append(Averagecontigs)
        LISTtotalcontigs.append(Numbercontigs)
        LISTshortestcontigs.append(shortestcontigs)
        LISTlongestcontigs.append(longestcontigs)
        LISTN50.append(N50contigs)
        LISTBP300N50.append(contigs_above300N50)


    print("Total no. of contigs",LISTtotalcontigs)
    print("Average contigs:",LISTAveragecontigs)
    print("shortest contigs:", LISTshortestcontigs)
    print("longest contigs:",LISTlongestcontigs)
    print("N50 of all contigs:", LISTN50)
    print("N50 of all contigs longer than 300 BP:", LISTBP300N50)

