from Bio import SeqIO
import sys

#first input has to be config.fasta file
file1 = sys.argv[1]

archives = list(SeqIO.parse(file1, "fasta"))

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
    
print("Total no. of contigs",Numcontigs)
print("Average contigs:",Averagecontigs)