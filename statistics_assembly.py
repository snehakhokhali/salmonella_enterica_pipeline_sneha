from Bio import SeqIO
import sys

#first input has to be config.fasta file
file1 = sys.argv[1]

archives = list(SeqIO.parse(file1, "fasta"))

Averagecontigs= 0

for number in archives:
    Numbercontigs= len(number)
    Totallengthcontigs=0
    contigslength = []

    for record in number:
        contigslength.append(len(record))
        Totallengthcontigs += len(record)
    Averagecontigs = Totallengthcontigs/Numbercontigs
print("Total no. of contigs",Numcontigs)
print("Average contigs:",Averagecontigs)