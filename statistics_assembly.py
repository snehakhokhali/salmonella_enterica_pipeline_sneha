from Bio import SeqIO
from shutil import copyfile
import os

LISTAveragecontigs=[]
LISTtotalcontigs=[]
LISTshortestcontigs=[]
LISTlongestcontigs=[]
LISTN50 =[]
LISTBP300N50=[]
pathList= []
for inputpath in snakemake.input:
    pathList.append(inputpath)
    shortestcontigs = 100000000
    longestcontigs = 0
    Numbercontigs = 0
    Totallengthcontigs = 0
    Totallength300contigs = 0

    for record in SeqIO.parse(inputpath, "fasta"):
        Numbercontigs += 1
        num = len(record)
        Totallengthcontigs += num
        if (num > 300):
            Totallength300contigs += num
        if (num < shortestcontigs):
            shortestcontigs = num
        if (num > longestcontigs):
            longestcontigs = num

    Averagecontigs = Totallengthcontigs/Numbercontigs

    N50 = 0
    N50_sum = 0

    N50_BP300 = 0
    N50_sum_BP300 = 0

    for record in SeqIO.parse(inputpath, "fasta"):
        num = len(record)
        if (N50_sum < Totallengthcontigs/2):
            N50_sum += num
            if (N50_sum >= Totallengthcontigs/2):
                N50 = num
        if (N50_sum_BP300 < Totallength300contigs/2):
            N50_sum_BP300 += num
            if (N50_sum_BP300 >= Totallength300contigs/2):
                N50_BP300 = num




    LISTAveragecontigs.append(Averagecontigs)
    LISTtotalcontigs.append(Numbercontigs)
    LISTshortestcontigs.append(shortestcontigs)
    LISTlongestcontigs.append(longestcontigs)
    LISTN50.append(N50)
    LISTBP300N50.append(N50_BP300)


print("Total no. of contigs",LISTtotalcontigs)
print("Average contigs:",LISTAveragecontigs)
print("shortest contigs:", LISTshortestcontigs)
print("longest contigs:",LISTlongestcontigs)
print("N50 of all contigs:", LISTN50)
print("N50 of all contigs longer than 300 BP:", LISTBP300N50)

BESTN50 = lISTN50.index(max(LISTN50))
bestpath = pathList[BESTN50]

copypath = "/data/bestAssembly/" + str(snakemake.wildcards)

os.makedirs(copypath, exist_ok = True)
copypath = copypath + "/contigs.fasta"
copyfile(bestpath, copypath)

with open(snakemake.output[0], "w") as result:
    result.write("best assembly is in :" + bestpath)
    result.write("\nshortest contig:")
    for out in LISTshortestcontigs:
        result.write(str(out)+"\t")
    result.write("\nlongest contig:")
    for out in LISTlongestcontigs:
        result.write(str(out)+"\t")
    result.write("\nnumber of contigs: ")
    for out in LISTtotalcontigs:
        result.write(str(out) + "\t")
    result.write ("\navg contig length; ")
    for out in LISTAveragecontigs:
        result.write(str(out)+"\t")
    result.write("\nN50 of all contigs: ")
    for out in LISTN50:
        result.write(str(out)+"\t")
    result.write("\nN50 of contigs longer than 300 bp:")
    for out in LISTBP300N50:
        result.write(str(out)+"\t")