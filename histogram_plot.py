import matplotlib.pyplot as plt
from Bio import SeqIO
import pandas as pd
import seaborn as sns

lengthcontigs = []
with open(snakemake.input[0]) as handle:
    for archive in SEQIO.parse(handle,"fasta"):
        read = int(archive.id.split("_")[3])
        lengthcontigs.append(read)

plot=pd.DataFrame(lengthcontigs, columns = ["lengthcontigs"])
histog = sns.histplot(data=plot, x="lengthcontigs", bins=50)
histog.set(xlabel= "length of contigs")
histog.set(ylabel = "Number of contigs")
histog.set(title ="histogram of contig length")
plt.yscale("log")
plt.savefig(snakemake.output[0])