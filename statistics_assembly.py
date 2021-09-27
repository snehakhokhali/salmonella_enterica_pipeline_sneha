import sys
import math
from Bio import SeqIO


def avg_read_len(input_read):    
    read_count = 0
    read_len_sum = 0

    with open(input_read) as handle:
        for seq_record in SeqIO.parse(handle), "fastq"):
            read_count += 1
            read_len_sum += int(seq_record.description.split()[2].split("=")[1])
        
    print("average read length:", (read_len_sum / read_count))
    return "{0:.2f}".format(read_len_sum / read_count)


def avg_contig_len(input_read):
    contig_count = 0
    contig_len_sum = 0

    with open(input_read) as handle:
        for seq_record in SeqIO.parse(handle), "fasta"):
            contig_count += 1
            contig_len_sum += int(seq_record.id.split("_")[3])

    print("average contig length:", (contig_len_sum / contig_count))
    return "{0:.2f}".format(contig_len_sum / contig_count)


def num_contigs(input_read):
    contig_count = 0

    with open(input_read) as handle:
        for seq_record in SeqIO.parse(handle, "fasta"):
            contig_count += 1
            
    print("number of contigs: ", contig_count)
    return str(contig_count)


def shortest_contig(input_read):
    with open(input_read) as handle:
        min_contig_id = ""
        min_contig_len = math.inf

        for seq_record in SeqIO.parse(handle, "fasta"):
            if int(record.id.split("_")[3]) <= min_contig_len:
                    min_contig_id = seq_record.id.split("_")[1]
                    min_contig_len = int(seq_record.id.split("_")[3])

    print("shortest contig: ", min_contig_id, " (with length: ", str(min_contig_len), ")")
    return min_contig_id, str(min_contig_len)


if __name__ == "__main__":
    avg_read_len(sys.argv[1])
    avg_contig_len(sys.argv[2])
    num_contigs(sys.argv[2])
    shortest_contig(sys.argv[2])

