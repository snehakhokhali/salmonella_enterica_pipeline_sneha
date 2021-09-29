# Genomics pipeline for reproducible data analysis of Salmonella enterica genomes - Bioinformatics Practical Course
## Day 4: 
[TOC]
## Table of Contents
### Task 1: About Salmonella

Salmonella enterica is a highly diverse Gram negative bacterial species. A few S.enterica serovars includongy Typhi, Paratyphi A, B or C are highly adapted to the human host as their reservoir. They are the casuative agents of enteric fever(also known as typhoid fever or paratyphoid fever). Enteric fever is an invasive, life-threatening, systemic disease resulting more deaths. It is endemic in developing world in regions that lack water nad adequate sanitation, faciliating the spread of these pathogens via the faecal-oral route. 

Source: https://www.frontiersin.org/articles/10.3389/fmicb.2014.00391/full
### Task 2 and 3: Pipelining and Improving the assembly

#### creating the conda environment
 * environment.yml witth content:
``` 
channels:
        - bioconda
        - conda-forge
dependencies:
        - snakemake-minimal
        - sra-tools
        - SPAdes
```
 * conda environment 
```
conda env create -n salmonella-pipeline -f environment.yaml
```
 * activate the environment 
```
conda activate salmonella-pipeline
```
 * check installation:
```
snakemake --help
```

### Using sra-tools to download the sequencing data with the barcode (SRR1965341) from the SRA 
* first rule of sra-tools download
```
configfile: "config.yaml"
rule all:
    input:
      expand("/data/short-reads/{barcode}_1.fastq",barcode=config["barcode"]),
      expand("/data/short-reads/{barcode}_2.fastq",barcode=config["barcode"])
rule download_srr:
    output:
      "/data/short-reads/{barcode}_1.fastq",
      "/data/short-reads/{barcode}_2.fastq"
    threads: 6
    shell:
      """fastq-dump --split-files {wildcards.barcode} -O /data/short-reads/"""

```
* fasterq-dump helps to download the data from SRA and split the paired files automatically into forward and reverse reads. 
* configfile is "config.yaml" which contains the barcode (SRR1965341). And the config file is used in input for both {barcode}_1 and {barcode}_2 of rule all.
* to split and save the files into two separate files, --split-files parameter is used.
* Output files are: {barcode}_1.fastq, {barcode}_2.fastq 
* to specify the output directory, -O parameter is used.
* Threads: 6
 
## Day 5:

### **SPAdes assembly**
* config.yaml
```
barcode:["SRR1965341"]
kvalue:["31","55"]
```
* snakefile

```
configfile: "config.yaml"
rule all:
    input:
      #expand("/data/short-reads/{barcode}_1.fastq",barcode=config["barcode"]),
      #expand("/data/short-reads/{barcode}_2.fastq",barcode=config["barcode"])
      expand("/data/spades_assembled/{barcode}/k_{kvalue}/contigs.fasta",barcode=config["barcode"],kvalue=config["kvalue"])
rule download_srr:
    output:
      "/data/short-reads/{barcode}_1.fastq",
      "/data/short-reads/{barcode}_2.fastq"
    log: 
      "logs/short-reads/{barcode}.log"
    threads: 4
    shell:
      """fastq-dump --split-files  {wildcards.barcode} -O /data/short-reads/"""

rule Spades:
    input:
      forward_p = "/data/short-reads/{barcode}_1.fastq",
      reverse_p = "/data/short-reads/{barcode}_2.fastq"
    output:
      "/data/spades_assembled/{barcode}/k_{kvalue}/contigs.fasta"
    log: 
      "logs/Spades_assembled/{barcode}/k_{kvalue}/spades.log"
    params:
      kmer = "{kvalue}"
    threads: 4
    shell:
      """spades.py -k {params.kmer} -t {threads} --pe1-1 {input.forward_p} --pe1-2 {input.reverse_p} -o /data/spades_assembled/{wildcards.barcode}/k_{wildcards.kvalue}"""
```
* input files are from the output of rule download_srr, i.e SRR1965341_1. fastq and SRR1965341_2.fastq
* output files are /data/spades_assembled/{barcode}/k_31/contigs.fasta and /data/spades_assembled/{barcode}/k_55/contigs.fasta
* parameter are k-value ["31","55"]
* threads: 4
* parameter --pe1-1 PATH is for path to forward reads
* parameter --pe1-2 PATH is  for path to reverse reads
* -k INT is for k-values 
* -t is for number of threads
* -o is for output directory. 
* the assembly should be run with two different values for k (the k-mer size), therefore producing one contigs file per k value
* Note: -e parameter is not used in fastq-dump.(failure execution because it doesnot contain -e parameter for threads)
* ` Barcode: ["SRR1965341","SRR1968189","SRR7828287","SRR2075991","SRR5584993"]` : 4 new Barcodes were added
* And ` kvalue: ["auto","31","55"]` were added.
* And start the screen session:` screen -S salmonella-pipeline`
* activate conda environment:` conda activate salmonella-pipeline`
* run the pipeline:` snakemake --cores 4`
* interrupt the execution of pipeline: Ctrl+C
* unlock the screen: snakemake --unlock
* start with only one barcode 
* and k-value ["auto","31","55"] was used and only auto of k_value of 1 barcode was pushed to git. I delete the spades_assembled and short-reads of logs. Now only k_auto of 1 barcode was created??????
* attach (re-enter) session: screen -x or screen -r salmonella-pipeline
* run pipeline: snakemake --cores 4
* detach session when the analysis is finish: press Ctrl+a, d

## Day 6:Monday(27.09.2021)

### **Task 6: Adapter trimming**
* adapter. fasta was available on Teams
* adapter_1.fasta and adapter_2.fasta were created.
* add cutadapt in environment.yml for trimming:
```
channels:
        - bioconda
        - conda-forge
dependencies:
        - snakemake-minimal=6.8.0
        - sra-tools=2.11.0
        - SPAdes=3.15.3
        - cutadapt=3.4
```
* install cutadapt automatically: conda env update --file environment.yml --prune
* activate conda environment: conda activate salmonella-pipeline
* snakefile
```
expand("/data/trimmed/{barcode}_1.fastq",barcode=config["barcode"]),
expand("/data/trimmed/{barcode}_2.fastq",barcode=config["barcode"]),
expand("/data/spades_assembled/{barcode}_trimmed/k_{kvalue}/contigs.fasta",barcode=config["barcode"],kvalue=config["kvalue"])
```
* These are added in the rule of all to trim the barcode with the help pf adapter. 
 
```
rule adapter_trimming:
    input:
      "/data/short-reads/{barcode}_1.fastq",
      "/data/short-reads/{barcode}_2.fastq"
    output:
      "/data/trimmed/{barcode}_1.fastq", 
      "/data/trimmed/{barcode}_2.fastq" 
    threads: 4
    shell:
      "cutadapt -a file:adapter_1.fasta -A file:adapter_2.fasta -o /data/trimmed/{wildcards.barcode}_1.fastq -p /data/trimmed/{wildcards.barcode}_2.fastq "
      "/data/short-reads/{wildcards.barcode}_1.fastq /data/short-reads/{wildcards.barcode}_2.fastq -j 0"

rule spades_assembler:
    input:
      left="/data/trimmed/{barcode}_1.fastq",
      right="/data/trimmed/{barcode}_2.fastq"
    output:
      "/data/spades_assembled/{barcode}_trimmed/k_{kvalue}/contigs.fasta"
    threads: 4
    shell:
      "spades.py --pe1-1 {input.left} --pe1-2 {input.right} "
      "-o /data/spades_assembled/{wildcards.barcode}_trimmed/k_{wildcards.kvalue} -k {wildcards.kvalue} -t 4"

```
* rule adapter_trimming and spades_assembler were added beneath the rules of spades. Changed the barcode of rule of spades into {barcode}_untrimmed to specify the wildcards.
* `expand("/data/spades_assembled/{barcode}_trimmed/k_{kvalue}/contigs.fasta",barcode=config["barcode"],kvalue=config["kvalue"])` was added in rule of all. 
* input files from trimming: from /data/short-reads/{barcode}_1.fastq and from /data/short-reads/{barcode}_2.fastq
* output files: "/data/trimmed/{barcode}_1.fastq" and "/data/trimmed/{barcode}_2.fastq"
* parameter:` -a PATH`: adapter.fasta forward file and ` -A PATH`: adapter.fasta reverse file
* `-o`: output of forward fasta file and `-p`: output of reverse fasta file
* `-j 0` : automatically selects the available threads 

### **Task 7**
* ` Barcode: ["SRR1965341","SRR1968189","SRR7828287","SRR2075991","SRR5584993"]` : 4 new Barcodes were added in config.yaml
* And start the screen session:` screen -S salmonella-pipeline `
* run the pipeline:` snakemake --cores 4 `

## Day 7: Tuesday
* To check if all the assemblies of trimmed and untrimmed with different parameters were completed or not, this command has been used.
```
ls -la /data/spades_assembled/*/*/contigs.fasta
```
- if the contigs files exists, the assembly was completed. It was shown like below as a result. And every assembly of trimmed and untrimmed of every barcode with every parameters were existed.
```
-rw-r--r-- 1 root root 5086604 Sep 27 15:38 /data/spades_assembled/SRR1965341_trimmed/k_31/contigs.fasta 
-rw-r--r-- 1 root root 5067565 Sep 27 15:54 /data/spades_assembled/SRR1965341_trimmed/k_55/contigs.fasta 
-rw-r--r-- 1 root root 5067136 Sep 27 15:22 /data/spades_assembled/SRR1965341_trimmed/k_auto/contigs.fasta
-rw-r--r-- 1 root root 5153904 Sep 24 13:31 /data/spades_assembled/SRR1965341_untrimmed/k_31/contigs.fasta 
-rw-r--r-- 1 root root 5085537 Sep 24 14:16 /data/spades_assembled/SRR1965341_untrimmed/k_55/contigs.fasta 
-rw-r--r-- 1 root root 5085101 Sep 24 18:21 /data/spades_assembled/SRR1965341_untrimmed/k_auto/contigs.fasta 
```

###  **Task 4: Assembly statistics**
*  creating biopython script to collect the following statistics about the read data assemblies:
   * average read length
   * average contig length
   * total number of contigs
   * shortest contigs
   * longest contigs
   * N50 of all contigs
   * N50 of all contigs longer than 300 bp
* **installation of biopython in conda environment**
  * include `biopython=1.79` to ` environment.yml` file under dependencies.
  * update conda environment: `conda env update --file environment.yaml --prune`
  * check installation: `conda list`

**Calculation of assembly statistics**
* bioscript: `statistics_assembly.py`
```

m Bio import SeqIO
from shutil import copyfile
import os
import math

LISTAveragecontigs=[]  #average contig length
LISTtotalcontigs=[]    #total number of contigs
LISTshortestcontigs=[] #shortest contigs
LISTlongestcontigs=[]  #longest contigs
LISTN50 =[]            #N50 contigs
LISTBP300N50=[]        #N50 of all contigs longer than 300 bp
pathList= []           # list with input paths
for inputpath in snakemake.input:
    pathList.append(inputpath)
    shortestcontigs = 100000000
    longestcontigs = 0
    Numbercontigs = 0
    Totallengthcontigs = 0
    Totallength300contigs = 0
    #calculate average contig length, total number of contigs, shortest contigs and longest contigs
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
    #calculate N50
    for record in SeqIO.parse(inputpath, "fasta"):
        num = len(record)
        if (N50_sum < Totallengthcontigs/2):
            N50_sum += num
            if (N50_sum >= Totallengthcontigs/2):
                N50 = num
        if (N50_sum_BP300 < Totallength300contigs/2): #calculate N50 for contigs longer than 300 bp
            N50_sum_BP300 += num
            if (N50_sum_BP300 >= Totallength300contigs/2):
                N50_BP300 = num



    #append lists with calculated values
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

#calculate best assembly
BESTN50 = LISTN50.index(max(LISTN50))
bestpath = pathList[BESTN50]

#copy contig.fasta of best result into separate folder
copypath = "/data/bestAssembly/" + str(snakemake.wildcards)

os.makedirs(copypath, exist_ok = True)
copypath = copypath + "/contigs.fasta"
copyfile(bestpath, copypath)

#create new file with overview of results
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
```
* The biopython script was integrated by adding a new rule named rule statistics into the snakefile.
```
rule statistics:
    input:
      expand("/data/spades_assembled/{{barcode}}_trimmed/k_{kvalue}/contigs.fasta",kvalue=config["kvalue"]),
      expand("/data/spades_assembled/{{barcode}}_untrimmed/k_{kvalue}/contigs.fasta",kvalue=config["kvalue"])
    output:
      "statistics/statistics_{barcode}.txt",
      "/data/bestAssembly/{barcode}/contigs.fasta"
    script: 
      "statistics_assembly.py"
``` 
```
expand("statistics/statistics_{barcode}.txt",barcode=config["barcode"]),
expand("/data/bestAssembly/{barcode}/contigs.fasta",barcode=config["barcode"])
```
* These are added in the rule of all.
* 
![dag](dag.svg)

 
