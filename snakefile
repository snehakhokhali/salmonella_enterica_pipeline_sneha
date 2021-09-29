configfile: "config.yaml"
ADAPTER = ["untrimmed","trimmed"]

rule all:
    input:
      #expand("/data/short-reads/{barcode}_1.fastq",barcode=config["barcode"]),
      #expand("/data/short-reads/{barcode}_2.fastq",barcode=config["barcode"])
      #expand("/data/trimmed/{barcode}_1.fastq",barcode=config["barcode"]),
      #expand("/data/trimmed/{barcode}_2.fastq",barcode=config["barcode"]),
      #expand("/data/spades_assembled/{barcode}_untrimmed/k_{kvalue}/contigs.fasta",barcode=config["barcode"],kvalue=config["kvalue"]),
      #expand("/data/spades_assembled/{barcode}_trimmed/k_{kvalue}/contigs.fasta",barcode=config["barcode"],kvalue=config["kvalue"]),
      #expand("statistics/statistics_{barcode}.txt",barcode=config["barcode"]),
      #expand("/data/bestAssembly/{barcode}/contigs.fasta",barcode=config["barcode"]),
      expand("plots/histogram/{barcode}_{adapter}_k_{kvalue}.png",barcode=config["barcode"],kvalue=config["kvalue"],adapter=ADAPTER)
rule download_srr:
    output:
      "/data/short-reads/{barcode}_1.fastq",
      "/data/short-reads/{barcode}_2.fastq"
    log: 
      "logs/short-reads/{barcode}.log"
    threads: 4
    shell:
      """(fastq-dump --split-files  {wildcards.barcode} -O /data/short-reads/) > {log}"""

rule Spades:
    input:
      forward_p = "/data/short-reads/{barcode}_1.fastq",
      reverse_p = "/data/short-reads/{barcode}_2.fastq"
    output:
      "/data/spades_assembled/{barcode}_untrimmed/k_{kvalue}/contigs.fasta"
    log: 
      "logs/Spades_assembled/{barcode}_untrimmed/k_{kvalue}/spades.log"
    params:
      kmer = "{kvalue}"
    threads: 4
    shell:
      """(spades.py -k {params.kmer} -t {threads} --pe1-1 {input.forward_p} --pe1-2 {input.reverse_p} -o /data/spades_assembled/{wildcards.barcode}_untrimmed/k_{wildcards.kvalue}) > {log}"""

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
rule statistics:
    input:
      expand("/data/spades_assembled/{{barcode}}_trimmed/k_{kvalue}/contigs.fasta",kvalue=config["kvalue"]),
      expand("/data/spades_assembled/{{barcode}}_untrimmed/k_{kvalue}/contigs.fasta",kvalue=config["kvalue"])
    output:
      "statistics/statistics_{barcode}.txt",
      "/data/bestAssembly/{barcode}/contigs.fasta"
    script: 
      "statistics_assembly.py"

rule histogram_plot:
    input:
      "/data/spades_assembled/{barcode}_{adapter}/k_{kvalue}/contigs.fasta"
    output:
      "plots/histogram/{barcode}_{adapter}_k_{kvalue}.png"
    threads: 4
    script:
      "histogram_plot.py"


