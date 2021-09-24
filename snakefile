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
    threads: 6
    shell:
      """fastq-dump --split-files  {wildcards.barcode} -O /data/short-reads/ -e {threads}"""

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
    threads: 5
    shell:
      """spades.py -k {params.kmer} -t {threads} --pe1-1 {input.forward_p} --pe1-2 {input.reverse_p} -o /data/spades_assembled/{wildcards.barcode}/k_{wildcards.kvalue}"""
