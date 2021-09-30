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
      #expand("plots/histogram/{barcode}_{adapter}_k_{kvalue}.png",barcode=config["barcode"],kvalue=config["kvalue"],adapter=ADAPTER)
      #expand("/data/long-reads/{long_barcode}.fastq",long_barcode=config["long_barcode"])
      expand("/data/long_read_assembly/{long_barcode}.paf.gz",long_barcode=config["long_barcode"]),
      expand("/data/long_read_assembly/{long_barcode}.gfa",long_barcode=config["long_barcode"]),
      expand("/data/long_read_assembly/{long_barcode}/contigs.fasta",long_barcode=config["long_barcode"]),
      expand("assembly-statistik/short-reads/{barcode}.txt",barcode=config["barcode"]),
      expand("assembly-statistik/long-reads/{long_barcode}.txt",long_barcode=config["long_barcode"])
rule download_srr:
    output:
      "/data/short-reads/{barcode}_1.fastq",
      "/data/short-reads/{barcode}_2.fastq"
    log: 
      "logs/short-reads/{barcode}.log"
    threads: 4
    shell:
      """(fastq-dump --split-files  {wildcards.barcode} -O /data/short-reads/) > {log}"""
rule download_srr_long:
    output:
      "/data/long-reads/{long_barcode}.fastq"
    threads: 4
    log:
      "logs/long-reads/{long_barcode}.log"
    shell:
      """(fastq-dump {wildcards.long_barcode} -O /data/long-reads/) > {log}"""
rule minimap2:
    input:
      "/data/long-reads/{long_barcode}.fastq"
    output:
      "/data/long_read_assembly/{long_barcode}.paf.gz"
    log:
      "logs/minimap2/{long_barcode}.log"
    shell:
      """minimap2 -x ava-ont -t8 {input} {input} | gzip -1 > {output}"""
rule miniasm:
    input:
      paf="/data/long_read_assembly/{long_barcode}.paf.gz",
      fastq="/data/long-reads/{long_barcode}.fastq"
    output:
      "/data/long_read_assembly/{long_barcode}.gfa"
    log:
      "logs/miniasm/{long_barcode}.log"
    shell:
      """miniasm -f {input.fastq} {input.paf} > {output}"""
rule gfa_to_fasta:
    input:
      "/data/long_read_assembly/{long_barcode}.gfa"
    output:
      "/data/long_read_assembly/{long_barcode}/contigs.fasta"
    log:
      "logs/fasta/{long_barcode}.log"
    shell:
      """awk '/^S/{{print \">\"$2"\\n\"$3}}' {input} | fold > {output}"""
rule assembly_stats_short:
    input:
      short_in="/data/bestAssembly/{barcode}/contigs.fasta"
    output:
      short_out="assembly-statistik/short-reads/{barcode}.txt"
    shell:
      """assembly-stats {input.short_in} > {output.short_out}"""
rule assembly_stats_long:
    input: 
      long_in="/data/long_read_assembly/{long_barcode}/contigs.fasta"
    output:
      long_out="assembly-statistik/long-reads/{long_barcode}.txt"
    shell:
      """assembly-stats {input.long_in} > {output.long_out}"""

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


