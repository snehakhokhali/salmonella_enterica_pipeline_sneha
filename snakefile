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
      """fasterq-dump {wildcards.barcode} -O /data/short-reads/ -e {threads}"""

    
