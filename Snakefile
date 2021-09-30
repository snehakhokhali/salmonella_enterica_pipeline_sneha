configfile: "config.yaml"
rule all:
    input:
      expand("/data/prokka/short-reads/{barcode}.txt",barcode=config["barcode"]),
      expand("/data/prokka/short-reads/{barcode}.gff",barcode=config["barcode"])
rule prokka_short:
    input:
      "/data/bestAssembly/{barcode}/contigs.fasta"
    output:
      "/data/prokka/short-reads/{barcode}.txt",
      "/data/prokka/short-reads/{barcode}.gff"
    threads: 4
    shell:
      """prokka -O /data/prokka/short-reads/ --force --prefix {wildcards.barcode} --cpus {threads}  {input}"""
