rule all:
    input:
        "results/summary.txt"

rule quast:
    input:
        fasta="demo_data/assembly/EMBL1.spades-run-10.Chr1.fa",
        ref="demo_data/reference_genome/CEA10_Chr1.fasta",
    output:
        directory("results/quast/"),
    log:
        "logs/quast.log",
    params:
        extra="--no-check --no-plots --no-html --no-icarus --no-snps --no-gc --no-sv --no-read-stats --fungus -m 0",
    conda:
        "workflow/envs/quast.yaml",
    shell:
        """
        quast -r {input.ref} -o {output} {input.fasta} {params.extra}
        """

rule busco:
    input:
        fasta="demo_data/assembly/EMBL1.spades-run-10.Chr1.fa",
        busco_db_path="demo_data/busco_db/",
    params:
        busco_db="-l eurotiales_odb10 --miniprot",
    output:
        directory("results/busco/"),
    log:
        "logs/busco.log",
    conda:
        "workflow/envs/busco.yaml"
    threads: 4
    shell:
        """
        busco -m genome -i {input.fasta} -o busco --out_path {output}  --download_path {input.busco_db_path} -f -c {threads} {params.busco_db}
        """

rule summarise:
    input:
        quast="results/quast/report.tsv",
        busco="results/busco/busco/short_summary.specific.eurotiales_odb10.busco.txt",
    output:
        "results/summary.txt",
    log:
        "logs/summarise.log",
    shell:
        """
        LENGTH=$(cat {quast} | grep "Total length (>= 0 bp)" | awk '{print $NF}' )
        BUSCO=$(cat {busco} | grep "C:" | cut -d":" -f2 | cut -d"%" -f1)
        echo "${LENGTH}	${BUSCO}" > {output}    
        """



