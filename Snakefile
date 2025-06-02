rule quast:
    input:
        fasta="demo_data/EMBL1.spades-run-10.Chr1.fa",
        ref="demo_data/CEA10_Chr1.fasta",
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




