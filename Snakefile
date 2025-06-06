import os

GROUPED_DATA = [
    {
        "sample": "demo_data/assembly/EMBL1.spades-run-10.Chr1.fa",
        "reference": "demo_data/reference_genome/CEA10_Chr1.fasta",
        "busco_db": "eurotiales_odb10",
        "protein": "demo_data/protein/A1163.protein.faa",
    }
]



# Use basenames (without extensions) for output filenames
def strip_ext(path):
    return os.path.splitext(os.path.basename(path))[0]



rule all:
    input:
        "results/summary.txt"

rule quast:
    input:
        fasta="demo_data/assembly/EMBL1.spades-run-10.Chr1.fa",
        ref="demo_data/reference_genome/CEA10_Chr1.fasta",
    output:
        directory = directory("results/quast/"),
        report="results/quast/report.tsv",
    log:
        "logs/quast.log",
    params:
        extra="--no-check --no-plots --no-html --no-icarus --no-snps --no-gc --no-sv --no-read-stats --fungus -m 0",
    conda:
        "workflow/envs/quast.yaml",
    shell:
        """
        quast -r {input.ref} -o {output.directory} {input.fasta} {params.extra}
        """

rule busco:
    input:
        fasta=lambda wc: next(p["sample"] for p in GROUPED_DATA if strip_ext(p["sample"]) == wc.sample),
        busco_db_path="demo_data/busco_db/",
    params:
        busco_db_name=lambda wc: next(p["busco_db"] for p in GROUPED_DATA if strip_ext(p["sample"]) == wc.sample)
    output:
        directory = directory("results/samples/{sample}/{busco_db}"),
        report="results/samples/{sample}/{busco_db}/busco/short_summary.specific.{busco_db}.busco.txt",
    log:
        "logs/samples/{sample}/busco.{busco_db}.log",
    conda:
        "workflow/envs/busco.yaml"
    threads: 4
    shell:
        """
        busco -m genome -i {input.fasta} -o busco --out_path {output.directory} --download_path {input.busco_db_path} --miniprot -l {params.busco_db_name} -f -c {threads}
        """

rule miniprot:
    input:
        fasta=lambda wc: next(p["sample"] for p in GROUPED_DATA if strip_ext(p["sample"]) == wc.sample),
        protein=lambda wc: next(p["protein"] for p in GROUPED_DATA if strip_ext(p["sample"]) == wc.sample),
    output:
        "results/samples/{sample}/{sample}.miniprot.gff",
    log:
        "logs/samples/{sample}/miniprot.log",
    conda:
        "workflow/envs/busco.yaml"
    threads: 4
    shell:
        """
        miniprot -t {threads} --gff -I {input.fasta} {input.protein} > {output}
        """

rule miniprot_reference:
    input:
        fasta=lambda wc: next(p["reference"] for p in GROUPED_DATA if strip_ext(p["reference"]) == wc.reference),
        protein=lambda wc: next(p["protein"] for p in GROUPED_DATA if strip_ext(p["reference"]) == wc.reference)
    output:
        "results/references/{reference}/{reference}.miniprot.gff",
    log:
        "logs/references/{reference}/miniprot.log",
    conda:
        "workflow/envs/busco.yaml"
    threads: 4
    shell:
        """
        miniprot -t {threads} --gff -I {input.fasta} {input.protein} > {output}
        """
        
rule earlgrey:
    input:
        fasta=lambda wc: next(p["sample"] for p in GROUPED_DATA if strip_ext(p["sample"]) == wc.sample),
    output:
        report="results/samples/{sample}_EarlGrey/{sample}_summaryFiles/{sample}.highLevelCount.txt"
    log:
        "logs/samples/{sample}/earlgrey.log",
    conda:
        "workflow/envs/earlgrey.yaml"
    threads: 4
    shell:
        """
        earlGrey -g {input.fasta} -t {threads} -s {wildcards.sample} -o results/samples/       
        """


rule summarise:
    input:
        quast="results/quast/report.tsv",
        busco="results/busco/busco/short_summary.specific.eurotiales_odb10.busco.txt",
        miniprot="results/miniprot/miniprot.gff",
        earlgrey="results/earlgrey/test_EarlGrey/test_summaryFiles/test.highLevelCount.txt",
    output:
        "results/summary.txt",
    log:
        "logs/summarise.log",
    shell:
        """
        # QUAST     
        LENGTH=$(grep "Total length (>= 0 bp)" {input.quast} | cut -f2)
        # BUSCO
        BUSCO=$(grep "C:" {input.busco} | cut -d":" -f2 | cut -d"%" -f1)
        # Miniprot
        MINIPROT=$(cat {input.miniprot} | grep -c "mRNA" )
        # EarlGrey
        DNA=$(cat {input.earlgrey} | grep "DNA" |  cut -d$'\t' -f3 || echo "0" )
        LINE=$(cat {input.earlgrey} | grep "LINE" |  cut -d$'\t' -f3 || echo "0" )
        LTR=$(cat {input.earlgrey} | grep "LTR" |  cut -d$'\t' -f3 || echo "0" )
        OTHER=$(cat {input.earlgrey} | grep "Other (Simple Repeat, Microsatellite, RNA)" |  cut -d$'\t' -f3 || echo "0" )
        UNCLASS=$(cat {input.earlgrey} | grep "Unclassified" |  cut -d$'\t' -f3 || echo "0" )
        
        echo "${{LENGTH}}	${{BUSCO}}  ${{MINIPROT}}   ${{DNA}}    ${{LINE}}   ${{LTR}}    ${{OTHER}}  ${{UNCLASS}}" > {output}    
        """



