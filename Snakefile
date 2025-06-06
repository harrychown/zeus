# Declare samples



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
        fasta="demo_data/assembly/EMBL1.spades-run-10.Chr1.fa",
        busco_db_path="demo_data/busco_db/",
    params:
        busco_db="-l eurotiales_odb10 --miniprot",
    output:
        directory = directory("results/busco/"),
        report="results/busco/busco/short_summary.specific.eurotiales_odb10.busco.txt",
    log:
        "logs/busco.log",
    conda:
        "workflow/envs/busco.yaml"
    threads: 4
    shell:
        """
        busco -m genome -i {input.fasta} -o busco --out_path {output.directory}  --download_path {input.busco_db_path} -f -c {threads} {params.busco_db}
        """

rule miniprot:
    input:
        fasta="demo_data/assembly/EMBL1.spades-run-10.Chr1.fa",
        protein="demo_data/protein/A1163.protein.faa",
    output:
        "results/miniprot/miniprot.gff",
    log:
        "logs/miniprot.log",
    conda:
        "workflow/envs/busco.yaml"
    threads: 4
    shell:
        """
        miniprot -t {threads} --gff -I {input.fasta} {input.protein} > {output}
        """
        
rule earlgrey:
    input:
        fasta="demo_data/assembly/EMBL1.spades-run-10.Chr1.fa",
    output:
        directory = directory("results/earlgrey/"),
        report="results/earlgrey/test_EarlGrey/test_summaryFiles/test.highLevelCount.txt"
    log:
        "logs/earlgrey.log",
    conda:
        "workflow/envs/earlgrey.yaml"
    threads: 4
    shell:
        """
        earlGrey -g {input.fasta} -t {threads} -s test -o {output.directory}       
        """
rule dagchainer:
    input:
        ref_miniprot="demo_data/miniprot/CEA10_Chr1.gff",
        miniprot="results/miniprot/miniprot.gff",
    output:
        matchlist="results/dagchainer/test.matchList",
        filtered="results/dagchainer/test.matchList.filtered",
        coords="results/dagchainer/test.matchList.filtered.aligncoords",
    log:
        "logs/dagchainer.log",
    conda:
        "workflow/envs/dagchainer.yaml"
    shell:
        """
        python workflow/scripts/minicombiner.py -r={input.ref_miniprot} -d={miniprot} -o={output.matchlist}
        workflow/scripts/dagchainer/bin/accessory_scripts/filter_repetitive_matches.pl 5 < {output.matchlist} > {output.filtered}
        run_DAG_chainer.pl -i  {output.filtered}
        """
rule summarise:
    input:
        quast="results/quast/report.tsv",
        busco="results/busco/busco/short_summary.specific.eurotiales_odb10.busco.txt",
        miniprot="results/miniprot/miniprot.gff",
        earlgrey="results/earlgrey/test_EarlGrey/test_summaryFiles/test.highLevelCount.txt",
        dagchainer="results/dagchainer/test.matchList.filtered.aligncoords",
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
        # Dagchainer
        DAG_LENGTH=$(grep -v "#" {input.dagchainer} | awk '{{sum += ($7 > $8 ? $7 - $8 : $8 - $7)}} END {print sum}')        
        echo "${{LENGTH}}	${{BUSCO}}  ${{MINIPROT}}   ${{DNA}}    ${{LINE}}   ${{LTR}}    ${{OTHER}}  ${{UNCLASS}}    ${{DAG_LENGTH}}" > {output}    
        """



