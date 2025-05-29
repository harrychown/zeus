rule quast:
    input:
        fasta="demo_data/EMBL1.spades-run-10.Chr1.fa",
        ref="demo_data/CEA10_Chr1.fasta",
        #gff="annotations.gff",
        #pe1="reads_R1.fastq",
        #pe2="reads_R2.fastq",
        #pe12="reads.fastq",
        #mp1="matereads_R1.fastq",
        #mp2="matereads_R2.fastq",
        #mp12="matereads.fastq",
        #single="single.fastq",
        #pacbio="pacbio.fas",
        #nanopore="nanopore.fastq",
        #ref_bam="ref.bam",
        #ref_sam="ref.sam",
        #bam=["s1.bam","s2.bam"],
        #sam=["s1.sam","s2.sam"],
        #sv_bedpe="sv.bed",
    output:
        multiext("demo_data/report.", "html", "tex", "txt", "pdf", "tsv"),
        multiext("demo_data/transposed_report.", "tex", "txt", "tsv"),
        multiext(
            "demo_data/basic_stats/",
            "cumulative_plot.pdf",
            "GC_content_plot.pdf",
            "gc.icarus.txt",
            "genome_GC_content_plot.pdf",
            "NGx_plot.pdf",
            "Nx_plot.pdf",
        ),
        multiext(
            "demo_data/contigs_reports/",
            "all_alignments_genome.tsv",
            "contigs_report_genome.mis_contigs.info",
            "contigs_report_genome.stderr",
            "contigs_report_genome.stdout",
        ),
        "demo_data/contigs_reports/minimap_output/genome.coords_tmp",
        "demo_data/icarus.html",
        "demo_data/icarus_viewers/contig_size_viewer.html",
        "demo_data/quast.log",
    log:
        "logs/quast.log",
    params:
        extra="--min-contig 5 --min-identity 95.0",
    conda:
        "workflow/envs/quast.yaml",
    wrapper:
        "v3.0.4/bio/quast"

