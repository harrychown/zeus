#!/usr/bin/env nextflow
/*
* Pipeline parameters
*/

/*
=========================================================
nf-core modules
=========================================================
*/
include { QUAST } from "./modules/nf-core/quast/main.nf"
// Primary input
params.input_genome = "${projectDir}/demo_data/EMBL.spades-run-10.Chr1.fa"
// Optional input
params.reference = "${projectDir}/demo_data/CEA10_Chr1.fasta"
params.quast_output = 'quast_output'
params.threads = 4


workflow {
    // Create input genome channel
    genome_ch = Channel.fromPath(params.input_genome)
    // Create input reference channel
    reference_ch = Channel.fromPath(params.reference)
    // Run the QUAST process
   QUAST(
        genome_ch,
        reference_ch,
        [[],[]]
    )   
}