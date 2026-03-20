#!/usr/bin/env nextflow

nextflow.enable.dsl=2


// PARAMETERS

params.reads = 'data/*_{1,2}.fq.gz'
params.outdir = 'outputs'
params.adapters = 'adapters.fa'
params.genome = 'data/LG12.fasta'

log.info """
      LIST OF PARAMETERS
================================
Reads            : ${params.reads}
Output-folder    : ${params.outdir}
Adapters         : ${params.adapters}
"""


// CHANNELS

// Channel for paired-end reads

read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    .map { sample, reads -> tuple(sample, reads.collect{it.toAbsolutePath()}) }


// Channel for adapter file

adapter_ch = Channel.fromPath(params.adapters)

// channel for the genome file

genome_ch = Channel.fromPath("data/LG12.fasta*")

//Including the Modules

include { fastqc } from './modules/fastqc.nf'

include { trimmomatic } from './modules/trimmomatic.nf'

include { bwa_mem2 } from './modules/bwa_mem2.nf'


// WORKFLOW
          
workflow {
    
    read_pairs_ch.view()
    
    fastqc(read_pairs_ch)
       
    trimmomatic(read_pairs_ch, adapter_ch)
 
    bwa_mem2(trimmomatic.out.paired_reads, genome_ch)
}



