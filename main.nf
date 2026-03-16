#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = 'data/*_{1,2}.fq.gz'
params.outDir = 'outputs/'
params.adapters = 'adapters.fa'
log.info """
      LIST OF PARAMETERS
================================
Reads            : ${params.reads}
Output-folder    : ${params.outdir}
Adapters         : ${params.adapters}
"""

// Create read channel
read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true).map { sample, reads -> tuple(sample, reads.collect { it.toAbsolutePath() }) }
adapter_ch = Channel.fromPath(params.adapters)

// Define fastqc process
process fastqc {
    label 'fastqc'
    publishDir "${params.outdir}/quality-control-${sample}/", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(reads)

    output:
    path("*_fastqc.{zip,html}")
    script:
    """
    fastqc ${reads}
    """
}

// Process trimmomatic
process trimmomatic {
    label 'trimmomatic'
    publishDir "${params.outdir}/trimmed-reads-${sample}/", mode: 'copy'

    input:
    tuple val(sample), path(reads)
    path adapters_file

    output:
    file "${sample}_1.paired.fq.gz"
    file "${sample}_2.paired.fq.gz"

    script:
    """
trimmomatic PE -phred33 \
${reads[0]} ${reads[1]} \
${sample}_1.paired.fq.gz ${sample}_1.unpaired.fq.gz \
${sample}_2.paired.fq.gz ${sample}_2.unpaired.fq.gz \
ILLUMINACLIP:/Users/vinithanadar/miniconda/envs/nf_course/share/trimmomatic-0.40-0/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
"""
}

// Run the workflowww
workflow {
    read_pairs_ch.view()
    fastqc(read_pairs_ch)
    trimmomatic(read_pairs_ch, adapter_ch)
}




