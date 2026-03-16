#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// -------------------------
// PARAMETERS
// -------------------------
params.reads = 'data/*_{1,2}.fq.gz'       // paired-end reads
params.outdir = 'outputs'                // output folder for all results
params.adapters = 'adapters.fa'          // adapters file in repo

log.info """
      LIST OF PARAMETERS
================================
Reads            : ${params.reads}
Output-folder    : ${params.outdir}
Adapters         : ${params.adapters}
"""

// -------------------------
// CHANNELS
// -------------------------

// Channel for paired-end reads
read_pairs_ch = Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .map { sample, reads -> tuple(sample, reads) }

// Channel for adapter file
adapter_ch = Channel.fromPath(params.adapters)

// -------------------------
// PROCESS: FASTQC
// -------------------------
process fastqc {
    label 'fastqc'
    publishDir "${params.outdir}/quality-control", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(reads)

    output:
    path("*_fastqc.{zip,html}")

    script:
    """
    fastqc ${reads[0]} ${reads[1]}
    """
}

// -------------------------
// PROCESS: TRIMMOMATIC
// -------------------------
process trimmomatic {
    label 'trimmomatic'
    publishDir "${params.outdir}/trimmed-reads", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(reads)
    path adapters_file

    output:
    file "${sample}_1.paired.fq.gz"
    file "${sample}_1.unpaired.fq.gz"
    file "${sample}_2.paired.fq.gz"
    file "${sample}_2.unpaired.fq.gz"

    script:
    """
    trimmomatic PE -phred33 \
    ${reads[0]} ${reads[1]} \
    ${sample}_1.paired.fq.gz ${sample}_1.unpaired.fq.gz \
    ${sample}_2.paired.fq.gz ${sample}_2.unpaired.fq.gz \
    ILLUMINACLIP:${adapters_file}:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// -------------------------
// WORKFLOW
// -------------------------
workflow {

    // Show input samples
    read_pairs_ch.view()

    // Run FastQC
    fastqc(read_pairs_ch)

    // Run Trimmomatic
    trimmomatic(read_pairs_ch, adapter_ch)
}

