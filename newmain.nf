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

genome_ch = Channel.fromPath("data/LG12.fasta*").collect()
ref_prefix = file(params.genome).name


// PROCESS: FASTQC

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

// PROCESS: TRIMMOMATIC

process trimmomatic {
    label 'trimmomatic'
    publishDir "${params.outdir}/trimmed-reads-${sample}/", mode: 'copy', overwrite: true

    input: 
    tuple val(sample), path(reads)
    path adapters_file
    
    output:
    tuple val("${sample}"),path("${sample}*.paired.fq.gz"),emit:
    paired_reads
    tuple val("${sample}"),path("${sample}*.unpaired.fq.gz"),emit:
    unpaired_reads
    
    script:
    """   
    trimmomatic PE -phred33 \
    ${reads[0]} ${reads[1]} \
    ${sample}_1.paired.fq.gz ${sample}_1.unpaired.fq.gz \
    ${sample}_2.paired.fq.gz ${sample}_2.unpaired.fq.gz \
    ILLUMINACLIP:${adapters_file}:2:30:10
    """
}

//Process bwa_mem2
process bwa_mem2 {
    label 'bwa_mem2'
    publishDir "${params.outdir}/bwa_alignment/",mode: 'copy'

    input:
    tuple val(sample), path(paired_reads)
    path index_files
    val index_prefix

    output:
    tuple val(sample), path("${sample}.bam")

    script:
    """
    bwa-mem2 mem -t 2 {index_prefix} ${paired_reads[0]} ${paired_reads[1]} > ${sample}.sam
    """
}

// WORKFLOW
          
workflow {
    
    read_pairs_ch.view()
    
    fastqc(read_pairs_ch)
       
    trimmomatic(read_pairs_ch, adapter_ch)
 
    bwa_mem2(trimmomatic.out.paired_reads, genome_ch, ref_prefix)
}

