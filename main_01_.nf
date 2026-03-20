#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// -------------------------
// PARAMETERS
// -------------------------
params.reads = 'data/*_{1,2}.fq.gz'       // paired-end reads
params.outdir = 'outputs'                // output folder for all results
params.adapters = 'adapters.fa'          // adapters file in repo
params.genome = 'genome.fa'

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

// channel for the genome file
genome_ch = Channel.fromPath(params.genome)

// -------------------------
// PROCESS: FASTQC
// -------------------------
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
    trimmed_ch = trimmomatic(read_pairs_ch, adapter_ch)
      LIST OF PARAMETERS
================================
Reads            : data/*_{1,2}.fq.gz
Output-folder    : outputs
Adapters         : adapters.fa

[-        ] fastqc      -
[-        ] trimmomatic -
ERROR ~ No such variable: trimmed_ch

 -- Check script 'main_01_.nf' at line: 118
      LIST OF PARAMETERS
================================
Reads            : data/*_{1,2}.fq.gz
Output-folder    : outputs
Adapters         : adapters.fa

[-        ] fastqc      -
[-        ] trimmomatic -
ERROR ~ No such variable: trimmed_ch

 -- Check script 'main_01_.nf' at line: 118
      LIST OF PARAMETERS
================================
Reads            : data/*_{1,2}.fq.gz
Output-folder    : outputs
Adapters         : adapters.fa

[-        ] fastqc      -
[-        ] trimmomatic -
ERROR ~ No such variable: trimmed_ch

 -- Check script 'main_01_.nf' at line: 118
      LIST OF PARAMETERS
================================
Reads            : data/*_{1,2}.fq.gz
Output-folder    : outputs
Adapters         : adapters.fa

[-        ] fastqc      -
[-        ] trimmomatic -
ERROR ~ No such variable: trimmed_ch

 -- Check script 'main_01_.nf' at line: 118
      LIST OF PARAMETERS
================================
Reads            : data/*_{1,2}.fq.gz
Output-folder    : outputs
Adapters         : adapters.fa

[-        ] fastqc      -
[-        ] trimmomatic -
ERROR ~ No such variable: trimmed_ch

 -- Check script 'main_01_.nf' at line: 118
      LIST OF PARAMETERS
================================
Reads            : data/*_{1,2}.fq.gz
Output-folder    : outputs
Adapters         : adapters.fa

[-        ] fastqc      -
[-        ] trimmomatic -
ERROR ~ No such variable: trimmed_ch

 -- Check script 'main_01_.nf' at line: 118
      LIST OF PARAMETERS
================================
Reads            : data/*_{1,2}.fq.gz
Output-folder    : outputs
Adapters         : adapters.fa

[-        ] fastqc      -
[-        ] trimmomatic -
ERROR ~ No such variable: trimmed_ch

 -- Check script 'main_01_.nf' at line: 118
      LIST OF PARAMETERS
================================
Reads            : data/*_{1,2}.fq.gz
Output-folder    : outputs
Adapters         : adapters.fa

[-        ] fastqc      -
[-        ] trimmomatic -
ERROR ~ No such variable: trimmed_ch

 -- Check script 'main_01_.nf' at line: 118
      LIST OF PARAMETERS
================================
Reads            : data/*_{1,2}.fq.gz
Output-folder    : outputs
Adapters         : adapters.fa

[-        ] fastqc      -
[-        ] trimmomatic -
ERROR ~ No such variable: trimmed_ch

 -- Check script 'main_01_.nf' at line: 118

    script:
    """
    trimmomatic PE -phred33 \
    ${reads[0]} ${reads[1]} \
    ${sample}_1.paired.fq.gz ${sample}_1.unpaired.fq.gz \
    ${sample}_2.paired.fq.gz ${sample}_2.unpaired.fq.gz \
    ILLUMINACLIP:${adapters_file}:2:30:10
    """
}

//PROCESS: ALIGNMENT

process alignment {    label 'alignment'
    publishDir "${params.outdir}/aligned", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(reads)
    path genome

    output:
    file "${sample}.sorted.bam"

    script:
    """
    bwa-men2 mem -t 4 ${genome} ${reads[0]} ${reads[1]} | \
    samtools sort -@ 4 -o ${sample}.sorted.bam
    """
}

// PROCESS: 
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

    // Run Alignment
    paired_reads_ch = trimmed_ch.map { file ->
       def sample = file.name.split('_')[0..-2].join('_')
       tuple(sample, file)
    }.groupTuple()
   
    alignment(paired_reads_ch, genome_ch)
}

Trial Changes
