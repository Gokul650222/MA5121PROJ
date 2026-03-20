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

