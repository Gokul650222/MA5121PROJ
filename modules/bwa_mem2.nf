//Process bwa_mem2 
process bwa_mem2 {
    label 'bwa_mem2'
    publishDir "${params.outdir}/bwa_alignment/",mode: 'copy'

    input:
    tuple val(sample), path(paired_reads)
    path genome_files  
    

    output:
    tuple val(sample), path("${sample}.bam")
    
    script:
    """
    bwa-mem2 mem -t 2 ./data/LG12.fasta ${paired_reads[0]} ${paired_reads[1]} | samtools sort -@ 2 -o ${sample}.bam
    """
}

