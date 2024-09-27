#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.project = "EcoliRNASeq"

params.nbcpu = "8"

params.genome_name = "Ecoli_K12"


workflow {

    // Initialise les channels
    bio_conds = Channel.of([1,"WT_1"],[2,"WT_2"],[3,"MazF_1"],[4,"MazF_2"])
    ftp_ebis = Channel.of([1,"005/ERR2686025/ERR2686025"], [2,"006/ERR2686026/ERR2686026"], [3,"007/ERR2686027/ERR2686027"], [4,"008/ERR2686028/ERR2686028"])
    // channel combiné
    bio_cond_ftp_ebis = bio_conds.cross(ftp_ebis)

    // Telechargement du genome et indexation
    DL_GENOME | INDEX_GENOME

    // Telecharger les reads
    DL_FASTQ(bio_cond_ftp_ebis)
    
    // Trim des adapters
    CUTADAPT(bio_conds, DL_FASTQ.out)
    FASTQC(CUTADAPT.out)

    // Mapping et convert to bam
    MAPPING(INDEX_GENOME.out, CUTADAPT.out)
    //MAPPING.out.view()
    SAM_TO_BAM(MAPPING.out)

    // Comptage des reads sur chaque gene
    GENE_COUNT(SAM_TO_BAM.out, DL_GENOME.out.genome_gtf)

    // Final run avec multiqc
    FASTQC.out.collect().view()
    MULTIQC(FASTQC.out.collect(), GENE_COUNT.out.collect())

}

/*
 * Download genome files
 * We have now to download the genome file of Escherichis coli K-12. We download first the sequence (fasta file) and then the annotation (gff file).
 */
process DL_GENOME {
    publishDir params.genome

    output: 
    path "${params.genome_name}.fna", emit: genome_fna
    path "${params.genome_name}.gtf", emit: genome_gtf

    script:
    """
    fna_file=${params.genome_name}.fna
    gtf_file=${params.genome_name}.gtf
    echo "Downloaded fasta file" ${params.fna_ftp}
    wget -c -O \${fna_file}.gz ${params.fna_ftp}
    gunzip -d \${fna_file}.gz
    echo "Downloaded fasta file" ${params.gtf_ftp}   
    wget -c -O \${gtf_file}.gz ${params.gtf_ftp}
    gunzip -d \${gtf_file}.gz
    """
}

process INDEX_GENOME {
    publishDir params.genome

    input:
    path genome_fna
    path genome_gtf

    output:
    tuple val(params.genome_name), path("${params.genome_name}.*.bt2")
    
    script:
    """
    bowtie2-build --threads $params.nbcpu ${genome_fna} ${params.genome_name}
    """
}

process DL_FASTQ {
    publishDir params.rawdata

    input:
    tuple val(bio_cond), val(ftp_ebi)

    output:
    tuple path("${bio_cond[1]}_R1.fastq.gz"), path("${bio_cond[1]}_R2.fastq.gz")
    
    script:
    """
    data=${bio_cond[1]}_R1.fastq.gz
    if [ ! -f \$data ]; then    
        wget -O \$data ${params.ebi_fq}${ftp_ebi[1]}_1.fastq.gz
    fi
    data=${bio_cond[1]}_R2.fastq.gz
    if [ ! -f \$data ]; then    
        wget -O \$data ${params.ebi_fq}${ftp_ebi[1]}_2.fastq.gz
    fi
    """
}


process CUTADAPT {
    publishDir params.fastq

    input:
    val bio_cond
    tuple path(data_R1), path(data_R2)
    
    output:
    tuple val(bio_cond), path("${bio_cond[1]}_R1.fastq"), path("${bio_cond[1]}_R2.fastq")

    script:
    """
    cutadapt -j $params.nbcpu -b $params.read1Adapt -B $params.read2Adapt -o ${bio_cond[1]}_R1.fastq \
        -p ${bio_cond[1]}_R2.fastq ${data_R1} ${data_R2}
    """   
}

process FASTQC {
    publishDir params.fastqc

    input:
    tuple val(bio_cond), path(fastq_R1), path(fastq_R2)

    output:
    tuple path("${bio_cond[1]}_R1_fastqc*"), path("${bio_cond[1]}_R2_fastqc*")

    script:
    """
    fastqc -t $params.nbcpu ${fastq_R1} ${fastq_R2}
    """
}


process MAPPING {
    publishDir params.mapping

    input:
    tuple val(index_path), path(genome_indexes)
    tuple val(bio_cond), path(fastq_R1), path(fastq_R2)
    
    output:
    tuple val(bio_cond), path("${bio_cond[1]}.sam"), path("${bio_cond[1]}.log")

    script:
    """
    bowtie2 -p $params.nbcpu -x ${index_path} -1 ${fastq_R1} -2 ${fastq_R2} -S ${bio_cond[1]}.sam 2> ${bio_cond[1]}.log
    """
}

process SAM_TO_BAM {
    publishDir params.mapping

    input:
    tuple val(bio_cond), path(sam_file), path(log_file)

    output:
    tuple path("${bio_cond[1]}.bam"), path("${bio_cond[1]}.bam.bai")

    script:
    """
    samtools view --threads $params.nbcpu -b -q 1 ${sam_file} > ${bio_cond[1]}_raw.bam
    samtools sort --threads $params.nbcpu -o ${bio_cond[1]}.bam ${bio_cond[1]}_raw.bam
    samtools index -@ $params.nbcpu ${bio_cond[1]}.bam
    rm ${bio_cond[1]}_raw.bam
    """
}

process GENE_COUNT {
    publishDir params.count
    
    input:
    tuple path(bam_file), path(bai_file)
    path genome_gtf

    output:
    tuple path("${bam_file}.txt"), path("${bam_file}.txt.summary")

    script:
    """
    featureCounts -T $params.nbcpu -p -t "gene" -a ${genome_gtf} -o ${bam_file}.txt ${bam_file}
    """
}


process MULTIQC {
    publishDir params.multiqc

    input:
    tuple path(fastqc_R1), path(fastqc_R2)
    tuple path(gene_count), path(gene_count_summary)

    script:
    """
    multiqc --filename "MultiQC_${params.project}" ${projectDir}
    """
}
