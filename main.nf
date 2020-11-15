#!/usr/bin/env nextflow

params.project = "EcoliRNASeq"

params.nbcpu = "8"

params.bio_cond = ["WT_1","WT_2","MazF_1","MazF_2"]
params.ftp_ebi = ["005/ERR2686025/ERR2686025" , "006/ERR2686026/ERR2686026", "007/ERR2686027/ERR2686027", "008/ERR2686028/ERR2686028"]
params.ebi_fq = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR268/"

params.fna_ftp = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
params.gtf_ftp = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gtf.gz"
params.genome_name = "Ecoli_K12"


Channel
    .fromList(params.bio_cond)
    .set{bio_conds}
Channel
    .fromList(params.ftp_ebi)
    .set{ftp_ebis}


process init_folder {
    
    output:
    stdout into folderCreated

    script:
    """
    if [ ! -d ${params.rawdata} ]; then
        mkdir ${params.rawdata}
    fi
    if [ ! -d ${params.genome} ]; then
        mkdir ${params.genome}
    fi
    if [ ! -d ${params.fastq} ]; then
        mkdir ${params.fastq}
    fi
    if [ ! -d ${params.fastqc} ]; then
        mkdir ${params.fastqc}
    fi
    if [ ! -d ${params.mapping} ]; then
        mkdir ${params.mapping}
    fi
    if [ ! -d ${params.count} ]; then
        mkdir ${params.count}
    fi
    if [ ! -d ${params.multiqc} ]; then
        mkdir ${params.multiqc}
    fi
    """
}


/*
 * Download genome files
 * We have now to download the genome file of Escherichis coli K-12. We download first the sequence (fasta file) and then the annotation (gff file).
 */
process dlGenome {
    echo true

    input:
    stdin from folderCreated
    
    output: 
    val "${params.genome}/${params.genome_name}" into dl_genome
    val "${params.genome}/${params.genome_name}.gtf" into genome_annot

    script:
    """
    fna_file=${params.genome}/${params.genome_name}.fna
    gtf_file=${params.genome}/${params.genome_name}.gtf
    if [ ! -f \$fna_file ]; then
        echo "Downloaded fasta file" ${params.fna_ftp}
        wget -c -O \${fna_file}.gz ${params.fna_ftp}
        gunzip -d \${fna_file}.gz
    fi
    if [ ! -f \$gtf_file ]; then 
        echo "Downloaded fasta file" ${params.gtf_ftp}   
        wget -c -O \${gtf_file}.gz ${params.gtf_ftp}
        gunzip -d \${gtf_file}.gz
    fi
    """
}


process index_genome {
    input:
    val genome_path from dl_genome

    output:
    val "${genome_path}" into genome_index
    """
    data=${genome_path}.1.bt2
    if [ ! -f \$data ]; then
        bowtie2-build --threads $params.nbcpu ${genome_path}.fna ${genome_path}
    fi
    """
}


/*
 * Download fastq from EBI
 * Test for eachf ile if the dataset was not already downloaded 
 */
 process dlFastQ {
    echo true

    input:
    stdin from folderCreated
    val bio_cond from bio_conds
    val ftp_ebi from ftp_ebis

    output:
    val "${params.rawdata}/${bio_cond}_R1.fastq.gz" into data_R1
    val "${params.rawdata}/${bio_cond}_R2.fastq.gz" into data_R2
    
    script:
    """
    data=${params.rawdata}/${bio_cond}_R1.fastq.gz
    if [ ! -f \$data ]; then    
        wget -O \$data ${params.ebi_fq}${ftp_ebi}_1.fastq.gz
    fi
    data=${params.rawdata}/${bio_cond}_R2.fastq.gz
    if [ ! -f \$data ]; then    
        wget -O \$data ${params.ebi_fq}${ftp_ebi}_2.fastq.gz
    fi  
    """
}


/*
 * Run cutadapt for all pair of reads
 */
process cutadapt {
    input:
    val bio_cond from Channel.fromList(params.bio_cond)
    val fastq_R1 from data_R1
    val fastq_R2 from data_R2
    
    output:
    val "${bio_cond}_R1" into data_cut_R1
    val "${bio_cond}_R2" into data_cut_R2
    val "${bio_cond}_R1" into rawdata_cut_R1
    val "${bio_cond}_R2" into rawdata_cut_R2

    script:
    """
    data=${params.fastq}/${bio_cond}_R1.fastq
    if [ ! -f \$data ]; then    
        cutadapt -j $params.nbcpu -b $params.read1Adapt -B $params.read2Adapt -o ${params.fastq}/${bio_cond}_R1.fastq -p ${params.fastq}/${bio_cond}_R2.fastq ${fastq_R1} ${fastq_R2}
    fi
    """   
}


/*
 * Control quality of fastq files after cutadapt
 */
process fastQC {
    input:
    val data_fq_R1 from data_cut_R1
    val data_fq_R2 from data_cut_R2

    output: 
    stdout into fastqc

    script:
    """
    data=${params.fastqc}/${data_fq_R1}_fastqc.html
    if [ ! -f \$data ]; then
        fastqc -t $params.nbcpu -o ${params.fastqc} ${data_fq_R1}.fastq
    fi
    data=${params.fastqc}/${data_fq_R2}_fastqc.html
    if [ ! -f \$data ]; then
        fastqc -t $params.nbcpu -o ${params.fastqc} ${data_fq_R2}.fastq
    fi
    
    """
}


/*
 * Run mapping of reads with bowtie2
 */
process mapping {
    input:
    val genome from genome_index
    val data_R1 from rawdata_cut_R1
    val data_R1 from rawdata_cut_R2
    val biocond from Channel.fromList(params.bio_cond)

    output:
    val biocond into sam_files

    """
    data=${params.mapping}/${biocond}.sam
    if [ ! -f \$data ]; then
        bowtie2 -p $params.nbcpu -x $genome -1 ${params.fastq}/${data_R1}.fastq -2 ${params.fastq}/${data_R1}.fastq -S ${params.mapping}/${biocond}.sam > ${params.mapping}/${biocond}.log
    fi
    """
}

/*
 * Convert sam files to bam files
 */
process sam_to_bam {
    
    echo true
    input:
    val sam_file from sam_files

    output:
    val "${sam_file}" into bam_files 

    """
    data=${params.mapping}/${sam_file}.bam
    if [ ! -f \$data ]; then
        samtools view --threads $params.nbcpu -b -q 1 ${params.mapping}/${sam_file}.sam > ${params.mapping}/${sam_file}_raw.bam
        samtools sort --threads $params.nbcpu -o ${params.mapping}/${sam_file}.bam ${params.mapping}/${sam_file}_raw.bam
        samtools index -@ $params.nbcpu ${params.mapping}/${sam_file}.bam
        rm ${params.mapping}/${sam_file}_raw.bam
    fi
    """
}

/*
 * Count number of reads for every gene
 */ 
process gene_count {
    input:
    val bam_file from bam_files
    val genome from genome_annot

    output:
    stdout into mappFinish

    script:
    """
    data=${params.count}/${bam_file}.txt
    if [ ! -f \$data ]; then
        featureCounts -T $params.nbcpu -t "gene" -a ${genome} -o ${params.count}/${bam_file}.txt ${params.mapping}/${bam_file}.bam
    fi
    """
}

/*
 * Run final qc process
 */
process multiqc {
    input:
    val stdin from mappFinish
    val stdin2 from fastqc

    script:
    """
    data=${params.multiqc}/${params.project}.html
    if [ ! -f \$data ]; then
        multiqc -f -i $params.project -n $params.project -o $params.multiqc ${params.path}/.
    fi
    """

}