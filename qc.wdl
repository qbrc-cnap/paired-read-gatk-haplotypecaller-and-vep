workflow test_fastqc {
    File fastq

    call run_fastqc{
        input:
            fastq=fastq
    }

    output {
        File fqc_zip = run_fastqc.fastqc_zip
    }
}

task run_fastqc {
    File fastq

    String bn = basename(fastq)
    String target_zip = sub(bn, "\\.fastq\\.gz", "_fastqc.zip")
    Int disk_size = 100

    command {
        fastqc ${fastq} -o .
    }

    output {
        File fastqc_zip = "${target_zip}"
    }

    runtime {
        docker: "docker.io/blawney/star_rnaseq:v0.0.1"
        cpu: 2
        memory: "4 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}

task run_alignment_metrics {
    File input_bam
    File input_bam_index
    String sample_name = basename(input_bam, ".bam")

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    # Runtime parameters
    Int disk_dize = 100

    command {
        java -jar $PICARD_JAR \
            CollectAlignmentSummaryStatistics \
            R=${ref_fasta} \
            I=${input_bam} \
            O=${sample_name}.alignment_metrics.txt;
    }

    output {
        File alignment_metrics = "${sample_name}.alignment_metrics.txt"
    }

    runtime {
        docker: "docker.io/blawney/star_rnaseq:v0.0.1"
        cpu: 2
        memory: "4 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}