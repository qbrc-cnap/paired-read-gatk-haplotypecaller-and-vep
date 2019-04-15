workflow test_qc {
    File fastq

    call run_fastqc{
        input:
            fastq=fastq
    }

    output {
        File fqc_zip = run_fastqc.fastqc_zip
    }
}

task create_multi_qc {
    Array[File] alignment_metrics
    Array[File] dedup_metrics
    Array[File] r1_fastqc_zips
    Array[File] r2_fastqc_zips

    # Runtime parameters
    Int disk_size = 100

    command {
        multiqc .
    }

    output {
        File report = "multiqc_report.html"
    }

    runtime {
        docker: "docker.io/hsphqbrc/gatk-variant-detection-workflow-tools:1.1"
        cpu: 2
        memory: "4 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}