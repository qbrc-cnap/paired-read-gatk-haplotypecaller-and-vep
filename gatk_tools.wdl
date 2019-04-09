workflow gatk_test {
    File input_bam
    File input_bam_index
    String sample_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File dbsnp
    File dbsnp_index
    File known_indels
    File known_indels_index
    Array[String] contigs

    call base_recalibrator {
        input:
            input_bam = input_bam,
            input_bam_index = input_bam_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            dbsnp = dbsnp,
            dbsnp_index = dbsnp_index,
            known_indels = known_indels,
            known_indels_index = known_indels_index
    }

    call apply_recalibration {
        input:
            input_bam = input_bam,
            input_bam_index = input_bam_index,
            recalibration_report = base_recalibrator.recalibration_report,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict
    }

    scatter (scatter_interval in contigs) {
        # Identifies variants with GATK's HaplotypeCaller
        call haplotypecaller {
            input:
                input_bam = apply_recal.recalibrated_bam,
                input_bam_index = apply_recal.recalibrated_bam_index,
                interval = scatter_interval,
                gvcf_name = sample_name,
                ref_dict = ref_dict,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
        }
    }
    
    # Merges the scattered VCFs together
    call merge_vcf {
        input:
            input_vcfs = haplotypecaller.output_vcf,
            sample_name = sample_name,
            ref_dict = ref_dict,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index
    }

    output {
        File vcf = merge_vcf.output_vcf
    }
}

task base_recalibrator {
    File input_bam
    File input_bam_index
    String sample_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File dbsnp
    File dbsnp_index
    File known_indels
    File known_indels_index

    command {
        java -Xmx4000m -jar ${GATK} \
            BaseRecalibrator \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -known-sites ${dbsnp} \
            -known-sites ${known_indels} \
            -O recal_data.table;
    }
    
    output {
        File recalibration_report = "recal_data.table"
    }
}

task apply_recalibration {
    File input_bam
    File input_bam_index
    File recalibration_report
    String sample_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict

    command {
        java -Xmx4000m -jar ${GATK} \
            ApplyBQSR \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${sample_name}.bqsr.bam \
            -bqsr-recal-file ${recalibration_report};
        samtools index ${sample_name}.bqsr.bam;
    }

    output {
        File recalibrated_bam = "${sample_name}.bqsr.bam"
    }
}

task haplotypecaller {
    File input_bam
    File input_bam_index
    String sample_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    String interval

    command {
        java -Xmx8000m -jar ${GATK} \
            HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${sample_name}.vcf \
            -L ${interval};
    }

    output {
        File output_vcf = "${sample_name}.vcf"
    }
}

task merge_vcf {
    File input_vcfs
    String sample_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict

    command {
        java -Xmx3000m -jar ${PICARD_JAR} \
            MergeVcfs \
            INPUT=${sep=' INPUT=' input_vcfs} \
            OUTPUT=${sample_name}.vcf
    }

    output {
        File output_vcf = "${sample_name}.vcf"
    }
}