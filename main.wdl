import "single_sample_haplotypecaller.wdl" as single_sample_haplotypecaller
import "fastqc.wdl" as fastqc
import "report.wdl" as reporting

workflow PairedHaplotypecallerAndVepWorkflow {
    # This workflow is a 'super' workflow that parallelizes
    # HaplotypeCaller and VEP analysis over multiple samples.

    # Input files
    Array[File] r1_files
    Array[File] r2_files
    Boolean use_dedup

    Array[Pair[File, File]] fastq_pairs = zip(r1_files, r2_files)
    
    # Reference files
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_bwt
    File ref_sa
    File ref_amb
    File ref_ann
    File ref_pac

    # Inputs for GATK BQSR
    File dbsnp
    File dbsnp_index
    File known_indels
    File known_indels_index

    # Inputs for HaplotypeCaller scatter
    File contig_list
    Array[String] contigs = read_lines(contig_list)

    # Inputs for VEP
    String vep_species
    File vep_cache_tar
    
    # Other
    String output_zip_name
    String genome
    String git_repo_url
    String git_commit_hash

    scatter(item in fastq_pairs){
        call fastqc.run_fastqc as fastqc_for_read1 {
            input:
                fastq = item.left
        }
        
        call fastqc.run_fastqc as fastqc_for_read2 {
            input:
                fastq = item.right
        }

        call single_sample_haplotypecaller.SingleSampleHaplotypecallerWorkflow as single_sample_process {
            input:
                r1_fastq = item.left,
                r2_fastq = item.right,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                ref_bwt = ref_bwt,
                ref_sa = ref_sa,
                ref_amb = ref_amb,
                ref_ann = ref_ann,
                dbsnp = dbsnp,
                dbsnp_index = dbsnp_index,
                known_indels = known_indels,
                known_indels_index = known_indels_index,
                contigs = contigs
        }
    }

    call reporting.create_multi_qc as multiqc {
        input:
            alignment_metrics = single_sample_process.alignment_metrics,
            dedup_metrics = single_sample_process.deduplication_metrics,
            r1_fastqc_zips = fastqc_for_read1.fastqc_zip,
            r2_fastqc_zips = fastqc_for_read2.fastqc_zip
    }

    call reporting.generate_report as generate_report {
        input:
            r1_files = r1_files,
            r2_files = r2_files,
            genome = genome,
            git_commit_hash = git_commit_hash,
            git_repo_url = git_repo_url,
    }

    call zip_results {
        input:
            zip_name = output_zip_name,
            bam_files = single_sample_process.bam,
            vcf_files = single_sample_process.vcf,
            tsv_files = single_sample_process.annotated_vcf,
            vep_stats_files = single_sample_process.annotated_vcf_stats,
            multiqc_results = multiqc.report,
            analysis_report = generate_report.report
    }
}

task zip_results {
    String zip_name
    Array[File] bam_files
    Array[File] vcf_files
    Array[File] tsv_files
    Array[File] vep_stats_files
    File multiqc_results
    File analysis_report

    # runtime parameters
    Int disk_size = 1000

    command {
        mkdir alignments
        mv -t alignments ${sep=" " bam_files}
        zip -r "${zip_name}.alignments.zip" alignments
        
        mkdir output
        mkdir output/VCFs
        mkdir output/qc
        mkdir output/TSVs
        mkdir output/variant_stats
        mv ${multiqc_results} output
        mv ${analysis_report} output
        mv -t output/VCFs ${sep=" " vcf_files}
        mv -t output/TSVs ${sep=" " tsv_files}
        mv -t output/variant_stats ${sep=" " vep_stats_files}
        zip -r "${zip_name}.report_and_output.zip" output
    }

    output {
        File zip_out = "${zip_name}.report_and_output.zip"
        File alignments_zip = "${zip_name}.alignments.zip"
    }

    runtime {
        docker: "docker.io/hsphqbrc/gatk-variant-detection-workflow-tools:1.1"
        cpu: 2
        memory: "6 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}