import "single_sample_haplotypecaller.wdl" as single_sample_haplotypecaller
import "report.wdl" as reporting

workflow PairedHaplotypecallerAndVepWorkflow {
    # This workflow is a 'super' workflow that parallelizes
    # HaplotypeCaller and VEP analysis over multiple samples.

    # Input files
    Array[File] r1_files
    Array[File] r2_files
    Boolean use_dedup
    
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

        call qc.run_alignment_metrics as aln_metrics {
            input:
                input_bam = single_sample_process.bam,
                input_bam_index = single_sample_process.bam_index,
                sample_name = sample_name,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict
        }
    }

    call reporting.create_multi_qc as multiqc {
        input:
            alignment_metrics = aln_metrics.alignment_metrics,
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
}

task zip_results {
    String zip_name


    command {

    }

    output {

    }

    runtime {

    }
}