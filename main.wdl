import "single_sample_haplotypecaller.wdl" as single_sample_haplotypecaller

workflow PPairedHaplotypecallerAndVepWorkflow {
    # This workflow is a 'super' workflow that parallelizes
    # HaplotypeCaller and VEP analysis over multiple samples.
    Array[File] r1_files
    Array[File] r2_files
    # Inputs for alignment
    String genome
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_bwt
    File ref_sa
    File ref_amb
    File ref_ann
    File ref_pac
    File dbsnp
    File dbsnp_index
    File known_indels
    File known_indels_index
    # Inputs for HaplotypeCaller scatter
    File scattered_calling_intervals_list_file
    Array[String] scattered_calling_intervals = read_lines(scattered_calling_intervals_list_file)
    # Inputs for vep
    String vep_species
    File vep_cache_tar
    String output_zip_name
    String git_repo_url
    String git_commit_hash
    # Recommended sizes:
    # small_disk = 200
    # medium_disk = 300
    # large_disk = 400
    # preemptible_tries = 3
    Int small_disk
    Int medium_disk
    Int large_disk
    Int preemptible_tries

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
                
        }
    }
}