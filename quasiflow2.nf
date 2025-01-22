#!/usr/bin/env nextflow
/*
========================================================================================
                              Q U A S I F L O W  P I P E L I N E
========================================================================================
              
              A Nextflow pipeline for analysis of NGS-based HIV Drug Resistance data
               
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    ============================================================
    nedjoni/QuasiFlow2  ~  version ${params.version}
    ============================================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run QuasiFlow2 --reads <path to fastq files> --outdir <path to output directory>
    
    Optional arguments:
      --reads                        Path to input data (must be surrounded with quotes, default is the folder "fastq")
      --outdir                       Path to directory where results will be saved (default - results). 

    HyDRA arguments (optional):
      --reporting_threshold          Minimum mutation frequency percent to report.
      --consensus_pct                The minimum percentage of a base needs to be incorporated into the consensus sequence.
      --length_cutoff                Reads that fall short of the specified length will be filtered out.
      --score_cutoff                 Reads that have a median or mean quality score (depending on the score type specified) less than the score cutoff value will be filtered out.
      --min_variant_qual             Minimum quality for the variant to be considered later on in the pipeline.
      --min_dp                       Minimum required read depth for the variant to be considered later on in the pipeline.
      --min_ac                       The minimum required allele count for the variant to be considered later on in the pipeline.
      --min_freq                     The minimum required frequency for a mutation to be considered in the drug resistance report.

     Other arguments (optional):
      --overwrite                    Set to true to overwrite previous reports (default - false).
      --name                         Name of the run.
      --email                        Email to receive notification once the run is done.

    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Header log info
log.info "============================================================"
log.info " nedjoni/QuasiFlow2   ~  version ${params.version}"
log.info "============================================================"
log.info "  Use parameter - help for the full list of parameters"
log.info "************************************************************"
def summary = [:]
summary['Run Name'] = custom_runName ?: workflow.runName
summary['Reads'] = params.reads
summary['Output directory'] = params.outdir
summary['Reporting threshold'] = params.reporting_threshold
summary['Consensus percentage'] = params.consensus_pct
summary['Length cut off'] = params.length_cutoff
summary['Score cutoff'] = params.score_cutoff
summary['Minimum variant quality'] = params.min_variant_qual
summary['Minimum depth'] = params.min_dp
summary['Minimum allele count'] = params.min_ac
summary['Minimum mutation frequency'] = params.min_freq
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home'] = "$HOME"
summary['Current user'] = "$USER"
summary['Current path'] = "$PWD"
summary['Script dir'] = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(params.email) {
    summary['E-mail Address'] = params.email
}
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

process runFastQC {
    tag "${pairId}"
    publishDir "${params.outdir}/fastqc", mode: "copy", overwrite: false

    input:
    tuple val(pairId), path(in_fastq)

    output:
    tuple val(pairId), path("*.zip"), path("*.html")

    script:
    """
    # mkdir -p ${pairId}_fastqc
    fastqc ${in_fastq[0]} ${in_fastq[1]}
    """
}

process runMultiQC {
    publishDir "${params.outdir}", mode: "copy", overwrite: false

    input:
    path fastqc_results

    output:
    path("raw_reads_multiqc_report.html")

    script:
    """
    multiqc ${fastqc_results} -o .
    mv multiqc_report.html raw_reads_multiqc_report.html
    """
}

process runTrimGalore {
    tag "${pairId}"
    publishDir "${params.outdir}/adaptors-trimmed-reads", mode: "copy", overwrite: false

    input:
    tuple val(pairId), path(in_fastq)

    output:
    tuple val(pairId), path("*.fq") // Return the trimmed reads output

    script:
    """
    trim_galore --dont_gzip -q 30 --paired ${in_fastq[0]} ${in_fastq[1]} -o .
    """
}

process runHydra {
    tag "${pairId}"
    publishDir params.outdir, mode: "copy", overwrite: false

    input:
    tuple val(pairId), path(trimmed_reads)

    output:
    tuple val(pairId), path("consensus_${pairId}.fasta"), // Output JSON file
          path("dr_report_${pairId}.csv"),
          path("mutation_report_${pairId}.aavf"),
          path("filtered_${pairId}.fastq")

    script:
    """
    # Run Hydra
    quasitools hydra \
        ${trimmed_reads[0]} ${trimmed_reads[1]} \
        -o . \
        --generate_consensus \
        --reporting_threshold ${params.reporting_threshold} \
        --consensus_pct ${params.consensus_pct} \
        --length_cutoff ${params.length_cutoff} \
        --score_cutoff ${params.score_cutoff} \
        --min_variant_qual ${params.min_variant_qual} \
        --min_dp ${params.min_dp} \
        --min_ac ${params.min_ac} \
        --min_freq ${params.min_freq}

    mv consensus.fasta consensus_${pairId}.fasta
    mv dr_report.csv dr_report_${pairId}.csv
    mv mutation_report.aavf mutation_report_${pairId}.aavf
    mv filtered.fastq filtered_${pairId}.fastq
    """
}

// Main workflow block
workflow {
    // Define input channel for read pairs
    reads_channel = Channel.fromFilePairs(params.reads)

    // Check if the reads_channel is empty
    reads_channel.view { "Read pairs: $it" } // Optional, for debugging purposes
    reads_channel.count().subscribe { count ->
        if (count == 0) {
            error "Cannot find any reads matching: ${params.reads}"
        }
    }
    
    // Step 1: Run FastQC for each read pair
    fastqc_results_channel = runFastQC(reads_channel)

    // Step 2: Run MultiQC after collecting all FastQC results
    multiqc_input = fastqc_results_channel \
        .map { it[1] } \
        .collect()
    runMultiQC(multiqc_input)
    
    // Step 3: Run Trim Galore for adapter and quality trimming
    trimmed_reads_channel = runTrimGalore(reads_channel)
    trimmed_reads_channel.view { "Trimmed reads: $it" }
    
    // Step 4: Run Hydra for mapping, consensus sequences, and mutation and drug resistance reports
    hydra_results_channel = runHydra(trimmed_reads_channel)
    hydra_results_channel.view { "Hydra results: $it" }
}