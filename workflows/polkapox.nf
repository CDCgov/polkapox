/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowPolkapox.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.kraken_db, params.multiqc_config, params.fasta, params.fai]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { 
    ch_input = Channel.fromPath("${params.input}", type: 'file', checkIfExists: true) 
    }
else if (params.indir) { 
    ch_indir = file(params.indir)
    } 
else if (params.sra_ids) {
    ch_sra_id = file(params.sra_ids)
    }
//else { exit 1, 'Must specify samplesheet, input directory of fastq files, or sra id list!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK         } from '../subworkflows/local/input_check'
include { PREPARE_GENOME      } from '../subworkflows/local/prepare_genome'
include { SRA_TOOLS           } from '../subworkflows/nf-core/sra_tools'
include { CREATE_SAMPLESHEET  } from '../modules/local/create_samplesheet/create_samplesheet'
include { READ_FILTER         } from '../subworkflows/local/filter_reads'
include { DENOVO              } from '../subworkflows/local/denovo'
include { REFBASED            } from '../subworkflows/local/ref_based'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                                        } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                       } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                   } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { SUMMARIZE_QC                                  } from '../modules/local/summarize_qc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow POLKAPOX {

    ch_versions = Channel.empty() 

    if (params.indir) {
        //
        // MODULE: Create samplesheet based on directory of fastq files 
        //
        CREATE_SAMPLESHEET (
            ch_indir
        )
        ch_input = CREATE_SAMPLESHEET.out.samplesheet
        ch_versions = ch_versions.mix(CREATE_SAMPLESHEET.out.versions)
    }

    else if (params.sra_ids) {
        // This is modified from fetchngs https://github.com/nf-core/fetchngs/tree/master
        // get SRA ch for each accession
        ch_sra = Channel.fromPath(params.sra_ids)
            .splitCsv ( header: false )
            .flatten()
            .map { 
                accession -> [meta: accession, accession: accession] 
                }  // creating a map with id and accession
        //
        // Subworkflow: Create samplesheet from list of SRA Accessions 
        //
        SRA_TOOLS (
            ch_sra
        )
        ch_versions = ch_versions.mix(SRA_TOOLS.out.versions.first())
    
        ch_reads = SRA_TOOLS.out.forward.join(SRA_TOOLS.out.reverse)
        ch_input = ch_reads
            .map { meta, fastq1, fastq2 -> "$meta,$fastq1,$fastq2" }
            .collectFile(name: 'sra-samplesheet.csv', newLine: true, seed: 'sample,fastq_1,fastq_2')
    }
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)            

    //
    // SUBWORKFLOW: Prepare reference
    //
    PREPARE_GENOME ()
    ch_bwa_index = PREPARE_GENOME.out.bwa_index
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // SUBWORKFLOW: Only run Read Filter
    //

    if ( params.workflow == 'filter_reads'
        || params.workflow == 'ref_based'
        || params.workflow == 'denovo'
        || params.workflow == 'full' ) {
        READ_FILTER (
            INPUT_CHECK.out.reads,
        )
        ch_versions = ch_versions.mix(READ_FILTER.out.versions)
    }

    //
    // SUBWORKFLOW: Run Read Filter + Reference-based Assembly
    //

    if ( params.workflow == 'ref_based' || params.workflow == 'full' ) {
        REFBASED (
            READ_FILTER.out.trimmed_fastq,
            ch_bwa_index
        )
        ch_versions = ch_versions.mix(REFBASED.out.versions)

    }

    //
    // SUBWORKFLOW: Run Read Filter + Denovo Assembly
    //

    if ( params.workflow == 'denovo' || params.workflow == 'full' ) {
        DENOVO (
            READ_FILTER.out.trimmed_fastq,
        )
        ch_versions = ch_versions.mix(DENOVO.out.versions)
    }

    //
    // MODULE: Dump all software versions used
    //

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowPolkapox.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(READ_FILTER.out.json.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(READ_FILTER.out.kraken2_report.collect{it[1]}.ifEmpty([]))
    if ( params.workflow == 'ref_based' || params.workflow == 'full' ) {
        ch_multiqc_files = ch_multiqc_files.mix(REFBASED.out.flagstat.collect{it[1]}.ifEmpty([]))
    }
    if ( params.workflow == 'denovo' || params.workflow == 'full' ) {
        ch_multiqc_files = ch_multiqc_files.mix(DENOVO.out.quast_tsv.collect().ifEmpty([]))
    }

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)

    //
    // MODULE: summarize_qc.py
    //

    ch_summarizeqc_files = Channel.empty()
    ch_summarizeqc_files = ch_summarizeqc_files.mix(READ_FILTER.out.json.collect{it[1]}.ifEmpty([]))
    ch_summarizeqc_files = ch_summarizeqc_files.mix(READ_FILTER.out.kraken2_report.collect{it[1]}.ifEmpty([]))
    ch_summarizeqc_files = ch_summarizeqc_files.mix(READ_FILTER.out.classified_reads_assignment.collect{it[1]}.ifEmpty([]))
    ch_summarizeqc_files = ch_summarizeqc_files.mix(READ_FILTER.out.seqtk_reads.collect{it[1]}.ifEmpty([]))
    ch_summarizeqc_files = ch_summarizeqc_files.mix(MULTIQC.out.data.collect().ifEmpty([]))
    if ( params.workflow == 'ref_based' || params.workflow == 'full') {
        ch_summarizeqc_files = ch_summarizeqc_files.mix(REFBASED.out.depth_tsv.collect{it[1]}.ifEmpty([]))
        ch_summarizeqc_files = ch_summarizeqc_files.mix(REFBASED.out.tsv_vars.collect{it[1]}.ifEmpty([]))
        ch_summarizeqc_files = ch_summarizeqc_files.mix(REFBASED.out.ivar_tsv.collect{it[1]}.ifEmpty([]))
    }
    if ( params.workflow == 'denovo' || params.workflow == 'full' ) {
        ch_summarizeqc_files = ch_summarizeqc_files.mix(DENOVO.out.flagstat.collect{it[1]}.ifEmpty([]))
        ch_summarizeqc_files = ch_summarizeqc_files.mix(DENOVO.out.graph_recon_log.collect().ifEmpty([]))
        ch_summarizeqc_files = ch_summarizeqc_files.mix(DENOVO.out.gfa_assembly.collect{it[1]}.ifEmpty([]))
        ch_summarizeqc_files = ch_summarizeqc_files.mix(DENOVO.out.mummer_summary.collect{it[1]}.ifEmpty([]))
        ch_summarizeqc_files = ch_summarizeqc_files.mix(DENOVO.out.fasta.collect{it[1]}.ifEmpty([]))
        ch_summarizeqc_files = ch_summarizeqc_files.mix(DENOVO.out.mpileup.collect{it[1]}.ifEmpty([]))
    }

    SUMMARIZE_QC (
        ch_summarizeqc_files.collect(),
        ch_input,
   )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
