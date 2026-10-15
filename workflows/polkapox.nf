/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowPolkapox.initialise(params, log)

// Check input path parameters to see if they exist
//def checkPathParamList = [ params.kraken_db, params.multiqc_config, params.fasta, params.fai]
def checkPathParamList = [ params.multiqc_config, params.fasta, params.fai]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
ch_indir = Channel.empty()
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
ch_multiqc_custom_config = params.multiqc_config ? file(params.multiqc_config, checkIfExists: true) : []

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK         } from '../subworkflows/local/input_check/main'
include { SRA_TOOLS           } from '../subworkflows/local/sra_tools/main'
include { CREATE_SAMPLESHEET  } from '../modules/local/create_samplesheet/main'
include { READ_FILTER         } from '../subworkflows/local/filter_reads/main'
include { DENOVO              } from '../subworkflows/local/denovo/main'
include { REFBASED            } from '../subworkflows/local/ref_based/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                                       } from '../modules/nf-core/multiqc/main'
include { SUMMARIZE_QC                                  } from '../modules/local/summarize_qc/main'
include { BWA_INDEX                                     } from '../modules/nf-core/bwa/index/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Generate workflow version string
def getWorkflowVersion() {
    def version_string = "" as String
    if (workflow.manifest.version) {
        def prefix_v = workflow.manifest.version[0] != 'v' ? 'v' : ''
        version_string += "${prefix_v}${workflow.manifest.version}"
    }
    if (workflow.commitId) {
        def git_shortsha = workflow.commitId.substring(0, 7)
        version_string += "-g${git_shortsha}"
    }
    return version_string
}

// Parse a legacy versions.yml file's contents into a trimmed YAML string
def processVersionsFromYAML(yaml_file) {
    def yaml = new org.yaml.snakeyaml.Yaml()
    def versions = yaml.load(yaml_file).collectEntries { k, v -> [k.tokenize(':')[-1], v] }
    return yaml.dumpAsMap(versions).trim()
}

// Get workflow version block for the collated versions report
def workflowVersionToYAML() {
    return """
    Workflow:
        ${workflow.manifest.name}: ${getWorkflowVersion()}
        Nextflow: ${workflow.nextflow.version}
    """.stripIndent().trim()
}

// Get a channel of software versions (legacy versions.yml files) in YAML format
def softwareVersionsToYAML(ch_versions) {
    return ch_versions.unique().map { version -> processVersionsFromYAML(version) }.unique().mix(channel.of(workflowVersionToYAML()))
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow POLKAPOX {

    ch_versions = Channel.empty() 
    ch_versions = ch_versions.mix(Channel.topic('versions')) //pick up modules migrated to topic-based version tracking

    if (params.indir) {
        //
        // MODULE: Create samplesheet based on directory of fastq files 
        //
        CREATE_SAMPLESHEET (
            ch_indir
        )
        ch_input = CREATE_SAMPLESHEET.out.samplesheet
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
        ch_versions = ch_versions.mix(SRA_TOOLS.out.versions)
    
        ch_input = SRA_TOOLS.out.reads
            .map { meta, reads -> "$meta,${reads[0]},${reads.size() > 1 ? reads[1] : ''}" }
            .collectFile(name: 'sra-samplesheet.csv', newLine: true, seed: 'sample,fastq_1,fastq_2')
    }
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )

    //
    // SUBWORKFLOW: Prepare reference
    //
    BWA_INDEX ( [id: 'reference'],  file(params.fasta) )

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
    // SUBWORKFLOW: Run Reference-based Analyses
    //

    if ( params.workflow == 'ref_based' || params.workflow == 'full' ) {
        REFBASED (
            READ_FILTER.out.trimmed_fastq,
            BWA_INDEX.out.index,
        )
        ch_versions = ch_versions.mix(REFBASED.out.versions)

    }

    //
    // SUBWORKFLOW: Run Denovo Assembly
    //

    if ( params.workflow == 'denovo' || params.workflow == 'full' ) {
        DENOVO (
            READ_FILTER.out.trimmed_fastq,
        )
        ch_versions = ch_versions.mix(DENOVO.out.versions)
    }

    //
    // Collate software versions used (replaces CUSTOM_DUMPSOFTWAREVERSIONS; handles both
    // legacy versions.yml files and topic-channel version tuples)
    //
    ch_topic_versions = ch_versions
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    ch_topic_versions_string = ch_topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    ch_collated_versions = softwareVersionsToYAML(ch_topic_versions.versions_file)
        .mix(ch_topic_versions_string)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'collated_versions.yml', sort: true, newLine: true)
    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowPolkapox.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions.collect())
    ch_multiqc_files = ch_multiqc_files.mix(READ_FILTER.out.json.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(READ_FILTER.out.kraken2_report.collect{it[1]}.ifEmpty([]))
    if ( params.workflow == 'ref_based' || params.workflow == 'full' ) {
        ch_multiqc_files = ch_multiqc_files.mix(REFBASED.out.flagstat.collect{it[1]}.ifEmpty([]))
    }
    if ( params.workflow == 'denovo' || params.workflow == 'full' ) {
        ch_multiqc_files = ch_multiqc_files.mix(DENOVO.out.quast_tsv.collect().ifEmpty([]))
    }

    // MultiQC config: bundled default plus an optional user-supplied override (nf-core/modules multiqc dedicated `multiqc_config` input)
    ch_multiqc_config_files = ch_multiqc_custom_config ? [ ch_multiqc_config, ch_multiqc_custom_config ] : ch_multiqc_config

    MULTIQC (
        ch_multiqc_files.collect().map { files -> [ [id:'multiqc'], files, ch_multiqc_config_files, [], [], [] ] }
    )
    multiqc_report = MULTIQC.out.report.map { meta, report -> report }.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)

    //
    // MODULE: summarize_qc.py
    //

    ch_summarizeqc_files = Channel.empty()
    ch_summarizeqc_files = ch_summarizeqc_files.mix(READ_FILTER.out.json.collect{it[1]}.ifEmpty([]))
    ch_summarizeqc_files = ch_summarizeqc_files.mix(READ_FILTER.out.kraken2_report.collect{it[1]}.ifEmpty([]))
    ch_summarizeqc_files = ch_summarizeqc_files.mix(READ_FILTER.out.classified_reads_assignment.collect{it[1]}.ifEmpty([]))
    ch_summarizeqc_files = ch_summarizeqc_files.mix(READ_FILTER.out.orthopox_reads.collect{it[1]}.ifEmpty([]))
    ch_summarizeqc_files = ch_summarizeqc_files.mix(MULTIQC.out.data.map { meta, data -> data }.collect().ifEmpty([]))
    if ( params.workflow == 'ref_based' || params.workflow == 'full') {
        ch_summarizeqc_files = ch_summarizeqc_files.mix(REFBASED.out.depth_tsv.collect{it[1]}.ifEmpty([]))
        ch_summarizeqc_files = ch_summarizeqc_files.mix(REFBASED.out.tsv_vars.collect{it[1]}.ifEmpty([]))
        ch_summarizeqc_files = ch_summarizeqc_files.mix(REFBASED.out.ivar_tsv.collect{it[1]}.ifEmpty([]))
    }
    if ( params.workflow == 'denovo' || params.workflow == 'full' ) {
        ch_summarizeqc_files = ch_summarizeqc_files.mix(DENOVO.out.flagstat.collect{it[1]}.ifEmpty([]))
        ch_summarizeqc_files = ch_summarizeqc_files.mix(DENOVO.out.coverage.collect{it[1]}.ifEmpty([]))
        ch_summarizeqc_files = ch_summarizeqc_files.mix(DENOVO.out.graph_recon_log.collect().ifEmpty([]))
        ch_summarizeqc_files = ch_summarizeqc_files.mix(DENOVO.out.mummer_summary.collect{it[1]}.ifEmpty([]))
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
