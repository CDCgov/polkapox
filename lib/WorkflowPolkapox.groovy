//
// This file holds several functions specific to the workflow/mpxvassembly.nf in the mpxvassembly pipeline
//

class WorkflowPolkapox {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log, projectDir) {
        //ensure genome inpit fasta is defined and if not return error
        if (!params.fasta) {
            log.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
            System.exit(1)
        }

        // Set clade-specific reference defaults (command-line params override these)
        if (params.clade == "cladeI") {
            params.fasta     = params.fasta     ?: "${projectDir}/assets/NC_003310_clade_I.fasta"
            params.fai       = params.fai       ?: "${projectDir}/assets/NC_003310_clade_I.fasta.fai"
            params.gff       = params.gff       ?: "${projectDir}/assets/NC_003310.gff3"
            params.kraken_db = params.kraken_db ?: "s3://io-pe1-prod-ncezid-oamd-nextstrain/polkapox/orthopox_t2t_kdb/"
        } else if (params.clade == "cladeII") {
            params.fasta     = params.fasta     ?: "${projectDir}/assets/MPXV-UK_P2.noN.fasta"
            params.fai       = params.fai       ?: "${projectDir}/assets/MPXV-UK_P2.noN.fasta.fai"
            params.gff       = params.gff       ?: "${projectDir}/assets/UK-P2.noN.gff"
            // can also use older kraken db: s3://io-pe1-prod-ncezid-oamd-nextstrain/polkapox/mpx_human_kdb/
            params.kraken_db = params.kraken_db ?: "s3://io-pe1-prod-ncezid-oamd-nextstrain/polkapox/orthopox_t2t_kdb/"
        }
        //need to define a clade and if not return error
        else {
            log.error "Invalid clade specified: ${params.clade}. Must be one of: cladeI, cladeII"
            System.exit(1)
        }
    }

    //
    // Enforce per-task resource maximums
    //
    public static def checkMax(params, obj, type) {
        if (type == 'memory') {
            try {
                if (obj.compareTo(params.max_memory as nextflow.util.MemoryUnit) == 1)
                    return params.max_memory as nextflow.util.MemoryUnit
                else
                    return obj
            } catch (all) {
                println "   ### ERROR ###   Max memory '${params.max_memory}' is not valid! Using default value: $obj"
                return obj
            }
        } else if (type == 'time') {
            try {
                if (obj.compareTo(params.max_time as nextflow.util.Duration) == 1)
                    return params.max_time as nextflow.util.Duration
                else
                    return obj
            } catch (all) {
                println "   ### ERROR ###   Max time '${params.max_time}' is not valid! Using default value: $obj"
                return obj
            }
        } else if (type == 'cpus') {
            try {
                return Math.min( obj, params.max_cpus as int )
            } catch (all) {
                println "   ### ERROR ###   Max cpus '${params.max_cpus}' is not valid! Using default value: $obj"
                return obj
            }
        }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }
}
