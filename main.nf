#!/usr/bin/env nextflow

/*
 * RadarMeta: reconciled, ML-augmented metagenomics classification
 * https://github.com/Ifeanyi-omah/radarmeta
 */

nextflow.enable.dsl = 2

// ---------- params ----------
params.samplesheet     = null
params.outdir          = './results'
params.from_contigs    = false

// reference databases
params.kraken2_db      = null
params.bracken_readlen = 150
params.host_grch38     = null              // bowtie2 index prefix
params.rdrpcatch_db    = null              // dir with hmm_dbs/ and mmseqs_dbs/
params.diamond_nr      = null              // diamond_nr.dmnd
params.blast_nt        = null              // path to nt BLAST DB
params.esm2_embeddings = null              // .npz of {ids, vectors, labels} for kNN
params.pathogens_list  = "$projectDir/assets/pathogens_of_concern.tsv"

// behaviour
params.assembler           = 'megahit'     // 'megahit' or 'metaspades'
params.dehost_verify_thresh= 0.001         // fraction host reads allowed after depletion
params.min_contig_len      = 200
params.diamond_evalue      = '1e-5'
params.knn_neighbours      = 5
params.knn_min_cos_sim     = 0.65          // τ for T3 novelty acceptance

// ---------- imports ----------
include { VALIDATE_SAMPLESHEET } from './modules/local/validate'
include { FASTP                } from './modules/local/qc'
include { MULTIQC              } from './modules/local/qc'
include { DEHOST_HUMAN         } from './modules/local/dehost'
include { HOST_TAG_VERTEBRATE  } from './modules/local/dehost'
include { ASSEMBLE_MEGAHIT     } from './modules/local/assemble'
include { MAP_BACK             } from './modules/local/assemble'
include { KRAKEN2_BRACKEN      } from './modules/local/kraken'
include { RDRPCATCH            } from './modules/local/rdrpcatch'
include { DIAMOND_BLASTX       } from './modules/local/diamond'
include { BLAST_NT             } from './modules/local/blast'
include { ML_NOVELTY           } from './modules/local/ml_novelty'
include { RECONCILE            } from './modules/local/reconcile'
include { REPORT               } from './modules/local/report'

// ---------- workflow ----------
workflow {

    if (!params.samplesheet) {
        error "Provide --samplesheet"
    }

    ch_samples = VALIDATE_SAMPLESHEET( file(params.samplesheet) )
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id          : row.sample_id,
                host_type   : row.host_type,
                sample_class: row.sample_class
            ]
            if (params.from_contigs) {
                tuple(meta, file(row.contigs))
            } else {
                tuple(meta, [file(row.fastq_1), file(row.fastq_2)])
            }
        }

    // ----- branch by entry mode -----
    if (params.from_contigs) {
        ch_contigs = ch_samples
    } else {
        // QC
        FASTP( ch_samples )
        ch_qc = FASTP.out.reads

        // host-removal router
        ch_qc.branch {
            human       : it[0].host_type == 'human'
            vertebrate  : it[0].host_type == 'vertebrate'
            invertebrate: it[0].host_type == 'invertebrate'
            wastewater  : it[0].host_type == 'none'
        }.set { ch_routed }

        DEHOST_HUMAN(        ch_routed.human )
        HOST_TAG_VERTEBRATE( ch_routed.vertebrate )

        ch_for_assembly = DEHOST_HUMAN.out.reads
            .mix( HOST_TAG_VERTEBRATE.out.reads )
            .mix( ch_routed.invertebrate )
            .mix( ch_routed.wastewater )

        ASSEMBLE_MEGAHIT( ch_for_assembly )
        MAP_BACK( ASSEMBLE_MEGAHIT.out.contigs.join( ch_for_assembly ) )
        ch_contigs = ASSEMBLE_MEGAHIT.out.contigs

        // read-level classification (only when starting from FASTQ)
        KRAKEN2_BRACKEN( ch_for_assembly )
        ch_kraken = KRAKEN2_BRACKEN.out.report
    }

    // ----- contig-level classification (parallel) -----
    RDRPCATCH(      ch_contigs )
    DIAMOND_BLASTX( ch_contigs )
    BLAST_NT(       ch_contigs )

    // ML novelty: only run on contigs flagged by RdRpCATCH
    ML_NOVELTY( RDRPCATCH.out.hits.join( ch_contigs ) )

    // ----- reconcile -----
    ch_reconcile_in = ch_contigs
        .join( RDRPCATCH.out.hits )
        .join( DIAMOND_BLASTX.out.hits )
        .join( BLAST_NT.out.hits )
        .join( ML_NOVELTY.out.calls )

    RECONCILE( ch_reconcile_in, file(params.pathogens_list) )

    // ----- report -----
    REPORT( RECONCILE.out.calls.collect() )

    // ----- multiqc on QC outputs only when starting from FASTQ -----
    if (!params.from_contigs) {
        MULTIQC( FASTP.out.json.collect() )
    }
}

workflow.onComplete {
    log.info "RadarMeta finished: ${workflow.success ? 'OK' : 'FAILED'}"
    log.info "Results in: ${params.outdir}"
}
