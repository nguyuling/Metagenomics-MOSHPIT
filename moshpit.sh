#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/lib.sh"

run_pair() {
    sample_name="$1"
    r1="$2"
    r2="$3"

    # per-sample output directories
    Assembly="$ASSEMBLY/$sample_name"
    Binning="$BINNING/$sample_name"
    Dereplication="$DEREPLICATION/$sample_name"
    Taxonomic="$TAXONOMIC/$sample_name"

    log "=== Sample $sample_name ==="
    for stage in genome_assembly mag_binning mag_dereplication taxonomic_classification; do
        "$stage" || { log "Sample $sample_name stopped at $stage"; return 1; }
    done
}

list_pairs() {                              # to detect name, r1 and r2 of samples
    for r1 in "$READS"/*_R1*.fastq.gz; do
        [ -e "$r1" ] || continue
        r2="${r1/_R1/_R2}"
        [ -e "$r2" ] || { echo "WARNING: no mate for $r1" >&2; continue; }
        sample_name=$(basename "$r1" | sed -E 's/_R1.*//')
        printf '%s\t%s\t%s\n' "$sample_name" "$r1" "$r2"
    done
}

run_all_pairs() {                           # to loop run_pair() on all read pairs
    move_reads
    pairs=$(list_pairs)
    if [ -z "$pairs" ]; then
        log "ERROR: no paired reads (*_R1* and *_R2*) found in $READS."
        return 1
    fi

    failed=""
    while IFS=$'\t' read -r sample_name r1 r2 <&3; do
        run_pair "$sample_name" "$r1" "$r2" || failed="$failed $sample_name"
    done 3<<< "$pairs"

    if [ -n "$failed" ]; then
        log "Failed samples:$failed"
        return 1
    fi
}

move_reads() {                          # to move all raw reads from outdir to subdir
    mkdir -p "$READS"
    if ls "$OUT_DIR"/*.fastq* >/dev/null 2>&1; then  # if read files exist in outdir
        mv "$OUT_DIR"/*.fastq* "$READS" || return 1  # move them to reads folder
        log "Moved reads to $READS"
    else
        log "No .fastq* files in $OUT_DIR, skipping move."   # to handle subsequent runs
    fi
}

genome_assembly() {
    need "$r1" || return 1
    need "$r2" || return 1
    need "$quast_db" || return 1
    skip_if_done "genome_assembly" "$Assembly/contigs.qza" "$Assembly/contigs-qc-quast.qzv" && return 0
    mkdir -p "$Assembly"

    manifest="$Assembly/manifest.tsv"
    printf 'sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n%s\t%s\t%s\n' \
        "$sample_name" "$r1" "$r2" > "$manifest"
    
    # import reads
    run qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "$manifest" \
        --input-format PairedEndFastqManifestPhred33V2 \
        --output-path "$Assembly/reads.qza"
    # genome assembly
    run mosh assembly assemble-megahit \
        --i-reads "$Assembly/reads.qza" \
        --p-presets meta-sensitive \
        --p-num-cpu-threads "$threads" \
        --p-min-contig-len 500 \
        --o-contigs "$Assembly/contigs.qza" \
        --parallel-config "$parallelization" \
        --verbose || return 1
    # comprehensive contigs quality assessment
    run mosh assembly evaluate-quast \
        --i-contigs "$Assembly/contigs.qza" \
        --i-references "$quast_db" \
        --p-threads "$threads" \
        --p-min-contig 500 \
        --o-reference-genomes "$Assembly/ref-genomes.qza" \
        --o-results-table "$Assembly/quast-results.qza" \
        --o-visualization "$Assembly/contigs-qc-quast.qzv" \
        --verbose || return 1
}

mag_binning() {
    need "$Assembly/contigs.qza" genome_assembly || return 1
    need "$busco_db" || return 1
    skip_if_done "mag_binning" "$Binning/mags-filtered.qza" && return 0
    mkdir -p "$Binning"
    # contig indexing
    run mosh assembly index-contigs \
        --i-contigs "$Assembly/contigs.qza" \
        --p-threads "$threads" \
        --p-seed 100 \
        --o-index "$Binning/contigs-index.qza" \
        --parallel-config "$parallelization" \
        --verbose || return 1
    # read mapping
    run mosh assembly map-reads \
        --i-index "$Binning/contigs-index.qza" \
        --i-reads "$Assembly/reads.qza" \
        --p-threads "$threads" \
        --p-seed 100 \
        --o-alignment-maps "$Binning/reads-to-contigs-aln.qza" \
        --parallel-config "$parallelization" \
        --verbose || return 1
    # contig binning using MetaBAT 2
    run mosh mag bin-contigs-metabat \
        --i-contigs "$Assembly/contigs.qza" \
        --i-alignment-maps "$Binning/reads-to-contigs-aln.qza" \
        --p-num-threads "$threads" \
        --p-seed 100 \
        --o-mags "$Binning/mags.qza" \
        --o-contig-map "$Binning/contig-map.qza" \
        --o-unbinned-contigs "$Binning/unbinned-contigs.qza" \
        --verbose || return 1
    # bins quality assessment
    run mosh mag evaluate-busco \
        --i-mags "$Binning/mags.qza" \
        --i-db "$busco_db" \
        --i-unbinned-contigs "$Binning/unbinned-contigs.qza" \
        --p-lineage-dataset bacteria_odb12 \
        --p-cpu 2 \
        --o-results "$Binning/busco-results.qza" \
        --o-visualization "$Binning/mags.qzv" \
        --parallel-config "$parallelization" \
        --verbose || return 1
    # filter out low quality MAGs
    run mosh mag filter-mags \
        --i-mags "$Binning/mags.qza" \
        --m-metadata-file "$Binning/busco-results.qza" \
        --p-where "completeness>50 AND contamination<10" \
        --p-on "mag" \
        --o-filtered-mags "$Binning/mags-filtered.qza" \
        --verbose || return 1
}

mag_dereplication() {
    need "$Binning/mags-filtered.qza" mag_binning || return 1
    skip_if_done "mag_dereplication" "$Dereplication/mags-derep.qza" && return 0
    mkdir -p "$Dereplication" 
    # compute hash sketches of each MAG
    run qiime sourmash compute \
        --i-sequence-file "$Binning/mags-filtered.qza" \
        --p-ksizes 105 \
        --p-scaled 100 \
        --o-min-hash-signature "$Dereplication/min-hash.qza" \
        --verbose || return 1
    # compare the sketches to generate matrix of pairwise distance;
    run qiime sourmash compare \
        --i-min-hash-signature "$Dereplication/min-hash.qza" \
        --p-ksize 105 \
        --o-compare-output "$Dereplication/min-hash-compare.qza" \
        --verbose || return 1
    # dereplicate MAGs based on similarity and completeness
    run mosh mag dereplicate-mags \
        --i-mags "$Binning/mags-filtered.qza" \
        --i-distance-matrix "$Dereplication/min-hash-compare.qza" \
        --m-metadata-file "$Binning/busco-results.qza" \
        --p-metadata-column completeness \
        --p-threshold 0.9 \
        --o-dereplicated-mags "$Dereplication/mags-derep.qza" \
        --o-table "$Dereplication/table.qza" \
        --verbose || return 1
}

taxonomic_classification() {
    need "$Dereplication/mags-derep.qza" mag_dereplication || return 1
    need "$kraken_db" || return 1
    skip_if_done "taxonomic_classification" "$Taxonomic/mags-abundances.qza" "$Taxonomic/mags-taxa-barplot.qzv" && return 0
    mkdir -p "$Taxonomic"
    # classify the dereplicated MAGs
    run mosh annotate classify-kraken2 \
        --i-seqs "$Dereplication/mags-derep.qza" \
        --i-db "$kraken_db" \
        --p-threads "$threads" \
        --p-memory-mapping \
        --o-reports "$Taxonomic/kraken2-reports-mags.qza" \
        --o-outputs "$Taxonomic/kraken2-hits-mags.qza" \
        --parallel-config "$parallelization" \
        --verbose || return 1
    # convert artifacts
    run mosh annotate kraken2-to-mag-features \
        --i-reports "$Taxonomic/kraken2-reports-mags.qza" \
        --i-outputs "$Taxonomic/kraken2-hits-mags.qza" \
        --o-taxonomy "$Taxonomic/mags-taxonomy.qza" \
        --verbose || return 1
    # mag indexing
    run mosh assembly index-derep-mags \
        --i-mags "$Dereplication/mags-derep.qza" \
        --p-threads "$threads" \
        --p-seed 100 \
        --o-index "$Taxonomic/mags-derep-index.qza" \
        --verbose || return 1
    # read mapping
    run mosh assembly map-reads \
        --i-index "$Taxonomic/mags-derep-index.qza" \
        --i-reads "$Assembly/reads.qza" \
        --p-threads "$threads" \
        --p-seed 100 \
        --o-alignment-maps "$Taxonomic/reads-to-mags-aln.qza" \
        --parallel-config "$parallelization" \
        --verbose || return 1
    # mag length estimation
    run mosh mag get-feature-lengths \
        --i-features "$Dereplication/mags-derep.qza" \
        --o-lengths "$Taxonomic/mags-derep-lengths.qza" \
        --verbose || return 1
    # abundance estimation
    run mosh mag estimate-abundance \
        --i-alignment-maps "$Taxonomic/reads-to-mags-aln.qza" \
        --i-feature-lengths "$Taxonomic/mags-derep-lengths.qza" \
        --p-metric tpm \
        --p-min-mapq 42 \
        --p-threads "$threads" \
        --o-abundances "$Taxonomic/mags-abundances.qza" \
        --verbose || return 1
    # taxonomic composition visualization
    run mosh taxa barplot \
        --i-table "$Taxonomic/mags-abundances.qza" \
        --i-taxonomy "$Taxonomic/mags-taxonomy.qza" \
        --o-visualization "$Taxonomic/mags-taxa-barplot.qzv" \
        --verbose || return 1
}

case "${1:-}" in
    all)      force=0; run_all_pairs ;;
    allforce) force=1; run_all_pairs ;;
    *)        echo "Usage: bash $0 [all|allforce]"; exit 1 ;;
esac