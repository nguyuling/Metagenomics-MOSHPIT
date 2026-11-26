#!/bin/bash
set -e

# SET
OUT_DIR=/T7/metagenomics
DB_DIR=/T7/databases
threads=6

# docker image
IMAGE="quay.io/qiime2/moshpit:2026.7"

# output directories
READS="$OUT_DIR/0_Reads"
ASSEMBLY="$OUT_DIR/1_Genome_Assembly"
BINNING="$OUT_DIR/2_MAG_Binning"
DEREPLICATION="$OUT_DIR/3_MAG_Dereplication"
TAXONOMIC="$OUT_DIR/4_Taxonomic_Classification"
LOG="$OUT_DIR/logs"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# database directories
quast_db="$DB_DIR/quast/reference-genomes.qza"
busco_db="$DB_DIR/busco/busco-db-prokaryota.qza"
kraken_db="$DB_DIR/kraken/kraken2-db.qza"

# parallelization 
parallelization="$SCRIPT_DIR/parallelization/local.config.toml"

mkdir -p "$LOG"
log_file="$LOG/pipeline_$(date '+%Y%m%d_%H%M%S').log"

log() {                                     # log datetime and message to terminal and log file
    msg="[$(date '+%Y-%m-%d %H:%M:%S')] $*"
    echo "$msg"
    echo "$msg" >> "$log_file"
}

run() {                                     # usage: run <command> <parameter>
    log "START: $*"                         # log the starting step's command
    start=$(date +%s)                       # store the start time
    if "$@" >> "$log_file" 2>&1; then       # run the command, write stdout / stderr to log file
        ran_count=$(( ${ran_count:-0} + 1 ))    # add 1 for each successful step run
        log "OK ($(( $(date +%s) - start ))s)"  # exit 0 : log run time (current time - start time)
    else
        rc=$?                               # exit != 0 : hold exit code of the last command
        log "FAILED (exit $rc): $*"         # log the exit code and the failed command
        return "$rc"                        # return the exit code
    fi
}

need() {                                    # usage: need <file or dir> <step number that creates it>
    if [ ! -e "$1" ]; then                  # if the prerequisite file is absent
        log "ERROR: missing $1. Run step $2 first."  && return 1   # log reminder message about running the upstream step first
    fi
}

skip_if_done() {                            # usage: skip_if_done <step name> <output>...
    name=$1; shift                          # remove first arguments (step name), leaving only the output files
    [ "${force:-0}" = 1 ] && { log "RUN $name: force=1"; return 1; }
    for f in "$@"; do                       # for each file
        if [ ! -e "$f" ]; then              # if the file is missing
            log "RUN $name: missing $f"     # log message to run the step
            return 1                        # stop the function
        fi
    done
    log "SKIPPED $name: outputs already exist"  # if file(s) exist, log skipping step message
    return 0
}