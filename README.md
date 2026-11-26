# Metagenomics-MOSHPIT
Metagenomics analysis using MOSHPIT (MOdular SHotgun metagenome Pipelines with Integrated provenance Tracking).

### Databases
0. Install MOSHPIT [here](https://moshpit.qiime2.org/en/stable/chapters/howtos/setup/)

1. Quast database
```bash
    mkdir -p <DB_PATH>/quast
    wget -O reference-genomes.qza https://polybox.ethz.ch/index.php/s/SDczx2T6sxnNL9x/download
```

2. Busco database
```bash
qiime mag fetch-busco-db --p-lineages bacteria_odb12 --o-db <DB_PATH>/busco/busco-db-prokaryota.qza --verbose
```

3. Kraken database
```bash
qiime annotate build-kraken-db --p-collection standard8 --o-kraken2-db <DB_PATH>/kraken/kraken2-db.qza --o-bracken-db <DB_PATH>/kraken/bracken-db.qza --verbose
```

### Guide to Run the Pipeline
0. Set output and databases directories path in `lib.sh`, locate all paired-end reads at output directory.

1. Activate `conda env` OR run `docker` 
- Windows / Linux
```bash
    conda activate <ENV_NAME>
``` 
-  MacOS (update output directory path)
```bash
    docker run --rm -it \
        -v /Volumes/T7:/T7 \
        -v "$HOME/Metagenomics-MOSHPIT":/work \
        -w /work \
        quay.io/qiime2/moshpit:2026.7
```

2. Run the pipeline
```bash
    bash moshpit.sh allforce        # run all unfinished steps, OR
    bash moshpit.sh all             # run all steps regardlessly
```