# Metagenomics-MOSHPIT
Metagenomics analysis using MOSHPIT (MOdular SHotgun metagenome Pipelines with Integrated provenance Tracking).

### Databases
1. Quast database
```bash
    cd <DB_PATH>
    wget -O reference-genomes.qza https://polybox.ethz.ch/index.php/s/SDczx2T6sxnNL9x/download
```

2. Busco database
```bash
qiime mag fetch-busco-db --p-lineages bacteria_odb12 --o-db <DB_PATH>/busco-db-prokaryota.qza --verbose
```

3. Kraken database
```bash
qiime annotate build-kraken-db --p-collection standard8 --o-kraken2-db <DB_PATH>/kraken2-db.qza --o-bracken-db <DB_PATH>/bracken-db.qza --verbose
```
