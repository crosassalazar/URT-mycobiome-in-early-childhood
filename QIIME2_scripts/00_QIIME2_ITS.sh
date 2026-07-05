#########################################
#              ITS ANALYSIS             #
#########################################

# Data location:
# /data/das_lab/raw_data/vanderbilt/7540

# Provided by BAS

#############################################################
# -------------------- REFERENCE DATABASE ------------------#
#############################################################

# Download UNITE database

wget https://files.plutof.ut.ee/public/orig/C5/54/C5547B97AAA979E45F79DC4C8C4B12113389343D7588716B5AD330F8BDB300C9.tgz

tar zxvf C5547B97AAA979E45F79DC4C8C4B12113389343D7588716B5AD330F8BDB300C9.tgz

# UNITE release:
# sh_refs_qiime_ver8_dynamic_10.05.2021.fasta

#############################################################
# -------------------- QIIME2 INSTALL ----------------------#
#############################################################

wget https://data.qiime2.org/distro/core/qiime2-2022.2-py38-osx-conda.yml

conda env create \
  -n qiime2-2022 \
  --file qiime2-2022.2-py38-osx-conda.yml

rm qiime2-2022.2-py38-osx-conda.yml

conda activate qiime2-2022

#############################################################
# -------------------- SAMPLE SELECTION --------------------#
#############################################################

# Manifest file:
# 6224_manifest_500cutoff.csv

# Prior to QIIME2 import:
# - samples with <500 reads removed
# - negative controls removed
# - positive controls removed

#############################################################
# -------------------- IMPORT FASTQ FILES ------------------#
#############################################################

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path 6224_manifest_500cutoff.csv \
  --output-path input_seqs.qza \
  --input-format PairedEndFastqManifestPhred33

qiime demux summarize \
  --i-data input_seqs.qza \
  --o-visualization seqs_summary.qzv

qiime tools view seqs_summary.qzv

#############################################################
# -------------------- DADA2 DENOISING ---------------------#
#############################################################

# Paired-end

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs input_seqs.qza \
  --p-trunc-len-f 239 \
  --p-trunc-len-r 239 \
  --p-n-threads 8 \
  --o-table paired-table.qza \
  --o-representative-sequences paired-rep_seqs.qza \
  --o-denoising-stats paired-denoise_stats.qza

# Single-end

qiime dada2 denoise-single \
  --i-demultiplexed-seqs input_seqs.qza \
  --p-trunc-len 249 \
  --p-trim-left 0 \
  --p-n-threads 8 \
  --o-table single-table.qza \
  --o-representative-sequences single-rep_seqs.qza \
  --o-denoising-stats single-denoise_stats.qza

#############################################################
# -------------------- INSPECT OUTPUTS ---------------------#
#############################################################

qiime feature-table summarize \
  --i-table single-table.qza \
  --o-visualization single-table-summary.qzv

qiime tools view single-table-summary.qzv

qiime feature-table tabulate-seqs \
  --i-data single-rep_seqs.qza \
  --o-visualization single-rep-seqs-summary.qzv

qiime tools view single-rep-seqs-summary.qzv

qiime metadata tabulate \
  --m-input-file single-denoise_stats.qza \
  --o-visualization single-denoise-stats.qzv

qiime tools view single-denoise-stats.qzv

#############################################################
# -------------------- IMPORT UNITE ------------------------#
#############################################################

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path sh_refs_qiime_ver8_dynamic_10.05.2021.fasta \
  --output-path unite.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path sh_taxonomy_qiime_ver8_dynamic_10.05.2021.txt \
  --output-path unite-taxonomy.qza

#############################################################
# -------------------- TAXONOMIC ASSIGNMENT ----------------#
#############################################################

qiime feature-classifier classify-consensus-vsearch \
  --i-query single-rep_seqs.qza \
  --i-reference-reads unite.qza \
  --i-reference-taxonomy unite-taxonomy.qza \
  --p-threads 8 \
  --o-classification single-taxonomy.qza

qiime metadata tabulate \
  --m-input-file single-taxonomy.qza \
  --o-visualization single-taxonomy.qzv

qiime tools view single-taxonomy.qzv

#############################################################
# -------------------- TAXONOMY BARPLOTS -------------------#
#############################################################

qiime taxa barplot \
  --i-table single-table.qza \
  --i-taxonomy single-taxonomy.qza \
  --m-metadata-file metadata_qualitymetafilter.txt \
  --o-visualization taxa-barplot.qzv

qiime tools view taxa-barplot.qzv

#############################################################
# -------------------- EXPORT FOR R ------------------------#
#############################################################

qiime tools export \
  --input-path single-table.qza \
  --output-path exported_single

qiime tools export \
  --input-path single-taxonomy.qza \
  --output-path exported_single

# Edit taxonomy.tsv header:
# OTUID    taxonomy    confidence

biom add-metadata \
  -i exported_single/feature-table.biom \
  -o exported_single/table-with-taxonomy.biom \
  --observation-metadata-fp exported_single/taxonomy.tsv \
  --sc-separated taxonomy
