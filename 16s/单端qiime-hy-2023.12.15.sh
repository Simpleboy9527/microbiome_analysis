module load miniconda/4.9.2
source activate
conda activate qiime2-2022.2

#fastp filter
for i in ./clean_data/*.fastq.gz; do
    fastp -i "$i" -o cleandata_fastp/"$(basename "$i" .effective.fastq.gz)".fastq
done

#manifest
echo -e sample-id'\t'absolute-filepath > manifest
for i in cleandata_fastp/*fastq; do path=$(readlink -f $i); id=$(basename $i | cut -d '.' -f 1); echo -e $id'\t'$path >>manifest;done

#import data 
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path manifest \
  --output-path single-end-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2 

# demux-view
qiime demux summarize \
  --i-data single-end-demux.qza \
  --o-visualization demux.qzv

#deblur
qiime quality-filter q-score \
  --i-demux single-end-demux.qza \
  --o-filtered-sequences demux-filtered.qza \
  --o-filter-stats demux-filter-stats.qza

qiime metadata tabulate \
  --m-input-file demux-filter-stats.qza \
  --o-visualization demux-filter-stats.qzv

qiime demux summarize \
  --i-data demux-filtered.qza \
  --o-visualization demux-filtered.qzv

qiime deblur denoise-16S \
  --i-demultiplexed-seqs demux-filtered.qza \
  --p-trim-length 180 \
  --p-jobs-to-start 16 \
  --p-no-hashed-feature-ids \
  --o-representative-sequences rep-seqs-deblur.qza \
  --o-table table-deblur.qza \
  --p-sample-stats \
  --o-stats deblur-stats.qza

qiime deblur visualize-stats \
  --i-deblur-stats deblur-stats.qza \
  --o-visualization deblur-stats.qzv

mv rep-seqs-deblur.qza rep-seqs.qza

mv table-deblur.qza table.qza

#FeatureTable and FeatureData summaries
qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv

qiime feature-table tabulate-seqs \
 --i-data rep-seqs.qza \
 --o-visualization rep-seqs.qzv

#Generate a tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --p-n-threads 28 \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

#Alpha and beta diversity analysis
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-steps 15 \
  --p-max-depth 50000 \
  --o-visualization alpha-rarefaction-50000.qzv

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 64425 \
  --m-metadata-file metadata_*.txt \
  --output-dir core-metrics-results

#faith_pd指数的显著性检验

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file metadata_*.txt \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file metadata_*.txt \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/observed_features_vector.qza \
  --m-metadata-file metadata_*.txt \
  --o-visualization core-metrics-results/observed-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file metadata_*.txt \
  --o-visualization core-metrics-results/shannon-group-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata_*.txt \
  --m-metadata-column group \
  --o-visualization core-metrics-results/unweighted-unifrac-group-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata_*.txt \
  --m-metadata-column group \
  --o-visualization core-metrics-results/weighted-unifrac-group-significance.qzv \
  --p-pairwise


qiime feature-classifier classify-sklearn \
  --i-classifier /beegfs/songying/2022_QUCHI/gg-13-8-99-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --p-n-jobs 28 \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file metadata_*.txt \
  --o-visualization taxa-bar-plots.qzv


qiime composition add-pseudocount \
  --i-table table.qza \
  --o-composition-table comp-table.qza

#对属进行ancom
qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table table-l7.qza

qiime composition add-pseudocount \
  --i-table table-l7.qza \
  --o-composition-table comp-table-l7.qza

qiime composition ancom \
  --i-table comp-table-l7.qza \
  --m-metadata-file metadata_*.txt \
  --m-metadata-column group \
  --o-visualization l7-ancom-group.qzv

#ancom
#抽平
qiime feature-table rarefy \
 --i-table table.qza \
 --p-sampling-depth 58783\
 --o-rarefied-table raref-table.qza
#去0
qiime composition add-pseudocount \
  --i-table raref-table.qza \
  --o-composition-table compraref-table.qza
#ancom
qiime composition ancom \
  --i-table compraref-table.qza \
  --m-metadata-file metadata_*.txt \
  --m-metadata-column group \
  --o-visualization ancom-group.qzv

qiime feature-table filter-features \
  --i-table filtered-table.qza \
  --p-min-samples 3  \
  --p-min-frequency 20 \
  --o-filtered-table filtered-sample-table.qza

qiime feature-table summarize \
  --i-table filtered-sample-table.qza \
  --o-visualization filtered-sample-table.qzv

qiime feature-table rarefy \
 --i-table table.qza \
 --p-sampling-depth 64425\
 --o-rarefied-table filtered-table.qza
 #去0
qiime composition add-pseudocount \
  --i-table raref-filtered-table.qza \
  --o-composition-table compraref-filtered-table.qza
#ancom
qiime composition ancom \
  --i-table compraref-filtered-table.qza \
  --m-metadata-file metadata_*.txt \
  --m-metadata-column group \
  --o-visualization ancom-filtered-group.qzv




## bugbase
qiime feature-table merge \
    --i-tables raref-filtered-table.qza raref-table.qza \
    --o-merged-table merged-table.qza
#raref-table.qza is fmt .raref-filtered-table.qza is c_pp

qiime feature-table rarefy \
 --i-table merged-table.qza \
 --p-sampling-depth 58783\
 --o-rarefied-table merged-table-58783.qza


#picrust2


