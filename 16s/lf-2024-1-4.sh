module load miniconda/4.9.2
source activate
conda activate qiime2-2022.2

#fastp
awk -F'\t' 'NR > 1 { system("fastp -i ./rawdata/00.CleanData/" $1 "/" $1 ".effective.fastq.gz -o ./cleandata/" $1 ".fastq") }' metadata.txt

#manifest
echo -e sample-id'\t'absolute-filepath > manifest
for i in cleandata/*fastq; do path=$(readlink -f $i); id=$(basename $i | cut -d '.' -f 1); echo -e $id'\t'$path >>manifest;done

#import data 
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path manifest \
  --output-path single-end-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2 


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
  --p-max-depth 150000 \
  --o-visualization alpha-rarefaction-150000.qzv

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 38275 \
  --m-metadata-file metadata.txt \
  --output-dir core-metrics-results

#faith_pd指数的显著性检验

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file metadata.txt \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file metadata.txt \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/observed_features_vector.qza \
  --m-metadata-file metadata.txt \
  --o-visualization core-metrics-results/observed-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file metadata.txt \
  --o-visualization core-metrics-results/shannon-group-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.txt \
  --m-metadata-column group\
  --o-visualization core-metrics-results/unweighted-unifrac-group-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.txt \
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
  --m-metadata-file metadata.txt \
  --o-visualization taxa-bar-plots.qzv

qiime feature-table rarefy \
 --i-table table.qza \
 --p-sampling-depth 38275\
 --o-rarefied-table raref-table.qza



 #Linear mixed effect models
qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata.txt \
  --m-metadata-file /home/xuhuan/qiime-analysis/20231225-16s-lf/core-metrics-results/shannon_vector.qza \
  --p-metric shannon_entropy \
  --p-group-columns group \
  --p-random-effects volunteer-id \
  --p-formula shannon_entropy~group+time \
  --p-state-column time \
  --p-individual-id-column volunteer-id \
  --o-visualization linear-mixed-effects.qzv

qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata.txt \
  --m-metadata-file /home/xuhuan/qiime-analysis/20231225-16s-lf/core-metrics-results/shannon_vector.qza \
  --p-metric shannon_entropy \
  --p-group-columns group \
  --p-random-effects volunteer-id \
  --p-formula pielou_evenness~group+time \
  --p-state-column time \
  --p-individual-id-column volunteer-id \
  --o-visualization linear-mixed-effects_shannon.qzv

qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata.txt \
  --m-metadata-file /home/xuhuan/qiime-analysis/20231225-16s-lf/core-metrics-results/faith_pd_vector.qza \
  --p-metric faith_pd \
  --p-group-columns group,time \
  --p-random-effects volunteer-id \
  --p-state-column time \
  --p-individual-id-column volunteer-id \
  --o-visualization linear-mixed-effects_faith_pd.qzv

qiime longitudinal linear-mixed-effects \
  --m-metadata-file metadata.txt \
  --m-metadata-file /home/xuhuan/qiime-analysis/20231225-16s-lf/core-metrics-results/observed_features_vector.qza \
  --p-metric observed_features \
  --p-group-columns group,time \
  --p-random-effects volunteer-id \
  --p-state-column time \
  --p-individual-id-column volunteer-id \
  --o-visualization linear-mixed-effects_observed.qzv



qiime longitudinal volatility \
  --m-metadata-file metadata.txt \
  --m-metadata-file /home/xuhuan/qiime-analysis/20231225-16s-lf/core-metrics-results/shannon_vector.qza \
  --p-default-metric shannon_entropy \
  --p-default-group-column group \
  --p-state-column time \
  --p-individual-id-column volunteer-id \
  --o-visualization volatility.qzv

qiime longitudinal volatility \
  --m-metadata-file metadata.txt \
  --m-metadata-file /home/xuhuan/qiime-analysis/20231225-16s-lf/core-metrics-results/evenness_vector.qza \
  --p-default-metric pielou_evenness \
  --p-default-group-column group \
  --p-state-column time \
  --p-individual-id-column volunteer-id \
  --o-visualization volatility_evenness.qzv

qiime longitudinal volatility \
  --m-metadata-file metadata.txt \
  --m-metadata-file /home/xuhuan/qiime-analysis/20231225-16s-lf/core-metrics-results/faith_pd_vector.qza \
  --p-default-metric faith_pd \
  --p-default-group-column group \
  --p-state-column time \
  --p-individual-id-column volunteer-id \
  --o-visualization volatility_faithpd.qzv

  qiime longitudinal volatility \
  --m-metadata-file metadata.txt \
  --m-metadata-file /home/xuhuan/qiime-analysis/20231225-16s-lf/core-metrics-results/observed_features_vector.qza \
  --p-default-metric observed_features \
  --p-default-group-column group \
  --p-state-column time \
  --p-individual-id-column volunteer-id \
  --o-visualization volatility_observed.qzv



#依据time提取子集
qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file metadata.txt \
  --p-where '[time]="0"' \
  --o-filtered-table table_d0.qza

qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file metadata.txt \
  --p-where '[time]="14"' \
  --o-filtered-table table_d14.qza

qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file metadata.txt \
  --p-where '[time]="21"' \
  --o-filtered-table table_d21.qza

qiime diversity-lib unweighted-unifrac \
  --i-table table_d0.qza \
  --i-phylogeny rooted-tree.qza \
  --o-distance-matrix unweight_distance_d0.qza 

qiime diversity-lib unweighted-unifrac \
  --i-table table_d14.qza \
  --i-phylogeny rooted-tree.qza \
  --o-distance-matrix unweight_distance_d14.qza 

qiime diversity-lib unweighted-unifrac \
  --i-table table_d21.qza \
  --i-phylogeny rooted-tree.qza \
  --o-distance-matrix unweight_distance_d21.qza 

qiime diversity beta-group-significance \
  --i-distance-matrix unweight_distance_d0.qza \
  --m-metadata-file metadata.txt \
  --m-metadata-column group \
  --p-pairwise \
  --o-visualization unweight_distance_d0_significance.qzv

qiime diversity adonis \
  --i-distance-matrix unweight_distance_d0.qza \
  --m-metadata-file metadata.txt \
  --p-formula group \
  --o-visualization unweight_distance_d0_adonis.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix unweight_distance_d14.qza \
  --m-metadata-file metadata.txt \
  --m-metadata-column group \
  --p-pairwise \
  --o-visualization unweight_distance_d14_significance.qzv

qiime diversity adonis \
  --i-distance-matrix unweight_distance_d14.qza \
  --m-metadata-file metadata.txt \
  --p-formula group \
  --o-visualization unweight_distance_d14_adonis.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix unweight_distance_d21.qza \
  --m-metadata-file metadata.txt \
  --m-metadata-column group \
  --p-pairwise \
  --o-visualization unweight_distance_d21_significance.qzv

qiime diversity adonis \
  --i-distance-matrix unweight_distance_d21.qza \
  --m-metadata-file metadata.txt \
  --p-formula group \
  --o-visualization unweight_distance_d21_adonis.qzv

#pcoa 图
qiime emperor plot \
  --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file metadata.txt \
  --p-custom-axes time \
  --o-visualization core-metrics-results/unweighted-unifrac-emperor-time.qzv

qiime emperor plot \
  --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
  --m-metadata-file metadata.txt \
  --p-custom-axes time \
  --o-visualization core-metrics-results/bray-curtis-emperor-time.qzv








  #absolute qiime2 analysis

  #import absolute features table
  biom convert -i feature_absolute.txt -o feature_absolute.biom --table-type="OTU table" --to-hdf5  (qiime2-2022.2) 
   
  qiime tools import \
   --input-path feature_absolute.biom \
   --type 'FeatureTable[Frequency]' \
   --input-format BIOMV210Format \
   --output-path feature-table.qza

qiime diversity beta-phylogenetic 
 --i-table feature-table.qza \
 --i-phylogeny rooted-tree.qza\ 
 --p-metric 'weighted_unifrac'\
 --o-distance-matrix weighted_unifrac_distance.qza

qiime diversity beta-phylogenetic \
--i-table feature-table.qza  \
--i-phylogeny rooted-tree.qza \
--p-metric 'unweighted_unifrac' \
--o-distance-matrix unweighted_unifrac_distance.qza

qiime emperor plot --i-pcoa weighted_unifrac_distance.qza --m-metadata-file metadata.txt --o-visualization weighted_unifrac_emperor.qzv

qiime emperor plot --i-pcoa unweighted_unifrac_distance.qza --m-metadata-file metadata.txt --o-visualization unweighted_unifrac_emperor.qzv


#qiime2 混合线性模型
