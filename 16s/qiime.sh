# sra to fasta
mkdir rawdata
awk -F ',' 'NR >1 {system("fastq-dump --split-3 sra/"$1"/"$1".sra --outdir rawdata")}' SraRunTable.txt

# clean data
mkdir clean_data
awk -F ',' 'NR >1 {system("fastp -i rawdata/"$1"_1.fastq -I rawdata/"$1"_2.fastq -o clean_data/"$1"_1.fastq -O clean_data/"$1"_2.fastq")}' SraRunTable.txt

# 生成manifest文件, 双端数据
ls clean_data/*_1* | awk -F '/' '{print $2}' | awk -F '_' 'NR ==1 {print "id"} NR>1 {print $1}' > sample.txt
  
# 数据导入qiime2，格式为双端33格式
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

#序列切除引物和质量控制
#切除引物
#双端数据去除引物
qiime cutadapt trim-paired \
        --i-demultiplexed-sequences paired-end-demux.qza \
        --p-cores 8 \
        --p-front-f GTGCCAGCMGCCGCGGTAA \
		--p-front-r GGACTACHVGGGTWTCTAAT \
        --o-trimmed-sequences trimmed-seqs.qza \
        --verbose

# 过滤
qiime quality-filter q-score \
  --i-demux trimmed-seqs.qza \
  --o-filtered-sequences demux-filtered.qza \
  --o-filter-stats demux-filter-stats.qza

# 查看过滤后的数据
qiime demux summarize \
  --i-data demux-filtered.qza \
  --o-visualization demux-filtered.qzv

# deblur
qiime deblur denoise-16S \
  --i-demultiplexed-seqs demux-filtered.qza \
  --p-trim-length 150 \
  --o-representative-sequences rep-seqs-deblur.qza \
  --o-table table-deblur.qza \
  --p-sample-stats \
  --p-jobs-to-start 16 \
  --o-stats deblur-stats.qza
  
qiime deblur visualize-stats \
  --i-deblur-stats deblur-stats.qza \
  --o-visualization deblur-stats.qzv


# 导入主要分析流程
cp table-deblur.qza table.qza
cp rep-seqs-deblur.qza rep-seqs.qza

# close reference before
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv

# close reference
qiime vsearch cluster-features-closed-reference \
    --i-sequences rep-seqs.qza \
    --i-table table.qza \
    --i-reference-sequences /home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/data/99_otus.qza \
    --p-perc-identity 0.97 --p-threads 12 \
    --o-clustered-table clustered-table \
    --o-clustered-sequences clustered-sequences \
    --o-unmatched-sequences unmatched-sequences
    
# close reference after
qiime feature-table summarize \
  --i-table clustered-table.qza \
  --o-visualization clustered-table.qzv

# 导出丰度表
qiime tools export \
  --input-path table.qza \
  --output-path feature
