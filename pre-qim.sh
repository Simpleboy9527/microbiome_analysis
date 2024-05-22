module load miniconda/4.9.2
source activate
conda activate qiime2-2022.2

mkdir rawdata
mkdir clean_data
mkdir sra
#ncbi download

nohup prefetch -O . $(<list.txt) & 


#fastq
nohup awk -F ',' 'NR {system("fastq-dump --split-3 sra/"$1"/"$1".sra* --outdir rawdata")}' list.txt &

#fastp

# clean data
mkdir clean_data

nohup awk -F '\r$' '{system("fastp -i rawdata/"$1"_R1.fastq.gz -I rawdata/"$1"_R2.fastq.gz -o clean_data/"$1"_R1.fastq.gz -O clean_data/"$1"_R2.fastq.gz -f 1")}' *.txt &

nohup awk -F '\r$' '{system("fastp -i rawdata/"$1".fastq  -o clean_data/"$1".fastq -w 8")}' sra/*.txt &


#manifest

单端：

  echo -e sample-id'\t'absolute-filepath > manifest
  for i in clean_data/*fastq; do path=$(readlink -f $i); id=$(basename $i | cut -d '.' -f 1); echo -e $id'\t'$path >>manifest;done


echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > manifest && for f in clean_data/*_1.fastq; do echo -e "$(basename $f _1.fastq)\t$(realpath $f)\t$(realpath ${f/_1/_2})"; done >> manifest




#导入

#双端
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2 

#单端
nohup qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path manifest \
  --output-path single-end-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2 &



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
nohup qiime quality-filter q-score \
  --i-demux *.qza \
  --o-filtered-sequences demux-filtered.qza \
  --o-filter-stats demux-filter-stats.qza &

# 查看过滤后的数据
nohup qiime demux summarize \
  --i-data demux-filtered.qza \
  --o-visualization demux-filtered.qzv &

  # deblur
qiime deblur denoise-16S \
  --i-demultiplexed-seqs demux-filtered.qza \
  --p-left-trim-len 0 \
  --p-trim-length 249 \
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
    --i-reference-sequences /home/xuhuan/qiime-analysis/reference-databases/gg_13_5_otus/rep_set/99_otus.qza\
    --p-perc-identity 0.99 --p-threads 12 \
    --o-clustered-table clustered-table \
    --o-clustered-sequences clustered-sequences \
    --o-unmatched-sequences unmatched-sequences
    
# close reference after
qiime feature-table summarize \
  --i-table clustered-table.qza \
  --o-visualization clustered-table.qzv


qiime feature-table merge 
  --i-tables PRJEB13051/clustered-table.qza PRJEB13619/clustered-table.qza PRJEB13679/clustered-table.qza PRJEB13680/clustered-table.qza PRJEB13895/clustered-table.qza PRJEB14674/clustered-table.qza PRJEB19825/clustered-table.qza PRJEB23009/clustered-table.qza PRJEB6518/clustered-table.qza PRJNA308319/clustered-table.qza PRJNA317429/clustered-table.qza PRJNA418765/clustered-table.qza PRJNA436359/clustered-table.qza
  --o-merged-table merge_table.qza

qiime feature-table filter-samples \
  --i-table clustered-table.qza \
  --p-min-frequency 3000 \
  --o-filtered-table clustered-filtered-table.qza

qiime tools export --input-path sample-frequency-filtered-table.qza --output-path feature-filtered

qiime tools export --input-path merge_table.qza --output-path feature











