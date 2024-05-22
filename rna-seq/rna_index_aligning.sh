#!/bin/bash 

#download ref
cd /home/xuhuan/rna-seq/ref

#ref gene
nohup wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz &
#ref gene Annotation 
nohup wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gtf.gz &

gunzip *gz

#build index
cd /home/xuhuan/rna-seq/liping-preservative-RNA-seq-data/

subread-buildindex -o my_index /home/xuhuan/rna-seq/ref/Mus-musculus/*fna

#create a sample id file to alignment

find ./2.cleandata -maxdepth 1 -type d | while read -r folder; do echo "${folder:14}"; done > samplename.txt

for i in `less /home/xuhuan/rna-seq/Pseudomonas\ aeruginosa/thd/sample.txt`; do
    subread-align -t 0 -i my_index \
                  -r ./clean_data/${i}_1.clean.fq.gz \
                  -R ./clean_data/${i}_2.clean.fq.gz \
                  -o ./result_mapping/${i}_result.bam \
                  -T 14

done

# Read summarization
featureCounts -a /home/xuhuan/rna-seq/Pseudomonas_aeruginosa/ref/Pseudomonas_aeruginosa_ref/*gtf -p --countReadPairs \
              -F \
              -t  -g gene_id \
              -o /home/xuhuan/rna-seq/Pseudomonas_aeruginosa/analysis/counts.txt \
              /home/xuhuan/rna-seq/Pseudomonas_aeruginosa/analysis/result/*bam

#-C chimeric fragments(那些两端与不同染色体对齐的片段)将不被计算在内。此选项仅适用于fragments count (read pairs)
#-t exon feature类型为外显子水平
#-p：指定输入数据包含paired-end reads



#featurecount 
GeneID
featureCounts -a /home/xuhuan/rna-seq/Pseudomonas_aeruginosa/ref/Pseudomonas_aeruginosa_ref/my.gtf\
              -p --countReadPairs \
              -t "CDS,exon"\
              -g gene_name \
              -o /home/xuhuan/rna-seq/Pseudomonas_aeruginosa/analysis/counts.txt \
              /home/xuhuan/rna-seq/Pseudomonas_aeruginosa/analysis/result/*bam



