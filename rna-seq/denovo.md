### 1.pre-assembly quanlity control and filtering
 awk -F '\r$' '{system("fastp -i rawdata/"$1"_R1_1.fastq.gz -I rawdata/"$1"_R2_1.fastq.gz -o clean_data/"$1"_R1_1.fastq -O clean_data/"$1"_R2_1.fastq -f 1")}' sample.txt &

awk -F '\r$' '{system("fastp -i rawdata/"$1"_R1_2.fastq.gz -I rawdata/"$1"_R2_2.fastq.gz -o clean_data/"$1"_R1_2.fastq -O clean_data/"$1"_R2_2.fastq -f 1")}' sample.txt &  

awk -F '\r$' '{system("fastp -i rawdata/"$1"_R1_3.fastq.gz -I rawdata/"$1"_R2_3.fastq.gz -o clean_data/"$1"_R1_3.fastq -O clean_data/"$1"_R2_3.fastq -f 1")}' sample.txt &

### 2.Indexing a reference genome (mapping)
bowtie2-build /home/xuhuan/rna-seq/Pseudomonas_aeruginosa/ref/Pseudomonas_aeruginosa_ref/Pseudomonas_aeruginosa_genomic.fna Pseudomonas_aeruginosa

nohup for i in `less /home/xuhuan/rna-seq/Pseudomonas_aeruginosa/analysis/sample.txt`; do
	bowtie2 -x Pseudomonas_aeruginosa -1 ./clean_data/${i}_R1_1.fastq -2 ./clean_data/${i}_R2_1.fastq -S ./result/${i}_1.sam
	bowtie2 -x Pseudomonas_aeruginosa -1 ./clean_data/${i}_R1_2.fastq -2 ./clean_data/${i}_R2_2.fastq -S ./result/${i}_2.sam
	bowtie2 -x Pseudomonas_aeruginosa -1 ./clean_data/${i}_R1_3.fastq -2 ./clean_data/${i}_R2_3.fastq -S ./result/${i}_3.sam
done &

#sam > bam
for file in *;do
samtools sort "$file" > "${file%.sam}.bam"
done

#StringTie
for file in /home/xuhuan/rna-seq/Pseudomonas_aeruginosa/analysis/result/*.bam;do
	StringTie $file \
		  -p 8\
		  -o /home/xuhuan/rna-seq/Pseudomonas_aeruginosa/analysis/stringtie_result/"${file%.bam}.gtf"\
		  -G /home/xuhuan/rna-seq/Pseudomonas_aeruginosa/ref/Pseudomonas_aeruginosa_ref/*gff\
		  -A "/home/xuhuan/rna-seq/Pseudomonas_aeruginosa/analysis/stringtie_result/${file%.bam}.tab"\
		  -C "/home/xuhuan/rna-seq/Pseudomonas_aeruginosa/analysis/stringtie_result/${file%.bam}_cov_refs.gtf"
done

for i in `less /home/xuhuan/rna-seq/Pseudomonas_aeruginosa/analysis/stringtie_result/sample_list.txt`; do
    data=$(grep -Ei 'ExoU|ExoS|ExoT|ExoY|LasI|lasR|PCNbla IMP-1|blaVIM-2|blaKPC-2|blaGES|bla OXA-10|gyrB|bla TEM|OprD' "$i".gtf)
    id=$(basename "$i")
    echo -e "$id\t$data" >> gene_result.txt
done


### 2.De novo assembly  （de novo assembly）
#
Trinity --seqType fq \
	--max_memory 20G  \
    --left ./clean_data/*R1.fastq.gz \ 
	--right ./clean_data/*R2.fastq.gz \ 
	--CPU 6  \
	--output /home/xuhuan/rna-seq/Pseudomonas_aeruginosa/analysis/trinity \
	--min_contig_length 200 \
    --no_normalize_reads



