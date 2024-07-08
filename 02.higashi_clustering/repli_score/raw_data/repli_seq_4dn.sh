# https://data.4dnucleome.org/resources/data-analysis/repli-seq-processing-pipeline
cutadapt="/mnt/software/anaconda3/envs/Bioinfo/bin/cutadapt"
bwa="/mnt/software/anaconda3/envs/Bioinfo/bin/bwa"
samtools="/mnt/software/anaconda3/envs/Bioinfo/bin/samtools"
bedtools="/mnt/software/anaconda3/envs/Bioinfo/bin/bedtools"
bedGraphToBigWig='/research/zhangchen/software/bedGraphToBigWig'

sample="$1"
thread=4
bin=1000000
ref="/research/xieyeming1/db/genome/hg19/Sequence/bwa_index_hg19mt/genome.fa"
chrom_size="/research/xieyeming1/db/genome/hg19/Sequence/bwa_index_hg19mt/hg19.chrom.sizes.txt"
raw_fq="/research/xieyeming1/proj_2023/hichew_paper_20230710/repli_seq/hek293/raw_data/${sample}.fastq.gz"

date

mkdir -p ${sample}
cd ${sample}
#  ${cutadapt} -q 0 -O 1 -m 0 -a AGATCGGAAGAGCACACGTCTG ${raw_fq} > ${sample}.clean.fq
#
#  ${bwa} mem ${ref} ${sample}.clean.fq > ${sample}.sam
#  ${samtools} view -@ ${thread} -bS ${sample}.sam > ${sample}.bam
#  ${samtools} sort -@ ${thread} ${sample}.bam -o ${sample}.sorted.bam
#  ${samtools} index ${sample}.sorted.bam
#  ${samtools} markdup -r -s ${sample}.sorted.bam ${sample}.dedup.bam
#
${bedtools} makewindows -g ${chrom_size} -w ${bin} > genome_${bin}.txt
${bedtools} coverage -counts -sorted -a genome_${bin}.txt -b ${sample}.dedup.bam|sort -k1,1 -k2,2n > ${sample}_${bin}.bedgraph

gzip -c ${sample}_${bin}.bedgraph > ${sample}_${bin}.bedgraph.gz

${bedGraphToBigWig} ${sample}_${bin}.bedgraph ${chrom_size} ${sample}_${bin}.bw

date
