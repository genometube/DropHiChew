cooler="/research/zhangchen/software/anaconda3/envs/liwenqing/bin/cooler" # v0.8.2
pairix="/research/xieyeming1/software/pairix/pairix-master/bin/pairix" # v0.3.7
bgzip="/research/xieyeming1/software/pairix/pairix-master/bin/bgzip" # v0.3.7
bedtools="/research/xieyeming1/software/Miniconda/envs/xieyeming1/bin/bedtools" # v2.31.1

date
proj_id="tanchen_scHiC10x293t_20240505"

ncores=16
sample="${proj_id}"
ref_fai="/research/xieyeming1/db/hic_pro_ref/mm10/mm10_clean.fa.fai"

proj_id="tanchen_scHiC10x293t_20240505"
cell_num="7000" # number of demultiplexed cells
proj_dir="/research/xieyeming1/proj_2023/${proj_id}" # demultiplexed library path
lib_id_raw=($(ls -d ${proj_dir}/raw_data/*M*))
lib_id=(${lib_id_raw[@]##*/})

fq_dir="${proj_dir}/demultiplex/${lib_id}/top${cell_num}"
batchfile=($(find ${fq_dir} -name "*_R1.fq.gz" -type f))
sample_id_tmp=(${batchfile[@]##*/})
sample_id=(${sample_id_tmp[@]%_R1.fq.gz})

for i in ${sample_id[@]}; do
  # echo ${i}
  cat ${proj_dir}/hic_pro/${lib_id}/${i}/outdir/hic_results/data/${i}/${i}.allValidPairs >> ${sample}.allValidPairs # hic pro single cell non dup valid pairs ".allValidPairs" file
done

mkdir -p tmp
cat ${sample}.rmdvalidPairs|awk '{if ($2 > $5 || $2 == $5 && $3 > $6){print i++"\t"$5"\t"$6"\t"$2"\t"$3"\t"$7"\t"$4}else {print i++"\t"$2"\t"$3"\t"$5"\t"$6"\t"$4"\t"$7}}'|\
sort -T tmp -k2,2 -k4,4 -k3,3n -k5,5n > ${sample}.pairs

${bgzip} ${sample}.pairs
${pairix} -p pairs -f ${sample}.pairs.gz

cat ${ref_fai}|awk '{print $1"\t"$2}'|grep -v chrM > chrom_size.txt
cat /research/xieyeming1/db/hic_pro_ref/hg19/chrom_hg19.sizes|grep -v chrM > chrom_size.txt
pairs_file=${sample}.pairs.gz  # bgzipped, with .px2 index
chrsize_file=chrom_size.txt
bin_size=25000
out_prefix=${sample}
max_split=2 # (default 2)

${cooler} makebins -o $chrsize_file:$bin_size $chrsize_file $bin_size
${cooler} cload pairix -p $ncores -s $max_split $chrsize_file:$bin_size $pairs_file $out_prefix.cool
chunksize=10000000
${cooler} zoomify --balance --balance-args '--convergence-policy store_nan' -n $ncores -o $out_prefix.mcool -c $chunksize $RES_ARG $out_prefix.cool

hicFindTADs="/mnt/software/anaconda3/envs/hicexplorer/bin/hicFindTADs"
${hicFindTADs} -m ${out_prefix}.mcool::resolutions/25000 \
--outPrefix ${out_prefix} \
--minDepth 100000 \
--maxDepth 750000 \
--step 50000 \
--thresholdComparisons 0.05 \
--delta 0.01 \
--correctForMultipleTesting fdr \
-p 20

cat ${out_prefix}_score.bedgraph|awk '$4>0' > tad_region.bed
${bedtools} merge -i tad_region.bed > tad_region_merged.bed

cat tad_region_merged.bed|awk '{print $1"\t"$2"\t"$2+50000}' > tad_region_merged_left_50k.bed
cat tad_region_merged.bed|awk '{print $1"\t"$3-50000"\t"$3}' > tad_region_merged_right_50k.bed

mkdir -p tmp
for i in ${sample_id[@]}; do
  echo ${i}
  cat ${rmdvalidPairs_dir}/${i}.rmdvalidPairs|awk '{print $2"\t"$3-1"\t"$3"\n"$4"\t"$5-1"\t"$5}' > tmp/${i}_1D.bed
  ${bedtools} intersect -a tad_region_merged_left_50k.bed -b tmp/${i}_1D.bed -wa -c |awk -v s="${i}" '{print $0"\t"s"\t"1+i++}' >> ${sample}_tad_left_50k.bed
  ${bedtools} intersect -a tad_region_merged_right_50k.bed -b tmp/${i}_1D.bed -wa -c |awk -v s="${i}" '{print $0"\t"s"\t"1+i++}'>> ${sample}_tad_right_50k.bed
done

#  rm chrom_size* $out_prefix.cool

date