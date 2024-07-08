HIC_pro="/research/xieyeming1/software/hic_pro/hic_pro_3.1.0/HiC-Pro_3.1.0/bin/HiC-Pro"
fq_dir="${proj_dir}/demultiplex/${lib_id}/top${cell_num}"
sample="$1"

echo ==========START at `date`==========

mkdir -p ${proj_dir}/hic_pro/${lib_id}/${sample}/rawdata/${sample}
ln -s ${fq_dir}/${sample}_R1.fq.gz ${proj_dir}/hic_pro/${lib_id}/${sample}/rawdata/${sample}/${sample}_R1.fastq
ln -s ${fq_dir}/${sample}_R2.fq.gz ${proj_dir}/hic_pro/${lib_id}/${sample}/rawdata/${sample}/${sample}_R2.fastq

${HIC_pro} -i ${proj_dir}/hic_pro/${lib_id}/${sample}/rawdata \
-o ${proj_dir}/hic_pro/${lib_id}/${sample}/outdir \
-c ${proj_dir}/hic_pro/${lane_id}_script/config-hicpro.txt

rm -rf outdir/bowtie_results/

echo ==========END at `date`==========