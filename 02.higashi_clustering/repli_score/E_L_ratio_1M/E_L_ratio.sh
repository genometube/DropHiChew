# https://data.4dnucleome.org/biosamples/4DNBSR3CBAW3/#expsets-table
bin=1000000
# early replicate 1 4DNFIIVKU149
early_bg="/research/xieyeming1/proj_2023/hichew_paper_20230710/repli_seq/hek293/clean_align_dedup/4DNFIIVKU149/4DNFIIVKU149_${bin}.bedgraph"

# late replicate 1 4DNFI2J9Z6YE
late_bg="/research/xieyeming1/proj_2023/hichew_paper_20230710/repli_seq/hek293/clean_align_dedup/4DNFI2J9Z6YE/4DNFI2J9Z6YE_${bin}.bedgraph"

paste -d "\t" ${early_bg} ${late_bg} |awk '{if($8 != 0 && $4 != 0){print $1,$2,$3,log($4/$8)/log(2)}}' OFS='\t' > E_L_293_replicate_1_${bin}.bedgraph

echo -e "chr\tstart\tstop\t"`ls *.bedgraph`| sed 's/\ /\t/g' > merge_RT.txt.raw
cat *.bedgraph >> merge_RT.txt.raw

cat merge_RT.txt.raw|grep -v chrM|grep -v chrY > merge_RT.txt

Rscript="/mnt/software/anaconda3/envs/Bioinfo/bin/Rscript"
${Rscript} repli_score.R E_L_293_replicate_1_${bin}.bedgraph

