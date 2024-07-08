fq_dir="/research/xieyeming1/proj_2023/hichew_paper_20230710/repli_seq/hek293/raw_data"
batchfile=($(find ${fq_dir} -name "*.fastq.gz" -type f))
sample_id_tmp=(${batchfile[@]##*/})
sample_id=(${sample_id_tmp[@]%.fastq.gz})
echo ${sample_id[@]}

for i in ${sample_id[@]}; do
#  if test ! -d "${i}";then
  echo ${i}
  sh repli_seq_4dn.sh ${i}
#  fi
done
