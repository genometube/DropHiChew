date
python="/research/xieyeming1/software/Miniconda/envs/xieyeming1/bin/python"
bedtools="/research/xieyeming1/software/Miniconda/envs/xieyeming1/bin/bedtools"

proj_id="tanchen_scHiC10x293t_20240505"
cell_num="7000" # number of demultiplexed cells
proj_dir="/research/xieyeming1/proj_2023/${proj_id}" # demultiplexed library path
lib_id_raw=($(ls -d ${proj_dir}/raw_data/*M*))
lib_id=(${lib_id_raw[@]##*/})
bin=1000000
raw_dir="${proj_dir}/hic_pro/${lib_id}/${i}/outdir/hic_results/matrix/${i}/raw/${bin}"
ice_dir="${proj_dir}/hic_pro/${lib_id}/${i}/outdir/hic_results/matrix/${i}/iced/${bin}"

i="$1"
mkdir -p data_tmp label_info_tmp

echo ${i}
${python} iced_input.py -i ${ice_dir}/${i}_${bin}_iced.matrix -b ${raw_dir}/${i}_${bin}_abs.bed -o data_tmp/${i}_data.txt
cat data_tmp/${i}_data.txt|awk -v b="${bin}" '{print $1"\t"$2"\t"$2+b"\t"$5"\n"$3"\t"$4"\t"$4+b"\t"$5}' > label_info_tmp/${i}.bed
early_ratio=($(${bedtools} intersect -a label_info_tmp/${i}.bed -b E_L_293_replicate_1_${bin}.bedgraph -wao|awk '{print $4 * $8}'|awk '{ total += $1 } END { print total/NR }' ))

echo "${early_ratio}"|tr [:blank:] \\t > label_info_tmp/${i}.early_ratio

date

