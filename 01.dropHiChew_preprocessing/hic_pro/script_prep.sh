proj_id="tanchen_scHiC10x293t_20240505"
cell_num="7000" # number of demultiplexed cells
proj_dir="/research/xieyeming1/proj_2023/${proj_id}" # demultiplexed library path
lib_id_raw=($(ls -d ${proj_dir}/raw_data/*M*))

lib_id=(${lib_id_raw[@]##*/})

for i in ${lib_id[@]}; do
  echo ${i}

  lane_id=($(echo ${lib_id[@]}|cut -d "_" -f3,4))
  echo ${lane_id[@]}

  mkdir -p ${lane_id}_script

  cp template_script/*txt ${lane_id}_script

  echo "proj_dir=${proj_dir} lane_id=${lane_id} lib_id=${i} proj_id=${proj_id} cell_num=${cell_num}"|tr [:blank:] \\t > ${lane_id}_script/batch-generator.sh
  cat template_script/batch-generator.sh >> ${lane_id}_script/batch-generator.sh

  echo "proj_dir=${proj_dir} lane_id=${lane_id} lib_id=${i} proj_id=${proj_id} cell_num=${cell_num}"|tr [:blank:] \\t > ${lane_id}_script/hic-pro-single.sh
  cat template_script/hic-pro-single.sh >> ${lane_id}_script/hic-pro-single.sh
done
