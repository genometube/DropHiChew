proj_id="tanchen_scHiC10x293t_20240505"
proj_path="/research/xieyeming1/proj_2023/${proj_id}"

lib_id_raw=($(ls -d ${proj_path}/raw_data/*M*_*_L*))
lib_id=(${lib_id_raw[@]##*/})
lane_id=($(echo ${lib_id[@]}|cut -d "_" -f3,4))
raw_R1_fq_raw="${proj_path}/raw_data/${lib_id}/${lane_id}_read_1.fq.gz" # PATH to raw MGI R1 fastq.gz
raw_R2_fq_raw="${proj_path}/raw_data/${lib_id}/${lane_id}_read_2.fq.gz" # PATH to raw MGI R2 fastq.gz

top_cell_num="7000" # number of cells to demultiplex, estimated by wet-lab experiment
sampling_read_num="10000000" # number of reads to sample from raw fastq
BUFFER_SIZE="100000000" # lines of fastq to store in ram
THREAD="10" # number of thread for multiprocessing
_10x_barcode="737K-cratac-v1.txt.rc" # reverse complementary of 737K-cratac-v1 10x barcode
python='/research/zhangchen/software/anaconda3/bin/python' # path to python3, tested on Python 3.8.6
TOP_CELL_BARCODE_LIST="top_${top_cell_num}_cells.list"

################################################################################################################
echo "get top cells ranked by read numbers"
echo ==========START at `date`==========

${python} demultiplex_dropHiChew_top_cell.py -b ${_10x_barcode} -R1 ${raw_R1_fq_raw} -R2 ${raw_R2_fq_raw} -d ./ -s ${sampling_read_num} -t ${top_cell_num}
echo ==========END at `date`==========

################################################################################################################
echo "demultiplex"|awk '{print "\n"$0}'
echo ==========START at `date`==========

mkdir -p chunk
${python} demultiplex_dropHiChew_multicore_gz.py -b ${_10x_barcode} -d chunk -R1 ${raw_R1_fq_raw} -R2 ${raw_R2_fq_raw} -B ${BUFFER_SIZE} -t ${THREAD} -l ${TOP_CELL_BARCODE_LIST}
echo ==========END at `date`==========

################################################################################################################
#  echo "optional: cp to target dir"|awk '{print "\n"$0}'
#  echo ==========START at `date`==========
#  target_dir="${proj_path}/demultiplex/${lib_id}/top${top_cell_num}"
#
#  mkdir -p ${target_dir}
#  cp chunk/* ${target_dir}/
#  echo ==========END at `date`==========
################################################################################################################