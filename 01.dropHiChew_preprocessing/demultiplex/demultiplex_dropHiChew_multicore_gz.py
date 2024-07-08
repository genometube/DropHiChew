import argparse
from operator import itemgetter
import os
import gzip
from multiprocessing import Process, Queue
import time

parser = argparse.ArgumentParser(description='DropHiChew raw fastq demultiplexing. Multi-thread functions are adapted from scHicDemultiplex.py in scHiCExplorer. Must be running on a high I/O speed SSD for optimal performance.')
parser.add_argument('-R1', '--R1', required=True, help = 'PATH to raw MGI R1 fastq.gz')
parser.add_argument('-R2', '--R2', required=True, help = 'PATH to raw MGI R2 fastq.gz')
parser.add_argument('-b', '--barcode_10x', required=True, help = 'reverse complementary of 737K-cratac-v1 10x barcode')
parser.add_argument('-d', '--out_dir', required=True, help = 'specify output path')
parser.add_argument('-B', '--buffer_size', required=True, help = 'lines of fastq to store in ram')
parser.add_argument('-t', '--thread', required=True, help = 'number of thread for multiprocessing')
parser.add_argument('-l', '--top_cell_barcode_list', required=True, help = 'cells to demultiplex')

args = parser.parse_args()

def writeFile(pFileName, pReadArray):
    with gzip.open(pFileName, 'at') as f:
        for read in pReadArray:
            f.write(read)

def writeSampleFiles(pFileNameList, pBuffer, pQueue):
    for i, cells in enumerate(pFileNameList):
        for j, cell_ in enumerate(cells):
            writeFile(cell_, pBuffer[i][j])
    pQueue.put('Done')

def handleCompressingMulticore(pFileNameList, pBuffer, pThreads):
    filesPerThread = len(pFileNameList) // pThreads
    if filesPerThread == 0:
        writeSampleFiles(pFileNameList, pBuffer, None)
    else:
        queue = [None] * pThreads
        process = [None] * pThreads
        file_list_sample = [None] * pThreads
        buffer_sample = [None] * pThreads

        all_data_collected = False

        for i in range(pThreads):
            if i < pThreads - 1:
                file_list_sample = pFileNameList[i * filesPerThread:(i + 1) * filesPerThread]
                buffer_sample = pBuffer[i * filesPerThread:(i + 1) * filesPerThread]
            else:
                file_list_sample = pFileNameList[i * filesPerThread:]
                buffer_sample = pBuffer[i * filesPerThread:]

            queue[i] = Queue()
            process[i] = Process(target=writeSampleFiles, kwargs=dict(pFileNameList=file_list_sample,pBuffer=buffer_sample,pQueue=queue[i]))
            process[i].start()

        while not all_data_collected:
            for i in range(pThreads):
                if queue[i] is not None and not queue[i].empty():
                    _ = queue[i].get()
                    process[i].join()
                    process[i].terminate()
                    process[i] = None

            all_data_collected = True

            for i in range(pThreads):
                if process[i] is not None:
                    all_data_collected = False
            time.sleep(1)

def rev_comp(DNAstr):
    cis = 'ATCGNatcgn'
    trans = 'TAGCNtagcn'
    rc = str.maketrans(cis,trans)
    rc_DNAstr = DNAstr.translate(rc)[::-1]
    return rc_DNAstr

demultiplex_count = {}
demultiplex_count['scHiC10x'] = 0

with open(args.top_cell_barcode_list, 'r') as f:
    line = f.read().strip().split('\n')

barcode_set= set()
for i in range(len(line)):
    tmp = line[i].split()
    barcode_set.add(tmp[0])

# read and split fastq
count,line_count,cell_counter = 0,0,0
cell_index = {}

with gzip.open(args.R1,'rt') as f1,gzip.open(args.R2,'rt') as f2:
    lines_out_buffer,file_writer=[],[]
    pBufferSize, pThreads = int(args.buffer_size), int(args.thread)

    for line1,line2 in zip(f1,f2):
        if count % 4 == 0:
            if line1.strip()[1:29] != line2.strip()[1:29]:
                print('readID from R1 not equal to readID from R2')
            read_id = line1.strip()[1:29]
        if count % 4 == 1:
            seq = line1.strip()+line2.strip()
            barcode_10x = seq[200:216]
        if count % 4 == 3:
            qual = line1.strip()+line2.strip()
            R1_insert_qual, R2_insert_qual, R1_insert, R2_insert = qual[:100], qual[100:200], seq[0:100], seq[100:200]
            read_id_w_tag = '@' + read_id + '_' + barcode_10x
            barcode_tag = barcode_10x
            if barcode_10x in barcode_set:
                if barcode_tag in cell_index:
                    lines_out_buffer[cell_index[barcode_tag]][0].extend([read_id_w_tag+ '\n', R1_insert+ '\n', '+\n', R1_insert_qual+ '\n'])
                    lines_out_buffer[cell_index[barcode_tag]][1].extend([read_id_w_tag+ '\n', R2_insert+ '\n', '+\n', R2_insert_qual+ '\n'])
                else:
                    cell_index[barcode_tag] = cell_counter
                    file_writer.append([os.path.abspath(args.out_dir) + '/' + barcode_tag + '_R1.fq.gz',
                                        os.path.abspath(args.out_dir) + '/' + barcode_tag + '_R2.fq.gz'])
                    lines_out_buffer.append([[read_id_w_tag+ '\n', R1_insert+ '\n', '+\n', R1_insert_qual+ '\n'],
                                             [read_id_w_tag+ '\n', R2_insert+ '\n', '+\n', R2_insert_qual+ '\n']])
                    cell_counter += 1

                demultiplex_count['scHiC10x'] += 1
        count += 1
        line_count += 4

        if line_count % pBufferSize == 0 or line1 is False:
            buffered_elements = 0
            for lines_buffered in lines_out_buffer:
                buffered_elements += len(lines_buffered[0])
                buffered_elements += len(lines_buffered[1])
            if buffered_elements > pBufferSize or line1 is False:
                handleCompressingMulticore(file_writer, lines_out_buffer, pThreads)
                lines_out_buffer, file_writer = [], []
                line_count = 0
                cell_index = {}
                cell_counter = 0

    if lines_out_buffer is not None and len(lines_out_buffer) > 0:
        handleCompressingMulticore(file_writer, lines_out_buffer, pThreads)
        lines_out_buffer, file_writer = [], []
        cell_counter = 0
        cell_index = {}
        
total_demultiplex_count = 0
for k in demultiplex_count:
    print(k)
    print(demultiplex_count[k])
    total_demultiplex_count += demultiplex_count[k]

print("total demultiplex count")
print(total_demultiplex_count)
print("total read count")
print(int(count/4))
