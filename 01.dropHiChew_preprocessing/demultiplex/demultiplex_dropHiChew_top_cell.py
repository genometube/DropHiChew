import argparse
from operator import itemgetter
import os
import gzip

parser = argparse.ArgumentParser(description='get top cells ranked by read numbers')
parser.add_argument('-R1', '--R1', required=True, help = 'raw MGI R1 fastq.gz')
parser.add_argument('-R2', '--R2', required=True, help = 'raw MGI R2 fastq.gz')
parser.add_argument('-b', '--barcode_10x', required=True, help = 'reverse complementary of 737K-cratac-v1 10x barcode')
parser.add_argument('-d', '--out_dir', required=True, help = 'specify output path for cell list')
parser.add_argument('-s', '--sampling_read_number', required=True, help = 'number of reads to sample from raw fastq')
parser.add_argument('-t', '--top_cell_threshold', required=True, help = 'number of cells to demultiplex, estimated by wet-lab experiment')

args = parser.parse_args()

def rev_comp(DNAstr):
    cis = 'ATCGNatcgn'
    trans = 'TAGCNtagcn'
    rc = str.maketrans(cis,trans)
    rc_DNAstr = DNAstr.translate(rc)[::-1]
    return rc_DNAstr

barcode_10x_set=set()
with open(args.barcode_10x, 'r') as f:
    line = f.read().strip().split('\n')
    for i in range(len(line)):
        barcode_10x_set.add(line[i])

read_id_w_tag_set= set()

barcode_counter_dict = {}
line_count = 0
demultiplex_count = {}
demultiplex_count['scHiC10x'] = 0

with gzip.open(args.R1,'rt') as f1,gzip.open(args.R2,'rt') as f2:
    for line1,line2 in zip(f1,f2):
        if line_count >= int(args.sampling_read_number) * 4:
            break
        if line_count % 4 == 0:
            read_id = line1.strip()[1:29]
        if line_count % 4 == 1:
            seq = line1.strip()+line2.strip()
            barcode_10x = seq[200:216]
        if line_count % 4 == 3:
            qual = line1.strip()+line2.strip()
            R1_insert_qual, R2_insert_qual, R1_insert, R2_insert = qual[:100], qual[100:200], seq[0:100], seq[100:200]
            if barcode_10x in barcode_10x_set:
                if barcode_10x in barcode_counter_dict:
                    barcode_counter_dict[barcode_10x] += 1
                else:
                    barcode_counter_dict[barcode_10x] = 1

                demultiplex_count['scHiC10x'] += 1
        line_count += 1

# Sort the dictionary by values (counts) in descending order and retrieve the top XXX cells
top_sequences = sorted(barcode_counter_dict.items(), key=lambda item: item[1], reverse=True)[:int(args.top_cell_threshold)]

# Print the sequences along with their counts
with open(os.path.abspath(args.out_dir)+ '/top_' + args.top_cell_threshold + '_cells.list','w') as f:
    for sequence, count in top_sequences:
        f.write('\t'.join([sequence,str(count)]) + '\n')

total_demultiplex_count = 0
for k in demultiplex_count:
    print("sample_id: "+k)
    print(demultiplex_count[k])
    total_demultiplex_count += demultiplex_count[k]

print("total demultiplex read count")
print(total_demultiplex_count)
print("total read count")
print(int(line_count/4))

