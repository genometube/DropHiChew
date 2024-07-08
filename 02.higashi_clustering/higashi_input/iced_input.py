import argparse
from operator import itemgetter
# cell_name         cell_id  chrom1  pos1     chrom2  pos2       count
# CATGATGTGATCAGCT  0        chr1    3547027  chr1    4325341    1
# CATGATGTGATCAGCT  0        chr1    3755233  chr1    22841329   1
parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--input', required=True, help='')
parser.add_argument('-b', '--bed', required=True, help='')
# parser.add_argument('-c', '--cell_name', required=True, help='')
# parser.add_argument('-n', '--cell_id', required=True, help='')
parser.add_argument('-o', '--output', required=True, help='')
args = parser.parse_args()

if __name__ == '__main__':
    o = open(args.output, 'w')
    # cell_name,cell_id = args.cell_name,args.cell_id
    # cell_name = args.cell_name
    locus_dict = {}
    with open(args.bed, 'r') as f:
        line = f.read().strip().split('\n')
        for i in range(len(line)):
            tmp = line[i].split('\t')
            locus_dict[tmp[3]] = tmp[0] + '\t' + tmp[1]

    with open(args.input, 'r') as f:
        line = f.read().strip().split('\n')

        for i in range(len(line)):
            tmp = line[i].split('\t')
            if tmp[0] in locus_dict and tmp[1] in locus_dict:
                locus_1,locus_2 = locus_dict[tmp[0]],locus_dict[tmp[1]]

            o.write('\t'.join([locus_1,locus_2,tmp[2]])+'\n')
    o.close()


