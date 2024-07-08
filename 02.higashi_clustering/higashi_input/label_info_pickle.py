import pickle
import pandas as pd
import numpy as np

import argparse
from operator import itemgetter
parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--input', required=True, help = '')
parser.add_argument('-o', '--output', required=True, help = '')
args = parser.parse_args()

column_names = ['cell_name', 'cell_id', 'cell_type']

df = pd.read_csv(args.input, delimiter='\t', names=column_names, header=None)

output_label_file = open(args.output, "wb")
#df pandas df
label_info = {k:np.asarray(df[k]) for k in df.columns}
pickle.dump(label_info, output_label_file)

