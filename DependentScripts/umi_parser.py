from tqdm import tqdm
import argparse
from multiprocessing import Pool
import multiprocessing
import pandas as pd
description = '''
A tool for removing redundant UMI & location alignments from a SAM file. 

The output file will appear in the same directory, with the same input name, and \'_unique_only.sam\' appended to its end. 

NOTE: The tool also will remove the initial '@SQ' lines. 
'''
parser = argparse.ArgumentParser(description=description)
parser.add_argument('--sam_file', dest='f_loc', required=True, help='The location of the sam file to be parsed.')
args = parser.parse_args()
f_loc = args.f_loc
ignore_tags = ['@SQ', '@PG', '@CO', '@HD']
i = 0
with open(f_loc, 'r') as f:
    line = f.readline()
    while any(x in line for x in ignore_tags):
        i += 1
        line = f.readline()
header_amt = i
df = pd.read_csv(f_loc, delimiter='\t', header=None, skiprows=range(header_amt))
df.columns = range(15)
def get_umi(line):
    return line.split('UMI:')[1].split('_')[0]
df['unique_id'] = df[0].map(get_umi) + '_' + df[2] + '_' + df[3].map(str)
df_filtered = df.drop_duplicates(subset='unique_id', keep='first')
df_filtered = df_filtered.drop(columns=['unique_id'])
output_file = f_loc.split('.sam')[0] + '_unique_only.sam'
df_filtered.to_csv(output_file, header=False, index=False, sep='\t')
