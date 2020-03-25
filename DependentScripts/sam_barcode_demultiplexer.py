from tqdm import tqdm
import argparse
description = '''
A quick tool for splitting a sam file into 5 separate sam files based on barcodes. 

The output files will appear in the same directory, with the same input name, but have the barcode integer (1-5) appended to the end of the name.
'''
parser = argparse.ArgumentParser(description=description)
parser.add_argument('--sam_file', dest='f_loc', required=True, help='The location of the sam file to be split into 5 different barcodes.')
args = parser.parse_args()
f_loc = args.f_loc
#f_loc = '/home/hedwar/projects/def-awheeler/shared/G1_demultiplexed_starAligned.sortedByCoord.out.sam'

with open(f_loc, 'r') as f:
    lines = f.readlines()

for i in tqdm(range(1,6)):
    new_lines = []
    print('\nprocessing:\n',i)
    for line in tqdm(lines):
        if '_CB:{}_'.format(i) in line:
            new_lines.append(line)
    print('\nwriting file:\n',i)
    f_name = f_loc.split('.sam')[0] + str(i) + '.sam'
    with open(f_name, 'w') as f:
        f.writelines(new_lines)
