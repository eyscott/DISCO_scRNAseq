from multiprocessing import Pool
import multiprocessing
import Bio
from functools import partial
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import sys, time, argparse, yaml

short_cutoff = 20 
umi_length = 12
write_long_reads = False
description = '''
	A tool for parallel demultiplexing of Illumina RNA-seq reads with random variable front regions.

	Expects two fastq files. Left-side reads (R1) in the format:
	|NNN|CELL_BARCODE|UNIQUE_MOLECULAR_IDENTIFIER|
	Where the initial |NNN| bases are of variable length and sequence.

	Right side reads (R2) are expected to only contain cDNA sequence information.

	This file will output two fastqs which have the cell barcode and UMI information
	embedded in the fastq read_id lines in both R1 and R2. 

	example R1.fastq output:
	@@M06283:12:000000000-C9WHH:1:1102:11937:1107_CB:1_UMI:TGTTATGGT_SCORE:120.0_FULL_R1:GAGCTGTGCGAGTGTTATGGT
	GAGCTGTGCGAGTGTTATGGT
	+
	CCCCCGGGD+CFFF@FGF@EF,EFFE

	The read id in R2.fastq will be identical. 
	'''
parser = argparse.ArgumentParser(description=description, 
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--fq1', dest='fq1_loc', required=True, help='The location of the R1 fastq')
parser.add_argument('--fq2', dest='fq2_loc', required=True, help='The location of the R2 fastq')
parser.add_argument('--barcodes', dest='barcode_loc', required=True,
	help='A dictionary of the barcodes used')
parser.add_argument('--umi_length', dest='umi_length', 
	help='The length of the unique molecular identifier sequence used.', default=umi_length)
parser.add_argument('--short_read_cutoff', dest='short_cutoff', 
	help='The minimum length R1 reads to be included. All reads under this cutoff are removed',
	default=short_cutoff)
parser.add_argument('--write_long_reads', dest='write_long_reads', type=bool,
	help='Whether or not to output the intermediate fastqs with short reads removed.',
	default=write_long_reads)
args = parser.parse_args()
if not len(sys.argv) > 2:
	print(description)
	print('\n\nAdd \'-h\' for more help.\n')
	sys.exit()
fq1_loc = args.fq1_loc
fq2_loc = args.fq2_loc
barcode_loc = args.barcode_loc
cell_barcodes = yaml.load(open(barcode_loc))
barcode_length = len(list(cell_barcodes)[0])
write_long_reads = args.write_long_reads

# these are the optional intermediate fixed fastqs that will have all R1 reads with less than 20bp removed
fq1_loc_fixed = '{}_short_removed.fastq'.format(fq1_loc.split('.fastq')[0])
fq2_loc_fixed = '{}_short_removed.fastq'.format(fq2_loc.split('.fastq')[0])

# this will be the final R2 fastq with the new IDs
processed_r1_loc = '{}_demultiplexed.fastq'.format(fq1_loc.split('.fastq')[0])
processed_r2_loc = '{}_demultiplexed.fastq'.format(fq2_loc.split('.fastq')[0])

f1_lines = open(fq1_loc,'r').readlines()
f2_lines = open(fq2_loc,'r').readlines()

print('checking input fastqs for formatting and read id match...')
i = 0
f1_lines_fixed = f1_lines.copy()
f2_lines_fixed = f2_lines.copy()
for i in range(0,len(f1_lines), 4):
	read_id1 = f1_lines[i].split(' ')[0]
	read_id2 = f2_lines[i].split(' ')[0]
	if read_id1 != read_id2: 
		print('fastqs failed consistency check!')            
		print('info: line {} read ids do not match:'.format(i))
		print('R1: {} \t R2 {}'.format(read_id1, read_id2))
		print('exiting...')
		sys.exit()

print('removing short reads (<{}bp)...'.format(args.short_cutoff))
for i in range(0,len(f1_lines), 4):
	read_sequence = f1_lines[i+1]
	if len(read_sequence) < args.short_cutoff:        
		f1_lines_fixed[i:i+4] = [None]*4
		f2_lines_fixed[i:i+4] = [None]*4

f1_lines_fixed = [x for x in f1_lines_fixed if x != None]
f2_lines_fixed = [x for x in f2_lines_fixed if x != None]	

print('checking remaining reads for consistency...')
for i in range(0,len(f1_lines_fixed), 4):
	if f1_lines_fixed[i] == None: 
		print(i)
	read_id1 = f1_lines_fixed[i].split(' ')[0]
	read_id2 = f2_lines_fixed[i].split(' ')[0]	
	if read_id1 != read_id2: 
		print('fastqs failed consistency check!')            
		print('this is my fault and im sorry')         
		print('exiting...')
		sys.exit()

f1_lines_fixed_final = [x for x in f1_lines_fixed if x != None]
f2_lines_fixed_final = [x for x in f2_lines_fixed if x != None]

if write_long_reads:
	print('writing R1 long read file...')
	with open(fq1_loc_fixed, 'w') as f:
		f.writelines(f1_lines_fixed_final)
	print('writing R2 long read file...')
	with open(fq2_loc_fixed, 'w') as f:
		f.writelines(f2_lines_fixed_final)    

def process_alignment(read_sequence): 
#     we want to find the alignment with the best score, then we want
#     to find the offset of that alignment, clip the sequence, and 
#     then return the cell barcode, UMI, and score
	score = 0
	offset = 0
	for barcode_sequence, barcode_num in cell_barcodes.items():
		alignment = pairwise2.align.localms(read_sequence,
											barcode_sequence,
											10, #<-- 10 points for correct match
											-2, #<-- -2point for each non-match
											-2, #<-- -2 for each gap opened
											-1, #<-- -1 for extending the gap   
											one_alignment_only=True) 
		if len(alignment) == 0:
				# this is a case of alignment failure
			print('alignment failed! sequence:', read_sequence)
			return 'UNKNOWN', 'UNKNOWN', 'FAIL', read_sequence.split('\n')[0]
	
		alignment = alignment[0]
		if alignment[2] > score:
			score = alignment[2]
			offset = alignment[3]
			cell_barcode = cell_barcodes[barcode_sequence]
#     we now have the score, offset of the best alignment, and barcode number
	clipped_read_sequence = read_sequence[offset:].split('\n')[0]
#     now we have removed the variable front region
#     we want to now remove the first (umi_length) from the clipped sequence, and the 
# 	  remainder becomes the UMI 
	umi = clipped_read_sequence[barcode_length:] # this relies on the cell barcode integrity being maintained!!
	return cell_barcode, umi, score, clipped_read_sequence

def read_processor(line_number):
	if line_number % 1000 == 0: print('processing line...', line_number)
	read_id = f1_lines_fixed_final[line_number].split(' ')[0]
	read_sequence = f1_lines_fixed_final[line_number+1]
	cell_barcode, umi, score, clipped_read_sequence = process_alignment(read_sequence)
	return '{}_CB:{}_UMI:{}_SCORE:{}_FULL_R1:{}\n'.format(read_id, cell_barcode, umi, score, clipped_read_sequence)

print('demultiplexing longs reads (long running process)...')

p = Pool(multiprocessing.cpu_count())
result = p.map(read_processor, range(0,len(f2_lines_fixed_final), 4)) #
# now we just have to replace the id lines in the R2 file
print('WRITING FINAL IDENTIFIED FILES...')
f2_lines_final = f2_lines_fixed_final.copy()
f1_lines_final = f1_lines_fixed_final.copy()
f2_lines_final[0] = result[0]
f1_lines_final[0] = result[0]
f1_lines_final[1] = result[0].split('R1:')[1]
f1_lines_final[3] = f1_lines_final[3][-len(f1_lines_final[1]):]

for i in range(4, len(f2_lines_final), 4):
	f2_lines_final[i] = result[int(i/4)]
	f1_lines_final[i] = result[int(i/4)]
	f1_lines_final[i + 1] = result[int(i/4)].split('R1:')[1]
	f1_lines_final[i + 3] = f1_lines_final[i + 3][-len(f1_lines_final[i + 1]):]

with open(processed_r2_loc, 'w') as f:
	f.writelines(f2_lines_final)
with open(processed_r1_loc, 'w') as f:
	f.writelines(f1_lines_final)	
print('SUCCESS! \nExiting...')
time.sleep(5)
sys.exit()
