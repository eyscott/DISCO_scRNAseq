module load nixpkgs/16.09
module load python/3.7.0

source /home/eyscott/project/eyscott/eyscott/19.04.04_RG_LIVE/HarrisonEnviron/python_env/bin/activate
for f in $(ls *R2_001.fastq | awk -F'[-_.]' '{print $1 "_" $2 "_" $3}');
do python /home/eyscott/scratch/19.12.30_RedGreen/sequence_fixer.py --fq1 ${f}_R1_001.fastq --fq2 ${f}_R2_001.fastq --barcodes barcodes.yaml;
done 
