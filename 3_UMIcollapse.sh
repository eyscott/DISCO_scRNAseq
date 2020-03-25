##convert these .bam files into .sam files
for file in ./*.bam
do
    echo $file 
    samtools view -h $file > ${file/.bam/.sam}
done

####run UMI parser (Harrison's python script)####
#!/bin/sh
#SBATCH --job-name=run_umi.sh
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=12
#SBATCH -t 4:00:00
#SBATCH -o /home/eyscott/project/eyscott/eyscott/run_umi.sh.out
#SBATCH -e /home/eyscott/project/eyscott/eyscott/run_umi.sh.err
#SBATCH --mem 100000

module load nixpkgs/16.09
module load python/3.7.0


source /home/eyscott/project/eyscott/eyscott/19.04.04_RG_LIVE/HarrisonEnviron/python_env/bin/activate
for f in $(ls *sortedByCoord.out.sam | awk -F'[-_.]' '{print $1 "_" $2 "_" $3 "_" $4 "_" $5 "." $6 "." $7}');
do python /home/eyscott/scratch/19.12.30_RedGreen/Demult_fastqs/umi_parser.py --sam_file ${f}.sam; done
