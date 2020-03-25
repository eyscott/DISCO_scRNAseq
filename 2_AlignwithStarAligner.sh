module load nixpkgs/16.09
module load star/2.5.3a

for f in $(ls [G]*R2_001_demultiplexed.fastq | awk -F'[-_.]' '{print $1 "_" $2 "_" $3}');
do STAR --genomeDir /home/eyscott/project/eyscott/eyscott/RG_scRNA/human_genome \
--readFilesIn /home/eyscott/project/eyscott/eyscott/RG_scRNA/Wheeler_multi/Wheeler/${f}_R2_001_demultiplexed.fastq \
--outFileNamePrefix /home/eyscott/project/eyscott/eyscott/RG_scRNA/Wheeler_multi/Wheeler/${f}_human_star \
--outSAMtype BAM SortedByCoordinate --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --outFilterMismatchNmax 2; done
