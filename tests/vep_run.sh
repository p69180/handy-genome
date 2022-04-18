#!/usr/bin/bash
#SBATCH -c 2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o /dev/null

/home/users/pjh/scripts/conda_wrapper/vep_v104 \
	-i /home/users/pjh/practice/pipeline_test/vep/input.vcf.gz \
	-o ./vep_test_out.vcf.gz \
	--fasta /home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta \
	--cache \
	--species homo_sapiens \
	--assembly GRCh37 \
	--offline \
	--everything \
	--no_stats \
	--vcf \
	--dir /home/users/pjh/.vep