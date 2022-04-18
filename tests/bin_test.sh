#!/bin/bash
#211222
set -eu

fasta_path=/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta

vep_formatter() {
	vep_formatter=../bin/vep_formatter
	vcf_path=/home/users/pjh/practice/pipeline_test/vep/input.raw_vep.vcf.gz
	output_path=./vep_split.vcf.gz

	$vep_formatter -i $vcf_path -o $output_path -f $fasta_path --mode split
}


split_vcf() {
	vcf_path=/home/users/pjh/practice/pipeline_test/delly/delly_wrapper_test_v1.2_211122.vcf.gz
}


# main
vep_formatter
