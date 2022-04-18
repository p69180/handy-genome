#!/home/users/pjh/conda_bin/python

from pprint import pprint
"""
import pysam

from julib.variantplus.vcfplus import VcfPlus
from julib.variantplus.variantplus import VariantPlus
from julib.variantplus.vcfhandler import update_header_contigs
from julib.variantplus.vcfhandler import split_vcf
from julib.common import ChromDict
"""


#vcf_path = '/home/users/pjh/practice/pipeline_test/vcf_annotation/input.vcf.gz'
vcf_path = '/home/users/pjh/practice/pipeline_test/delly/delly_wrapper_test_v1.2_211122.vcf.gz'
#vcf_path = '/home/users/pjh/scripts/annotation/SV/vep_wrapper_for_SV/testinput_gremlin.vcf.gz'
#vcf_path = '/home/users/pjh/practice/pipeline_test/vep/input.raw_vep.vcf.gz'
#vcf_path = '/home/users/pjh/practice/pipeline_test/vep/zzzz.vcf.gz'
#vcf_path = '/home/users/pjh/practice/pipeline_test/example_SV_vcfs_by_caller/301d6ce3-4099-4c1d-8e50-c04b7ce91450.broad-snowman.20150918.somatic.sv.vcf.gz'
#vcf_path = '/home/users/pjh/practice/pipeline_test/vcf_annotation/input.vcf.gz'
#vcf_path = '/home/users/pjh/practice/pipeline_test/vep/input.vep_wrapper.vcf.gz'
#vcf_path = '/home/users/pjh/practice/pipeline_test/vep/input.raw_vep.vcf.gz'

VCF_PATH = vcf_path

FASTA_PATH = '/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta'
FASTA_PATH_HG38 = '/home/users/data/01_reference/human_g1k_v38/Homo_sapiens_assembly38.fasta'


def vep_formatter_test():
	fasta = pysam.FastaFile(fasta_path)
	chromdict = ChromDict(fasta_path = fasta_path)


def split_vcf_test():
	outdir = './split_vcfs'
	split_vcf(
			outdir,
			vcf_path = vcf_path,
			n_line = 1000,
			)


def slurmutils_test():
	import time
	from julib.workflow import Job

	jobscript_string = f'''\
#!/usr/bin/bash
#SBATCH -J TEST
set -eu

echo TEST JOB
'''
	job = Job(jobscript_string = jobscript_string)
	print('finished', job.check_finished())
	print('submitted', job.submitted)

	job.submit()
	print('jobid', job.jobid)

	while True:
		time.sleep(2)
		finished = job.check_finished()
		print('finished', finished)
		print('no_jobid', job.no_jobid)

		if job.no_jobid:
			break


def vep_test():
	import pysam
	import julib.vep as vep
	import julib.common as common
	from julib.variantplus.vcfhandler import create_pysamvr

	dic, can_dic = vep.run_vep_minimal('7', 55_086_971, 'A', 'C', fasta_path, species = 'homo_sapiens', assembly = 'GRCh37')
	pprint(can_dic)


def vcfhandler_test():
	import pysam
	import julib.variantplus.vcfplus as vcfplus
	import julib.variantplus.vcfhandler as vcfhandler

	fasta = pysam.FastaFile(FASTA_PATH)
	vcfp1 = vcfplus.VcfPlus(fasta, vcf_path = './test.vcf')
	vcfp2 = vcfplus.VcfPlus(fasta, vcf_path = './test2.vcf')
	merged_vcfp = vcfhandler.merge_vcfplus([vcfp1, vcfp2], isec = True)

	print(merged_vcfp.header)
	for vp in merged_vcfp.vp_list:
		print(vp.pysamvr)


def vep_test():
	import julib.bin.vep_annotate as vep_annotate

	infile_path = './input.vcf.gz'
	outfile_path = './input.vep.vcf.gz'
	vep_annotate.worker(infile_path, outfile_path, FASTA_PATH, 'homo_sapiens', 'GRCh37', sched = 'slurm', intv_check = 10)


def vpp_test():
	from julib.variantplus.vcfplus import VcfPlus
	import pysam
	fasta = pysam.FastaFile(FASTA_PATH)
	#infile_path = '/home/users/pjh/practice/pipeline_test/example_SV_vcfs_by_caller/301d6ce3-4099-4c1d-8e50-c04b7ce91450.broad-dRanger.20150918.somatic.sv.vcf.gz'
	infile_path = './mixed_SV.vcf'
	#infile_path = '/home/users/pjh/practice/pipeline_test/example_SV_vcfs_by_caller/301d6ce3-4099-4c1d-8e50-c04b7ce91450.broad-snowman.20150918.somatic.sv.vcf.gz'
	vcfp = VcfPlus(vcf_path = infile_path, fasta=fasta)

	for vpp in vcfp.vpp_list:
		print(vpp.vp1.pysamvr)
		print('is_bnd1', vpp.vp1.is_bnd1)
		print('is_SV', vpp.vp1.is_SV)
		if vpp.vp1.is_SV:
			print('HOMSEQ', vpp.vp1.HOMSEQ)
			vpp.vp1.bnds.show()

		if vpp.vp2 is not None:
			print(vpp.vp2.pysamvr)
			print('is_bnd1', vpp.vp2.is_bnd1)
			print('is_SV', vpp.vp2.is_SV)
			print('HOMSEQ', vpp.vp2.HOMSEQ)
			vpp.vp2.bnds.show()
		print('\n@@@@@@@@@@@@\n')


def vep_motif_annotation_test():
	from julib.common import ChromDict
	from julib.variantplus.vcfhandler import create_pysamvr
	from julib.veplib import run_vep_pysamvr
	from julib.veplib import get_vepannot_pysamvr
	from julib.ensembl_rest import overlap

	var = ('chr19', 7068400, 7071500, 'C', '<DEL>')
	chrom, pos, end, ref, alt = var

	chromdict = ChromDict(fasta_path = FASTA_PATH_HG38)
	pysamvr = create_pysamvr(chrom, pos, ref, alt, end, chromdict)
	pysamvr = run_vep_pysamvr(pysamvr, FASTA_PATH_HG38, 'homo_sapiens', 'GRCh38')
	vepannot = get_vepannot_pysamvr(pysamvr, 'CSQ')

#	print(var)
#	for feature in vepannot.features:
#		feature.show()
#	exit()

	rest_result = overlap(chrom, pos, end, hg19 = False)
	validated_motifs = list()
	unvalidated_motifs = list()
	for x in rest_result:
		if 'epigenomes_with_experimental_evidence' in x:
			validated_motifs.append(x['stable_id'])
		else:
			unvalidated_motifs.append(x['stable_id'])

	n_validated = 0
	n_unvalidated = 0
	for feature in vepannot.features:
		if feature.is_motif:
			if feature.raw_items['Feature'] in validated_motifs:
				n_validated += 1
			elif feature.raw_items['Feature'] in unvalidated_motifs:
				n_unvalidated += 1
			else:
				raise Exception()

	print(len(vepannot.features))
	print(n_validated)
	print(n_unvalidated)


def ensembl_biotype_listing():
	import requests

	url = 'https://rest.ensembl.org/info/biotypes/homo_sapiens?'
	result = requests.get(url, headers = {'Content-Type' : 'application/json'}).json()
	#pprint(result) ; exit()
	for dic in sorted(result, key = lambda x: x['biotype']):
		print(dic['biotype'])


def scratch():
	from julib.common import ChromDict
	from julib.variantplus.vcfhandler import create_pysamvr
	from julib.annotation.veplib import run_vep_pysamvr
	from julib.annotation.feature import get_feature_list_pysamvr
	from julib.annotation.ensembl_rest import overlap

	#var = ('chr12', 25_225_610, 'T', 'G')
	#var = ('chr2', 10_914_506, ':wA', 'G')
	var = ('chr9', 33_676_000, 'A', 'C')
	chrom, pos, ref, alt = var ; end = None

	#var = ('chr19', 7065000, 7093000, 'A', '<DEL>')
	#var = ('chr19', 44_905_600, 44_909_400, 'G', '<DEL>') 
	#var = ('chr19', 44_906_000, 44_906_200, 'C', '<DEL>')
	#chrom, pos, end, ref, alt = var

	chromdict = ChromDict(fasta_path = FASTA_PATH_HG38)
	pysamvr = create_pysamvr(chrom, pos, ref, alt, end, chromdict)
	pysamvr = run_vep_pysamvr(pysamvr, FASTA_PATH_HG38, 'homo_sapiens', 'GRCh38')
	feat_list = get_feature_list_pysamvr(pysamvr, 'CSQ', hg19 = False)

	print(var)
	for feat in feat_list:
		feat.show()


def scratch2():
	from julib.annotation.feature import parse_rest_regulatory
	from julib.annotation.feature import parse_rest_lookup_transcript
	from julib.annotation.feature import parse_rest_overlap
	from julib.annotation.feature import parse_rest_vep
	from julib.annotation.ensembl_rest import regulatory
	from julib.annotation.ensembl_rest import lookup_id
	from julib.annotation.ensembl_rest import overlap
	from julib.annotation.ensembl_rest import vep

	#result = overlap('chr17', 7476974, 7506975, hg19 = False)
	#parsed = parse_rest_overlap(result, hg19 = False)
	result = vep(chrom = 'chr17', pos = 7_657_780, ref = 'G', alt = 'C', hg19 = False)
	#pprint(result)
	parsed = parse_rest_vep(result, hg19 = False)
	pprint(parsed)


def main():
	scratch2()


main()
