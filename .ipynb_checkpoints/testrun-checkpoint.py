import sys

import pysam
import cyvcf2

import julib.common.funcs as common_funcs

refver = 'hg19'
chrom_dict = common_funcs.get_chrom_dict(refver=refver)

vcf_path = '/home/users/team_projects/Lung_Cancer_Panel_data/03_Data_from_YTKim/04_delly/14.tumor/14.tumor.delly.DEL.bcf'

#bam_path = '/home/users/team_projects/Lung_Cancer_Panel_data/03_Data_from_YTKim/02_BAM/14/LU-14.normal.bam'
bam_path = '/home/users/team_projects/Lung_Cancer_Panel_data/03_Data_from_YTKim/02_BAM/14/LU-14.tumor.bam'
#bam_path = '/home/users/team_projects/Lung_Cancer_Panel_data/03_Data_from_YTKim/02_BAM/14/LU-14.panel.bam'

fasta_path = '/home/users/pjh/References/reference_genome/GRCh37/hg1kv37/human_g1k_v37.fasta'


vcf = cyvcf2.VCF(vcf_path)
bam = pysam.AlignmentFile(bam_path)
fasta = pysam.FastaFile(fasta_path)



###


###

def readhandler_test():

	from julib.readplus.readhandler import get_fetch
	from julib.readplus.readhandler import filter_bad_read
	from julib.readplus.readhandler import get_template_range
	from julib.readplus.readhandler import get_pairorient
	from julib.readplus.readhandler import get_mate
	from julib.readplus.readhandler import get_threeprime_end
	from julib.readplus.readhandler import get_fiveprime_end

	CHROM = '2'
	POS = 29_446_540
	POS0 = POS - 1

	for read in get_fetch(bam, CHROM, POS0):
		print(filter_bad_read(read))

	print()

	for read in bam.fetch(CHROM, POS0, POS0+1):
		print(filter_bad_read(read))


####

def proper_pair_check():
	
	CHROM = '9'
	POS = 87_640_283
	POS0 = POS - 1
	
	readlist = list()
	for read in bam.fetch(CHROM, POS0, POS0+1):
		if read.query_name == 'NDX550168:116:HH35KBGXC:3:22412:13035:3958':
			readlist.append(read)
	
	for read in readlist:
		print('qname', read.query_name)
		print('is_reverse', read.is_reverse)
		print('is_read1', read.is_read1)
		print('start', read.reference_start)
		print('end', read.reference_end)
		print('proper_pair', read.is_proper_pair)
		print()

#### read.compare method test

def read_compare_method_test():

	CHROM = '1'
	POS = 27_107_330
	POS0 = POS - 1

	for read_new in bam.fetch(CHROM, POS0, POS0+1):
		if read_new.query_name == 'ST-E00130:582:HC7KJALXX:2:2101:11535:48371':
			read = read_new
			break
	
	print('original read', read.query_name)
	print()
	
	for read_new in bam.fetch(CHROM, POS0, POS0+1):
		print(read_new.query_name, read_new.compare(read))



####
def alleleclass_test():
	from julib.readplus.readplus import ReadPlus
	from julib.ria.rpedit import add_alleleclass

	CHROM = '1'
	POS = 1_397_636
	REF = 'AATGAGATGAAATG'
	ALT = 'A'

	read_list = list(bam.fetch(CHROM, POS0, POS0+1))
	
	for read in read_list:
		rp = ReadPlus(read)
		add_alleleclass(rp, CHROM, POS0, REF, ALT)
	
		print(rp.read.query_name)
		for k,v in rp.__dict__.items():
			print(k, v)
		print()
	
###

def myfun():
	import sys
	import pysam
	import julib.common.funcs as common_funcs

	refver = 'hg19'
	chrom_dict = common_funcs.get_chrom_dict(refver=refver)

	bam_path = '/home/users/team_projects/Lung_Cancer_Panel_data/03_Data_from_YTKim/02_BAM/14/LU-14.tumor.bam'
	bam = pysam.AlignmentFile(bam_path)

	CHROM = '9'
	POS = 87_642_408
	POS0 = POS - 1
	REF = 'T'
	ALT = 'C'

	import collections
	import julib.readplus.readpluspair_list as readpluspair_list
	from julib.ria.get_alleleclass_from_readpluspair import get_seq_REFrange_from_readpluspair
	from julib.ria.get_alleleclass_from_readpluspair import get_alleleclass_from_readpluspair

	rpplist = readpluspair_list.get_readpluspair_list(bam, CHROM, POS0, chrom_dict)
	for rpp in rpplist:
		seq_REFrange = get_seq_REFrange_from_readpluspair(rpp, CHROM, POS0, REF, ALT)
		alleleclass = get_alleleclass_from_readpluspair(rpp, CHROM, POS0, REF, ALT)
		print(seq_REFrange, alleleclass)
	
	#print(collections.Counter([rpp.alleleclass for rpp in rpplist]))


def myfun2():
	import sys
	import pysam
	import julib.varplus.repeats as repeats

	CHROM = '9'
	POS = 87_642_408
	POS0 = POS - 1
	REF = 'T'
	ALT = 'TTTT'

	repeat_list = repeats.get_repeat_list_from_seq(fasta.fetch(CHROM, POS0 - 10, POS0 + 40 + 1))
	for repeat in repeat_list:
		print(repeat.start_localcoord, repeat.end_localcoord, repeat.repeat_unit, repeat.repeat_count, sep='\t')
	return

	repeat_list = repeats.get_repeat_list(CHROM, POS0, fasta)
	relevant_repeat_list = repeats.get_relevant_repeat_list(repeat_list, CHROM, POS0, REF, ALT)

	for repeat in repeat_list:
		print(repeat.CHROM, repeat.start_localcoord, repeat.end_localcoord, repeat.start_refcoord, repeat.end_refcoord, repeat.repeat_unit, repeat.repeat_count, sep='\t')
	print()
	for repeat in relevant_repeat_list:
		print(repeat.CHROM, repeat.start_localcoord, repeat.end_localcoord, repeat.start_refcoord, repeat.end_refcoord, repeat.repeat_unit, repeat.repeat_count, sep='\t')

def myfun3():
	import julib.readplus.readhandler as readhandler

	CHROM = '9'
	POS = 87_645_193
	POS0 = POS - 1
	readlist = list(bam.fetch(CHROM, POS0, POS0+1))
	for read in readlist:
		pairs_dict = readhandler.get_pairs_dict(read)
		MM = readhandler.get_MM(read, pairs_dict)
		print(read.query_name, MM, sep='\t')

####

myfun3()
