#!/home/users/pjh/conda_bin/python

import re
pat_aln_parse = re.compile('^-*([^-](.*[^-])?)-*$')
from pprint import pprint


def show_attributes(obj):
	for attr in dir(obj):
		if not attr.startswith("_"):
			print(attr, getattr(obj, attr), sep = '\t')


def make_reads_into_fasta():
	import pysam
	from Bio.SeqRecord import SeqRecord
	from Bio.Seq import Seq
	from Bio import SeqIO

	NBAM_PATH = '/home/users/pjh/Projects/indel_realignment/changhyunnam/9177697-Blood.s.md.ir.br.bam'
	TBAM_PATH = '/home/users/pjh/Projects/indel_realignment/changhyunnam/9177697-Tumor.s.md.ir.br.bam'
	CHROM = '14'
	START = 67_144_260
	END = 67_144_570
	POS = 67_144_415

	NBAM_FASTA_PATH = './nbam_reads.fa'
	TBAM_FASTA_PATH = './tbam_reads.fa'

	seqrec_list = list()
	for read in pysam.AlignmentFile(NBAM_PATH).fetch(CHROM, POS - 1, POS):
		seqrec = SeqRecord(Seq(read.query_sequence), id = read.query_name)
		seqrec_list.append(seqrec)
	SeqIO.write(seqrec_list, NBAM_FASTA_PATH, "fasta")

	seqrec_list = list()
	for read in pysam.AlignmentFile(TBAM_PATH).fetch(CHROM, POS - 1, POS):
		seqrec = SeqRecord(Seq(read.query_sequence), id = read.query_name)
		seqrec_list.append(seqrec)
	SeqIO.write(seqrec_list, TBAM_FASTA_PATH, "fasta")


def get_aligner_1():
	import Bio.Align

	aligner = Bio.Align.PairwiseAligner(
			match_score = 1.0,
			mismatch_score = -2,
			query_internal_open_gap_score = -5,
			query_internal_extend_gap_score = -1,
			target_internal_open_gap_score = -5,
			target_internal_extend_gap_score = -1,
			)

	return aligner


def get_aligner_2():
	import Bio.Align

	def gap_score_function(start, length):
		return -5 - (length - 1) - (0.01*start)

	aligner = Bio.Align.PairwiseAligner()
	aligner.match_score = 1
	aligner.mismatch_score = -2
	aligner.gap_score = gap_score_function
	aligner.left_gap_score = 0
	aligner.right_gap_score = 0

	return aligner


# get_alignment_spec_1 #

def subtract_tup_pairs(tup_pair1, tup_pair2):
	'''
	returns tup_pair1 - tup_pair2
	tup_pair: ( (target_start, target_end), (query_start, query_end) )
	'''
	return ( ( tup_pair2[0][1], tup_pair1[0][0] ), ( tup_pair2[1][1], tup_pair1[1][0] ), )


#def get_idx_tuples(alignment):
#	idx_tuples = list()
#
#	target_len = len(alignment.target)
#	query_len = len(alignment.query)
#
#	first_tup_pair = ( (0,0), (0,0) )
#	last_tup_pair = ( (target_len, target_len), (query_len, query_len) )
#
#	match_idx_tuples = list(zip(*alignment.aligned))
#	
#	for idx in range(len(match_idx_tuples)):
#		if idx == 0:
#			if not (match_idx_tuples[idx][0][0] == 0 and match_idx_tuples[idx][1][0] == 0):
#				idx_tuples.append( subtract_tup_pairs(match_idx_tuples[idx], first_tup_pair) )
#			idx_tuples.append(match_idx_tuples[idx])
#
#		if idx > 0 and idx <= len(match_idx_tuples) - 1:
#			idx_tuples.append( subtract_tup_pairs(match_idx_tuples[idx], match_idx_tuples[idx-1]) )
#			idx_tuples.append( match_idx_tuples[idx] )
#
#		if idx == len(match_idx_tuples) - 1:
#			if not (match_idx_tuples[idx][0][1] == target_len and match_idx_tuples[idx][1][1] == query_len):
#				idx_tuples.append( subtract_tup_pairs(last_tup_pair, match_idx_tuples[idx]) )
#
#	return idx_tuples


def get_idx_tuples(alignment):
	idx_tuples = list()
	for idx in range(1, len(alignment.path)):
		idx_tuples.append( tuple(zip(alignment.path[idx-1], alignment.path[idx])) )

	return idx_tuples


def get_target_match_range(idx_tuples):
	candidates = list(filter(
			lambda x: x[0][1] - x[0][0] > 0 and x[1][1] - x[1][0] > 0,
			idx_tuples
			))
	return range(candidates[0][0][0], candidates[-1][0][1])


#def get_cigar_attributes(idx_tuples):
#	cigarstring = ''
#	cigarstats = [0]*9
#	for idx, tup_pair in enumerate(idx_tuples):
#		target_len = tup_pair[0][1] - tup_pair[0][0]
#		query_len = tup_pair[1][1] - tup_pair[1][0]
#
#		if target_len == 0 and query_len != 0: # insertion
#			cigarstring += str(query_len)
#			cigarstring += 'I'
#			cigarstats[1] += query_len
#		elif target_len != 0 and query_len == 0: # deletion
#			if idx != 0 and idx != len(idx_tuples)-1:
#				cigarstring += str(target_len)
#				cigarstring += 'D'
#				cigarstats[2] += target_len
#		else: # match
#			cigarstring += str(target_len)
#			cigarstring += 'M'
#			cigarstats[0] += target_len
#
#	return cigarstring, cigarstats


def get_cigartuples(idx_tuples):
	cigartuples_list = list()

	for idx, tup_pair in enumerate(idx_tuples):
		target_len = tup_pair[0][1] - tup_pair[0][0]
		query_len = tup_pair[1][1] - tup_pair[1][0]

		if target_len == 0 and query_len != 0: # insertion
			cigartuples_list.append((1, query_len))
		elif target_len != 0 and query_len == 0: # deletion
			if idx != 0 and idx != len(idx_tuples)-1:
				cigartuples_list.append((2, target_len))
		else: # match
			cigartuples_list.append((0, target_len))

	return tuple(cigartuples_list)


def get_alignment_spec_1(alignment):
	idx_tuples = get_idx_tuples(alignment)
	cigarstring, cigarstats = get_cigar_attributes(idx_tuples)

	alignment_spec = {
			'idx_tuples' : idx_tuples,
			'cigarstring' : cigarstring,
			'cigarstats' : cigarstats,
			}

	return alignment_spec

#####################################

def get_new_alignment(read, target_seq, aligner):
	alignments = aligner.align(target_seq, read.query_sequence)
	alignment_chosen = min(alignments, key =  lambda x: get_gap_start_score(x))
	return alignment_chosen


def get_query_reference_start_end(idx_tuples, target_reference_start):
	target_match_range = get_target_match_range(idx_tuples)
	query_reference_start = target_reference_start + target_match_range.start
	query_reference_end = target_reference_start + target_match_range.stop

	return query_reference_start, query_reference_end


def get_alignment_spec_3(alignment, target_reference_start):
	idx_tuples = get_idx_tuples(alignment)
	cigartuples = get_cigartuples(idx_tuples)
	query_reference_start, query_reference_end = get_query_reference_start_end(idx_tuples, target_reference_start)

	return {
			'idx_tuples' : idx_tuples,
			'cigartuples' : cigartuples,
			'query_reference_start' : query_reference_start,
			'query_reference_end' : query_reference_end,
			}
#	return (
#			idx_tuples,
#			cigartuples,
#			new_reference_start,
#			new_reference_end,
#			)


def get_new_read(read, cigartuples, reference_start, next_reference_start, template_length):
	import copy
	new_read = copy.deepcopy(read)
	new_read.cigartuples = cigartuples
	new_read.reference_start = reference_start
	new_read.next_reference_start = next_reference_start
	new_read.template_length = template_length

	return new_read



########################


# get_alignment_spec_2 #

def op_classifier(x):
	if x[1] == '|':
		return 0 # match
	elif x[1] == '.':
		return 1 # mismatch
	elif x[1] == '-':
		if x[0] == '-':
			return 2 # target gap
		elif x[2] == '-':
			return 3 # query gap


def get_alignment_spec_2(alignment):
	import itertools

	alignment_spec = {
			'target_left_gap' : 0,
			'query_left_gap' : 0,
			'target_internal_gap' : 0,
			'query_internal_gap' : 0,
			'target_right_gap' : 0,
			'query_right_gap' : 0,
			'match' : 0,
			'mismatch' : 0,
			}

	aln_str_list = alignment.format().strip().split('\n')
	op_class = map(op_classifier, zip(*aln_str_list))

	aln_str_list.append( range(len(aln_str_list[0])) )
	aln_str_list.append(op_class)

	aln_str_tuple_list = list( zip(*aln_str_list) ) # list of tuple (target, operation, query, index, op_class)
	aln_list_grouped = [ (k, list(v)) for k, v in itertools.groupby(aln_str_tuple_list, key = lambda x: x[4]) ]

	for idx, tup in enumerate(aln_list_grouped):
		k, v = tup
		if k == 0:
			alignment_spec['match'] += len(v)
		elif k == 1:
			alignment_spec['mismatch'] += len(v)
		elif k == 2:
			if idx == 0:
				alignment_spec['target_left_gap'] += len(v)
			elif idx == len(aln_list_grouped) - 1:
				alignment_spec['target_right_gap'] += len(v)
			else:
				alignment_spec['target_internal_gap'] += len(v)
		elif k == 3:
			if idx == 0:
				alignment_spec['query_left_gap'] += len(v)
			elif idx == len(aln_list_grouped) - 1:
				alignment_spec['query_right_gap'] += len(v)
			else:
				alignment_spec['query_internal_gap'] += len(v)

	return alignment_spec, aln_list_grouped

####################

def get_gap_start_score(alignment, idx_tuples = None):
	if idx_tuples == None:
		idx_tuples = get_idx_tuples(alignment)

	score = 0
	for idx, tup_pairs in enumerate(idx_tuples):
		if idx != 0 and idx != len(idx_tuples)-1:
			if tup_pairs[0][1] - tup_pairs[0][0] != tup_pairs[1][1] - tup_pairs[1][0]: # target or query gap
				score += tup_pairs[0][0]

	return score


####################

def fun1():
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	from Bio.Align import MultipleSeqAlignment
	from Bio import AlignIO
	
	align1 = MultipleSeqAlignment(
			[
				SeqRecord(Seq("ACTGCTAGCTAG"), id="Alpha"),
				SeqRecord(Seq("ACT-CTAGCTAG"), id="Beta"),
				SeqRecord(Seq("ACTGCTAGDTAG"), id="Gamma"),
			]
			)
	align2 = MultipleSeqAlignment(
			[
				SeqRecord(Seq("GTCAGC-AG"), id="Delta"),
				SeqRecord(Seq("GACAGCTAG"), id="Epsilon"),
				SeqRecord(Seq("GTCAGCTAG"), id="Zeta"),
			]
	)
	align3 = MultipleSeqAlignment(
			[
				SeqRecord(Seq("ACTAGTACAGCTG"), id="Eta"),
				SeqRecord(Seq("ACTAGTACAGCT-"), id="Theta"),
				SeqRecord(Seq("-CTACTACAGGTG"), id="Iota"),
			]
	)
	my_alignments = [align1, align2, align3]
	
	
	AlignIO.write(my_alignments, "my_example.phy", "phylip")


def fun2():
	from Bio import AlignIO
	count = AlignIO.convert("PF05356_seed.sth", "stockholm", "PF05356_seed.aln", "clustal")
	print("Conveted %i alignments" % count)


def fun3():
	from Bio import AlignIO
	aln = AlignIO.read("PF05356_seed.sth", "stockholm")
	print(aln + aln)


def fun4():
	from Bio import SeqIO
	rec = next(SeqIO.parse("example.fastq", "fastq"))
	left = rec[:20]
	right = rec[22:]
	right.id = 'hihi'
	new = left + right
	print(new)


def fun5():
	import numpy as np
	from Bio import AlignIO
	aln = AlignIO.read('PF05356_seed.sth', 'stockholm')
	rec = aln[0]
	print(list(rec))
	print(rec.seq)
	#aln_arr = np.array([list(rec) for rec in aln], np.character)
	#print(aln_arr)


def fun6():
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	from Bio.Align import MultipleSeqAlignment
	aln = MultipleSeqAlignment(
			[
				SeqRecord(Seq("ACTCCTA"), id='seq1'),
				SeqRecord(Seq("AAT-CTA"), id='seq2'),
				SeqRecord(Seq("CCTACT-"), id='seq3'),
				SeqRecord(Seq("TCTCCTC"), id='seq4'),
			]
	)
	print(type(aln.substitutions))


def fun7():
	from Bio import pairwise2
	from Bio import SeqIO
	from Bio.Align import substitution_matrices

	blosum62 = substitution_matrices.load("BLOSUM62")
	seq1 = SeqIO.read("alpha.faa", "fasta")
	seq2 = SeqIO.read("beta.faa", "fasta")

	#aln = pairwise2.align.globalds(seq1.seq, seq2.seq, blosum62, -10, -0.5)
	aln = pairwise2.align.localds("LSPADKTNVKAA", "PEEKSAV", blosum62, -10, -1)

	#print(blosum62)
	print(pairwise2.format_alignment(*aln[0], full_sequences = True))


def fun8():
	import Bio.Align

	aligner = Bio.Align.PairwiseAligner()
	#aligner.mode = "local"
	seq1 = 'AGAACTC'
	seq2 = 'GAACT'
	alignments = aligner.align(seq1, seq2)
	for alignment in alignments:
		print(alignment)


def fun9():
	import pysam
	import Bio.Seq
	import Bio.SeqRecord
	import Bio.Align
	import Bio.Align.AlignInfo

	from julib.readplus.readpluspairlist import Rpplist

	FASTA_PATH = '/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta'
	NBAM_PATH = '/home/users/pjh/Projects/indel_realignment/changhyunnam/9177697-Blood.s.md.ir.br.bam'
	TBAM_PATH = '/home/users/pjh/Projects/indel_realignment/changhyunnam/9177697-Tumor.s.md.ir.br.bam'
	TBAM_REALIGNED_PATH = '/home/users/pjh/scripts/python_genome_packages/tbam_realigned.bam'

	CHROM = '14'
	POS = 67_144_413
	FETCH_START = 67_144_405
	FETCH_END = 67_144_420
	START = 67_144_257
	END = 67_144_567

	POS0 = POS - 1

	fasta = pysam.FastaFile(FASTA_PATH)
	tbam = pysam.AlignmentFile(TBAM_PATH)
	nbam = pysam.AlignmentFile(NBAM_PATH)

	rpplist_tbam = Rpplist(tbam, fasta, CHROM, POS0)
	rpplist_nbam = Rpplist(nbam, fasta, CHROM, POS0)
	mate_search_range = range(
			min(rpplist_tbam.mate_search_range[0], rpplist_nbam.mate_search_range[0]), 
			max(rpplist_tbam.mate_search_range[-1]+1, rpplist_nbam.mate_search_range[-1]+1),
			)
	#rpplist = rpplist_tbam + rpplist_nbam
	target_reference_start, target_reference_end = mate_search_range[0], mate_search_range[-1]+1
	ref_seq = fasta.fetch(CHROM, target_reference_start, target_reference_end)

	#seqrec_list = list()
	#for rpp in rpplist_tbam.rpplist + rpplist_nbam.rpplist:
	#	seqrec_list.append( Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(rpp.rp1.read.query_sequence), id = rpp.rp1.read.query_name) )
	#	seqrec_list.append( Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(rpp.rp2.read.query_sequence), id = rpp.rp2.read.query_name) )

	#################################

	aligner = get_aligner_1()

	alignment_list_tumor = list()
	for rpp in rpplist_tbam.rpplist:
		for read in (rpp.rp1.read, rpp.rp2.read):
			alignments = aligner.align(ref_seq, read.query_sequence)
			alignment_chosen = min( 
					[ (x, get_alignment_spec_1(x), read.query_name,) for x in alignments ],
					key = lambda x: get_gap_start_score(x[1]['idx_tuples']),
					)
			alignment_list_tumor.append(alignment_chosen)

	alignment_list_normal = list()
	for rpp in rpplist_nbam.rpplist:
		for read in (rpp.rp1.read, rpp.rp2.read):
			alignments = aligner.align(ref_seq, read.query_sequence)
			alignment_chosen = min( 
					[ (x, get_alignment_spec_1(x), read.query_name,) for x in alignments ],
					key = lambda x: get_gap_start_score(x[1]['idx_tuples']),
					)
			alignment_list_normal.append(alignment_chosen)


	print(mate_search_range) ; print()	

	for alignment, spec, query_name in alignment_list_tumor:
		print('TUMOR', query_name)
		print(alignment)
		print(alignment.aligned)
		pprint(spec)
		print()

	for alignment, spec, query_name in alignment_list_normal:
		print('NORMAL', query_name)
		print(alignment)
		print(alignment.aligned)
		pprint(spec)
		print()


def fun10():
	import Bio.SeqIO
	import Bio.Align

	# AC--GTACCAT
	# ACGT-------

	aln = Bio.Align.PairwiseAlignment(
			target = 'ACGTACCAT',
			query = 'ACGT',
			path = (
				(0,0),
				(2,2),
				(2,2),
				),
			score = 0,
			)
	print(aln)


def fun11():
	import pysam
	import Bio.Seq
	import Bio.SeqRecord
	import Bio.Align
	import Bio.Align.AlignInfo

	from julib.readplus.readpluspairlist import Rpplist

	FASTA_PATH = '/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta'
	NBAM_PATH = '/home/users/pjh/Projects/indel_realignment/changhyunnam/9177697-Blood.s.md.ir.br.bam'
	TBAM_PATH = '/home/users/pjh/Projects/indel_realignment/changhyunnam/9177697-Tumor.s.md.ir.br.bam'
	TBAM_REALIGNED_PATH = '/home/users/pjh/scripts/python_genome_packages/tbam_realigned.bam'

	CHROM = '14'
	POS = 67_144_413
	FETCH_START = 67_144_405
	FETCH_END = 67_144_420
	START = 67_144_257
	END = 67_144_567

	POS0 = POS - 1

	fasta = pysam.FastaFile(FASTA_PATH)
	tbam = pysam.AlignmentFile(TBAM_PATH)
	nbam = pysam.AlignmentFile(NBAM_PATH)


	rpplist_tbam = Rpplist(tbam, fasta, CHROM, POS0)
	pad = 50
	target_reference_start, target_reference_end = rpplist_tbam.mate_search_range.start - pad, rpplist_tbam.mate_search_range.stop + pad
	ref_seq = fasta.fetch(CHROM, target_reference_start, target_reference_end)

	#################################

	aligner = get_aligner_1()

	with pysam.AlignmentFile(TBAM_REALIGNED_PATH, "wb", header = tbam.header) as tbam_realigned:
		for rpp in rpplist_tbam.rpplist:
			alignment_read1 = get_new_alignment(rpp.rp1.read, ref_seq, aligner)
			alignment_read2 = get_new_alignment(rpp.rp2.read, ref_seq, aligner)
			spec_read1 = get_alignment_spec_3(alignment_read1, target_reference_start)
			spec_read2 = get_alignment_spec_3(alignment_read2, target_reference_start)
			tlen = spec_read2['query_reference_end'] - spec_read1['query_reference_start']

			new_read1 = get_new_read(
					rpp.rp1.read, 
					spec_read1['cigartuples'], 
					spec_read1['query_reference_start'], 
					spec_read2['query_reference_start'], 
					tlen,
					)
			new_read2 = get_new_read(
					rpp.rp2.read, 
					spec_read2['cigartuples'], 
					spec_read2['query_reference_start'], 
					spec_read1['query_reference_start'], 
					-1*tlen,
					)

			tbam_realigned.write(new_read1)
			tbam_realigned.write(new_read2)


def fun12():
	import pysam

	TBAM_PATH = '/home/users/pjh/Projects/indel_realignment/changhyunnam/9177697-Tumor.s.md.ir.br.bam'
	tbam = pysam.AlignmentFile(TBAM_PATH)
	read = next(tbam)

	print(read.to_string())
	print(read.cigarstring)
	print(read.cigartuples)

	read.cigartuples = ((0, len(read.query_sequence)),)
	read.reference_start = 12345
	read.next_reference_start = 999
	read.template_length = 5555555

	print(read.to_string())
	print(read.cigarstring)
	print(read.cigartuples)


def fun13():
	TBAM_REALIGNED_PATH = '/home/users/pjh/scripts/python_genome_packages/tbam_realigned.bam'
	bam = pysam.AlignmentFile(TBAM_REALIGNED_PATH)


def main():
	fun11()


main()
