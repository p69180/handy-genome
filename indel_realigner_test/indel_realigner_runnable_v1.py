#!/home/users/pjh/conda_bin/python

import os
import re
import sys
import argparse

import pysam
import Bio.Align

import julib.common as common
from julib.readplus.readpluspairlist import Rpplist

from pprint import pprint


def printerr(*args):
	print(*args, file = sys.stdout, flush = True)


def get_aligner_blastn():
	aligner = Bio.Align.PairwiseAligner(
			match_score = 2,
			mismatch_score = -3,
			query_internal_open_gap_score = -7,
			query_internal_extend_gap_score = -2,
			target_internal_open_gap_score = -7,
			target_internal_extend_gap_score = -2,
			)

	return aligner


def get_aligner_1():
	aligner = Bio.Align.PairwiseAligner(
			match_score = 2,
			mismatch_score = -3,
			query_internal_open_gap_score = -7,
			query_internal_extend_gap_score = -2,
			target_internal_open_gap_score = -7,
			target_internal_extend_gap_score = -2,
			)

	return aligner


def get_aligner_2():
	aligner = Bio.Align.PairwiseAligner(
			match_score = 2,
			mismatch_score = -3,
			query_internal_open_gap_score = -7,
			query_internal_extend_gap_score = -1,
			target_internal_open_gap_score = -7,
			target_internal_extend_gap_score = -1,
			)

	return aligner


def get_gap_start_score(alignment, idx_tuples = None):
	if idx_tuples == None:
		idx_tuples = get_idx_tuples(alignment)

	score = 0
	for idx, tup_pairs in enumerate(idx_tuples):
		if idx != 0 and idx != len(idx_tuples)-1:
			if tup_pairs[0][1] - tup_pairs[0][0] != tup_pairs[1][1] - tup_pairs[1][0]: # target or query gap
				score += tup_pairs[0][0]

	return score


def get_idx_tuples(alignment):
	idx_tuples = list()
	for idx in range(1, len(alignment.path)):
		idx_tuples.append( tuple(zip(alignment.path[idx-1], alignment.path[idx])) )

	return idx_tuples


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


def get_target_match_range(idx_tuples):
	#printerr(idx_tuples)
	candidates = list(filter(
			lambda x: x[0][1] - x[0][0] > 0 and x[1][1] - x[1][0] > 0,
			idx_tuples
			))
	if len(candidates) > 0:
		return range(candidates[0][0][0], candidates[-1][0][1])
	else:
		return None


def get_query_reference_start_end(idx_tuples, target_reference_start):
	target_match_range = get_target_match_range(idx_tuples)
	if target_match_range == None:
		query_reference_start = None
		query_reference_end = None
	else:
		query_reference_start = target_reference_start + target_match_range.start
		query_reference_end = target_reference_start + target_match_range.stop

	return query_reference_start, query_reference_end


####################################################################


def show_all_alignments(read, target_seq, aligner):
	alignments = aligner.align(target_seq, read.query_sequence)
	for alignment in alignments:
		print(alignment)


def get_new_alignment(read, target_seq, aligner):
	alignments = aligner.align(target_seq, read.query_sequence)
	alignment_chosen = min(alignments, key =  lambda x: get_gap_start_score(x))
	return alignment_chosen


def get_alignment_spec(alignment, target_reference_start):
	idx_tuples = get_idx_tuples(alignment)
	cigartuples = get_cigartuples(idx_tuples)
	query_reference_start, query_reference_end = get_query_reference_start_end(idx_tuples, target_reference_start)

	return {
			'idx_tuples' : idx_tuples,
			'cigartuples' : cigartuples,
			'query_reference_start' : query_reference_start,
			'query_reference_end' : query_reference_end,
			}


def get_new_read(read, cigartuples, reference_start, next_reference_start, template_length):
	new_read = read.__copy__()
	new_read.cigartuples = cigartuples
	new_read.reference_start = reference_start
	new_read.next_reference_start = next_reference_start
	new_read.template_length = template_length

	return new_read


####################################################################

def sort_index_bam(bam_path, sorted_bam_path):
	import subprocess
	SAMTOOLS_PATH = '/home/users/pjh/conda_bin/samtools'

	p1 = subprocess.run([
		SAMTOOLS_PATH,
		'sort',
		'-o', sorted_bam_path,
		'-O', 'BAM',
		bam_path,
		])
	p2 = subprocess.run([
		SAMTOOLS_PATH,
		'index',
		sorted_bam_path,
		])


def realign_bam(
		bam_path, 
		realigned_bam_path,
		fasta,
		chrom,
		pos0,
		chromdict,
		target_pad,
		aligner,
		):
	bam = pysam.AlignmentFile(bam_path)

	rpplist = Rpplist(bam, fasta, chrom, pos0, chromdict)
	target_reference_start = rpplist.mate_search_range.start - target_pad
	target_reference_end = rpplist.mate_search_range.stop + target_pad
	ref_seq = fasta.fetch(chrom, target_reference_start, target_reference_end)

	with pysam.AlignmentFile(realigned_bam_path, 'wb', header = bam.header) as realigned_bam:
		for rpp in rpplist.rpplist:
			alignment_read1 = get_new_alignment(rpp.rp1.read, ref_seq, aligner)
			spec_read1 = get_alignment_spec(alignment_read1, target_reference_start)

			alignment_read2 = get_new_alignment(rpp.rp2.read, ref_seq, aligner)
			spec_read2 = get_alignment_spec(alignment_read2, target_reference_start)

			if \
			spec_read1['query_reference_start'] != None and \
			spec_read2['query_reference_start'] != None and \
			spec_read1['query_reference_end'] != None and \
			spec_read2['query_reference_end'] != None:
				tlen = spec_read2['query_reference_end'] - spec_read1['query_reference_start']
				new_read1 = get_new_read(
						rpp.rp1.read, 
						spec_read1['cigartuples'], 
						spec_read1['query_reference_start'], 
						spec_read2['query_reference_start'], 
						tlen,
						)
				realigned_bam.write(new_read1)
				new_read2 = get_new_read(
						rpp.rp2.read, 
						spec_read2['cigartuples'], 
						spec_read2['query_reference_start'], 
						spec_read1['query_reference_start'], 
						-1*tlen,
						)
				realigned_bam.write(new_read2)


def show_realignment(
		bam_path, 
		fasta,
		chrom,
		pos0,
		target_pad,
		aligner,
		):
	bam = pysam.AlignmentFile(bam_path)

	rpplist = Rpplist(bam, fasta, chrom, pos0)
	target_reference_start = rpplist.mate_search_range.start - target_pad
	target_reference_end = rpplist.mate_search_range.stop + target_pad
	ref_seq = fasta.fetch(chrom, target_reference_start, target_reference_end)
	
	rpplist.rpplist.sort(key = lambda rpp: rpp.rp1.read.reference_start)
	for rpp in rpplist.rpplist:
		alignment_read1 = get_new_alignment(rpp.rp1.read, ref_seq, aligner)
		spec_read1 = get_alignment_spec(alignment_read1, target_reference_start)

		alignment_read2 = get_new_alignment(rpp.rp2.read, ref_seq, aligner)
		spec_read2 = get_alignment_spec(alignment_read2, target_reference_start)

		if pos0 in range(spec_read1['query_reference_start'], spec_read1['query_reference_end']):
			print(rpp.rp1.read.query_name, 'read1')
			print(alignment_read1)
		if pos0 in range(spec_read2['query_reference_start'], spec_read2['query_reference_end']):
			print(rpp.rp2.read.query_name, 'read2')
			print(alignment_read2)


####################################################################


def get_pileup(
		bam,
		fasta,
		chrom,
		pos0,
		):
	pup = bam.pileup(
			chrom, pos0, pos0 + 1,
			truncate = True,
			stepper = 'nofilter',
			fastafile = fasta,
			)
	pupcol = next(pup)
	return pupcol


def get_alleles(
		bam_path,
		fasta,
		chrom, 
		pos0,
		):
	with pysam.AlignmentFile(bam_path) as bam:
		rpplist = Rpplist(bam, fasta, chrom, pos0)
		pup = bam.pileup(
				chrom, pos0, pos0 + 1,
				truncate = True,
				stepper = 'nofilter',
				fastafile = fasta,
				)
		pupcol = next(pup)

		print(pupcol.get_query_sequences(
			mark_matches = True,
			mark_ends = True,
			add_indels = True,
			))


####################################################################


class Params:
	target_pad = 50


def argument_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input', required = True, help = 'Input bam')
	parser.add_argument('-o', '--output', required = True, help = 'Output realigned bam')
	parser.add_argument('-f', '--fasta', required = True, help = 'Reference fasta')
	parser.add_argument('--chrom', required = True, help = 'CHROM of realignment position')
	parser.add_argument('--pos', required = True, type = int, help = 'POS of realignment position')

	args = parser.parse_args()
	
	return args


def main():
	args = argument_parser()

	chrom = args.chrom
	pos = args.pos
	pos0 = pos - 1
	fasta = pysam.FastaFile(args.fasta)

	chromdict = common.ChromDict(fasta_path = args.fasta)
	aligner = get_aligner_1()
	tmp_bam_path = common.get_tempfile_path(delete = False)

	realign_bam(
			args.input,
			tmp_bam_path,
			fasta,
			chrom,
			pos0,
			chromdict,
			Params.target_pad,
			aligner,
			)
	sort_index_bam(tmp_bam_path, args.output)

	os.remove(tmp_bam_path)


if __name__ == '__main__':
	main()
