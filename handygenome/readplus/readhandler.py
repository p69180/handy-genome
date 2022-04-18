import sys
import importlib
import itertools

import pysam

top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))


'''
class Params(common.Params):
	pass
'''


def filter_bad_read(read):
	# returns True if the read is to be excluded
	if \
	(not read.is_paired) or \
	read.is_unmapped or \
	read.mate_is_unmapped or \
	read.is_duplicate:
		return True
	else:
		return False

###

def get_fetch(
		bam, CHROM, 
		POS0 = None, 
		start = None, 
		end = None, 
		filter_fun = None,
		):
	"""
	- returns a generator
	- similar to bam.fetch iterator except read filtering
	"""

	# set default filtering function
	if filter_fun == None:
		filter_fun = filter_bad_read

	if POS0 != None:
		fetch = bam.fetch(CHROM, POS0, POS0+1)
	else:
		fetch = bam.fetch(CHROM, start, end)
	
	for read in fetch:
		if not filter_fun(read):
			yield read


def get_fetchresult_list(
		bam, CHROM, 
		POS0 = None, 
		start = None, 
		end = None, 
		filter_fun = None,
		):
	return list( get_fetch(bam, CHROM, POS0 = POS0, start = start, end = end, filter_fun = filter_fun) )


def get_fetchresult_dict(
		bam, CHROM, 
		POS0 = None, 
		start = None, 
		end = None, 
		filter_fun = None,
		):
	# groups the reads with queryname
	fetchresult_dict = dict()
	for read in get_fetch(bam, CHROM, POS0 = POS0, start = start, end = end, filter_fun = filter_fun):
		if read.query_name not in fetchresult_dict:
			fetchresult_dict[read.query_name] = list()
		fetchresult_dict[read.query_name].append(read)

	return fetchresult_dict

###

def get_pairs_dict(read):
	try:
		pairs_dict = dict()
		pairs_dict['raw'] = read.get_aligned_pairs(with_seq = True)
		pairs_dict['querypos0'], pairs_dict['refpos0'], pairs_dict['refseq'] = zip(*pairs_dict['raw'])
			# querypos0: tuple of query positions
			# refpos0: tuple of reference positions
			# refseq: tuple of reference sequences (substitutions in lowercase)
	except ValueError: # when MD tag is absent, get_aligned_pairs method fails with ValueError
		pairs_dict = None

	return pairs_dict


def get_MM(read, pairs_dict):
	if pairs_dict == None:
		MM = None
	else:
		cigarM_indices = [ 
				True if tup[0] != None and tup[1] != None else False 
				for tup in zip(pairs_dict['querypos0'], pairs_dict['refpos0'])
				]
		pairs_dict_seq_subset = itertools.compress(pairs_dict['refseq'], cigarM_indices)
		MM = len([ x for x in pairs_dict_seq_subset if x.islower() ])

	return MM


#####################################################


def check_cigar_clip_inside(read):
	if len(read.cigartuples) <= 2:
		return False
	else:
		cigartuples_inside = read.cigartuples[1:-1]
		cigarops_inside = [ x[0] for x in cigartuples_inside ]
		if 4 in cigarops_inside or 5 in cigarops_inside:
			return True
		else:
			return False


def check_cigar_DN_outside(read):
	# cigar I can be outside!
	cigartuples_rmclip = [ tup for tup in read.cigartuples if tup[0] not in (4,5) ]
	if \
	cigartuples_rmclip[0][0] in (2,3) or \
	cigartuples_rmclip[-1][0] in (2,3):
		return True
	else:
		return False


def check_unexpected_cigar_pattern(read):
	if check_cigar_DN_outside(read):
		return True
	elif check_cigar_clip_inside(read):
		return True
	elif 'B' in read.cigarstring:
		return True
	elif 'P' in read.cigarstring:
		return True
	else:
		return False

#############################################################

def get_padded_seqs(read, fasta):
	if check_unexpected_cigar_pattern(read):
		raise Exception(f'Unexpected cigar pattern. query name: {read.query_name} ; \
cigar string: {read.cigarstring} ; full read string: {read.to_string()}')

	ref_seq_padded = list() 
	read_seq_padded = list()

	ref_seq = list(fasta.fetch(
		read.reference_name, 
		read.reference_start, 
		read.reference_start + read.reference_length,
		))
	read_seq = list(read.query_sequence)

	for idx_cigar, tup in enumerate(read.cigartuples):
		op, l = tup # op: operation ; l: length
		if op in (0,7,8): # cigar M, =, X
			for i in range(l):
				ref_seq_padded.append(ref_seq.pop(0))
				read_seq_padded.append(read_seq.pop(0))
		elif op == 1: # cigar I
			ref_seq_padded.append(None)
			ins_seq_buffer = list()
			for i in range(l):
				ins_seq_buffer.append(read_seq.pop(0))
			read_seq_padded.append(''.join(ins_seq_buffer))
		elif op == 2: # cigar D
			read_seq_padded.append(None)
			del_seq_buffer = list()
			for i in range(l):
				del_seq_buffer.append(ref_seq.pop(0))
			ref_seq_padded.append(''.join(del_seq_buffer))
		elif op == 3: # cigar N
			del ref_seq[:l]
		elif op == 4: # cigar S
			del read_seq[:l]
		elif op == 5: # cigar H
			pass

	return ref_seq_padded, read_seq_padded


def get_MD(read, fasta = None, ref_seq_padded = None, read_seq_padded = None):
	if ref_seq_padded == None or read_seq_padded == None:
		ref_seq_padded, read_seq_padded = get_padded_seqs(read, fasta)

	MD_list = list()
	identical = 0
	for ref_base, read_base in zip(ref_seq_padded, read_seq_padded):
		if ref_base != None and read_base != None:
			if read_base == '=': # means a base identical to reference base
				identical += 1
			else:
				if ref_base == read_base: # identical match
					identical += 1
				else: # different match
					MD_list.append(str(identical))
					identical = 0
					MD_list.append(ref_base)
		elif ref_base == None: # cigar I
			pass
		elif read_base == None: # cigar D
			MD_list.append(str(identical))
			identical = 0
			MD_list.append('^'+ref_base)

	MD_list.append(str(identical))

	return ''.join(MD_list)


def get_NM(read, fasta = None, ref_seq_padded = None, read_seq_padded = None):
	if ref_seq_padded == None or read_seq_padded == None:
		ref_seq_padded, read_seq_padded = get_padded_seqs(read, fasta)

	NM = 0
	for ref_base, read_base in zip(ref_seq_padded, read_seq_padded):
		if ref_base != None and read_base != None:
			if read_base == '=': # means a base identical to reference base
				pass
			else:
				if ref_base == read_base: # identical match
					pass
				else: # different match
					NM += 1
		elif ref_base == None: # cigar I
			NM += len(read_base)
		elif read_base == None: # cigar D
			NM += len(ref_base)

	return NM


def get_MD_unused(read, fasta):
	if check_unexpected_cigar_pattern(read):
		raise Exception(f'Unexpected cigar pattern. query name: {read.query_name} ; cigar string: {read.cigarstring}')

	idx_SEQ = -1
	idx_refseq = -1
	refseq = fasta.fetch(read.reference_name, read.reference_start, read.reference_start + read.reference_length)
	match = 0
	MD_result_list = list()

	for idx_cigar, tup in enumerate(read.cigartuples):
		op, l = tup
		if op == 0:
			for i in range(l):
				idx_SEQ += 1
				idx_refseq += 1
				if read.query_sequence[idx_SEQ] == refseq[idx_refseq]:
					match += 1
				else:
					MD_result_list.append(str(match))
					match = 0
					MD_result_list.append(refseq[idx_refseq])
		elif op == 7: # cigar =
			for i in range(l):
				idx_SEQ += 1
				idx_refseq += 1
				match += 1
		elif op == 8: # cigar X
			for i in range(l):
				MD_result_list.append(str(match))
				match = 0

				idx_SEQ += 1
				idx_refseq += 1

				MD_result_list.append(refseq[idx_refseq])
		elif op == 1: # cigar I
			idx_SEQ += l
		elif op == 2: # cigar D
			MD_result_list.append(str(match))
			match = 0
			MD_result_list.append('^')
			for i in range(l):
				idx_refseq += 1
				MD_result_list.append(refseq[idx_refseq])
		elif op == 3: # cigar N
			idx_refseq += l
		elif op == 4: # cigar S
			idx_SEQ += l
		elif op == 5: # cigar H
			pass

	MD_result_list.append(str(match))

	return ''.join(MD_result_list)


###

def get_pairorient_substring(read, mate=False):
	if mate:
		orientation = 'R' if read.mate_is_reverse else 'F'
		read12 = '2' if read.is_read1 else '1'
	else:
		orientation = 'R' if read.is_reverse else 'F'
		read12 = '1' if read.is_read1 else '2'

	return orientation + read12


def get_pairorient(read):
	if read.mate_is_unmapped:
		pairorient = None
	else:
		if read.template_length == 0:
			pairorient = None
		else:
			if read.template_length > 0:
				substring1 = get_pairorient_substring(read, mate=False)
				substring2 = get_pairorient_substring(read, mate=True)
			elif read.template_length:
				substring1 = get_pairorient_substring(read, mate=True)
				substring2 = get_pairorient_substring(read, mate=False)

			pairorient = substring1 + substring2

	return pairorient


#####

def check_TRA(read):
	return read.reference_name != read.next_reference_name


def check_long_template(
		read, 
		threshold_template_length = common.THRESHOLD_TEMPLATE_LENGTH,
		):
	if abs(read.template_length) > threshold_template_length: # suggestive of SV
		return True
	else:
		return False


def check_primary_alignment(read):
	if \
	(not read.is_secondary) and \
	(not read.is_supplementary):
		return True
	else:
		return False

#####

def get_fiveprime_end(read):
	if read.is_reverse:
		fiveprime_end = read.reference_end - 1
	else:
		fiveprime_end = read.reference_start
	return fiveprime_end


def get_threeprime_end(read):
	if read.is_reverse:
		threeprime_end = read.reference_start
	else:
		threeprime_end = read.reference_end - 1
	return threeprime_end

#####

def get_template_range(read):
	'''
	Returns a range object representing the interval between 5' ends of the read pairs.
	It is in 0-based, half-open format (bed format).
	'''
	if read.template_length == 0:
		# TLEN == 0 cases include 
			# different chromosomes between read pairs
			# two reads are on different strands and have identical 5' end base position
		template_range = None
	else:
		if read.template_length > 0:
			start = get_fiveprime_end(read)
			end = start + read.template_length
		elif read.template_length < 0:
			end = get_fiveprime_end(read) + 1 # bed format
			start = end + read.template_length

		template_range = range(start, end)

	return template_range

#####

def check_mateless(read):
	if \
	(not read.is_paired) or \
	read.mate_is_unmapped:
		return True
	else:
		return False


def get_primary_mate_candidate(bam, read):
	primary_mate_candidate = list()
	for new_read in bam.fetch(
			read.next_reference_name, 
			read.next_reference_start, 
			read.next_reference_start + 1
			):
		if check_primary_alignment(new_read):
			if \
			new_read.query_name == read.query_name and \
			read.compare(new_read) != 0: # read.compare method returns 0 if the two reads are identical
				primary_mate_candidate.append(new_read)
		else:
			continue
	
	return primary_mate_candidate
		

def get_primary_mate(read, bam):
	"""
	Returns None if the input read is unpaired or the mate is unmapped.
	Else:
		1) Performs fetch with input bam in a position specified by RNEXT and PNEXT of the input read.
		2) From the fetch result, picks reads with the same QNAME as the input read.
		3) If there are only one picked read, returns it. Otherwise, raises an exception.
	"""

	if check_mateless(read):
		mate = None
	else:
		primary_mate_candidate = get_primary_mate_candidate(bam, read)

		if len(primary_mate_candidate) == 0:
			raise Exception(f'''\
No primary mate candidate found for this read:
{read.to_string()}''')
		elif len(primary_mate_candidate) == 1:
			mate = primary_mate_candidate[0]
		elif len(primary_mate_candidate) > 1:
			mate_candidate_strings = '\n'.join([read.to_string() for read in primary_mate_candidate])
			e_msg = f'''\
More than one primary mate candidates found for this read:
{read.to_string()}
Detected mate candidates for this read:
{mate_candidate_strings}'''
			raise Exception(e_msg)

	return mate

#####

def get_softclip_ends_range(read):
	if read.cigartuples[0][0] == 4: # first cigar is S
		start = read.reference_start - read.cigartuples[0][1]
	else:
		start = read.reference_start

	if read.cigartuples[-1][0] == 4: # last cigar is S
		end = read.reference_end + read.cigartuples[-1][1]
	else:
		end = read.reference_end

	return range(start, end)
