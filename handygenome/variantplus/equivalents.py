import itertools

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))


DEFAULT_FETCH_EXTEND_LENGTH = 10


def leftmost(vcfspec, fasta, fetch_extend_length = DEFAULT_FETCH_EXTEND_LENGTH):
	return indel_equivalents(vcfspec, fasta, fetch_extend_length)[0]


def rightmost(vcfspec, fasta, fetch_extend_length = DEFAULT_FETCH_EXTEND_LENGTH):
	return indel_equivalents(vcfspec, fasta, fetch_extend_length)[-1]


def check_equivalent(vcfspec1, vcfspec2, fasta):
	return leftmost(vcfspec1, fasta) == leftmost(vcfspec2, fasta)


def indel_equivalents(vcfspec, fasta, fetch_extend_length = DEFAULT_FETCH_EXTEND_LENGTH):
	"""
	Returns:
		[ vcfspec, vcfspec, ... ] 
			- each vcfspec is equivalent to the input vcfspec
			- sorted by pos (increasing order)

		Return value includes only input 'vcfspec' itself if one or more of the below is true:
			1) ref or alt sequence includes 'n' or 'N'
			2) input vcfspec is not a pure indel
	"""

	def progressive_fetch(chrom, pos, fasta, fetch_extend_length, diff_seq, is_ins, forward):
		extend_count = 0
		while True:
			extend_count += 1
			if forward: # for lookforward
				if is_ins:
					fetch_start = pos + ((extend_count-1) * fetch_extend_length)
					fetch_end =	pos + (extend_count * fetch_extend_length)
				else: # deletion
					fetch_start = pos + len(diff_seq) + ((extend_count-1) * fetch_extend_length)
					fetch_end =	pos + len(diff_seq) + (extend_count * fetch_extend_length)
				result = list(fasta.fetch(chrom, fetch_start, fetch_end))

			else: # for lookback
				fetch_start = pos - (extend_count * fetch_extend_length)
				fetch_end = pos - ((extend_count-1) * fetch_extend_length)
				result = reversed(list(fasta.fetch(chrom, fetch_start, fetch_end)))

			yield result

	def navigate(chrom, pos, fasta, diff_seq, new_pos_list, fetch_extend_length, is_ins, forward):
		"""
		Args:
			is_ins: True for insertion, False for deletion
			forward: True for forward naviation (to 3' side on the positive strand), False for backward navigation
		"""

		fetch_reservoir = list()
		fetched_seq = list()
		fetcher = progressive_fetch(chrom, pos, fasta, fetch_extend_length, diff_seq, is_ins, forward)
		diff_seq_compared = diff_seq if forward else diff_seq[::-1]

		#for start, end in fetch_coord_iter(candidate_residuals, diff_seq):
		for idx in itertools.count(0):
			while len(fetch_reservoir) <= idx + 1 + len(diff_seq):
				fetch_reservoir.extend(next(fetcher))

			idx_remainder = idx % len(diff_seq)

			fetched_seq.append(fetch_reservoir[idx])
			
			if fetched_seq[-1] != diff_seq_compared[idx_remainder]:
				break
			else:
				if is_ins:
					if forward:
						new_pos = pos + (idx + 1)
						new_ref = fetch_reservoir[idx]
						inserted_seq = diff_seq_compared[idx_remainder+1:] + diff_seq_compared[:idx_remainder+1]
					else:
						new_pos = pos - (idx + 1)
						new_ref = fetch_reservoir[idx+1]
						inserted_seq = list(reversed( diff_seq_compared[idx_remainder+1:] + diff_seq_compared[:idx_remainder+1] ))

					new_pos_list.append((new_pos, new_ref, inserted_seq))

				else: # deletion
					# overlapping deletion
					if idx < len(diff_seq):
						if forward:
							new_pos = pos + (idx + 1)
							deleted_seq = diff_seq_compared[idx+1:] + fetched_seq
							new_ref = diff_seq[idx] + ''.join(deleted_seq)
						else:
							new_pos = pos - (idx + 1)
							deleted_seq = list(reversed(diff_seq_compared[idx+1:] + fetched_seq))
							new_ref = fetch_reservoir[idx+1] + ''.join(deleted_seq)

						new_pos_list.append((new_pos, new_ref))

					# non-overlapping deletion
					beyond_seq = fetch_reservoir[ idx + 1 : idx + 1 + len(diff_seq) ]
					if (
							( diff_seq_compared[idx_remainder+1:] == beyond_seq[:-(idx_remainder+1)] ) and
							( diff_seq_compared[:idx_remainder+1] == beyond_seq[-(idx_remainder+1):] )
					   ):
						if forward:
							new_pos = pos + len(diff_seq) + (idx + 1)
							deleted_seq = beyond_seq
							new_ref = fetch_reservoir[idx] + ''.join(deleted_seq)
						else:
							new_pos = pos - (idx + 1) - len(diff_seq)
							deleted_seq = list(reversed(beyond_seq))
							new_ref = fetch_reservoir[idx + len(diff_seq) + 1] + ''.join(deleted_seq)

						new_pos_list.append((new_pos, new_ref))

	def get_new_pos_list(chrom, pos, fasta, diff_seq, fetch_extend_length, is_ins):
		new_pos_list = list()
		navigate(chrom, pos, fasta, diff_seq, new_pos_list, fetch_extend_length, is_ins, forward = True)
		navigate(chrom, pos, fasta, diff_seq, new_pos_list, fetch_extend_length, is_ins, forward = False)

		return new_pos_list

	def get_result(chrom, pos, ref, alt, fasta, is_ins, new_pos_list, diff_seq):
		result = list()
		result.append( common.Vcfspec(chrom, pos, ref, alt) )
		if is_ins:
			for new_pos, new_ref, inserted_seq in new_pos_list:
				new_alt = new_ref + ''.join(inserted_seq)
				result.append( common.Vcfspec(chrom, new_pos, new_ref, new_alt) )
		else: # deletion
			for new_pos, new_ref in new_pos_list:
				new_alt = new_ref[0]
				result.append( common.Vcfspec(chrom, new_pos, new_ref, new_alt) )

		result.sort(key = lambda x: x.pos)

		return result

	def main(vcfspec, fasta, fetch_extend_length):
		chrom, pos, ref, alt = vcfspec

		if set('nN').isdisjoint(ref) and set('nN').isdisjoint(alt):
			mttype = common.get_mttype(ref, alt)
			if mttype in ('ins', 'del'):
				is_ins = (mttype == 'ins')
				diff_seq = list(alt[1:]) if is_ins else list(ref[1:])
				new_pos_list = get_new_pos_list(chrom, pos, fasta, diff_seq, fetch_extend_length, is_ins)
				result = get_result(chrom, pos, ref, alt, fasta, is_ins, new_pos_list, diff_seq)
			else:
				result = [ vcfspec ]
		else:
			result = [ vcfspec ]

		return result
	
	return main(vcfspec, fasta, fetch_extend_length)


'''
def _check_commutative_for_proof(refseq, insseq):
	def cond1(refseq, insseq):
		n = 0
		while True:
			n += 1
			end = n*len(insseq)
			start = (n-1)*len(insseq)
			subseq = refseq[start:end]
			if subseq != insseq:
				return False

			if end in range(len(refseq) - len(insseq) + 1, len(refseq) + 1):
				break

		return end

	def cond2(refseq, insseq):
		return refseq[ -len(insseq) : ] == insseq

	def cond3(refseq, insseq, end):
		#assert end < len(refseq)
		idx = len(refseq) - end
		return refseq[ -idx : ] == insseq[ : idx ]

	def main(refseq, insseq):
		end = cond1(refseq, insseq)
		if end == len(refseq):
			return True
		else:
			if cond2(refseq, insseq):
				return cond3(refseq, insseq, end)
			else:
				return False

	return main(refseq, insseq)
'''

