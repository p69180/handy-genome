
def get_equivalents_ins(chrom, pos, ref, alt, fasta):
	def get_candidate_residuals(insseq):
		candidate_residuals = list()
		candidate_residuals.append(0)
		for idx in range(1, len(insseq)):
			if insseq[:idx] == insseq[-idx:]:
				candidate_residuals.append(idx)

		return candidate_residuals

	def new_pos_iter(candidate_residuals, insseq):
		"""infinite iterator"""
		cycle = 0
		first = True
		while True:
			for x in [ len(insseq)*cycle + i for i in candidate_residuals ]:
				if first:
					prev_val = x
					first = False
					continue
				current_val = x	
				yield (prev_val, current_val)
				prev_val = current_val

			cycle += 1

	def lookback(chrom, pos, fasta, insseq, candidate_residuals, fetch_range_step = 10):
		fetch_extend = 0
		fetched_seq = list()

		for new_pos in new_pos_iter(candidate_residuals, insseq):

			
			
				
			

		fetch_extend += 1
		fetched_seq.extend( 
				list(fasta.fetch(
						chrom, 
						pos - (fetch_extend * fetch_range_step), 
						pos - ((fetch_extend-1) * fetch_range_step)
						)
					).__reversed__() 
				)
			
			
			
			







def _get_equivalents_ins(chrom, pos, ref, alt, fasta):
	"""
	Args:
		pos: 1-based position (as VCF spec)
	"""

	def check_commutative(refseq, insseq):
		return refseq + insseq == insseq + refseq

	def lookback(chrom, pos, new_pos_list, insseq, fasta):
		new_pos = pos
		while True:
			new_pos -= 1
			refseq = fasta.fetch(chrom, new_pos, pos)
			"""
			fetch_start0 = (new_pos + 1) - 1
			fetch_end0 = pos
			"""

			print(new_pos)
			print(refseq)

			if check_commutative(refseq, insseq):
				new_pos_list.append(new_pos)
				continue
			else:
				break

	def lookforward(chrom, pos, new_pos_list, insseq, fasta):
		new_pos = pos
		while True:
			new_pos += 1
			refseq = fasta.fetch(chrom, pos, new_pos)
			"""
			fetch_start0 = (pos + 1) - 1
			fetch_end0 = (new_pos + 1) - 1
			"""

			print(new_pos)
			print(refseq)

			if check_commutative(refseq, insseq):
				new_pos_list.append(new_pos)
				continue
			else:
				break

	def get_new_pos_list(chrom, pos, fasta, insseq):
		new_pos_list = list()
		lookback(chrom, pos, new_pos_list, insseq, fasta)
		lookforward(chrom, pos, new_pos_list, insseq, fasta)

		return new_pos_list

	def get_result(chrom, pos, ref, alt, fasta, insseq, new_pos_list):
		result = list()
		result.append( (chrom, pos, ref, alt) )
		for new_pos in new_pos_list:
			new_ref = fasta.fetch(chrom, new_pos - 1, new_pos)
			new_alt = new_ref + insseq
			result.append( (chrom, new_pos, new_ref, new_alt) )

		return result

	def main(chrom, pos, ref, alt, fasta):
		insseq = alt[1:]
		new_pos_list = get_new_pos_list(chrom, pos, fasta, insseq)
		result = get_result(chrom, pos, ref, alt, fasta, insseq, new_pos_list)

		return result

	return main(chrom, pos, ref, alt, fasta)


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

