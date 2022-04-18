import re
import pprint
import itertools
import logging

import pysam
import Bio.Seq

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
workflow = importlib.import_module('.'.join([top_package_name, 'workflow']))
structvars = importlib.import_module('.'.join([top_package_name, 'svlib', 'structvars']))


LOGGER = workflow.get_logger(
		name = __name__,
		formatter = logging.Formatter(
			fmt = '[%(asctime)s  %(levelname)s] %(name)s - %(message)s',
			datefmt = workflow.DEFAULT_DATEFMT,
			),
		level = 'info', stderr = True, filename = None, append = False,
		)


class Breakends:
	"""
	Attributes:
		fasta : pysam.FastaFile instance
		chromdict : julib.common.ChromDict instance
		chrom1
		pos1
		pos1_endis5
		chrom2
		pos2
		pos2_endis5
		inserted_seq : list composed of characters of inserted sequence (e.g. ['A', 'A', 'G', 'C'])
		svtype
		score
		equivs: A list of Breakends objects with maximum reference coverage.
		pos1adv_form = None
		pos2adv_form = None
		homlen = None
		homseq = None
	"""

	BASIC_ATTRS = (
			'chrom1', 
			'pos1', 
			'pos1_endis5', 
			'chrom2', 
			'pos2', 
			'pos2_endis5', 
			'inserted_seq',
			'svtype',
			)

	def __init__(self, 
			chrom1,
			pos1,
			pos1_endis5,
			chrom2,
			pos2,
			pos2_endis5,
			fasta,
			inserted_seq = list(),
			chromdict = None,
			svtype = None,
			score = None,
			#set_bnds_equivs = True,
			):
		# basic attributes
		self.chrom1 = chrom1
		self.pos1 = pos1
		self.pos1_endis5 = pos1_endis5
		self.chrom2 = chrom2
		self.pos2 = pos2
		self.pos2_endis5 = pos2_endis5
		self.inserted_seq = inserted_seq
		self.fasta = fasta

		if chromdict is None:
			self.chromdict = common.ChromDict(fasta = fasta)
		else:
			self.chromdict = chromdict

		# svtype
		if svtype is None:
			self._set_svtype()
		else:
			self.svtype = svtype

		# score
		if score is None:
			self.score = 0
		else:
			self.score = score

		# bnds_equivs
		self.equivs = None
		self.pos1adv_form = None
		self.pos2adv_form = None
		self.homlen = None
		self.homseq = None

	def __repr__(self):
		def func1(self):
			result = ['<Breakends>']
			result.extend([f'{key}: {getattr(self, key)}' 
			               for key in ('chrom1', 'pos1', 'pos1_endis5',
			                           'inserted_seq',
			                           'chrom2', 'pos2', 'pos2_endis5')])
			return '\n'.join(result)

		def func2(self):
			brkt1 = '[' if self.pos1_endis5 else ']'
			brkt2 = '[' if self.pos2_endis5 else ']'
			str1 = f'{brkt1}{self.chrom1}:{self.pos1}{brkt1}'
			str2 = f'{brkt2}{self.chrom2}:{self.pos2}{brkt2}'
			insseq = ''.join(self.inserted_seq)

			return f'<Breakends> {str1} {insseq} {str2}'

		return func1(self)

	def __hash__(self):
		return hash((
					self.chrom1,
					self.pos1,
					self.pos1_endis5,
					self.chrom2,
					self.pos2,
					self.pos2_endis5,
					tuple(self.inserted_seq),
					))

#	def __hash__(self):
#		pos1adv_form = get_pos1adv_form(self)
#		return hash((
#					pos1adv_form.chrom1,
#					pos1adv_form.pos1,
#					pos1adv_form.pos1_endis5,
#					pos1adv_form.chrom2,
#					pos1adv_form.pos2,
#					pos1adv_form.pos2_endis5,
#					pos1adv_form.inserted_seq,
#					))

	def __eq__(self, other):
		#assert isinstance(other, Breakends)
		return hash(self) == hash(other)

	###########################################

	def get_vcfspec_pos1(self):
		chrom = self.chrom1
		pos = self.pos1
		ref = self.fasta.fetch(self.chrom1, self.pos1 - 1, self.pos1)

		if self.pos1_endis5:
			t = ''.join(self.inserted_seq) + ref
		else:
			t = ref + ''.join(self.inserted_seq)
		bracket = '[' if self.pos2_endis5 else ']'
		alt_matestring = f'{bracket}{self.chrom2}:{self.pos2}{bracket}'
		if self.pos1_endis5:
			alt = alt_matestring + t
		else:
			alt = t + alt_matestring

		return common.Vcfspec(chrom, pos, ref, alt)

	def get_vcfspec_pos2(self):
		chrom = self.chrom2
		pos = self.pos2
		ref = self.fasta.fetch(self.chrom2, self.pos2 - 1, self.pos2)

		if self.pos2_endis5:
			t = self._convert_seq_between_pos(''.join(self.inserted_seq)) + ref
		else:
			t = ref + self._convert_seq_between_pos(''.join(self.inserted_seq))
		bracket = '[' if self.pos1_endis5 else ']'
		alt_matestring = f'{bracket}{self.chrom1}:{self.pos1}{bracket}'
		if self.pos2_endis5:
			alt = alt_matestring + t
		else:
			alt = t + alt_matestring

		return common.Vcfspec(chrom, pos, ref, alt)

	def get_simplesv(self):
		if self.chrom1 != self.chrom2:
			raise Exception(f'Translocation cannot be converted to a SimpleStructuralVariant object.')
		else:
			chrom = self.chrom1
			if self.svtype == 'DEL':
				start1 = self.pos1 + 1
				end1 = self.pos2 - 1
				simplesv = structvars.Deletion(chrom = chrom, start1 = start1, end1 = end1, fasta = self.fasta)
			elif self.svtype == 'INV':
				if self.pos1_endis5:
					start1 = self.pos1
					end1 = self.pos2 - 1
				else:
					start1 = self.pos1 + 1
					end1 = self.pos2
				simplesv = structvars.Inversion(chrom = chrom, start1 = start1, end1 = end1, fasta = self.fasta)
			elif self.svtype == 'DUP':
				start1 = self.pos1
				end1 = self.pos2
				simplesv = structvars.TandemDuplication(chrom = chrom, start1 = start1, end1 = end1, fasta = self.fasta)

			return simplesv

	def get_hgvsg(self):
		return self.get_simplesv().get_hgvsg()

	def get_id_pos1(self):
		return f'{self.chrom1}_{self.pos1}_{self.chrom2}_{self.pos2}_1'

	def get_id_pos2(self):
		return f'{self.chrom1}_{self.pos1}_{self.chrom2}_{self.pos2}_2'

	###########################################

	def get_equivs(self):
		if self.equivs is None:
			self.equivs = get_bnds_equivalents(self)
		return self.equivs

	def get_pos1adv_form(self):
		if self.pos1adv_form is None:
			self.pos1adv_form = self.get_equivs()[0]
		return self.pos1adv_form

	def get_pos2adv_form(self):
		if self.pos2adv_form is None:
			self.pos2adv_form = self.get_equivs()[-1]
		return self.pos2adv_form

	def get_homlen(self):
		if self.homlen is None:
			self.homlen, self.homseq = get_microhomology_spec(self.get_equivs())
		return self.homlen

	def get_homseq(self):
		if self.homseq is None:
			self.homlen, self.homseq = get_microhomology_spec(self.get_equivs())
		return self.homseq
	
	###################################################################

	def copy(self):
		return Breakends(
				fasta = self.fasta,
				chromdict = self.chromdict,
				chrom1 = self.chrom1,
				pos1 = self.pos1,
				pos1_endis5 = self.pos1_endis5,
				chrom2 = self.chrom2,
				pos2 = self.pos2,
				pos2_endis5 = self.pos2_endis5,
				inserted_seq = self.inserted_seq.copy(),
				svtype = self.svtype,
				score = self.score,
				)

	def identical(self, other):
		return hash(self) == hash(other)

	def sameseq1(self, other):
		return get_pos1adv_form(self) == get_pos1adv_form(other)

	def sameseq2(self, other):
		for key in ('chrom1', 'pos1_endis5', 'chrom2', 'pos2_endis5'):
			if getattr(self, key) != getattr(other, key):
				return False

		target = self.copy()
		query = other.copy()

		if target.pos1_endis5:
			goal_pos1 = max(target.pos1, query.pos1)
		else:
			goal_pos1 = min(target.pos1, query.pos1)

		if target.pos2_endis5:
			goal_pos2 = max(target.pos2, query.pos2)
		else:
			goal_pos2 = min(target.pos2, query.pos2)

		for bnds in (target, query):
			while True:
				if bnds.pos1 == goal_pos1:
					break
				else:
					bnds.retract_pos1()
					continue
			while True:
				if bnds.pos2 == goal_pos2:
					break
				else:
					bnds.retract_pos2()
					continue

		return target.equal(query)

	###################################################################

	def advance_pos1(self):
		def get_newpos12(self):
			newpos1 = self._get_advanced_pos1()
			if len(self.inserted_seq) > 0:
				newpos2 = self.pos2
			else:
				newpos2 = self._get_retracted_pos2()

			return newpos1, newpos2

		def base_match_check(self, newpos1):
			newbase = self.fasta.fetch(self.chrom1, newpos1 - 1, newpos1)
			if len(self.inserted_seq) > 0:
				oldbase = self.inserted_seq[-1] if self.pos1_endis5 else self.inserted_seq[0]
			else:
				oldbase = self._convert_seq_between_pos( self.fasta.fetch(self.chrom2, self.pos2 - 1, self.pos2) )

			return (newbase == oldbase)

		def action(self, newpos1, newpos2):
			self.pos1 = newpos1
			self.pos2 = newpos2
			if len(self.inserted_seq) > 0:
				self.score += 1
				if self.pos1_endis5:
					del self.inserted_seq[-1]
				else:
					del self.inserted_seq[0]

		# main
		newpos1, newpos2 = get_newpos12(self)
		if self._newpos_range_check(newpos1, newpos2):
			if base_match_check(self, newpos1):
				action(self, newpos1, newpos2)
				success = True
			else:
				success = False
		else:
			success = False

		return success

	def advance_pos2(self):
		def get_newpos12(self):
			newpos2 = self._get_advanced_pos2()
			if len(self.inserted_seq) > 0:
				newpos1 = self.pos1
			else:
				newpos1 = self._get_retracted_pos1()

			return newpos1, newpos2

		def base_match_check(self, newpos2):
			newbase = self._convert_seq_between_pos( self.fasta.fetch(self.chrom2, newpos2 - 1, newpos2) )
			if len(self.inserted_seq) > 0:
				oldbase = self.inserted_seq[0] if self.pos1_endis5 else self.inserted_seq[-1]
			else:
				oldbase = self.fasta.fetch(self.chrom1, self.pos1 - 1, self.pos1)

			return (newbase == oldbase)

		def action(self, newpos1, newpos2):
			self.pos1 = newpos1
			self.pos2 = newpos2
			if len(self.inserted_seq) > 0:
				self.score += 1
				if self.pos1_endis5:
					del self.inserted_seq[0]
				else:
					del self.inserted_seq[-1]

		# main
		newpos1, newpos2 = get_newpos12(self)
		if self._newpos_range_check(newpos1, newpos2):
			if base_match_check(self, newpos2):
				action(self, newpos1, newpos2)
				success = True
			else:
				success = False
		else:
			success = False

		return success

	def retract_pos1(self):
		added_base = self.fasta.fetch(self.chrom1, self.pos1 - 1, self.pos1)
		if self.pos1_endis5:
			self.inserted_seq.append(added_base)
		else:
			self.inserted_seq.insert(0, added_base)

		self.pos1 = self._get_retracted_pos1()
		self.score -= 1

	def retract_pos2(self):
		added_base = self._convert_seq_between_pos(self.fasta.fetch(self.chrom2, self.pos2 - 1, self.pos2))
		if self.pos1_endis5:
			self.inserted_seq.insert(0, added_base)
		else:
			self.inserted_seq.append(added_base)

		self.pos2 = self._get_retracted_pos2()
		self.score -= 1

	#######################################################

	def _set_svtype(self):
		if self.chrom1 != self.chrom2:
			self.svtype = 'TRA'
		else:
			if (not self.pos1_endis5) and self.pos2_endis5:
				self.svtype = 'DEL'
			elif self.pos1_endis5 and (not self.pos2_endis5):
				self.svtype = 'DUP'
			elif self.pos1_endis5 == self.pos2_endis5:
				self.svtype = 'INV'

	def _newpos_range_check(self, newpos1, newpos2):
		return (
				newpos1 >= 1 and 
				newpos1 <= self.chromdict[self.chrom1] and
				newpos2 >= 1 and 
				newpos2 <= self.chromdict[self.chrom2]
			   )

	def _convert_seq_between_pos(self, seq):
		return convert_seq_between_pos(seq, self.pos1_endis5, self.pos2_endis5)

	def _get_advanced_pos1(self):
		if self.pos1_endis5:
			return self.pos1 - 1
		else:
			return self.pos1 + 1

	def _get_retracted_pos1(self):
		if self.pos1_endis5:
			return self.pos1 + 1
		else:
			return self.pos1 - 1

	def _get_advanced_pos2(self):
		if self.pos2_endis5:
			return self.pos2 - 1
		else:
			return self.pos2 + 1

	def _get_retracted_pos2(self):
		if self.pos2_endis5:
			return self.pos2 + 1
		else:
			return self.pos2 - 1


class NonStandardSVAlt(Exception):
	pass


class NonStandardSVRecord(Exception):
	pass


#######################################################


def convert_seq_between_pos(seq, pos1_endis5, pos2_endis5):
	if pos1_endis5 == pos2_endis5:
		return Bio.Seq.reverse_complement(seq)
	else:
		return seq


def get_bnds_equivalents(bnds):
	"""
	Returns: A list of Breakends objects, all equivalent to the input object, sorted such that the first item is the most advanced form with respect to pos1, and the last item is the most advanced form with respect to pos2.
	"""

	input_copy = bnds.copy()

	bnds_list_pos1adv = list()
	bnds_list_pos1adv.append(input_copy)
	while True:
		new = bnds_list_pos1adv[-1].copy()
		success = new.advance_pos1()
		if success:
			bnds_list_pos1adv.append(new)
			continue
		else:
			break

	bnds_list_pos2adv = list()
	bnds_list_pos2adv.append(input_copy)
	while True:
		new = bnds_list_pos2adv[-1].copy()
		success = new.advance_pos2()
		if success:
			bnds_list_pos2adv.append(new)
			continue
		else:
			break

	bnds_list_pos1adv.reverse()
	bnds_list = bnds_list_pos1adv[:-1] + [input_copy] + bnds_list_pos2adv[1:]
	max_score = max(x.score for x in bnds_list)

	bnds_equivs = [ x for x in bnds_list if x.score == max_score ]

	return bnds_equivs


def get_pos1adv_form(bnds):
	return get_bnds_equivalents(bnds)[0]


def get_pos2adv_form(bnds):
	return get_bnds_equivalents(bnds)[-1]


def get_microhomology_spec(bnds_equivs):
	homlen = len(bnds_equivs) - 1
	if homlen == 0:
		homseq = ''
	else:
		bnds = bnds_equivs[0]
		if bnds.pos1_endis5:
			pos_list = [x.pos1 for x in sorted(bnds_equivs, key = lambda x: x.pos1)[:-1]]
		else:
			pos_list = [x.pos1 for x in sorted(bnds_equivs, key = lambda x: x.pos1)[1:]]
		homseq = bnds.fasta.fetch(bnds.chrom1, pos_list[0]-1, pos_list[-1])

	return homlen, homseq


#######################################################


def get_bnds_from_vr(vr, fasta, chromdict):
	"""
	Raises:
		If input vr does not conform to a known SV variant record format (including a valid non-SV variant record)
	"""

	assert len(vr.alts) == 1, f'Multiallelic variant record is not allowed:\n{vr}'

	vr_svinfo = get_vr_svinfo_standard_vr(vr, fasta, chromdict)
	bnds = get_bnds_from_vr_svinfo(vr, vr_svinfo, fasta, chromdict)

	#return bnds, vr_svinfo
	return bnds


########################################################


def get_bnds_from_vr_svinfo(vr, vr_svinfo, fasta, chromdict):
	"""
	inserted_seq: As seen in the viewpoint where pos1-side sequence is on the plus(Crick) strand
	"""

	def warn():
		LOGGER.warning(f'"t" portion of SV ALT string is not an extension of REF string for this variant record:\n{vr}')

	if vr_svinfo['is_bnd1']:
		chrom1 = vr.contig
		chrom2 = vr_svinfo['chrom_mate']
		pos1_endis5 = vr_svinfo['current_endis5']
		pos2_endis5 = vr_svinfo['mate_endis5']

		if pos1_endis5:
			if vr_svinfo['t'][-1] == vr_svinfo['ref']:
				pos1 = vr.pos
				inserted_seq = list( vr_svinfo['t'][:-1] )
			else:
				warn()
				pos1 = vr.pos + 1
				inserted_seq = list( vr_svinfo['t'] )
		else:
			if vr_svinfo['t'][0] == vr_svinfo['ref']:
				pos1 = vr.pos
				inserted_seq = list( vr_svinfo['t'][1:] )
			else:
				warn()
				pos1 = vr.pos - 1
				inserted_seq = list( vr_svinfo['t'] )
			
		pos2 = vr_svinfo['pos_mate']

	else:
		chrom1 = vr_svinfo['chrom_mate']
		chrom2 = vr.contig
		pos1_endis5 = vr_svinfo['mate_endis5']
		pos2_endis5 = vr_svinfo['current_endis5']

		pos1 = vr_svinfo['pos_mate']

		if pos2_endis5:
			if vr_svinfo['t'][-1] == vr_svinfo['ref']:
				pos2 = vr.pos
				inserted_seq = list(convert_seq_between_pos(
							vr_svinfo['t'][:-1], pos1_endis5, pos2_endis5
							))
			else:
				warn()
				pos2 = vr.pos + 1
				inserted_seq = list(convert_seq_between_pos(
							vr_svinfo['t'], pos1_endis5, pos2_endis5
							))
		else:
			if vr_svinfo['t'][0] == vr_svinfo['ref']:
				pos2 = vr.pos
				inserted_seq = list(convert_seq_between_pos(
							vr_svinfo['t'][1:], pos1_endis5, pos2_endis5
							))
			else:
				warn()
				pos2 = vr.pos - 1
				inserted_seq = list(convert_seq_between_pos(
							vr_svinfo['t'], pos1_endis5, pos2_endis5
							))

	bnds = Breakends(
			chrom1 = chrom1,
			pos1 = pos1,
			pos1_endis5 = pos1_endis5,
			chrom2 = chrom2,
			pos2 = pos2,
			pos2_endis5 = pos2_endis5,
			inserted_seq = inserted_seq,
			fasta = fasta,
			chromdict = chromdict,
			)

	return bnds


def get_vr_svinfo_standard_vr(vr, fasta, chromdict):
	#assert len(vr.alts) == 1, f'Multiallelic variant record is not allowed:\n{vr}'

	#if vr.ref == 'N':
	#	raise NonStandardSVRecord(f'Input variant record REF is "N":\n{vr}')

	vr_svinfo = dict()
	vr_svinfo['ref'] = vr.ref

	try:
		vr_svinfo['t'], \
		vr_svinfo['chrom_mate'], \
		vr_svinfo['pos_mate'], \
		vr_svinfo['current_endis5'], \
		vr_svinfo['mate_endis5'] = parse_sv_altstring(vr.alts[0])
	except NonStandardSVAlt as e:
		e_msg = f'{str(e)}\nInput variant record:\n{vr}'
		raise NonStandardSVRecord(e_msg)

	vr_svinfo['is_bnd1'] = get_is_bnd1(vr, vr_svinfo, chromdict)

	return vr_svinfo


def parse_sv_altstring(sv_altstring):
	mats = [ common.RE_PATS['alt_bndstring_1'].match(sv_altstring), common.RE_PATS['alt_bndstring_2'].match(sv_altstring) ]
	mats_isNotNone = [ (x is not None) for x in mats ]

	nTrue = mats_isNotNone.count(True)
	if nTrue == 0: # not a bnd string
		raise NonStandardSVAlt(f'ALT string "{sv_altstring}" does not match the standard SV string pattern.') 

	elif nTrue == 1:
		mat = next(itertools.compress(mats, mats_isNotNone))
		t = mat.group('t') # t : according to VCF spec documentation
		chrom_mate = mat.group('matechrom')
		pos_mate = int(mat.group('matepos'))

		if sv_altstring.startswith('[') or sv_altstring.startswith(']'):
			endtype_current_is5 = True
		else:
			endtype_current_is5 = False
		
		if mat.group('bracket1') == '[':
			endtype_mate_is5 = True
		else:
			endtype_mate_is5 = False

		sv_altstring_parsed = (t, chrom_mate, pos_mate, endtype_current_is5, endtype_mate_is5)

	elif nTrue == 2: # not a valid bnd string
		raise NonStandardSVAlt(f'ALT string "{sv_altstring}" matches both pat1 and pat2.') 

	return sv_altstring_parsed


def get_is_bnd1(vr, vr_svinfo, chromdict):
	order = common.get_order(vr.contig, vr.pos, vr_svinfo['chrom_mate'], vr_svinfo['pos_mate'], chromdict)
	if order < 0:
		is_bnd1 = True
	elif order > 0:
		is_bnd1 = False
	elif order == 0:
		raise Exception(f'Current chrom/pos and mate chrom/pos are identical for vcf record :\n{vr}')

	return is_bnd1


