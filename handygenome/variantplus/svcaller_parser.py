# svcparser: SV caller parser
import re
#import warnings

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
breakends_module = importlib.import_module('.'.join([top_package_name, 'variantplus', 'breakends']))
headerhandler = importlib.import_module('.'.join([top_package_name, 'variantplus', 'headerhandler']))


INFO_END_PAT = re.compile('^(.*;)?END=([^;]+);.*$')


def write_standard_vcf(
		infile_path, 
		outfile_path, 
		fasta,
		mode_bcftools = 'z', 
		mode_pysam = None, 
		):
	def sanitycheck(infile_path, outfile_path):
		common.check_infile_validity(infile_path)
		common.check_outfile_validity(outfile_path, must_not_exist = False)

	sanitycheck(infile_path, outfile_path)
	mode_pysam = common.write_mode_arghandler(mode_bcftools, mode_pysam)
	chromdict = common.ChromDict(fasta = fasta)

	in_vcf = pysam.VariantFile(infile_path)
	headerhandler.add_MATEID(in_vcf.header)
	out_vcf = pysam.VariantFile(outfile_path, mode = mode_pysam, header = in_vcf.header)

	# set vrdict
	input_vrdict = dict()
	for vr in in_vcf.fetch():
		bnds = get_bnds_from_caller_vr(vr, fasta, chromdict)
		if bnds is None:
			continue
		else:
			if bnds in input_vrdict:
				continue
			else:
				input_vrdict[bnds] = vr

	# set output vr list
	output_vrlist = list()
	ID_prefix_set = set()
	for bnds, vr in input_vrdict.items():
		vcfspec_pos1 = bnds.get_vcfspec_pos1()
		vcfspec_pos2 = bnds.get_vcfspec_pos2()

		ID_prefix = '_'.join([
				vcfspec_pos1.chrom,
				str(vcfspec_pos1.pos),
				vcfspec_pos2.chrom,
				str(vcfspec_pos2.pos),
				])
		assert ID_prefix not in ID_prefix_set
		ID_prefix_set.add(ID_prefix)
				
		ID_pos1 = ID_prefix + '_1'
		ID_pos2 = ID_prefix + '_2'

		vr_pos1 = vr.copy()
		vr_pos1.contig = vcfspec_pos1.chrom
		vr_pos1.pos = vcfspec_pos1.pos
		vr_pos1.alleles = ( vcfspec_pos1.ref, vcfspec_pos1.alt )
		vr_pos1.id = ID_pos1
		vr_pos1.info['MATEID'] = ID_pos2

		vr_pos2 = in_vcf.header.new_record()
		vr_pos2.contig = vcfspec_pos2.chrom
		vr_pos2.pos = vcfspec_pos2.pos
		vr_pos2.alleles = ( vcfspec_pos2.ref, vcfspec_pos2.alt )
		vr_pos2.id = ID_pos2
		vr_pos2.info['MATEID'] = ID_pos1

		output_vrlist.append(vr_pos1)
		output_vrlist.append(vr_pos2)

	# sort output vrs according to chrom1/pos1
	output_vrlist.sort(
			key = lambda vr: common.coord_sortkey(vr.contig, vr.pos, chromdict)
			)

	for vr in output_vrlist:
		out_vcf.write(vr)

	in_vcf.close()
	out_vcf.close()


#######################################################


def get_info_end(vr):
	candidates = list()

	infosplit = str(vr).split('\t')[7].strip().split(';')
	for item in infosplit:
		itemsp = item.split('=')
		if itemsp[0] == 'END':
			candidates.append(int(itemsp[1]))
			break

	if len(candidates) == 0:
		raise Exception(f'INFO does not include END in the input variant record:\n{vr}')

	return candidates[0]


#######################################################


def get_bnds_from_caller_vr(vr, fasta, chromdict):
	"""
	Raises:
		All raises are from get_vr_svinfo_caller_vr function.
	"""

	assert len(vr.alts) == 1, f'Multiallelic variant record is not allowed:\n{vr}'

	vr_svinfo = get_vr_svinfo_caller_vr(vr, fasta, chromdict) # this function does not return None ; it rather raises an Exception
	if vr_svinfo is None:
		return None
	else:
		bnds = breakends_module.get_bnds_from_vr_svinfo(vr, vr_svinfo, fasta, chromdict)
		bnds_pos1_adv = breakends_module.get_bnds_equivalents(bnds)[0]

		#return bnds, vr_svinfo
		return bnds_pos1_adv


def get_vr_svinfo_caller_vr(vr, fasta, chromdict):
	"""
	Raises:
		If input vr does not conform to a known SV variant record format (including a valid non-SV variant record)
	"""

	if check_vr_format_delly(vr):
		vr_svinfo = get_vr_svinfo_delly(vr, fasta, chromdict)

	elif check_vr_format_manta(vr):
		vr_svinfo = get_vr_svinfo_manta(vr, fasta, chromdict)

	else: # Assumed that other SV callers follow the SV vcf standard
		vr_svinfo = breakends_module.get_vr_svinfo_standard_vr(vr, fasta, chromdict)
		modify_vr_svinfo_for_dRanger(vr_svinfo, fasta)

	return vr_svinfo


########################################################


def modify_vr_svinfo_for_dRanger(vr_svinfo, fasta):
	if vr_svinfo['chrom_mate'] not in fasta.references:
		assert vr_svinfo['chrom_mate'].startswith('chr')
		vr_svinfo['chrom_mate'] = re.sub('^chr', '', vr_svinfo['chrom_mate'])


def get_vr_svinfo_manta(vr, fasta, chromdict):
	if vr.alts[0].startswith('<') and vr.alts[0].endswith('>'): # non-TRA
		vr_svinfo = dict()

		vr_svinfo['chrom_mate'] = vr.contig
		vr_svinfo['pos_mate'] = get_info_end(vr)

		assert vr_svinfo['pos_mate'] > vr.pos, f'INFO/END is smaller than POS for this presumed MANTA non-BND variant record:\n{vr}'
		assert vr.alts[0] in ('<DUP:TANDEM>', '<DEL>'), f'Unexpected ALT string for this presumed MANTA non-BND variant record:\n{vr}'

		if vr.alts[0] == '<DUP:TANDEM>':
			vr_svinfo['current_endis5'] = True
			vr_svinfo['mate_endis5'] = False
		elif vr.alts[0] == '<DEL>':
			vr_svinfo['current_endis5'] = False
			vr_svinfo['mate_endis5'] = True

		vr_svinfo['ref'] = fasta.fetch(vr.contig, vr.pos - 1, vr.pos)
		vr_svinfo['t'] = vr_svinfo['ref']

		vr_svinfo['is_bnd1'] = True
	
	else: # standard bnd string representation is assumed
		vr_svinfo = breakends_module.get_vr_svinfo_standard_vr(vr, fasta, chromdict)

	return vr_svinfo


def get_vr_svinfo_delly(vr, fasta, chromdict):
	#assert vr.ref == 'N', f'REF is not "N" for the input delly variant record:\n{vr}'

	if vr.alts[0] == '<INS>':
		#warnings.warn(f'ALT is <INS> for input delly variant record.')
		return None
	else:
		vr_svinfo = dict()

		vr_svinfo['chrom_mate'] = vr.info['CHR2']
		#vr_svinfo['pos_mate'] = vr.stop 
		vr_svinfo['pos_mate'] = get_info_end(vr)
			# delly INFO/END value seems to be 1-based
			# 220323: vr.stop does not return INFO/END value for some delly variant records

		CTsplit = vr.info['CT'].split('to')
		vr_svinfo['current_endis5'] = (CTsplit[0] == '5')
		vr_svinfo['mate_endis5'] = (CTsplit[1] == '5')

		vr_svinfo['ref'] = fasta.fetch(vr.contig, vr.pos - 1, vr.pos)
		vr_svinfo['t'] = vr_svinfo['ref']

		vr_svinfo['is_bnd1'] = breakends_module.get_is_bnd1(vr, vr_svinfo, chromdict)

		return vr_svinfo


######################################################################################


def check_vr_format_delly(vr):
	return (
			('CHR2' in vr.info) and 
			('CHR1' not in vr.info) and 
			('CT' in vr.info) and 
			(vr.ref == 'N')
		   )


def check_vr_format_manta(vr):
	# For BND records, ALT is VCF breakends string
	# For non-BND, ALT is like <DUP:TANDEM>
	# ID starts with 'MANTA'

	if ( vr.id is not None and vr.id.startswith('Manta') ):
		return True
	else:
		return False


def check_vr_format_brass(vr):
	existing_infokey_set = set(('TRDS'))
	if existing_infokey_set.issubset(set(vr.info.keys())):
		return True
	else:
		return False


def check_vr_format_dRanger(vr):
	existing_infokey_set = set(('SVCLASS', 'DRQUAL'))
	if existing_infokey_set.issubset(set(vr.info.keys())):
		return True
	else:
		return False


def check_vr_format_snowman(vr):
	existing_infokey_set = set(('EVDNC', 'SCTG'))
	if existing_infokey_set.issubset(set(vr.info.keys())):
		return True
	else:
		return False
