import itertools
import importlib

top_package_name = __name__.split('.')[0]
readhandler = importlib.import_module('.'.join([top_package_name, 'readplus', 'readhandler']))
ReadPlus = importlib.import_module('.'.join([top_package_name, 'readplus', 'readplus'])).ReadPlus
ReadPlusPair = importlib.import_module('.'.join([top_package_name, 'readplus', 'readpluspair'])).ReadPlusPair


class Rpplist:
	def __init__(self, bam, fasta, CHROM, POS0, chromdict, fetch_width = 100, searchrange_width = 0):
		"""
		read filtering done by readhandler.filter_bad_read
		mate-unmapped reads are filtered out
		"""

		self.chromdict = chromdict

		# qname_list, mate_search_range
		qname_list, mate_search_range = get_qnamelist_matesearchrange(bam, CHROM, POS0, fetch_width, searchrange_width)
			# both None if no reads found
		if qname_list is None:
			self.qname_list = None
			self.mate_search_range = None
			self.rpplist = None
			self.start = None
			self.end = None

		else:
			self.qname_list = qname_list
			self.mate_search_range = mate_search_range
	
			# fetch
			fetchresult_dict = get_fetchresult_dict(bam, CHROM, mate_search_range, qname_list)
	
			# assemble into rpplist
			self.rpplist = list()
			for readlist in fetchresult_dict.values():
				rpp = get_readpluspair(readlist, bam, fasta, chromdict)
				if rpp is not None:
					self.rpplist.append(rpp)
	
			# set start/end
			self.start = min([rpp.rp1.read.reference_start for rpp in self.rpplist])
			self.end = max([rpp.rp2.read.reference_end for rpp in self.rpplist])


def append_start_end_lists(read, start_list, end_list):
	if not readhandler.check_long_template(read):
		template_range = readhandler.get_template_range(read) 
		# TRA reads are removed here since template_length is zero and template_range is None
		if template_range is not None:
			start_list.append(template_range.start) # template_range is a range object
			end_list.append(template_range.stop)


def get_qnamelist_matesearchrange(bam, CHROM, POS0, fetch_width, searchrange_width):
	"""
	Read filtering by readhandler.filter_bad_read
	Additionally excluded reads (different from readhandler.filter_bad_read):
	- too long TLEN
	- TLEN == 0
	- TRA reads
	"""

	qname_list = list()
	start_list = list()
	end_list = list()

	# fetch start/end is padded with "width" to include reads which overlap POS0 with soft-clipped bases on IGV
	for read in readhandler.get_fetch(
			bam, CHROM, 
			start = POS0 - fetch_width,
			end = POS0 + fetch_width,
			filter_fun = readhandler.filter_bad_read,
			):
		#softclip_ends_range = readhandler.get_softclip_ends_range(read)
		#if POS0 not in softclip_ends_range:
		#	continue

		qname_list.append(read.query_name)
		append_start_end_lists(read, start_list, end_list)

	if len(qname_list) == 0:
		qname_list = None
		mate_search_range = None
	else:
		qname_list = list(set(qname_list))
		mate_search_range = range(min(start_list) - searchrange_width, max(end_list) + searchrange_width)

	return qname_list, mate_search_range


def get_matesearchrange(bam, CHROM, POS0):
	'''
Read filtering by filter_bad_read
Additionally excluded reads (different from readhandler.filter_bad_read):
- too long TLEN
- TLEN == 0
- TRA reads
'''
	qname_list = list()
	start_list = list()
	end_list = list()

	for read in readhandler.get_fetch(
			bam, CHROM, 
			POS0 = POS0, 
			filter_fun = readhandler.filter_bad_read,
			):
		qname_list.append(read.query_name)
		append_start_end_lists(read, start_list, end_list)

	qname_list = list(set(qname_list))
	mate_search_range = range(min(start_list), max(end_list))

	return qname_list, mate_search_range


def get_fetchresult_dict(bam, CHROM, mate_search_range, qname_list):
	'''
	Read filtering by filter_bad_read
	'''
	fetchresult_dict = dict()

	for qname in qname_list:
		fetchresult_dict[qname] = list()

	for read in readhandler.get_fetch(
			bam, CHROM, 
			start = mate_search_range.start, end = mate_search_range.stop,
			filter_fun = readhandler.filter_bad_read
			):
		try:
			fetchresult_dict[read.query_name].append(read)
		except KeyError:
			pass

	return fetchresult_dict


def split_readlist(readlist):
	is_primary = list(map(readhandler.check_primary_alignment, readlist))
	is_not_primary = [not x for x in is_primary]
	readlist_primary = list(itertools.compress(readlist, is_primary))
	readlist_nonprimary = list(itertools.compress(readlist, is_not_primary))

	return readlist_primary, readlist_nonprimary


def get_readpluspair(readlist, bam, fasta, chromdict):
	readlist_primary, readlist_nonprimary = split_readlist(readlist)

	if len(readlist_primary) > 2:
		e_message = '\n'.join( ['More than two primary reads found:'] + [ read.to_string() for read in readlist_primary ] )
		raise Exception(e_message)
	elif len(readlist_primary) == 0: # found reads are all non-primary
		rpp = None
	else:
		rplist_nonprimary = [ ReadPlus(x, fasta) for x in readlist_nonprimary ]

		if len(readlist_primary) == 1: # mate read is somewhere far away
			mate = readhandler.get_primary_mate(readlist_primary[0], bam)
			if mate is None:
				raise Exception('Mate read not found.')
			rplist_primary = [ ReadPlus(readlist_primary[0], fasta), ReadPlus(mate, fasta) ]
		else: # len(readlist_primary) == 2
			rplist_primary = [ ReadPlus(x, fasta) for x in readlist_primary ]

		rpp = ReadPlusPair(rplist_primary, rplist_nonprimary, chromdict)

	return rpp
