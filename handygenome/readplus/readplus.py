import importlib
top_package_name = __name__.split('.')[0] 

readhandler = importlib.import_module(top_package_name + '.readplus.readhandler')


class ReadPlus:
	def __init__(self, read, fasta = None):
		'''
		These attributes are set regardless of whether ria functions are used.
		'''
		self.read = read
		if fasta != None:
			set_NMMD(self, fasta)
		set_basic_attributes(self)
		set_irrelevant(self)
			

def set_basic_attributes(rp):
	rp.pairs_dict = readhandler.get_pairs_dict(rp.read)
	rp.fiveprime_end = readhandler.get_fiveprime_end(rp.read)
	rp.threeprime_end = readhandler.get_threeprime_end(rp.read)
	rp.cigarstats = rp.read.get_cigar_stats()[0] 
		# length 11 array which contains value for M, I, D, N, S, ...
	rp.margins = ( rp.read.reference_start, rp.read.reference_end )
	rp.MQ = rp.read.mapping_quality


def set_NMMD(rp, fasta):
	if (not rp.read.has_tag('NM')) or (not rp.read.has_tag('MD')):
		ref_seq_padded, read_seq_padded = readhandler.get_padded_seqs(rp.read, fasta)
		if not rp.read.has_tag('NM'):
			rp.read.set_tag(
					tag = 'NM', value_type = 'i',
					value = readhandler.get_NM(
						rp.read, 
						ref_seq_padded = ref_seq_padded, 
						read_seq_padded = read_seq_padded,
						), 
					)
		if not rp.read.has_tag('MD'):
			rp.read.set_tag(
					tag = 'MD', value_type = 'Z',
					value = readhandler.get_MD(
						rp.read, 
						ref_seq_padded = ref_seq_padded, 
						read_seq_padded = read_seq_padded,
						), 
					)


def set_irrelevant(rp):
	rp.irrelevant_cause = list()
	rp.irrelevant = False
