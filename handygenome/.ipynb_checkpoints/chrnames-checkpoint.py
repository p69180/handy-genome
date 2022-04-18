import itertools
import pprint

import collections

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))


URL_GRCH38 = 'https://rest.ensembl.org/info/assembly/homo_sapiens?synonyms=1'
URL_GRCH37 = 'https://grch37.rest.ensembl.org/info/assembly/homo_sapiens?synonyms=1'

INVALID_NAMES = [ 'NT_187507.1' ]


class ChromNameMap:
	def __init__(self, data):
		"""
		Args:
			data: { 'ucsc' : ['chr1', 'chr2', ... ], 'refseq' : [ 'NC_000001.11', 'NC_000002.12', ... ], ... }
		"""

		self.data = data
		self.dicts = dict()
		for key1, key2 in itertools.permutations(data.keys(), 2):
			dic = dict()
			for val1, val2 in zip(self.data[key1], self.data[key2]):
				if val1 is None:
					continue
				dic[val1] = val2

			self.dicts[(key1,key2)] = dic

	def __repr__(self):

		tmp_keys = list()
		tmp_vals = list()
		for key, val in self.data.items():
			tmp_keys.append(key)
			tmp_vals.append(val)

		show_result = list()
		for tup in zip(*tmp_vals):
			show_result.append(dict(zip(tmp_keys, tup)))

		return pprint.pformat(show_result)


def parse_raw_result(raw_result):
	def handle_na(val):
		if val == 'na':
			return None
		else:
			return val

	def check_exception(curr_ucsc, curr_refseq, curr_genbank):
		return (
				len(curr_ucsc) > 1 or 
				len(curr_refseq) > 1 or
				len(curr_genbank) > 1
			   )

	data = {
		'ensembl_basic' : list(),
		'ucsc' : list(),
		'genbank' : list(),
		'refseq' : list(),
	}

	for dic in raw_result['top_level_region']:
		assert dic['name'] != 'na'

		data['ensembl_basic'].append(dic['name'])

		curr_ucsc = list()
		curr_refseq = list()
		curr_genbank = list()

		for subdic in dic['synonyms']:
			if subdic['dbname'] == 'UCSC':
				curr_ucsc.append(subdic['name'])

			elif subdic['dbname'] == 'RefSeq_genomic':
				curr_refseq.append(subdic['name'])

			elif subdic['dbname'] == 'INSDC':
				curr_genbank.append(subdic['name'])
			else:
				raise Exception(f'Unexpected subdic dbname: {subdic["dbname"]}')

		# exception handling
		if check_exception(curr_ucsc, curr_refseq, curr_genbank):
			if (
					dic['name'] == 'KI270752.1' and
					curr_ucsc == [ 'chrUn_KI270752v1' ] and
					set(curr_refseq) == { 'NT_187507.1', 'na' } and
					len(curr_genbank) == 0
			   ):
				curr_refseq = [ None ]
				curr_genbank = [ dic['name'] ]
			else:
				raise Exception(f'Unexpected exception: {dic}')

		# non-exception
		else:
			if len(curr_ucsc) == 0:
				curr_ucsc.append(None)
			if len(curr_refseq) == 0:
				curr_refseq.append(None)
			if len(curr_genbank) == 0:
				curr_genbank.append(dic['name'])

		data['ucsc'].extend(curr_ucsc)
		data['refseq'].extend(curr_refseq)
		data['genbank'].extend(curr_genbank)

	return data


##############################################################################


NAMEMAP_GRCH37 = ChromNameMap( parse_raw_result( common.get_url_contents(URL_GRCH37) ) )
NAMEMAP_GRCH38 = ChromNameMap( parse_raw_result( common.get_url_contents(URL_GRCH38) ) )
