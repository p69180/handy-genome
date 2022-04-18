import re
from pprint import pprint

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
infoformat = importlib.import_module('.'.join([top_package_name, 'variantplus', 'infoformat']))
ensembl_rest = importlib.import_module('.'.join([top_package_name, 'annotation', 'ensembl_rest']))
ensembl_parser = importlib.import_module('.'.join([top_package_name, 'annotation', 'ensembl_parser']))


def rest_lookup_transcript(vr, ID, refver = None):
	if refver is None:
		refver = common.infer_refver_vr(vr)
	assert refver in ('hg19', 'hg38')
	hg19 = (refver == 'hg19')

	raw_result = ensembl_rest.lookup_id(ID, hg19 = hg19, expand = True)
	ensembl_result = \
		ensembl_parser.parse_rest_lookup_transcript(raw_result, hg19 = hg19, set_gene_name = True)
	


	
		

