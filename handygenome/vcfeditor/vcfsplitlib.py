import itertools 

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))


def get_lineno_vcf_path(vcf_path):
	with pysam.VariantFile(vcf_path, 'r') as vcf:
		lineno = get_lineno_vcf(vcf)
	return lineno


def get_lineno_vcf(vcf):
	lineno = 0
	for pysamvr in vcf.fetch():
		lineno += 1
	return lineno


def get_output_lineno_list(total_lineno, n_file = None, n_line = None):
	"""
	Args:
		total_lineno: Total number of variant records of input vcf.
		n_file: Number of output files. \
If greater than the total line number, reduced to the total line number.
		n_line: Number of variant records per one output file. \
If greater than the total line number, reduced to the total line number.
		* One and only one of n_file or n_line argument must be set.

	Returns:
		A list composed of the number of variant records per output file.
	"""
	
	common.check_num_None(1, (n_file, n_line), ('n_file', 'n_line'))
	assert total_lineno > 0, f'"total_lineno" must be a positive integer.'

	if n_file is not None:
		if n_file > total_lineno:
			warnings.warn(f'"n_file" is greater than "total_lineno". It will be changed to be the same with "total_lineno".')
		n_file = min(n_file, total_lineno)
		q, r = divmod(total_lineno, n_file)
		result = list(itertools.repeat(q+1, r)) + list(itertools.repeat(q, n_file - r))

	elif n_line is not None:
		if n_line > total_lineno:
			warnings.warn(f'"n_line" is greater than "total_lineno". It will be changed to be the same with "total_lineno".')
		n_line = min(n_line, total_lineno)
		q, r = divmod(total_lineno, n_line)
		result = list(itertools.repeat(n_line, q)) + list(itertools.repeat(r, 1))

	return result
