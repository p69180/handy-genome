import re
import os
import argparse

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
workflow = importlib.import_module('.'.join([top_package_name, 'workflow']))
concat_module = importlib.import_module('.'.join([top_package_name, 'vcfeditor', 'concat']))


def argument_parser(cmdargs):
	def sanity_check(args):
		pass

	parser_dict = workflow.init_parser()
	workflow.add_infilelist_arg(parser_dict['required'], required = True)
	workflow.add_outfile_arg(parser_dict['required'], required = True, must_not_exist = 'ask')
	workflow.add_outfmt_arg(parser_dict['optional'], required = False)

	args = parser_dict['main'].parse_args(cmdargs)
	sanity_check(args)

	return args


def main(cmdargs):
	args = argument_parser(cmdargs)
	concat_module.main(
			infile_path_list = args.infile_path_list,
			outfile_path = args.outfile_path,
			mode_pysam = args.mode_pysam, 
			outfile_must_not_exist = 'no',
			)
