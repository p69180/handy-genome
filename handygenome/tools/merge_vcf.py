import re
import os
import argparse
import textwrap

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
workflow = importlib.import_module('.'.join([top_package_name, 'workflow']))
merge_module = importlib.import_module('.'.join([top_package_name, 'vcfeditor', 'merge']))


def argument_parser(cmdargs):
    def sanity_check(args):
        e_msg = ('Invalid --isec-indices value. Please refer to the help '
                 'message.')
        if args.isec_indices is not None:
            if len(args.isec_indices) != len(args.infile_path_list):
                raise Exception(e_msg)
            elif not set(isec_indices).issubset({'0', '1'}):
                raise Exception(e_msg)

    def postprocess(args):
        if args.isec_indices is not None:
            args.isec_indices = [int(x) for x in args.isec_indices]

    parser_dict = workflow.init_parser()

    workflow.add_infilelist_arg(parser_dict['required'], required=True)
    workflow.add_outfile_arg(parser_dict['required'], required=True, 
                             must_not_exist='ask')
    workflow.add_fasta_arg(parser_dict['required'], required=True)
    parser_dict['required'].add_argument(
        '--mode', required=True, choices=('isec', 'union'),
        help=f'Must be "isec" (which means intersection) or "union".')

    parser_dict['optional'].add_argument(
        '--isec-indices', dest='isec_indices', required=False,
        help=('(Only applied for intersection) A string composed of 0 or 1, ' 
              'with the same length as the number of input files. Files '
              'marked with 0 are excluded and those with 1 are included.'))
    workflow.add_logging_args(parser_dict)
    workflow.add_outfmt_arg(parser_dict['optional'], required=False)

    args = parser_dict['main'].parse_args(cmdargs)
    sanity_check(args)
    postprocess(args)

    return args


def main(cmdargs):
    args = argument_parser(cmdargs)
    logger = workflow.get_logger(name='merge_vcf', 
                                 stderr=(not args.silent),
                                 filename=args.log, append=False)

    merge_module.main_file(
        infile_path_list=args.infile_path_list,
        outfile_path=args.outfile_path,
        fasta_path=args.fasta_path,
        isec=(args.mode == 'isec'),
        isec_indices=args.isec_indices, 
        mode_pysam=args.mode_pysam,
        outfile_must_not_exist='no',
        logger=logger,
        )
