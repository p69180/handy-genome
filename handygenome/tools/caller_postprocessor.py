import re
import os

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
workflow = importlib.import_module('.'.join([top_package_name, 'workflow']))
svcaller_parser = importlib.import_module('.'.join([top_package_name, 'variantplus', 'svcaller_parser']))
equivalents = importlib.import_module('.'.join([top_package_name, 'variantplus', 'equivalents']))
varianthandler = importlib.import_module('.'.join([top_package_name, 'variantplus', 'varianthandler']))
headerhandler = importlib.import_module('.'.join([top_package_name, 'variantplus', 'headerhandler']))
index = importlib.import_module('.'.join([top_package_name, 'vcfeditor', 'index']))


def argument_parser(cmdargs):
    def sanity_check(args):
        pass

    parser_dict = workflow.init_parser()

    workflow.add_infile_arg(parser_dict['required'], required=True)
    workflow.add_outfile_arg(parser_dict['required'], required=True, 
                             must_not_exist='ask')
    workflow.add_fasta_arg(parser_dict['required'], required=True)
    workflow.add_outfmt_arg(parser_dict['optional'], required=False)

    workflow.add_logging_args(parser_dict)

    args = parser_dict['main'].parse_args(cmdargs)
    sanity_check(args)

    return args


def process_nonsv(vr, out_vr_list, nonsv_set, fasta, logger):
    vcfspec_original = common.Vcfspec(vr.contig, vr.pos, vr.ref, vr.alts[0])
    vcfspec_left = equivalents.leftmost(vcfspec_original, fasta)
    ID = vcfspec_left.get_id()
    if ID in nonsv_set:
        logger.info(f'Input variant record is discarded since vcfspec is '
                    f'duplicated:\n{vr}')
        return
    else:
        nonsv_set.add(ID)
        vr.pos = vcfspec_left.pos
        vr.alleles = ( vcfspec_left.ref, vcfspec_left.alt )
        vr.id = ID
        out_vr_list.append(vr)


def process_sv(vr, out_vr_list, bnds_set, fasta, chromdict, logger):
    # pos1adv form
    bnds = svcaller_parser.get_bnds_from_caller_vr(vr, fasta, chromdict) 
    if bnds is None:
        raise Exception(f'Error while parsing this SV variant record:\n{vr}')

    if bnds in bnds_set:
        logger.info(f'Input variant record is discarded since breakends '
                    f'spec is duplicated:\n{vr}')
        return
    else:
        bnds_set.add(bnds)
        vcfspec_pos1 = bnds.get_vcfspec_pos1()
        vcfspec_pos2 = bnds.get_vcfspec_pos2()
        ID_pos1 = bnds.get_id_pos1()
        ID_pos2 = bnds.get_id_pos2()

        vr.contig = vcfspec_pos1.chrom
        vr.pos = vcfspec_pos1.pos
        vr.alleles = (vcfspec_pos1.ref, vcfspec_pos1.alt)
        vr.id = ID_pos1
        vr.info['MATEID'] = ID_pos2

        vr_pos2 = vr.header.new_record()
        vr_pos2.contig = vcfspec_pos2.chrom
        vr_pos2.pos = vcfspec_pos2.pos
        vr_pos2.alleles = (vcfspec_pos2.ref, vcfspec_pos2.alt)
        vr_pos2.id = ID_pos2
        vr_pos2.info['MATEID'] = ID_pos1

        out_vr_list.append(vr)
        out_vr_list.append(vr_pos2)


def process_monoallelic_vr(vr, out_vr_list, nonsv_set, bnds_set, fasta, 
                           chromdict, logger):
    # skip delly <INS> records
    if svcaller_parser.check_vr_format_delly(vr):
        if vr.alts[0] == '<INS>':
            return

    if varianthandler.check_SV(vr):
        process_sv(vr, out_vr_list, bnds_set, fasta, chromdict, logger)
    else:
        process_nonsv(vr, out_vr_list, nonsv_set, fasta, logger)


def main(cmdargs):
    args = argument_parser(cmdargs)
    logger = workflow.get_logger(
            name = 'caller_postprocessor',
            stderr = (not args.silent),
            filename = args.log,
            append = False,
            )
    logger.info('BEGINNING')

    fasta = pysam.FastaFile(args.fasta_path)
    chromdict = common.ChromDict(fasta = fasta)

    nonsv_set = set()
    bnds_set = set()
    out_vr_list = list()
    
    with pysam.VariantFile(args.infile_path) as in_vcf:
        headerhandler.add_MATEID(in_vcf.header)
        with pysam.VariantFile(args.outfile_path, mode = args.mode_pysam, 
                               header = in_vcf.header) as out_vcf:
            for vr in in_vcf.fetch():
                if len(vr.alts) == 1:
                    process_monoallelic_vr(vr, out_vr_list, nonsv_set, 
                                           bnds_set, fasta, chromdict, logger)
                else:
                    for alt in vr.alts:
                        new_vr = vr.copy()
                        new_vr.alts = [ alt ]
                        process_monoallelic_vr(new_vr, out_vr_list, nonsv_set, 
                                               bnds_set, fasta, chromdict, 
                                               logger)

            out_vr_list.sort(key=common.get_vr_sortkey(chromdict))

            for vr in out_vr_list:
                out_vcf.write(vr)

    index.index_vcf(args.outfile_path, force = True)

    logger.info('ALL SUCCESSFULLY FINISHED')

