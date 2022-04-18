import os
import itertools
import warnings
import textwrap

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
workflow = importlib.import_module('.'.join([top_package_name, 'workflow']))
varianthandler = importlib.import_module('.'.join([top_package_name, 'variantplus', 'varianthandler']))
headerhandler = importlib.import_module('.'.join([top_package_name, 'variantplus', 'headerhandler']))


def sanity_check_args(infile_path_list, outfile_path, isec, isec_indices,
                      outfile_must_not_exist):
    infile_path_list = workflow.arghandler_infile_list(infile_path_list)

    arghandler = workflow.get_arghandler_outfile(outfile_must_not_exist)
    outfile_path = arghandler(outfile_path)

    sanity_check_isec_indices(isec_indices, len(infile_path_list))

    return infile_path_list, outfile_path


def sanity_check_isec_indices(isec_indices, n_vcf):
    e_msg = textwrap.dedent(f"""\
        'isec_indices' must be a list composed of 0 or 1, \
        with the length the same as the number of vcf inputs.""")
    if isec_indices is not None:
        if len(isec_indices) != n_vcf:
            raise Exception(e_msg)
        elif not set(isec_indices).issubset({0, 1}):
            raise Exception(e_msg)


############################################


def isec_indices_to_bool(isec_indices, n_vcf):
    if isec_indices is None:
        isec_indices_bool = list(itertools.repeat(True, n_vcf))
    else:
        isec_indices_bool = [bool(x) for x in isec_indices]

    return isec_indices_bool


def load_vcf_data(vcf_list):
    """
    Returns: 
        [
            { vcfspec, vcfspec, ... },
            { vcfspec, vcfspec, ... },
            ...
        ] (same length as len(infile_path_list)) 
    """

    raw_vcfspecs = list()
    vr_dict = dict()
    for in_vcf in vcf_list:
        subset = set()
        for vr in in_vcf.fetch():
            vcfspec = varianthandler.get_vcfspec(vr)
            subset.add(vcfspec)
            if vcfspec not in vr_dict:
                vr_dict[vcfspec] = list()
            vr_dict[vcfspec].append(vr)

        raw_vcfspecs.append(subset)

    return raw_vcfspecs, vr_dict


def load_vcf_data_vcfpath(infile_path_list):
    vcf_list = [pysam.VariantFile(x) for x in infile_path_list]
    raw_vcfspecs, vr_dict = load_vcf_data(vcf_list)
    for vcf in vcf_list:
        vcf.close()

    return raw_vcfspecs, vr_dict


########################################


def get_merged_header(infile_path_list):
    """Any conflicting INFO or FORMAT keys 
    (with regard to any of Type, Number, or Description) are discarded."""

    header_list = list()
    for vcf_path in infile_path_list:
        with pysam.VariantFile(vcf_path, 'r') as in_vcf:
            header_list.append(in_vcf.header)
    #compatible = headerhandler.check_header_compatibility(header_list)
    merged_header = headerhandler.merge_vcfheaders(header_list)

    return merged_header


def get_output_vcfspecs_isec(raw_vcfspecs, isec_indices_bool, chromdict):
    vcfspecs_included = set.intersection( 
        *itertools.compress(raw_vcfspecs, isec_indices_bool))
    for subset in itertools.compress(raw_vcfspecs, 
                                     [not x for x in isec_indices_bool]):
        vcfspecs_included.difference_update(subset)
    output_vcfspecs = sorted(vcfspecs_included, 
                             key=common.get_vcfspec_sortkey(chromdict))

    return output_vcfspecs


def get_output_vcfspecs_union(raw_vcfspecs, chromdict):
    vcfspec_set = set.union(*raw_vcfspecs)
    output_vcfspecs = sorted(vcfspec_set, 
                             key=common.get_vcfspec_sortkey(chromdict))

    return output_vcfspecs


def reform_vr_dict(vr_dict, merged_header, output_vcfspecs):
    merged_vr_dict = dict()
    for vcfspec in output_vcfspecs:
        merged_vr_dict[vcfspec] = varianthandler.merge(
            vr_dict[vcfspec], merged_header)
        del vr_dict[vcfspec]

    return merged_vr_dict


def write(merged_vr_dict, output_vcfspecs, outfile_path, mode_pysam, 
          merged_header):
    with pysam.VariantFile(outfile_path, mode=mode_pysam, 
                           header=merged_header) as out_vcf:
        for vcfspec in output_vcfspecs:
            out_vcf.write(merged_vr_dict[vcfspec])


#########################################


def main_file(
        infile_path_list, 
        outfile_path, 
        fasta_path, 
        isec=False, 
        isec_indices=None, 
        mode_bcftools='z', 
        mode_pysam=None,
        outfile_must_not_exist='ask',
        logger=None,
        ):
    if logger is None:
        logger = workflow.get_logger(name='merge_vcf', stderr=True)

    logger.info('Beginning')
    infile_path_list, outfile_path = sanity_check_args(
        infile_path_list, outfile_path, isec, isec_indices, 
        outfile_must_not_exist)
    isec_indices_bool = isec_indices_to_bool(isec_indices, 
                                             len(infile_path_list))
    mode_pysam = common.write_mode_arghandler(mode_bcftools, mode_pysam)
    merged_header = get_merged_header(infile_path_list)
    chromdict = common.ChromDict(fasta_path=fasta_path)

    with pysam.FastaFile(fasta_path) as fasta:
        vcf_list = [pysam.VariantFile(x) for x in infile_path_list]

        logger.info('Loading input vcf data')
        raw_vcfspecs, vr_dict = load_vcf_data(vcf_list)

        logger.info('Getting union/intersection of vcfspecs')
        if isec:
            output_vcfspecs = get_output_vcfspecs_isec(
                raw_vcfspecs, isec_indices_bool, chromdict)
        else:
            output_vcfspecs = get_output_vcfspecs_union(
                raw_vcfspecs, chromdict)

        logger.info('Merging overlapping variant records')
        merged_vr_dict = reform_vr_dict(vr_dict, merged_header, 
                                        output_vcfspecs)

        logger.info('Writing final result')
        write(merged_vr_dict, output_vcfspecs, outfile_path, mode_pysam, 
              merged_header)

        for vcf in vcf_list:
            vcf.close()


#def main_vcfplus(vcfplus_list, isec=False, isec_indices=None, merge=True):
#    #assert all(isinstance(x, vcfplus.VcfPlus) for x in vcfplus_list)
#    def get_vp_list(output_vcfspecs, pysamvr_dict, fasta, chromdict):
#        vp_list = list()
#        for vcfspec in output_vcfspecs:
#            for pysamvr in pysamvr_dict[vcfspec]:
#                vp = variantplus.VariantPlus(pysamvr, fasta, chromdict)
#                vp_list.append(vp)
#
#        return vp_list
#
#    sanity_check_isec_indices(isec_indices, len(vcfplus_list))
#    chromdict = vcfplus_list[0].chromdict
#    fasta = vcfplus_list[0].fasta
#    isec_indices_bool = isec_indices_to_bool(isec_indices, len(vcfplus_list))
#
#    merged_header = merge_pysamhdr(vcfp.vcf.header for vcfp in vcfplus_list)
#    raw_vcfspecs = load_vcf_data(vcf_list)
#    output_vcfspecs = vcfmergelib.get_output_vcfspecs(raw_vcfspecs, isec, isec_indices_bool, chromdict)
#    pysamvr_dict = vcfmergelib.get_pysamvr_dict(vcfplus_list, output_vcfspecs, merged_header, isec_indices_bool, merge)
#    vp_list = get_vp_list(output_vcfspecs, pysamvr_dict, fasta, chromdict)
#
#    vcfp = vcfplus.VcfPlus(fasta = fasta)
#    vcfp.set_header(header = merged_header)
#    vcfp.init_vp(vp_list = vp_list)
#
#    return vcfp
