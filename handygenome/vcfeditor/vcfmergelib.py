import itertools 

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
variantplus = importlib.import_module('.'.join([top_package_name, 'variantplus', 'variantplus']))
varianthandler = importlib.import_module('.'.join([top_package_name, 'variantplus', 'varianthandler']))


def sanity_check_isec_indices(isec_indices, n_vcf):
    if isec_indices is not None:
        if (
                (len(isec_indices) != n_vcf) or 
                (not set(isec_indices).issubset({0,1}))):
            raise Exception("'isec_indices' must be a list composed of 0 or 1, with length same as the number of vcf inputs.")


def edit_isec_indices(isec_indices, n_vcf):
    if isec_indices is None:
        isec_indices_edit = list(itertools.repeat(True, n_vcf))
    else:
        isec_indices_edit = [ bool(x) for x in isec_indices ]

    return isec_indices_edit


def get_raw_coordIDs(vcfplus_list):
    """
    Returns: A list composed of sublists corresponding to each input vcf file. Each sublist consist of tuples (contig, pos, alleles).
    """

    raw_coordIDs = list()
    for vcfplus in vcfplus_list:
        sublist = list()
        for pysamvr in vcfplus.vcf.fetch():
            ID = varianthandler.get_coordID_from_pysamvr(pysamvr)
            if ID in sublist:
                pass
                #raise Exception(f'Duplicate CHROM/POS/REF/ALT {ID} in vcf file {vcfp.vcf_path}')
            sublist.append(ID)
        raw_coordIDs.append(sublist)

    return raw_coordIDs


def get_output_coordIDs(raw_coordIDs, isec, isec_indices_edit, chromdict):
    IDset = set()
    if isec:
        for idx, sublist in enumerate( itertools.compress(raw_coordIDs, isec_indices_edit) ):
            if idx == 0:
                IDset.update(sublist)
            else:
                IDset.intersection_update(sublist)
        for sublist in itertools.compress(raw_coordIDs, [not x for x in isec_indices_edit]):
            IDset.difference_update(sublist)
    else:
        for sublist in raw_coordIDs:
            IDset.update(sublist)

    output_coordIDs = sorted(IDset, key = lambda x: common.coord_sortkey(x[0], x[1], chromdict))

    return output_coordIDs


def get_pysamvr_dict(vcfplus_list, output_vcfspecs, merged_header, 
                     isec_indices_edit, merge):
    """
    Args:
        vcfplus_list: A list of vcfplus object. Each must have 'vp_list' or 
            'vcf' attribute set. 'vp_list' has precedence.
    Returns:
        A dictionary
            - keys: vcfspec
            - values: a list of pysamvr corresponding to the vcfspec 
                among all vcfs 
    """

    # init
    pysamvr_dict = dict()
    for vcfspec in output_vcfspecs:
        pysamvr_dict[vcfspec] = list()

    # gather raw pysamvr's
    for vcfp in itertools.compress(vcfplus_list, isec_indices_edit):
        if vcfp.vp_list is None:
            if vcfp.vcf is None:
                raise Exception('One of input vcfplus object does not have both "vp_list" or "vcf" attributes.')
            pysamvr_iter = vcfp.vcf.fetch()
        else:
            pysamvr_iter = (x.pysamvr for x in vcfp.vp_list)
            
        for pysamvr in pysamvr_iter:
            vcfspec = varianthandler.get_vcfspec_from_pysamvr(pysamvr)
            if vcfspec in output_vcfspecs:
                pysamvr_dict[vcfspec].append(pysamvr)

    # reform raw pysamvr's
    for vcfspec in output_vcfspecs:
        if merge:
            if len(pysamvr_dict[vcfspec]) == 1:
                pysamvr_dict[vcfspec] = [ varianthandler.reform_samples(pysamvr_dict[vcfspec][0], pysamhdr = merged_header) ]
            else:
                pysamvr_dict[vcfspec] = [ varianthandler.merge(pysamvr_dict[vcfspec], pysamhdr = merged_header, clear = False) ]
        else:
            pysamvr_dict[vcfspec] = [ 
                varianthandler.reform_samples(x, pysamhdr = merged_header) 
                for x in pysamvr_dict[vcfspec] 
                ]

    return pysamvr_dict
