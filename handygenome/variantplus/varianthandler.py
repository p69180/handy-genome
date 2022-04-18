import re
import warnings
import statistics
import textwrap

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
infoformat = importlib.import_module('.'.join([top_package_name, 'variantplus', 'infoformat']))
headerhandler = importlib.import_module('.'.join([top_package_name, 'variantplus', 'headerhandler']))
initvcf = importlib.import_module('.'.join([top_package_name, 'vcfeditor', 'initvcf']))


def get_vcfspec(vr):
    return common.Vcfspec(vr.contig, vr.pos, vr.ref, vr.alts[0])


def check_SV(vr):
    return any( 
            common.RE_PATS['nucleobases'].fullmatch(x) is None
            for x in vr.alts
            )


def check_is_bnd1(pysamvr, bndspec):
    matches_bnd1 = any(
            (pysamvr.contig == form.chrom1 and pysamvr.pos == form.pos1)
            for form in bndspec.maxscore_forms
            )
    matches_bnd2 = any(
            (pysamvr.contig == form.chrom2 and pysamvr.pos == form.pos2)
            for form in bndspec.maxscore_forms
            )

    if matches_bnd1 and matches_bnd2:
        raise Exception('pysamvr matches both bnd1 and bnd2')
    elif (not matches_bnd1) and (not matches_bnd2):
        raise Exception('pysamvr matches neither bnd1 nor bnd2')
    else:
        return matches_bnd1


# sorting helpers

def vr_sortkey(vr, chromdict):
    return common.coord_sortkey(vr.contig, vr.pos, chromdict)


# applying vcfspec to vr

def apply_vcfspec(vr, vcfspec):
    vr.contig = vcfspec.chrom
    vr.pos = vcfspec.pos
    vr.alleles = (vcfspec.ref, vcfspec.alt)


# pysamvr reformatting

#def merge(pysamvr_list, pysamhdr = None, clear = False):
#    def check_coord_identity(pysamvr_list):
#        for key in ('contig', 'pos', 'alleles'):
#            vals = ( getattr(x, key) for x in pysamvr_list )
#            if len(set(vals)) != 1:
#                raise Exception(f"'{key}' attributes of input variant records are different.")
#
#    def check_sampleid_identity(pysamvr_list, pysamhdr):
#        if pysamhdr is not None:
#            samples_pysamvr = set()
#            for pysamvr in pysamvr_list:
#                samples_pysamvr.update(pysamvr.samples.keys())
#
#            samples_pysamhdr = set(pysamhdr.samples)
#
#            if not samples_pysamvr.issubset(samples_pysamhdr):
#                raise Exception(f"Sample IDs of input 'pysamhdr' are not included among those of input 'pysamvr_list'.")
#
#    def get_merged_header(pysamhdr):
#        if pysamhdr is None:
#            #if not headerhandler.check_header_compatibility(x.header for x in pysamvr_list):
#                #warnings.warn('Headers of input variant records do not coincide with regards to INFO or FORMAT metadata Type or Number.')
#            merged_header = headerhandler.merge_pysamhdr([x.header for x in pysamvr_list])
#        else:
#            merged_header = pysamhdr
#
#        return merged_header
#
#    def get_merged_info(pysamvr_list):
#        """
#        If there is an overlapping field between multiple variant records, only the first one is represented.
#        """
#
#        merged_info = dict()
#        if not clear:
#            for pysamvr in pysamvr_list:
#                for key, val in pysamvr.info.items():
#                    if key not in merged_info:
#                        merged_info[key] = val
#
#        return merged_info
#
#    def get_merged_format(pysamvr_list, merged_header):
#        """
#        If there is an overlapping field between multiple variant records, only the first one is represented.
#        """
#
#        merged_format_dict = dict()
#        for sampleid in merged_header.samples:
#            merged_format_dict[sampleid] = dict()
#
#        for pysamvr in pysamvr_list:
#            for sampleid, annots in pysamvr.samples.items():
#                for key, val in annots.items():
#                    if key == 'GT':
#                        continue
#                    if key not in merged_format_dict[sampleid]:
#                        merged_format_dict[sampleid][key] = val
#
#        merged_format = list()
#        for sampleid in merged_header.samples:
#            merged_format.append(merged_format_dict[sampleid])
#
#        return merged_format
#
#    def init_new_pysamvr(merged_header, clear):
#        if clear:
#            new_pysamvr = merged_header.new_record()
#        else:
#            merged_info = get_merged_info(pysamvr_list)
#            merged_format = get_merged_format(pysamvr_list, merged_header)
#            new_pysamvr = merged_header.new_record(info = merged_info, samples = merged_format)
#
#        return new_pysamvr
#
#    def add_other_attrs(new_pysamvr, pysamvr_list):
#        new_pysamvr.contig = pysamvr_list[0].contig
#        new_pysamvr.pos = pysamvr_list[0].pos
#        new_pysamvr.alleles = pysamvr_list[0].alleles
#    
#    def set_qual(new_pysamvr, pysamvr_list):
#        qual_list = list(x.qual for x in pysamvr_list if x.qual is not None)
#        if len(qual_list) > 0:
#            new_pysamvr.qual = int(statistics.mean(x.qual for x in pysamvr_list))
#
#    def add_genotype(new_pysamvr, pysamvr_list):
#        used_sampleids = set()
#        for pysamvr in pysamvr_list:
#            for sampleid, annots in pysamvr.samples.items():
#                if sampleid not in used_sampleids:
#                    used_sampleids.add(sampleid)
#                    new_pysamvr.samples[sampleid].allele_indices = annots.allele_indices
#                    new_pysamvr.samples[sampleid].phased = annots.phased
#
#    def add_filters(new_pysamvr, pysamvr_list):
#        if all(x.filter.keys() == ['PASS'] for x in pysamvr_list):
#            new_pysamvr.filter.add('PASS')
#        else:
#            for pysamvr in pysamvr_list:
#                for key in pysamvr.filter.keys():
#                    if key != 'PASS':
#                        new_pysamvr.filter.add(key)
#
#    def main(pysamvr_list, pysamhdr, clear):
#        check_coord_identity(pysamvr_list)
#        check_sampleid_identity(pysamvr_list, pysamhdr)
#        merged_header = get_merged_header(pysamhdr)
#        new_pysamvr = init_new_pysamvr(merged_header, clear)
#        add_other_attrs(new_pysamvr, pysamvr_list)
#        set_qual(new_pysamvr, pysamvr_list)
#        if not clear:
#            add_genotype(new_pysamvr, pysamvr_list)
#            add_filters(new_pysamvr, pysamvr_list)
#
#        return new_pysamvr
#
#    return main(pysamvr_list, pysamhdr, clear)


# helpers
def _copy_basic_attrs(old_vr, new_vr):
    for key in ('contig', 'pos', 'id', 'qual'):
        setattr(new_vr, key, getattr(old_vr, key))
    if old_vr.alleles is not None:
        # assigning None to vr.alleles raises an Exception
        new_vr.alleles = old_vr.alleles
    if old_vr.stop is not None:
        new_vr.stop = old_vr.stop


def _copy_filter(old_vr, new_vr, skip_pass=True):
    for filter_name in old_vr.filter:
        # check if new vr header defines the key
        if filter_name in new_vr.header.filters:
            if skip_pass:
                # copy filter only if it is not PASS
                if filter_name != 'PASS':
                    new_vr.filter.add(filter_name)
            else:
                new_vr.filter.add(filter_name)


def _copy_info(old_vr, new_vr, dont_copy_if_exist=True):
    for info_key in old_vr.info:
        # check if new vr header defines the key
        if info_key in new_vr.header.info:  
            if dont_copy_if_exist:
                # copy the value only if the value is not already set
                if info_key not in new_vr.info:
                    infoformat.set_info(new_vr, info_key, 
                                        infoformat.get_info(old_vr, info_key))
            else:
                infoformat.set_info(new_vr, info_key, 
                                    infoformat.get_info(old_vr, info_key))


def _copy_format(old_vr, new_vr, dont_copy_if_exist=True):
    def assign_val(old_vr, new_vr, key, sampleid):
        if key == 'GT':
            new_vr.samples[sampleid].allele_indices = \
                old_vr.samples[sampleid].allele_indices
            new_vr.samples[sampleid].phased = old_vr.samples[sampleid].phased
        else:
            infoformat.set_format(new_vr, sampleid, key, 
                                  infoformat.get_format(old_vr, sampleid, key))

    for sampleid in old_vr.samples:
        # check if the new vr header defines the sample id of input vr
        if sampleid in tuple(new_vr.samples):
            for key in old_vr.samples[sampleid]:
                # check if new vr header defines the key
                if key in new_vr.header.formats:
                    if dont_copy_if_exist:
                        # copy the value only if the value is not already set
                        if key not in new_vr.samples[sampleid]:
                            assign_val(old_vr, new_vr, key, sampleid)
                    else:
                        assign_val(old_vr, new_vr, key, sampleid)


########


def merge(vr_list, vcfheader=None):
    """Intended to be run with variant records with identical vcfspecs.
    chrom, pos, id, ref, alt, qual, INFO/END is derived from the first vr.
    FILTER values are union'ed except PASS value.
    INFO and FORMAT values are chosen from the first variant record which
        contains the key.
    """

    if vcfheader is None:
        vcfheader = headerhandler.merge_vcfheaders([x.header for x in vr_list])

    new_vr = vcfheader.new_record()
    _copy_basic_attrs(vr_list[0], new_vr)

    for idx, vr in enumerate(vr_list):
        _copy_filter(vr, new_vr, skip_pass=True)
        _copy_info(vr, new_vr, dont_copy_if_exist=True)
        _copy_format(vr, new_vr, dont_copy_if_exist=True)

    return new_vr


def reheader(vr, vcfheader):
    """Return a new variant record created using the template
    header, filled with the contents of the input variant record.

    INFO or FORMAT keys or samples not included in the template header 
    are discarded.
    """

    new_vr = vcfheader.new_record()
    _copy_basic_attrs(vr, new_vr)
    _copy_filter(old_vr, new_vr, skip_pass=False)
    _copy_info(old_vr, new_vr, dont_copy_if_exist=False)
    _copy_format(old_vr, new_vr, dont_copy_if_exist=False)

    return new_vr


def rename(vr, samples):
    if len(samples) != len(vr.header.samples):
        raise Exception(
            f'The number of samples must be the same.')

    chromdict = common.ChromDict(vcfheader=vr.header)
    new_header = initvcf.create_header(chromdict=chromdict, samples=samples, 
                                       vcfheader=vr.header)
    new_vr = new_header.new_record()
    _copy_basic_attrs(vr, new_vr)
    _copy_filter(old_vr, new_vr, skip_pass=False)
    _copy_info(old_vr, new_vr, dont_copy_if_exist=False)

    for old_sampleid, new_sampleid in zip(vr.header.samples, samples):
        for key in vr.samples[old_sampleid]:
            if key == 'GT':
                new_vr.samples[new_sampleid].allele_indices = \
                    vr.samples[old_sampleid].allele_indices
                new_vr.samples[new_sampleid].phased = \
                    vr.samples[old_sampleid].phased
            else:
                infoformat.set_format(
                    new_vr, key, new_sampleid,
                    infoformat.get_format(vr, old_sampleid, key))

    return new_vr


#def translate(pysamvr, samples = None, pysamhdr = None, rename = False):
#    return reform_samples(pysamvr, samples = None, pysamhdr = None, rename = False)
#
#
#def reform_samples(pysamvr, samples = None, pysamhdr = None, rename = False):
#    """
#    Args:
#        pysamvr: pysam.VariantRecord object
#        samples: An iterable containing new sample names
#        pysamhdr: pysam.VariantHeader object
#        rename: If True, only sample names are changed without editing format entries
#
#    Returns:
#        A new pysam.VariantRecord object with format columns in the order represented by input argument 'samples' or 'pysamhdr'.
#        
#        The output is created with 'new_record' method of pysam.VariantHeader object.
#        The pysam.VariantHeader object is 'pysamhdr' argument if it is given, otherwise created using 'header' attribute of 'pysamvr' argument and 'samples' argument.
#    """
#
#    def sanity_check(pysamvr, samples, pysamhdr, rename):
#        common.check_num_notNone(1, (samples, pysamhdr), ('samples', 'pysamhdr'))
#        if rename:
#            e_msg = f'Numbers of samples are different between original pysam.VariantRecord object and new sample list.'
#            if samples is not None:
#                if len(samples) != len(pysamvr.samples.keys()):
#                    raise ValueError(e_msg)
#            elif pysamhdr is not None:
#                if len(pysamhdr.samples) != len(pysamvr.samples.keys()):
#                    raise ValueError(e_msg)
#
#    def check_identical_samples(pysamvr, samples, pysamhdr):
#        if samples is not None:
#            return tuple(pysamvr.samples.keys()) == tuple(samples)
#        elif pysamhdr is not None:
#            return tuple(pysamvr.samples.keys()) == tuple(pysamhdr.samples)
#
#    def get_params(pysamvr, samples, pysamhdr):
#        pysamvr_copy = pysamvr.copy()
#        infoformat.refine_vr_InfoFormatValue(pysamvr_copy)
#        if pysamhdr is not None:
#            new_pysamhdr = pysamhdr
#            new_samples = tuple(pysamhdr.samples)
#        elif samples is not None:
#            new_pysamhdr = initvcf.create_header(samples = samples, header = pysamvr.header)
#            new_samples = tuple(samples)
#
#        return pysamvr_copy, new_pysamhdr, new_samples
#
#    def get_format_items(pysamvr_copy, new_samples):
#        format_items = list()
#        for sampleid in new_samples:
#            if sampleid in pysamvr_copy.samples.keys():
#                items = dict( pysamvr_copy.samples[sampleid] )
#                if 'GT' in items:
#                    del items['GT']
#            else:
#                items = make_empty_format(pysamvr_copy)
#    
#            format_items.append(items)
#
#        return format_items
#
#    def get_format_items_rename(pysamvr_copy):
#        format_items = list()
#        for val in pysamvr_copy.samples.values():
#            items = dict( val )
#            if 'GT' in items:
#                del items['GT']
#            format_items.append(items)
#
#        return format_items
#
#    def init_new_pysamvr(pysamvr_copy, new_pysamhdr, format_items):
#        new_pysamvr = new_pysamhdr.new_record(filter = pysamvr_copy.filter, info = dict(pysamvr_copy.info), samples = format_items)
#
#        return new_pysamvr
#
#    def add_genotype(new_pysamvr, pysamvr_copy, new_samples):
#        # genotype information
#        new_pysamvr.alleles = pysamvr_copy.alleles
#        for sampleid in new_samples:
#            if sampleid in pysamvr_copy.samples.keys():
#                new_pysamvr.samples[sampleid].allele_indices = pysamvr_copy.samples[sampleid].allele_indices
#                new_pysamvr.samples[sampleid].phased = pysamvr_copy.samples[sampleid].phased
#            else:
#                new_pysamvr.samples[sampleid].allele_indices = None
#                new_pysamvr.samples[sampleid].phased = False
#
#    def add_genotype_rename(new_pysamvr, pysamvr_copy):
#        # genotype information
#        new_pysamvr.alleles = pysamvr_copy.alleles
#        for val1, val2 in zip(new_pysamvr.samples.values(), pysamvr_copy.samples.values()):
#            val1.allele_indices = val2.allele_indices
#            val1.phased = val2.phased
#    
#    def add_other_attrs(new_pysamvr, pysamvr_copy):
#        new_pysamvr.contig = pysamvr_copy.contig
#        new_pysamvr.pos = pysamvr_copy.pos
#        new_pysamvr.id = pysamvr_copy.id
#        new_pysamvr.qual = pysamvr_copy.qual
#
#    def main(pysamvr, samples, pysamhdr, rename):
#        sanity_check(pysamvr, samples, pysamhdr, rename)
#
#        if check_identical_samples(pysamvr, samples, pysamhdr):
#            return pysamvr
#
#        else:
#            pysamvr_copy, new_pysamhdr, new_samples = get_params(pysamvr, samples, pysamhdr)
#            format_items = ( 
#                    get_format_items_rename(pysamvr_copy)
#                    if rename else
#                    get_format_items(pysamvr_copy, new_samples)
#                    )
#            new_pysamvr = init_new_pysamvr(pysamvr_copy, new_pysamhdr, format_items)
#            if rename:
#                add_genotype_rename(new_pysamvr, pysamvr_copy)
#            else:
#                add_genotype(new_pysamvr, pysamvr_copy, new_samples)
#            add_other_attrs(new_pysamvr, pysamvr_copy)
#    
#            return new_pysamvr
#
#    return main(pysamvr, samples, pysamhdr, rename)
#
#
#def make_empty_format(pysamvr):
#    # Attempting to create a FORMAT field with Type=Flag results in this message:
#        # [E::vcf_parse_format] The format type 0 is currently not supported
#    result = dict()
#    for key, val in pysamvr.format.items():
#        if val.type == 'Flag':
#            raise Exception(f"FORMAT field '{key}' has 'Flag' type.")
#        if key == 'GT':
#            continue # GT value is assigned in a separate way
#        result[key] = infoformat.modify_InfoFormatValue(
#                key = key,
#                val = None,
#                transl_number = infoformat.translate_metadata_number(pysamvr, val.number),
#                keytype = val.type,
#                for_write = True,
#                typeconv = False,
#                sanitycheck = False,
#                )
#
#    return result


# functions dealing with list of pysamvr
    
def check_has_id(pysamvr_list):
    """
    Args:
        pysamvr_list: A list of pysam.VariantRecord objects
    """

    result = True
    for pysamvr in pysamvr_list:
        if pysamvr.id is None:
            result = False
            break

    return result


def check_has_duplicate_pos_alleles(pysamvr_list):
    """
    Args:
        pysamvr_list: A list of pysam.VariantRecord objects
    """

    result = False
    IDset = set()
    for pysamvr in pysamvr_list:
        ID = (pysamvr.contig, pysamvr.pos, pysamvr.alleles)
        if ID in IDset:
            result = True
            break
        else:
            IDset.add(ID)

    return result


def check_has_duplicate_id(pysamvr_list):
    """
    Args:
        pysamvr_list: A list of pysam.VariantRecord objects

    Returns:
        True when two or more pysamvr have no ID.
    """

    has_duplicate_id = False
    IDset = set()
    for pysamvr in pysamvr_list:
        if pysamvr.id in IDset:
            has_duplicate_id = True
            break
        else:
            IDset.add(pysamvr.id)

    return has_duplicate_id


def check_has_mateid(pysamvr_list):
    """
    Args:
        pysamvr_list: A list of pysam.VariantRecord objects
    """

    has_mateid = True
    for pysamvr in pysamvr_list:
        if 'MATEID' not in pysamvr.info.keys():
            has_mateid = False
            break
    
    return has_mateid

