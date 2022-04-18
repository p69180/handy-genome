import textwrap
import math

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))


NA_VALUES = ('', '.', None)


'''
class Params(common.Params):
    NA_values = ('', '.', None)
    '''


def show_info(vr):
    for key in vr.info:
        print(
                key, 
                get_value_info(vr, key), 
                sep = '\t'
                )


def show_format(vr):
    for sampleid in vr.header.samples:
        print('\n' + sampleid)
        for key in vr.samples[sampleid]:
            print(
                    key, 
                    get_value_format(vr, sampleid, key), 
                    sep = '\t'
                    )


def get_info(vr, key, collapse_tuple = True):
    """
    Returns:
        NA-equivalent values ("", ".", None) are converted into None
        An empty tuple is converted into None
        If the original value is a length-1 tuple, returns its element, with NA conversion.
        If the original value is a tuple longer than 1, return type is tuple.
        If the original value is a tuple longer than 1 composed of only NA-equivalents, returns a tuple with the same length composed of None.
    """

    if key not in vr.header.info.keys():
        raise KeyError(f'Input key "{key}" is not included in the input variant record INFO metadata.')

    if key in vr.info:
        transl_number = get_transl_number_info(vr, key)
        keytype = vr.header.info[key].type
        return modify_InfoFormatValue(
            key,
            vr.info[key], 
            transl_number,
            keytype,
            for_write = False,
            typeconv = False,
            collapse_tuple = collapse_tuple,
            )
    else:
        return None

get_value_info = get_info


def get_format(vr, sampleid, key, collapse_tuple = True):
    """
    Returns:
        NA-equivalent values ("", ".", None) are converted into None
        An empty tuple is converted into None
        If the original value is a length-1 tuple, returns its element, with NA conversion.
        If the original value is a tuple longer than 1, return type is tuple.
        If the original value is a tuple longer than 1 composed of only NA-equivalents, returns a tuple with the same length composed of None.
    """

    if key == 'GT':
        raise Exception(f'"GT" is not treated with this function.')

    if key not in vr.header.formats.keys():
        raise KeyError(f'Input key "{key}" is not included in the input variant record FORMAT metadata.')

    if sampleid not in vr.header.samples:
        raise KeyError(f'Input sample ID "{sampleid}" is not included in the input variant record sample name list.')

    if key in vr.samples[sampleid]:
        transl_number = get_transl_number_format(vr, key)
        keytype = vr.header.formats[key].type
        return modify_InfoFormatValue( 
            key,
            vr.samples[sampleid][key], 
            transl_number,
            keytype,
            for_write = False,
            typeconv = False,
            collapse_tuple = collapse_tuple,
            )
    else:
        return None

get_value_format = get_format


def get_genotype(vr, sampleid):
    return vr.samples[sampleid].allele_indices, vr.samples[sampleid].phased


def check_NA_info(vr, key):
    if key in vr.info:
        return _check_InfoFormatValue_isNA(vr.info[key])
    else:
        return True


def check_NA_format(vr, sampleid, key):
    if key in vr.samples[sampleid]:
        return _check_InfoFormatValue_isNA(vr.samples[sampleid][key])
    else:
        return True


def _check_InfoFormatValue_isNA(val):
    """
    Args:
        val: A raw value obtained from vr.info[key] or vr.samples[sampleid][key]
    """

    if isinstance(val, tuple):
        if len(val) == 0:
            """
            Value of FOO1 can be an empty tuple when:
                - Number > 1
                AND
                - INFO column looks like: FOO1;FOO2=1;FOO3=10
            """
            return True
        else:
            return set(val).issubset(set(NA_VALUES))
    else:
        return val in NA_VALUES


#############################################


def set_info(vr, key, val, typeconv = True):
    if key not in vr.header.info:
        raise Exception(f'Key {key} is absent from INFO metadata header.')
    else:
        transl_number = get_transl_number_info(vr, key)
        keytype = vr.header.info[key].type
        modified_val = modify_InfoFormatValue(
            key, 
            val, 
            transl_number, 
            keytype, 
            for_write = True,
            typeconv = typeconv,
            collapse_tuple = True,
            )
        vr.info[key] = modified_val

set_value_info = set_info


def set_NA_info(vr, key):
    set_value_info(vr, key, None)


def set_format(vr, sampleid, key, val, typeconv = True):
    if key == 'GT':
        raise Exception(f'"GT" is not treated with this function.')

    if key not in vr.header.formats:
        raise Exception(f'Key {key} is absent from FORMAT metadata header.')
    else:
        transl_number = get_transl_number_format(vr, key)
        keytype = vr.header.formats[key].type
        modified_val = modify_InfoFormatValue(
            key,
            val, 
            transl_number, 
            keytype, 
            for_write = True,
            typeconv = typeconv,
            collapse_tuple = True,
            )
        vr.samples[sampleid][key] = modified_val

set_value_format = set_format


def set_genotype(vr, sampleid, gt, phased):
    """
    Args:
        gt: A tuple composed of integers
        phased: True or False
    """

    vr.samples[sampleid].allele_indices = gt
    vr.samples[sampleid].phased = phased


def set_NA_format(vr, sampleid, key):
    set_value_format(vr, sampleid, key, None)


#############################################

"""
def make_vr_for_write(vr):
    new_vr = vr.copy()
    # INFO
    for key, old_val in new_vr.info.items():
        set_value_info(new_vr, key, old_val)
    # FORMAT
    for sampleid in new_vr.header.samples:
        for key, old_val in new_vr.samples[sampleid].items():
            set_value_format(new_vr, sampleid, key, old_val)

    return new_vr
"""


def refine_vr_InfoFormatValue(vr):
    """
    Modifies each INFO or FORMAT value into
    vcf writing-compatible form. (e.g. For Type=String 
    (None,None) into ('.','.'))
    """

    # INFO
    for key, old_val in vr.info.items():
        set_value_info(vr, key, old_val)
    # FORMAT
    for sampleid in vr.header.samples:
        for key, old_val in vr.samples[sampleid].items():
            set_value_format(vr, sampleid, key, old_val)


#############################################

# not used
def get_genotype_count(vr):
    if 'GT' in vr.format.keys():
        GT_count_set = set(len(x['GT']) for x in vr.samples.values())
        if len(GT_count_set) != 1:
            raise Exception(f'Genotype counts are different between samples in this variant record:\n{vr}')
        genotype_count = GT_count_set.pop()
    else:
        genotype_count = None

    return genotype_count


def translate_metadata_number(vr, metadata_number):
    if metadata_number == 'A':
        transl_number = len(vr.alts)
    elif metadata_number == 'R':
        transl_number = len(vr.alts) + 1
    elif metadata_number == 'G':
        # combinations with repetition H(len(vr.alleles), 2)
        transl_number = math.comb(len(vr.alleles) + 2 - 1, 2)
    else:
        transl_number = metadata_number

    return transl_number


def get_transl_number_info(vr, key):
    if key in vr.header.info.keys():
        metadata_number = vr.header.info[key].number
        return translate_metadata_number(vr, metadata_number)
    else:
        return None


def get_transl_number_format(vr, key):
    if key in vr.header.formats.keys():
        metadata_number = vr.header.formats[key].number
        return translate_metadata_number(vr, metadata_number)
    else:
        return None


#############################################


def modify_InfoFormatValue(
        key, 
        val, 
        transl_number, 
        keytype, 
        for_write = False, 
        typeconv = False, 
        collapse_tuple = True,
        ):
    """
    If input is an empty tuple: returns None
    Otherwise:
        If input is a length-1 tuple, it is converted into its element.
        Then,
        for each atomic value, if it is one of '' or '.',
            it is converted into None,
            then the result is returned.
    Boolean input is returned unchanged.

    Special cases - value is unchanged for these ones:
        INFO/AS_SB_TABLE (Mutect2)
        FORMAT/GT
        Type=Flag

    Args:
        val: A raw value obtained from vr.info[key] or vr.samples[sampleid][key]
        for_write: If True, modification goes for writing into a vcf file.
        keytype: VCF metadata 'Type' value. Only relevant when for_write == True.
        collapse_tuple: If True, length-1 tuple is converted to its element
    """

    def handle_special_cases(key, val, keytype):
        if key == 'AS_SB_TABLE':
            # Exception for Mutect2 INFO annotation AS_SB_TABLE
            # Number == 1 but its value looks like '20,14|4,4'
            stop = True
            return_val = ','.join(val)
        elif keytype == 'Flag' or key == 'GT':
            # Flag or GT should not be modified
            stop = True
            return_val = val
        else:
            stop = False
            return_val = None

        return stop, return_val

    def get_params(val, for_write, keytype):
        if isinstance(val, tuple):
            len_val = 1 if len(val) == 0 else len(val)
        else:
            len_val = 1

        isNA = _check_InfoFormatValue_isNA(val)
            # True with zero-length tuple
            # True with a non-length-zero tuple composed of NA equivalents

        if for_write:
            unifiedNA = '.' if keytype == 'String' else None
        else:
            unifiedNA = None

        return len_val, isNA, unifiedNA

    def sanity_check(
            key, val, transl_number, keytype, len_val, isNA,
            ):
        if (isinstance(transl_number, int) and
            (len_val != transl_number) and
            (not isNA)):
            raise Exception(
                (f'"Number" metadata and actual number of the value '
                 f'does not match.\n'
                 f'key = {key}, val = {val}, '
                 f'transl_number = {transl_number}, keytype = {keytype}'))

    def get_modified_val(
            val,
            transl_number,
            isNA,
            len_val,
            unifiedNA,
            collapse_tuple,
            ):
        if isNA:
            if transl_number == '.':
                if len_val == 1:
                    modified_val = unifiedNA
                else:
                    modified_val = tuple([unifiedNA]) * len_val
            elif transl_number == 1:
                modified_val = unifiedNA
            else:
                modified_val = tuple([unifiedNA]) * transl_number

        else:
            # modify 1-length tuple
            if collapse_tuple:
                if isinstance(val, tuple) and (len(val) == 1):
                    modified_val = val[0]
                else:
                    modified_val = val
            else:
                modified_val = val

            # turn NAvalues into unifiedNA
            if isinstance(modified_val, tuple):
                modified_val = tuple( 
                        unifiedNA if x in NA_VALUES else x
                        for x in modified_val 
                        )
            else:
                if modified_val in NA_VALUES:
                    modified_val = unifiedNA

        return modified_val

    def do_typeconv(modified_val, keytype):
        if keytype == 'Flag':
            return bool(modified_val)

        else:
            if keytype == 'String' or keytype == 'Character':
                converter = str
            elif keytype == 'Integer':
                converter = int
            elif keytype == 'Float':
                converter = float

            def subfun(val):
                return val if val in NA_VALUES else converter(val)

            if isinstance(modified_val, tuple):
                new_val = tuple( subfun(val) for val in modified_val )
            else:
                new_val = subfun(modified_val)

            return new_val

    def main(key, val, transl_number, keytype, for_write, typeconv):
        assert transl_number == '.' or isinstance(transl_number, int)

        stop, return_val = handle_special_cases(key, val, keytype)
        if stop:
            return return_val

        len_val, isNA, unifiedNA = get_params(val, for_write, keytype)
        sanity_check(key, val, transl_number, keytype, len_val, isNA)
        modified_val = get_modified_val(val, transl_number, isNA, len_val, 
                                        unifiedNA, collapse_tuple)
        if typeconv:
            modified_val = do_typeconv(modified_val, keytype)

        return modified_val

    return main(key, val, transl_number, keytype, for_write, typeconv)


