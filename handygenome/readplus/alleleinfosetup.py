import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
readplus = importlib.import_module('.'.join([top_package_name, 'readplus', 'readplus']))


def make_alleleinfoitem_readpluspair(vcfspec, aiitem_rp1, aiitem_rp2):
    result = dict()

    if aiitem_rp1['spans'] and aiitem_rp2['spans']:
        if aiitem_rp1['alleleclass'] == aiitem_rp2['alleleclass']:
            result['alleleclass'] = aiitem_rp1['alleleclass']
        else:
            result['alleleclass'] = None
    else:
        if (not aiitem_rp1['spans']) and (not aiitem_rp2['spans']):
            result['alleleclass'] = None
        else:
            if aiitem_rp1['spans']:
                result['alleleclass'] = aiitem_rp1['alleleclass']
            elif aiitem_rp2['spans']:
                result['alleleclass'] = aiitem_rp2['alleleclass']

    return result


def make_alleleinfoitem_readplus(vcfspec, rp, 
                                 flanklen=readplus.DEFAULT_FLANKLEN):
    assert flanklen >= 1, f'"flanklen" argument must be at least 1.'

    preflank_range0, postflank_range0 = get_flank_ranges(vcfspec, flanklen)
    spans = check_read_spans_allele(rp.read, vcfspec, preflank_range0, 
                                    postflank_range0)

    if spans:
        (pairs_indexes_pre, 
         pairs_indexes_post) = get_pairs_indexes(rp, preflank_range0, 
                                                 postflank_range0)
        allele_seq = get_allele_seq(rp, pairs_indexes_pre, pairs_indexes_post)
        alleleclass = get_alleleclass(vcfspec, rp, allele_seq, 
                                      pairs_indexes_pre, pairs_indexes_post)
    else:
        allele_seq = None
        alleleclass = None

    result = dict()
    result['spans'] = spans
    result['allele_seq'] = allele_seq
    result['alleleclass'] = alleleclass

    return result


def get_flank_ranges(vcfspec, flanklen):
    preflank_range0_candidates = set()
    for idx in range(len(vcfspec.alts)):
        preflank_range0_candidates.add(
            vcfspec.get_preflank_range0(idx=idx, flanklen=flanklen))

    if len(preflank_range0_candidates) != 1:
        raise Exception(
            f'There are more than one possible pre-flanking ranges '
            f'for the input vcfspec: {vcfspec}')
    else:
        preflank_range0 = preflank_range0_candidates.pop()

    postflank_range0 = vcfspec.get_postflank_range0(flanklen=flanklen)

    return preflank_range0, postflank_range0


def get_pairs_indexes(rp, preflank_range0, postflank_range0):
    pairs_indexes_pre = [(rp.pairs_dict['refpos0'].index(x)
                          if x in rp.pairs_dict['refpos0'] else
                          None) 
                         for x in preflank_range0]
    pairs_indexes_post = [(rp.pairs_dict['refpos0'].index(x)
                           if x in rp.pairs_dict['refpos0'] else
                           None) 
                          for x in postflank_range0]

    return pairs_indexes_pre, pairs_indexes_post


def check_read_spans_allele(read, vcfspec, preflank_range0, postflank_range0):
    """Returns:
        True if the read spans the flanking bases around the variable region
    """

    return (read.reference_name == vcfspec.chrom and
            read.reference_start <= preflank_range0.start and
            read.reference_end >= postflank_range0.stop)


def check_rp_matches_flanks(rp, pairs_indexes_pre, pairs_indexes_post):
    """Returns:
        True if the read spans the base immediately before & after the allele
            and matches with identical base on those positions
    """

    def check_matches(rp, pairs_indexes):
        for pairs_idx in pairs_indexes:
            querypos0 = rp.pairs_dict['querypos0'][pairs_idx]
            refseq = rp.pairs_dict['refseq'][pairs_idx]
            matches = (querypos0 is not None) and refseq.isupper()
            if not matches:
                return False
        return True

    matches_pre = check_matches(rp, pairs_indexes_pre)
    matches_post = check_matches(rp, pairs_indexes_post)

    return matches_pre and matches_post


def get_allele_seq(rp, pairs_indexes_pre, pairs_indexes_post):
    slice_start = pairs_indexes_pre[-1] + 1
    slice_end = pairs_indexes_post[0]

    allele_querypos0_list = [
        x for x in rp.pairs_dict['querypos0'][slice_start:slice_end]
        if x is not None]

    if len(allele_querypos0_list) == 0:
        allele_seq = ''
    else:
        idx_start = min(allele_querypos0_list)
        idx_end = max(allele_querypos0_list) + 1
        allele_seq = rp.read.query_sequence[idx_start:idx_end]

    return allele_seq


def get_alleleclass(vcfspec, rp, allele_seq, pairs_indexes_pre, 
                    pairs_indexes_post):
    """Returns:
        None: not informative
        -1: other than REF and ALTS
        0: REF
        1: 1st ALT
        2: 2nd ALT
        ...
    """

    if check_rp_matches_flanks(rp, pairs_indexes_pre, pairs_indexes_post):
        alleles = (vcfspec.ref,) + vcfspec.alts

        if allele_seq in alleles:
            alleleclass = alleles.index(allele_seq)
        else:
            alleleclass = -1
    else:
        alleleclass = -1

    return alleleclass
