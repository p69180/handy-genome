import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
readplus = importlib.import_module('.'.join([top_package_name, 'readplus', 'readplus']))


DEFAULT_FLANKLEN = 1


def make_alleleinfoitem_readpluspair(vcfspec, aiitem_rp1, aiitem_rp2):
    alleleinfoitem = dict()

    if aiitem_rp1['spans'] and aiitem_rp2['spans']:
        if aiitem_rp1['alleleclass'] == aiitem_rp2['alleleclass']:
            alleleinfoitem['alleleclass'] = aiitem_rp1['alleleclass']
        else:
            alleleinfoitem['alleleclass'] = None
    else:
        if (not aiitem_rp1['spans']) and (not aiitem_rp2['spans']):
            alleleinfoitem['alleleclass'] = None
        else:
            if aiitem_rp1['spans']:
                alleleinfoitem['alleleclass'] = aiitem_rp1['alleleclass']
            elif aiitem_rp2['spans']:
                alleleinfoitem['alleleclass'] = aiitem_rp2['alleleclass']

    return alleleinfoitem


def make_alleleinfoitem_readplus(vcfspec, rp, flanklen=DEFAULT_FLANKLEN):
    assert flanklen >= 1, f'"flanklen" argument must be at least 1.'

    preflank_range0, postflank_range0 = get_flank_ranges(vcfspec, flanklen)
    spans = check_read_spans_allele(rp.read, vcfspec, preflank_range0, 
                                    postflank_range0)

    if spans:
        pairs_indexes_pre = rp.get_pairs_indexes(preflank_range0)
        pairs_indexes_post = rp.get_pairs_indexes(postflank_range0)
        pairs_indexes_site = rp.get_pairs_indexes(vcfspec.get_REF_range0())

        allele_seq = rp.get_seq_from_pairs_indexes(pairs_indexes_site)
        alleleclass = get_alleleclass(vcfspec, rp, allele_seq, 
                                      pairs_indexes_pre, pairs_indexes_post)
    else:
        allele_seq = None
        alleleclass = None

    alleleinfoitem = dict()
    alleleinfoitem['spans'] = spans
    alleleinfoitem['allele_seq'] = allele_seq
    alleleinfoitem['alleleclass'] = alleleclass

    return alleleinfoitem


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


def check_matches(rp, pairs_indexes):
    for pairs_idx in pairs_indexes:
        querypos0 = rp.pairs_dict['querypos0'][pairs_idx]
        refseq = rp.pairs_dict['refseq'][pairs_idx]
        matches = ((querypos0 is not None) and 
                   (refseq is not None) and
                   refseq.isupper())
        if not matches:
            return False
    return True


def check_read_spans_allele(read, vcfspec, preflank_range0, postflank_range0):
    """Returns:
        True if the read spans the flanking bases around the variable region
    """

    return (read.reference_name == vcfspec.chrom and
            read.reference_start <= preflank_range0.start and
            read.reference_end >= postflank_range0.stop)


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

    if (
            rp.check_matches(pairs_indexes_pre) and
            rp.check_matches(pairs_indexes_post)):
        alleles = (vcfspec.ref,) + vcfspec.alts

        if allele_seq in alleles:
            alleleclass = alleles.index(allele_seq)
        else:
            alleleclass = -1
    else:
        alleleclass = -1

    return alleleclass
