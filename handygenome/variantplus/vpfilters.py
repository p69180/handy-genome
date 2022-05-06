import re

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))


def get_transcript_subtype_filter(key):
    def vpfilter(vp):
        return any(feature['transcript_subtype_flags'][key]
                   for feature in vp.annotdb.transcript_canon_ovlp.values())

    return vpfilter

CODING_GENE_INVOLVED = get_transcript_subtype_filter('is_coding')


def get_consequence_filter(key):
    def vpfilter(vp):
        return any(feature['consequence_flags'][key]
                   for feature in vp.annotdb.transcript_canon_ovlp.values())

    return vpfilter

PROTEIN_ALTERING = get_consequence_filter('is_protein_altering')
NOT_PROTEIN_ALTERING = get_consequence_filter('is_not_protein_altering')
    

def get_genename_filter(gene_list, canonical_only=True):
    gene_set = set(gene_list)
    if canonical_only:
        def vpfilter(vp):
            vp_genes = set(vp.get_gene_names(canonical_only=True))
            return not vp_genes.isdisjoint(gene_set)
    else:
        def vpfilter(vp):
            vp_genes = set(vp.get_gene_names(canonical_only=False))
            return not vp_genes.isdisjoint(gene_set)

    return vpfilter
