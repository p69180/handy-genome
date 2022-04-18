import os
import warnings
import re

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
infoformat = importlib.import_module('.'.join([top_package_name, 'variantplus', 'infoformat']))
equivalents = importlib.import_module('.'.join([top_package_name, 'variantplus', 'equivalents']))
varianthandler = importlib.import_module('.'.join([top_package_name, 'variantplus', 'varianthandler']))
annotationdb = importlib.import_module('.'.join([top_package_name, 'annotation', 'annotationdb']))
ensembl_parser = importlib.import_module('.'.join([top_package_name, 'annotation', 'ensembl_parser']))


DEFAULT_FETCHWIDTH = 10

PAT_GENESET_ID = re.compile('ID=(?P<type>[^:]+):(?P<id>[^;]+)')
    # type: CDS, gene, transcript (obtained from: zcat geneset_gff3_sorted.gz | awk '{if ($9 ~ /ID=([^:]+):/) {print gensub(/^.*ID=([^:]+):.*$/, "\\1", "g", $9)}}' | sort | uniq)
PAT_GENESET_PARENT = re.compile('Parent=[^:]+:(?P<id>[^;]+)')

PAT_REGULATORY_ID = re.compile('(.+;)?id=([^;]+)(;.+)?')


# VCF


class NoVcfIndexError(Exception):
    pass


def fetch_relevant_vr(vcfspec, vcf, fasta = None, search_equivs = True):
    """
    Args:
        fasta: Not required when 'search_equivs' is False
        search_equivs: If True, all equivalent forms are fetched. If False, 
            only vrs with identical chrom, pos, ref, alt are fetched.
    """

    assert not ((fasta is None) and search_equivs), \
        f'If "search_equivs" is True, "fasta" must be given.'

    def search_equivalents(vcfspec, vcf, fasta):
        mttype = common.get_mttype(vcfspec.ref, vcfspec.alt)

        if mttype in ('ins', 'del'):
            equivs = equivalents.indel_equivalents(vcfspec, fasta)
            poslist = [x.pos for x in equivs]
            try:
                fetchresult = vcf.fetch(contig = vcfspec.chrom, 
                                        start = poslist[0] - 1,
                                        end = poslist[-1])
            except ValueError as e:
                warnings.warn(f'{__name__}.fetch_relevant_vr: {str(e)}')
                relevant_vr_candidates = list()
            else:
                relevant_vr_candidates = list(filter(
                    lambda vr: varianthandler.get_vcfspec(vr) in equivs,
                    fetchresult))
        else:
            relevant_vr_candidates = search_identical(vcfspec, vcf)

        return relevant_vr_candidates

    def search_identical(vcfspec, vcf):
        try:
            fetchresult = vcf.fetch(contig = vcfspec.chrom,
                                    start = vcfspec.pos - 1,
                                    end = vcfspec.pos)
        except ValueError as e:
            warnings.warn(f'{__name__}.fetch_relevant_vr: {str(e)}')
            relevant_vr_candidates = list()
        else:
            relevant_vr_candidates = list(filter(
                lambda vr: varianthandler.get_vcfspec(vr) == vcfspec,
                fetchresult))
    
        return relevant_vr_candidates

    def get_relevant_vr(relevant_vr_candidates):
        if len(relevant_vr_candidates) == 0:
            relevant_vr = None
        elif len(relevant_vr_candidates) == 1:
            relevant_vr = relevant_vr_candidates[0]
        else:
            e_msg = list()
            e_msg.append(f'There are more than one relevant variant records:')
            for vr in relevant_vr_candidates:
                e_msg.append(vr.to_string())
            raise Exception('\n'.join(e_msg))

        return relevant_vr

    def main(vcfspec, vcf, fasta, search_equivs):
        if vcf.index is None:
            raise NoVcfIndexError(f'Input vcf is not indexed.')

        if search_equivs:
            relevant_vr_candidates = search_equivalents(vcfspec, vcf, fasta)
        else:
            relevant_vr_candidates = search_identical(vcfspec, vcf)

        relevant_vr = get_relevant_vr(relevant_vr_candidates)

        return relevant_vr

    return main(vcfspec, vcf, fasta, search_equivs)


def fetch_relevant_vr_vcfpath(vcfspec, vcf_path, fasta=None, 
                              search_equivs=True):
    with pysam.VariantFile(vcf_path, 'r') as vcf:
        try:
            relevant_vr = fetch_relevant_vr(vcfspec, vcf, fasta = fasta, 
                                            search_equivs = search_equivs)
        except NoVcfIndexError:
            print(f'Input vcf path: {vcf_path}')
            raise

    return relevant_vr


# cosmic


# geneset gff3

def parse_transcript_tabixline(tabixline):
    annotitem = annotationdb.AnnotItem()

    attrs_parsed = parse_gff3_attrs(tabixline.attributes)
    annotitem['id'] = attrs_parsed['ID'].split(':')[1]

    ensembl_parser.set_feature_type(annotitem, 'transcript')
    annotitem['biotype'] = attrs_parsed['biotype']
    ensembl_parser.set_transcript_subtypes(annotitem, annotitem['biotype'], annotitem['id'])

    if tabixline.strand == '+':
        annotitem['is_forward'] = True
    elif tabixline.strand == '-':
        annotitem['is_forward'] = False
    else:
        raise Exception(f'Unexpected transcript tabixline strand value: "{tabixline.strand}"')

    annotitem['chrom'] = tabixline.contig
    annotitem['start0'] = tabixline.start
    annotitem['start1'] = annotitem['start0'] + 1
    annotitem['end0'] = tabixline.end
    annotitem['end1'] = annotitem['end0']

    annotitem['transcript_name'] = attrs_parsed['Name']
    annotitem['gene_id'] = attrs_parsed['Parent'].split(':')[1]

    return annotitem


def fetch_transcript(chrom, start0, end0, tabixfile_geneset):
    transcript_dict = annotationdb.AnnotItemDict()
    for tabixline in tabixfile_geneset.fetch(chrom, start0, end0):
        if check_tabixline_type_transcript(tabixline):
            annotitem = parse_transcript_tabixline(tabixline)
            transcript_dict[annotitem['id']] = annotitem

    return transcript_dict


def fetch_transcript_tabixline(chrom, start0, end0, transcript_id_list, tabixfile_geneset):
    candidates = list()
    fetched = tabixfile_geneset.fetch(chrom, start0, end0)
    for tabixline in fetched:
        mat = PAT_GENESET_ID.search(tabixline.attributes)
        if mat is not None:
            if mat.group('type') == 'transcript':
                if mat.group('id') in transcript_id_list:
                    candidates.append(tabixline)

    if len(candidates) != len(transcript_id_list):
        raise Exception(f'Unsuccessful transcript fetch: transcript_id_list: {transcript_id_list}, coords: ({chrom}, {start0}, {end0})')

    return candidates


def check_geneset_tabixline_type(tabixline, type):
    assert type in ('gene', 'transcript', 'CDS', 'exon', 'UTR'), \
        f"Allowed 'type' values are: 'gene', 'transcript', 'CDS', 'exon', 'UTR'"

    if type == 'exon':
        return tabixline.feature == 'exon'
    elif type == 'UTR':
        return tabixline.feature in ('five_prime_UTR', 'three_prime_UTR')
    else:
        mat = PAT_GENESET_ID.search(tabixline.attributes)
        if mat is None:
            return False
        else:
            return (mat.group('type') == type)


def check_tabixline_type_transcript(tabixline):
    mat = PAT_GENESET_ID.search(tabixline.attributes)
    if mat is None:
        return False
    else:
        return (mat.group('type') == 'transcript')


def get_parent(tabixline):
    mat = PAT_GENESET_PARENT.search(tabixline.attributes)
    if mat is None:
        return None
    else:
        return mat.group('id')


# repeats bed

def parse_repeat_tabixline(tabixline):
    annotitem = annotationdb.AnnotItem()

    annotitem['chrom'] = tabixline.contig
    annotitem['start0'] = tabixline.start
    annotitem['end0'] = tabixline.end
    annotitem['start1'] = annotitem['start0'] + 1
    annotitem['end1'] = annotitem['end0']

    ensembl_parser.set_feature_type(annotitem, 'repeat')

    raw_split = tabixline.name.split(':')
    clsfml = raw_split[0].split('/')
    if len(clsfml) == 1:
        annotitem['class'] = clsfml[0]
        annotitem['family'] = clsfml[0]
    elif len(clsfml) == 2:
        annotitem['class'] = clsfml[0]
        annotitem['family'] = clsfml[1]
    else:
        raise Exception(f'Unexpected repeat name pattern: {tabixline.name}')
    annotitem['name'] = raw_split[1]

    annotitem['id'] = f'{annotitem["name"]}_{tabixline.contig}_{tabixline.start}_{tabixline.end}'

    return annotitem


def fetch_repeat(chrom, start0, end0, tabixfile_repeats):
    repeat_dict = annotationdb.AnnotItemDict()
    for tabixline in tabixfile_repeats.fetch(chrom, start0, end0):
        annotitem = parse_repeat_tabixline(tabixline)
        repeat_dict[annotitem['id']] = annotitem

    return repeat_dict


# regulatory gff3

def parse_regulatory_tabixline(tabixline):
    annotitem = annotationdb.AnnotItem()

    annotitem['chrom'] = tabixline.contig
    annotitem['start0'] = tabixline.start
    annotitem['end0'] = tabixline.end
    annotitem['start1'] = annotitem['start0'] + 1
    annotitem['end1'] = annotitem['end0']

    ensembl_parser.set_feature_type(annotitem, 'regulatory')
    ensembl_parser.set_regulatory_subtypes(annotitem, tabixline.feature, 'vep')

    attrs_parsed = parse_gff3_attrs(tabixline.attributes)
    annotitem['id'] = attrs_parsed['id']
    annotitem['activity'] = dict(
        x.split(':') for x in attrs_parsed['activity'].split(','))

    annotitem['bound_start1'] = int(attrs_parsed['bound_start'])
    annotitem['bound_end1'] = int(attrs_parsed['bound_end'])
    annotitem['bound_start0'] = annotitem['bound_start1'] - 1
    annotitem['bound_end0'] = annotitem['bound_end1']

    return annotitem


def fetch_regulatory_tabixline(chrom, start0, end0, regulatory_id_list, tabixfile_regulatory):
    candidates = list()
    fetched = tabixfile_regulatory.fetch(chrom, start0, end0)
    for tabixline in fetched:
        ID = get_ID_regulatory_tabixline(tabixline)
        if ID in regulatory_id_list:
            candidates.append(tabixline)

    if len(candidates) != len(regulatory_id_list):
        raise Exception(f'Unsuccessful regulatory fetch: regulatory_id_list: {regulatory_id_list}, coords: ({chrom}, {start0}, {end0})')

    return candidates


def fetch_regulatory(chrom, start0, end0, tabixfile_regulatory):
    regulatory_dict = annotationdb.AnnotItemDict()
    for tabixline in tabixfile_regulatory.fetch(chrom, start0, end0):
        annotitem = parse_regulatory_tabixline(tabixline)
        regulatory_dict[annotitem['id']] = annotitem

    return regulatory_dict


def get_ID_regulatory_tabixline(tabixline):
    return PAT_REGULATORY_ID.search(tabixline.attributes).group(2)

# misc.

def parse_gff3_attrs(raw_string):
    return dict(x.split('=') for x in raw_string.split(';'))

