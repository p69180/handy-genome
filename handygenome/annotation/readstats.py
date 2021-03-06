import os
import warnings
import re
import logging
import copy

import pysam
import numpy as np
import scipy.stats

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
workflow = importlib.import_module('.'.join([top_package_name, 'workflow']))
annotationdb = importlib.import_module('.'.join([top_package_name, 'annotation', 'annotationdb']))
readplus = importlib.import_module('.'.join([top_package_name, 'readplus', 'readplus']))
alleleinfosetup = importlib.import_module('.'.join([top_package_name, 'readplus', 'alleleinfosetup']))


READSTATS_FORMAT_KEY = 'readstats'


class ReadStats(annotationdb.AnnotItem):
    meta = {'ID': 'readstats', 'Number': 1, 'Type': 'String',
            'Description': 'Read-based statistics for the variant.'}

    def get_alleleindexes_mean(self, key, alleleindexes):
        numerator = sum(self[key][x] * self['rppcounts'][x]
                        for x in alleleindexes)
        denominator = sum(self['rppcounts'][x] for x in alleleindexes)

        if denominator == 0:
            return np.nan
        else:
            return numerator / denominator

    def get_total_rppcount(self):
        return sum(val for (key, val) in self['rppcounts'].items()
                   if isinstance(key, int))

    def postprocess(self):
        if 'rppcounts' not in self:
            tmp = dict()
            tmp[0] = self['ref_rppcount']
            for idx, val in enumerate(self['alt_rppcounts']):
                tmp[idx + 1] = val
            tmp[-1] = self['other_rppcount']
            tmp[None] = self['none_rppcount']
            tmp['softclip_overlap'] = self['clipoverlap_rppcount']
            self['rppcounts'] = tmp

        if 'mean_BQs' not in self:
            tmp = dict()
            tmp[0] = self['ref_mean_BQ']
            for idx, val in enumerate(self['alt_mean_BQs']):
                tmp[idx + 1] = val
            self['mean_BQs'] = tmp

        if 'mean_MQs' not in self:
            tmp = dict()
            tmp[0] = self['ref_mean_MQ']
            for idx, val in enumerate(self['alt_mean_MQs']):
                tmp[idx + 1] = val
            self['mean_MQs'] = tmp

        if 'mean_cliplens' not in self:
            tmp = dict()
            tmp[0] = self['ref_mean_cliplen']
            for idx, val in enumerate(self['alt_mean_cliplens']):
                tmp[idx + 1] = val
            self['mean_cliplens'] = tmp

        self['total_rppcount'] = sum(self['rppcounts'].values())


def get_readstats(vcfspec, bam, fasta, chromdict, no_matesearch=True, 
                  **kwargs):
    readstats_data = get_readstats_data(vcfspec, bam, fasta, chromdict, 
                                        no_matesearch=no_matesearch,
                                        **kwargs)
    result = summarize_readstats_data(readstats_data)

    return result


def rpplist_to_readstats_data(rpplist, vcfspec, 
                              flanklen=alleleinfosetup.DEFAULT_FLANKLEN):
    def add_MQ_rpp(rpp, data, alleleclass_rpp):
        if rpp.rp2 is None:
            MQ = rpp.rp1.read.mapping_quality
        else:
            MQ = (rpp.rp1.read.mapping_quality 
                  + rpp.rp2.read.mapping_quality) / 2

        data['MQ'][alleleclass_rpp].append(MQ)

    def add_BQ_rp(rp, data, alleleclass_rpp, vcfspec):
        data['BQ'][alleleclass_rpp].extend(rp.get_BQlist(vcfspec))

    def add_cliplen_rpp(rpp, data, alleleclass_rpp):
        if rpp.rp2 is None:
            cliplen = rpp.rp1.cigarstats[4]
        else:
            cliplen = rpp.rp1.cigarstats[4] + rpp.rp2.cigarstats[4]

        data['cliplen'][alleleclass_rpp].append(cliplen)

    def add_var_pos0s_rp(rp, data, alleleclass_rpp, vcfspec):
        var_querypos0s = rp.get_querypos0_of_range0_allmodes(
            vcfspec.REF_range0)
        if var_querypos0s['left'] is not None:
            data['pos0_left'][alleleclass_rpp].append(
                var_querypos0s['left'])
            data['pos0_right'][alleleclass_rpp].append(
                var_querypos0s['right'])
            data['pos0_5prime'][alleleclass_rpp].append(
                var_querypos0s['5prime'])
            data['pos0_3prime'][alleleclass_rpp].append(
                var_querypos0s['3prime'])

            data['pos0_left_fraction'][alleleclass_rpp].append(
                var_querypos0s['left_fraction'])
            data['pos0_right_fraction'][alleleclass_rpp].append(
                var_querypos0s['right_fraction'])
            data['pos0_5prime_fraction'][alleleclass_rpp].append(
                var_querypos0s['5prime_fraction'])
            data['pos0_3prime_fraction'][alleleclass_rpp].append(
                var_querypos0s['3prime_fraction'])

    def add_pairorient_rpp(rpp, data, alleleclass_rpp):
        data['pairorient'][alleleclass_rpp].append(rpp.pairorient)

    # initialize
    alleleclass_keys = (None,) + tuple(range(-1, len(vcfspec.alts) + 1))
    data = dict()
    # fields initialized as integer
    data['count'] = {x: 0 for x in alleleclass_keys}
    data['count']['softclip_overlap'] = 0
    # fields initialized as list
    for key in ('MQ', 'BQ', 'cliplen', 'pairorient',
                'pos0_left', 'pos0_right', 'pos0_5prime', 'pos0_3prime',
                'pos0_left_fraction', 'pos0_right_fraction', 
                'pos0_5prime_fraction', 'pos0_3prime_fraction'):
        data[key] = {x: list() for x in alleleclass_keys}
    # add data
    for rpp in rpplist:
        if vcfspec not in rpp.alleleinfo:
            rpp.update_alleleinfo(vcfspec=vcfspec, flanklen=flanklen)
        alleleclass_rpp = rpp.alleleinfo[vcfspec]['alleleclass']
        # count
        data['count'][alleleclass_rpp] += 1
        if alleleclass_rpp is None:
            if rpp.check_softclip_overlaps_vcfspec(vcfspec):
                data['count']['softclip_overlap'] += 1
        #if alleleclass_rpp is not None:
        # MQ
        add_MQ_rpp(rpp, data, alleleclass_rpp)
        # BQ
        add_BQ_rp(rpp.rp1, data, alleleclass_rpp, vcfspec)
        if rpp.rp2 is not None:
            add_BQ_rp(rpp.rp2, data, alleleclass_rpp, vcfspec)
        # cliplen
        add_cliplen_rpp(rpp, data, alleleclass_rpp)
        # var_pos0s
        add_var_pos0s_rp(rpp.rp1, data, alleleclass_rpp, vcfspec)
        if rpp.rp2 is not None:
            add_var_pos0s_rp(rpp.rp2, data, alleleclass_rpp, vcfspec)
        # pairorient
        add_pairorient_rpp(rpp, data, alleleclass_rpp)
    # turn into numpy arrays
    for key in data.keys():
        if key != 'count':
            data[key] = {subkey: np.array(subval)
                         for subkey, subval in data[key].items()}

    return data


def get_readstats_data(
        vcfspec, bam, fasta, chromdict, 
        no_matesearch=True,
        fetch_padding_common=readplus.FETCH_PADDING_COMMON,
        #fetch_padding_sv=readplus.FETCH_PADDING_SV,
        new_fetch_padding=readplus.NEW_FETCH_PADDING,
        long_insert_threshold=readplus.LONG_INSERT_THRESHOLD,
        flanklen=alleleinfosetup.DEFAULT_FLANKLEN,
        #flanklen_parside=alleleinfosetup_sv.DEFAULT_FLANKLEN_PARSIDE,
        #flanklen_bndside=alleleinfosetup_sv.DEFAULT_FLANKLEN_BNDSIDE,
        ):
    """Only for non-sv cases"""

    rpplist = readplus.get_rpplist_nonsv(
        bam=bam, fasta=fasta, chromdict=chromdict, 
        chrom=vcfspec.chrom, start0=vcfspec.pos0, 
        end0=vcfspec.end0, view=False, 
        no_matesearch=no_matesearch,
        fetch_padding_common=fetch_padding_common,
        new_fetch_padding=new_fetch_padding,
        long_insert_threshold=long_insert_threshold)
    rpplist.update_alleleinfo(
        vcfspec=vcfspec, flanklen=flanklen)
    readstats_data = rpplist_to_readstats_data(rpplist, vcfspec)

    return readstats_data


def summarize_readstats_data(readstats_data):
    def ref_mean(data):
        return np.mean(data[0])

    def alt_means(data, alt_keys):
        return [np.mean(data[x]) for x in alt_keys]

    def allele_means(data):
        return dict((alleleclass, np.mean(array)) 
                    for (alleleclass, array) in data.items())

    def allele_medians(data):
        return dict((alleleclass, np.median(array)) 
                    for (alleleclass, array) in data.items())

    def varpos_kstest(data):
        summary = dict()
        for alleleclass, array in data.items():
            if len(array) == 0:
                pval = np.nan
            else:
                pval = scipy.stats.kstest(array, 'uniform').pvalue
            summary[alleleclass] = pval

        return summary

    def setup_old(readstats, readstats_data):
        alt_keys = sorted(x for x in readstats_data['count'].keys()
                          if isinstance(x, int) and (x > 0))

        readstats['ref_rppcount'] = readstats_data['count'][0]
        readstats['alt_rppcounts'] = [readstats_data['count'][x] 
                                      for x in alt_keys]
        readstats['other_rppcount'] = readstats_data['count'][-1]
        readstats['none_rppcount'] = readstats_data['count'][None]
        readstats['clipoverlap_rppcount'] = (
            readstats_data['count']['softclip_overlap'])

        readstats['ref_mean_BQ'] = ref_mean(readstats_data['BQ'])
        readstats['alt_mean_BQs'] = alt_means(readstats_data['BQ'], alt_keys)

        readstats['ref_mean_MQ'] = ref_mean(readstats_data['MQ'])
        readstats['alt_mean_MQs'] = alt_means(readstats_data['MQ'], alt_keys)

        readstats['ref_mean_cliplen'] = ref_mean(readstats_data['cliplen'])
        readstats['alt_mean_cliplens'] = alt_means(readstats_data['cliplen'], 
                                                   alt_keys)

    def setup_new(readstats, readstats_data):
        readstats['rppcounts'] = readstats_data['count'].copy()

        readstats['mean_BQs'] = allele_means(readstats_data['BQ'])
        readstats['median_BQs'] = allele_medians(readstats_data['BQ'])

        readstats['mean_MQs'] = allele_means(readstats_data['MQ'])
        readstats['median_MQs'] = allele_medians(readstats_data['MQ'])

        readstats['mean_cliplens'] = allele_means(readstats_data['cliplen'])
        readstats['median_cliplens'] = allele_medians(readstats_data['cliplen'])

        readstats['varpos_uniform_pvalues'] = varpos_kstest(readstats_data['pos0_left_fraction'])
        readstats['mean_varpos_fractions'] = allele_means(readstats_data['pos0_left_fraction'])
        readstats['median_varpos_fractions'] = allele_medians(readstats_data['pos0_left_fraction'])

    # main
    readstats = ReadStats()
    setup_new(readstats, readstats_data)

    return readstats


def add_meta(vcfheader):
    vcfheader.add_meta(
        key='FORMAT',
        items=[('ID', READSTATS_FORMAT_KEY),
               ('Type', 'String'),
               ('Number', 1),
               ('Description', 
                'Read information statistics for the sample.')])
    

