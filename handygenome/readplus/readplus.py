import itertools
import pprint

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
workflow = importlib.import_module('.'.join([top_package_name, 'workflow']))
readhandler = importlib.import_module('.'.join([top_package_name, 'readplus', 'readhandler']))
alleleinfosetup = importlib.import_module('.'.join([top_package_name, 'readplus', 'alleleinfosetup']))


RPPLIST_MAX_SHOW_LEN = 30
THRESHOLD_TEMPLATE_LENGTH = 1000

FETCH_PADDING_COMMON = 0
FETCH_PADDING_SV = 850
FETCH_PADDING_VIEW = 1000
NEW_FETCH_PADDING = 0
LONG_INSERT_THRESHOLD = 1000

DEFAULT_FLANKLEN = 1

LOGGER_RPPLIST = workflow.get_logger(name='ReadPlusPairList', 
                                     level='warning')


class ReadPlus:

    def __init__(self, read, fasta=None):
        """
        These attributes are set regardless of whether ria functions are used.
        """

        self.read = read
        self.fasta = (common.infer_refver(bamheader=self.read.header)
                      if fasta is None else fasta)

        self._set_NMMD()
        self._set_SAlist()

        #self.aligned_pairs = self.read.get_aligned_pairs(with_seq=True)
        self.pairs_dict = readhandler.get_pairs_dict(self.read, self.fasta)
            # keys: 'querypos0', 'refpos0','refseq'

        self.fiveprime_end = readhandler.get_fiveprime_end(self.read)
        self.threeprime_end = readhandler.get_threeprime_end(self.read)
        #self.cigarstats = self.read.get_cigar_stats()[0] 
            # length 11 array which contains value for M, I, D, N, S, ...
        self.ref_range0 = range(self.read.reference_start, 
                                self.read.reference_end)
        #self.MQ = self.read.mapping_quality

        #self.irrelevant_cause = list()
        #self.irrelevant = False

        self.alleleinfo = dict()

    def update_alleleinfo(self, vcfspec, flanklen=DEFAULT_FLANKLEN):
        self.alleleinfo[vcfspec] = (
            alleleinfosetup.make_alleleinfoitem_readplus(vcfspec, self, 
                                                         flanklen))

    def _set_NMMD(self):
        if (not self.read.has_tag('NM')) or (not self.read.has_tag('MD')):
            (ref_seq_padded, read_seq_padded) = readhandler.get_padded_seqs(
                self.read, self.fasta)
            if not self.read.has_tag('NM'):
                self.read.set_tag(
                    tag='NM', value_type='i',
                    value=readhandler.get_NM(
                        self.read, ref_seq_padded=ref_seq_padded, 
                        read_seq_padded=read_seq_padded))
            if not self.read.has_tag('MD'):
                self.read.set_tag(
                    tag='MD', value_type='Z',
                    value=readhandler.get_MD(
                        self.read, ref_seq_padded=ref_seq_padded, 
                        read_seq_padded=read_seq_padded))

    def _set_SAlist(self):
        """cigartuples pattern check is done:
            asserts that one end is softclip and the other end is match.
        """

        def SA_cigartuples_sanitycheck(cigartuples, read):
            if not (
                    (cigartuples[0][0] == 4 and cigartuples[-1][0] == 0) or
                    (cigartuples[-1][0] == 4 and cigartuples[0][0] == 0)):
                raise Exception(
                    f'Unexpected SA cigarstring pattern:\n{read.to_string()}')

        if self.read.has_tag('SA'):
            self.SAlist = list()
            SAtag_split = self.read.get_tag('SA').strip(';').split(';')
            for field in SAtag_split:
                field_sp = field.split(',')

                # 
                MQ = int(field_sp[4])
                #if MQ == 0:
                #    continue

                chrom = field_sp[0]
                pos = int(field_sp[1])

                if field_sp[2] == '+':
                    is_forward = True
                elif field_sp[2] == '-':
                    is_forward = False
                else:
                    raise Exception(
                        f'Unexpected strand string from SA tag. read:\n'
                        f'{self.read.to_string()}')

                cigarstring = field_sp[3]
                cigartuples = readhandler.get_cigartuples(cigarstring)
                #SA_cigartuples_sanitycheck(cigartuples, self.read)

                NM = int(field_sp[5])

                SAitem = {'chrom': chrom, 'pos': pos, 
                          'is_forward': is_forward,
                          'MQ': MQ, 'NM': NM,
                          'cigarstring': cigarstring, 
                          'cigartuples': cigartuples}
                self.SAlist.append(SAitem)

            if len(self.SAlist) == 0:
                self.SAlist = None
        else:
            self.SAlist = None

#    def get_SA_seqs(self):
#        if self.SAlist is None:
#            SA_seqs = None
#        else:
#            SA_seqs = list()
#            for SAitem in self.SAlist:
#                first_cigarop, first_cigarlen = self.read.cigartuples[0]
#                first_is_candidate = (first_cigarop == 4 and
#                                      first_cigarlen == SAitem['seqlen'])
#                last_cigarop, last_cigarlen = self.read.cigartuples[0]
#                last_is_candidate = (last_cigarop == 4 and
#                                     last_cigarlen == SAitem['seqlen'])
#
#                if first_is_candidate and last_is_candidate:
#                    raise Exception(
#                        f'Both of the two softclips have the same length with '
#                        f'the SA tag in this read:\n'
#                        f'{self.read.to_string()}')
#                else:
#                    if first_is_candidate:
#                        SA_seq = self.read.query_sequence[:first_cigarlen]
#                    elif last_is_candidate:
#                        SA_seq = self.read.query_sequence[-last_cigarlen:]
#                    else:
#                        SA_seq = None
#                    SA_seqs.append(SA_seq)
#
#        return SA_seqs


    def _get_infostring(self):
        return (f'qname: {self.read.query_name}, '
                f'chrom: {self.read.reference_name}, '
                f'pos: {self.read.reference_start + 1:,}, '
                f'alleleinfo: {self.alleleinfo}')

    def __repr__(self):
        return f'<ReadPlus object ({self._get_infostring()})>'


class ReadPlusPair:
    """
    Requires two readplus objects which are mates of each other and all of 
        primary alignment.
    Non-primary reads with the same queryname is stored in 
        'rplist_nonprimary' attribute.
    """

    def __init__(self, rplist_primary, rplist_nonprimary, chromdict,
                 threshold_tlen=THRESHOLD_TEMPLATE_LENGTH):
        assert len(rplist_primary) == 2

        self._set_rp1_rp2(rplist_primary, chromdict)
        self.rplist_nonprimary = rplist_nonprimary

        self.query_name = self.rp1.read.query_name

        self.mate_chroms_differ = (
            self.rp1.read.reference_name != self.rp2.read.reference_name)
        self.is_TRA = self.mate_chroms_differ # alias

        self.pairorient = readhandler.get_pairorient(self.rp1.read) 
            # may be None when: mate unmapped, TLEN == 0
        self.tlen = (None 
                     if self.mate_chroms_differ else 
                     abs(self.rp1.read.template_length))
        self.template_length = self.tlen  # alias

        self.alleleinfo = dict()

        #self._set_is_proper_pair()
        #self._set_sv_supporting()
        #self.irrelevant = (self.rp1.irrelevant or self.rp2.irrelevant)

    def update_alleleinfo(self, vcfspec, flanklen=DEFAULT_FLANKLEN):
        # rp1
        if vcfspec not in self.rp1.alleleinfo:
            self.rp1.update_alleleinfo(vcfspec, flanklen=flanklen)
        aiitem_rp1 = self.rp1.alleleinfo[vcfspec]
        # rp2
        if vcfspec not in self.rp2.alleleinfo:
            self.rp2.update_alleleinfo(vcfspec, flanklen=flanklen)
        aiitem_rp2 = self.rp2.alleleinfo[vcfspec]

        self.alleleinfo[vcfspec] = (
            alleleinfosetup.make_alleleinfoitem_readpluspair(
                vcfspec, aiitem_rp1, aiitem_rp2))

    def __repr__(self):
        qname = self.rp1.read.query_name
        rp1_chrom = self.rp1.read.reference_name
        rp1_pos = self.rp1.read.reference_start + 1
        rp2_chrom = self.rp2.read.reference_name
        rp2_pos = self.rp2.read.reference_start + 1

        return (f'<ReadPlusPair object (qname: {qname}; '
                f'rp1_pos: {rp1_chrom}:{rp1_pos:,}; '
                f'rp2_pos: {rp2_chrom}:{rp2_pos:,}; '
                f'alleleinfo: {self.alleleinfo})>')

    def _set_rp1_rp2(self, rplist_primary, chromdict):
        order = common.compare_coords(
            rplist_primary[0].read.reference_name, 
            rplist_primary[0].fiveprime_end,
            rplist_primary[1].read.reference_name, 
            rplist_primary[1].fiveprime_end,
            chromdict)
        if order <= 0:
            self.rp1 = rplist_primary[0]
            self.rp2 = rplist_primary[1]
        else:
            self.rp1 = rplist_primary[1]
            self.rp2 = rplist_primary[0]

#    def _set_is_proper_pair(self):
#        if self.pairorient is None:
#            self.is_proper_pair = False
#        else:
#            self.is_proper_pair = (self.pairorient[0] == 'F' and 
#                                   self.pairorient[2] == 'R')

#    def _set_sv_supporting(self):
#        if self.mate_chroms_differ:
#            self.SV_supporting = True
#        else:
#            if self.tlen == 0:
#                self.SV_supporting = False
#            else:
#                if not self.is_proper_pair:
#                    self.SV_supporting = True
#                else:
#                    self.SV_supporting = (
#                        self.tlen > THRESHOLD_TEMPLATE_LENGTH)


# rpplist classes and functions

class ReadPlusPairList(list):

    def __init__(self):
        pass

    def update_alleleinfo(self, vcfspec, flanklen=DEFAULT_FLANKLEN):
        for rpp in self:
            rpp.update_alleleinfo(vcfspec, flanklen=flanklen)


def get_rpplist_nonsv(bam, fasta, chromdict, chrom, start0, end0, view=False,
                      fetch_padding_common=FETCH_PADDING_COMMON,
                      fetch_padding_view=FETCH_PADDING_VIEW,
                      new_fetch_padding=NEW_FETCH_PADDING,
                      long_insert_threshold=LONG_INSERT_THRESHOLD):
    LOGGER_RPPLIST.info('Beginning initial fetch')
    (relevant_qname_list, new_fetch_range) = initial_fetch_nonsv(
        bam, chrom, start0, end0, view,
        fetch_padding_common, fetch_padding_view,
        new_fetch_padding, long_insert_threshold)

    if relevant_qname_list is None:
        rpplist = None
    else:
        LOGGER_RPPLIST.info('Beginning refined fetch')
        fetchresult_dict = refined_fetch_nonsv(bam, chrom, new_fetch_range, 
                                               relevant_qname_list)

        LOGGER_RPPLIST.info('Beginning assembly into readpluspair')
        rpplist = ReadPlusPairList()
        for readlist in fetchresult_dict.values():
            rpp = get_rpp_from_refinedfetch(readlist, bam, fasta, chromdict)
            if rpp is not None:
                rpplist.append(rpp)

        # sort
        rpplist.sort(key=(lambda rpp: min(rpp.rp1.read.reference_start,
                                          rpp.rp2.read.reference_start)))

    return rpplist


def get_rpplist_sv(bam, fasta, chromdict, bnds, view=False,
                   fetch_padding_common=FETCH_PADDING_COMMON,
                   fetch_padding_sv=FETCH_PADDING_SV,
                   fetch_padding_view=FETCH_PADDING_VIEW,
                   new_fetch_padding=NEW_FETCH_PADDING,
                   long_insert_threshold=LONG_INSERT_THRESHOLD):
    LOGGER_RPPLIST.info('Beginning initial fetch')
    (relevant_qname_list_bnd1, 
     new_fetch_range_bnd1,
     relevant_qname_list_bnd2, 
     new_fetch_range_bnd2,
     relevant_qname_list_union) = initial_fetch_sv(bam, bnds, view, 
                                                   fetch_padding_common,
                                                   fetch_padding_sv, 
                                                   fetch_padding_view,
                                                   new_fetch_padding, 
                                                   long_insert_threshold)

    if relevant_qname_list_union is None:
        rpplist = None
    else:
        LOGGER_RPPLIST.info('Beginning refined fetch')
        fetchresult_dict = refined_fetch_sv(bam, bnds, new_fetch_range_bnd1, 
                                            new_fetch_range_bnd2,
                                            relevant_qname_list_union)

        LOGGER_RPPLIST.info('Beginning assembly into readpluspair')
        rpplist_bnd1 = ReadPlusPairList()
        rpplist_bnd2 = ReadPlusPairList()
        for readlist in fetchresult_dict.values():
            rpp = get_rpp_from_refinedfetch(readlist, bam, fasta, chromdict)
            if rpp is not None:
                if relevant_qname_list_bnd1 is not None:
                    if rpp.query_name in relevant_qname_list_bnd1:
                        rpplist_bnd1.append(rpp)
                if relevant_qname_list_bnd2 is not None:
                    if rpp.query_name in relevant_qname_list_bnd2:
                        rpplist_bnd2.append(rpp)

        # sort
        rpplist_sortkey = lambda rpp: min(rpp.rp1.read.reference_start,
                                          rpp.rp2.read.reference_start)
        rpplist_bnd1.sort(key=rpplist_sortkey)
        rpplist_bnd2.sort(key=rpplist_sortkey)

    return rpplist_bnd1, rpplist_bnd2


#### initial fetch ####

# readfilter

def readfilter_qname_nonsv(read, start0, end0):
    """For use in "initial_fetch" function.
    The read is selected if True.
    For decision of whether the read query_name will be included in the
        "relevant_qname_list", for non-SV.
    """

    #softclip_range0 = readhandler.get_softclip_ends_range0(read)
    #return (softclip_range0.stop > start0 and
    #        softclip_range0.start < end0)

    return True


def readfilter_qname_sv(read, start0, end0, endis5):
    """For use in "initial_fetch" function.
    The read is selected if True.
    For decision of whether the read query_name will be included in the
        "relevant_qname_list", for SV.
    """

    cond_overlapping = (read.reference_start < end0 and 
                        read.reference_end > start0)
    if endis5:
        cond_distant = (read.is_reverse and read.reference_start >= end0)
    else:
        cond_distant = (read.is_forward and read.reference_end <= start0)

    result = (cond_overlapping or cond_distant)

    return result


def readfilter_matesearch(read, long_insert_threshold):
    """For use in "initial_fetch" function.
    The read is selected if True.
    For decision of whether mate start position will be included in the 
        "new_fetch_range".
    """

    return (read.reference_name == read.next_reference_name and
            abs(read.template_length) < long_insert_threshold)


# get_initial_fetcher

def get_initial_fetcher_nonsv(bam, chrom, start0, end0, fetch_padding_common):
    fetcher = readhandler.get_fetch(
        bam, chrom, 
        start=(start0 - fetch_padding_common),
        end=(end0 + fetch_padding_common),
        filter_fun=readhandler.check_bad_read)

    return fetcher


def get_initial_fetcher_nonsv_view(bam, chrom, start0, end0, 
                                   fetch_padding_view):
    fetcher = readhandler.get_fetch(
        bam, chrom, 
        start=(start0 - fetch_padding_view),
        end=(end0 + fetch_padding_view),
        filter_fun=readhandler.check_bad_read)

    return fetcher


def get_initial_fetcher_sv(bam, bnds, fetch_padding_sv, fetch_padding_common):
    def subfun(bam, chrom, pos_range0, endis5, fetch_padding_common, 
               fetch_padding_sv):
        if endis5:
            fetcher = readhandler.get_fetch(
                bam, chrom,
                start=(pos_range0.start 
                       - fetch_padding_common),
                end=(pos_range0.stop 
                     + fetch_padding_common
                     + fetch_padding_sv),
                filter_fun=readhandler.check_bad_read)
        else:
            fetcher = readhandler.get_fetch(
                bam, chrom, 
                start=(pos_range0.start
                       - fetch_padding_common
                       - fetch_padding_sv),
                end=(pos_range0.stop 
                     + fetch_padding_common),
                filter_fun=readhandler.check_bad_read)

        return fetcher

    fetcher_bnd1 = subfun(bam, bnds.chrom_bnd1, bnds.pos_range0_bnd1,
                          bnds.endis5_bnd1,
                          fetch_padding_common, fetch_padding_sv)
    fetcher_bnd2 = subfun(bam, bnds.chrom_bnd2, bnds.pos_range0_bnd2,
                          bnds.endis5_bnd2,
                          fetch_padding_common, fetch_padding_sv)

    return fetcher_bnd1, fetcher_bnd2


def get_initial_fetcher_sv_view(bam, bnds, fetch_padding_view):
    fetcher_bnd1 = readhandler.get_fetch(
        bam, chrom,
        start=(bnds.pos_range0_bnd1.start - fetch_padding_view),
        end=(bnds.pos_range0_bnd1.stop + fetch_padding_view),
        filter_fun=readhandler.check_bad_read)
    fetcher_bnd2 = readhandler.get_fetch(
        bam, chrom,
        start=(bnds.pos_range0_bnd2.start - fetch_padding_view),
        end=(bnds.pos_range0_bnd2.stop + fetch_padding_view),
        filter_fun=readhandler.check_bad_read)

    return fetcher_bnd1, fetcher_bnd2


# initial_fetch

def initial_fetch_nonsv(bam, chrom, start0, end0, view,
                        fetch_padding_common, fetch_padding_view,
                        new_fetch_padding, long_insert_threshold):
    """Read filtering done by readhandler.check_bad_read"""

    relevant_qname_list = list()
    start0_list = list()

    if view:
        fetcher = get_initial_fetcher_nonsv_view(bam, chrom, start0, end0, 
                                                 fetch_padding_view)
    else:
        fetcher = get_initial_fetcher_nonsv(bam, chrom, start0, end0, 
                                            fetch_padding_common)

    for read in fetcher:
        if view:
            read_is_relevant = True
        else:
            read_is_relevant = readfilter_qname_nonsv(read, start0, end0)

        if read_is_relevant:
            relevant_qname_list.append(read.query_name)
            start0_list.append(read.reference_start)
            if readfilter_matesearch(read, long_insert_threshold):
                start0_list.append(read.next_reference_start)

    if len(relevant_qname_list) == 0:
        relevant_qname_list = None
        new_fetch_range = None
    else:
        relevant_qname_list = list(set(relevant_qname_list))
        new_fetch_range = range(min(start0_list) - new_fetch_padding, 
                                max(start0_list) + 1 + new_fetch_padding)

    return relevant_qname_list, new_fetch_range


def initial_fetch_sv(bam, bnds, view, fetch_padding_common,
                     fetch_padding_sv, fetch_padding_view,
                     new_fetch_padding, long_insert_threshold):
    """Read filtering done by readhandler.check_bad_read"""

    def fetcher_routine(fetcher, view, pos_range0, endis5, 
                        relevant_qname_list, start0_list, 
                        long_insert_threshold):
        """Adds info to 'relevant_qname_list' and 'start0_list'"""

        for read in fetcher:
            if view:
                read_is_relevant = True
            else:
                read_is_relevant = readfilter_qname_sv(
                    read, pos_range0.start, pos_range0.stop, endis5)

            if read_is_relevant:
                relevant_qname_list.append(read.query_name)
                start0_list.append(read.reference_start)
                if readfilter_matesearch(read, long_insert_threshold):
                    start0_list.append(read.next_reference_start)

    def postprocess(relevant_qname_list, start0_list, new_fetch_padding):
        if len(relevant_qname_list) == 0:
            relevant_qname_list = None
            new_fetch_range = None
        else:
            relevant_qname_list = list(set(relevant_qname_list))
            new_fetch_range = range(
                min(start0_list) - new_fetch_padding, 
                max(start0_list) + 1 + new_fetch_padding)

        return relevant_qname_list, new_fetch_range

    # set bnds.pos_range0s
    if (bnds.pos_range0_bnd1 is None) or (bnds.pos_range0_bnd2 is None):
        bnds.set_pos_range0s()

    # init
    relevant_qname_list_bnd1 = list()
    relevant_qname_list_bnd2 = list()
    start0_list_bnd1 = list()
    start0_list_bnd2 = list()

    # get fetchers
    if view:
        fetcher_bnd1, fetcher_bnd2 = get_initial_fetcher_sv_view(
            bam, bnds, fetch_padding_view)
    else:
        fetcher_bnd1, fetcher_bnd2 = get_initial_fetcher_sv(
            bam, bnds, fetch_padding_sv, fetch_padding_common)

    # extract information from fetchers
    fetcher_routine(fetcher_bnd1, view, bnds.pos_range0_bnd1, 
                    bnds.endis5_bnd1, relevant_qname_list_bnd1, 
                    start0_list_bnd1, long_insert_threshold)
    fetcher_routine(fetcher_bnd2, view, bnds.pos_range0_bnd2,
                    bnds.endis5_bnd2, relevant_qname_list_bnd2, 
                    start0_list_bnd2, long_insert_threshold)

    # postprocess
    relevant_qname_list_bnd1, new_fetch_range_bnd1 = postprocess(
        relevant_qname_list_bnd1, start0_list_bnd1, new_fetch_padding)
    relevant_qname_list_bnd2, new_fetch_range_bnd2 = postprocess(
        relevant_qname_list_bnd2, start0_list_bnd2, new_fetch_padding)

    # qname list union
    relevant_qname_list_union = set(
        (list() if relevant_qname_list_bnd1 is None else 
         relevant_qname_list_bnd1)
        + (list() if relevant_qname_list_bnd2 is None else 
           relevant_qname_list_bnd2))
    if len(relevant_qname_list_union) == 0:
        relevant_qname_list_union = None

    return (relevant_qname_list_bnd1, new_fetch_range_bnd1,
            relevant_qname_list_bnd2, new_fetch_range_bnd2,
            relevant_qname_list_union)


#### refined fetch ####

def refined_fetch_nonsv(bam, chrom, new_fetch_range, relevant_qname_list):
    """Read filtering is done by check_bad_read.
    Store the fetched read only if its qname is included in the
        "relevant_qname_list".
    """

    fetchresult_dict = dict()

    for qname in relevant_qname_list:
        fetchresult_dict[qname] = list()

    for read in readhandler.get_fetch(
            bam, chrom, start=new_fetch_range.start, 
            end=new_fetch_range.stop, 
            filter_fun=readhandler.check_bad_read):
        if read.query_name in fetchresult_dict:
            fetchresult_dict[read.query_name].append(read)

    return fetchresult_dict


def refined_fetch_sv(bam, bnds, new_fetch_range_bnd1, new_fetch_range_bnd2,
                     relevant_qname_list_union):
    """Read filtering is done by check_bad_read.
    Store the fetched read only if its qname is included in the
        "relevant_qname_list".
    """

    fetchresult_dict = dict()

    for qname in relevant_qname_list_union:
        fetchresult_dict[qname] = list()

    if new_fetch_range_bnd1 is not None:
        for read in readhandler.get_fetch(
                bam, bnds.chrom_bnd1,
                start=new_fetch_range_bnd1.start,
                end=new_fetch_range_bnd1.stop,
                filter_fun=readhandler.check_bad_read):
            if read.query_name in fetchresult_dict:
                fetchresult_dict[read.query_name].append(read)

    if new_fetch_range_bnd2 is not None:
        for read in readhandler.get_fetch(
                bam, bnds.chrom_bnd2,
                start=new_fetch_range_bnd2.start,
                end=new_fetch_range_bnd2.stop,
                filter_fun=readhandler.check_bad_read):
            if read.query_name in fetchresult_dict:
                fetchresult_dict[read.query_name].append(read)

    return fetchresult_dict


def get_rpp_from_refinedfetch(readlist, bam, fasta, chromdict):
    # classify reads into primary and non-primary ones
    readlist_primary = list()
    readlist_nonprimary = list()
    for read in readlist:
        if readhandler.check_primary_alignment(read):
            readlist_primary.append(read)
        else:
            readlist_nonprimary.append(read)

    if len(readlist_primary) > 2:
        e_msg_list = list()
        e_msg_list.append('More than two primary reads found:')
        for read in readlist_primary:
            e_msg_list.append(read.to_string())
        raise Exception('\n'.join(e_msg_list))
    elif len(readlist_primary) == 0: 
        # found reads are all non-primary
        # supplementary reads overlapping pos0 can be missed
        rpp = None
    else:
        rplist_nonprimary = [ReadPlus(x, fasta)
                             for x in readlist_nonprimary]
        if len(readlist_primary) == 1: # mate read is somewhere far away
            mate = readhandler.get_primary_mate(readlist_primary[0], bam)
            if mate is None:
                raise Exception(
                    f'Mate not found for this read:\n'
                    f'{readlist_primary[0].to_string()}')
            rplist_primary = [
                ReadPlus(readlist_primary[0], fasta), 
                ReadPlus(mate, fasta)]

        else: # len(readlist_primary) == 2
            rplist_primary = [ReadPlus(x, fasta) 
                              for x in readlist_primary]

        rpp = ReadPlusPair(rplist_primary, rplist_nonprimary, chromdict)

    return rpp



#class ReadPlusPairList_old(list):
#
#    @common.get_deco_arg_choices({'mode': ('sv', 'nonsv', 'view')})
#    def __init__(self, bam, chrom, start0, end0, fasta, chromdict, 
#                 mode, endis5=None,
#                 fetch_padding_common=FETCH_PADDING_COMMON,
#                 fetch_padding_sv=FETCH_PADDING_SV,
#                 long_insert_threshold=LONG_INSERT_THRESHOLD,
#                 new_fetch_padding=NEW_FETCH_PADDING,
#                 verbose=False):
#        """Read filtering is done by "readhandler.check_bad_read".
#        Mate-unmapped reads are filtered out.
#        Args:
#            mode: must be one of "sv", "nonsv", "view"
#                sv: for SV breakends
#                nonsv: for non-SV
#                view: for creating a new bam with tags added
#        """
#
#        # sanitycheck
#        if mode == 'sv':
#            if endis5 is None:
#                raise Exception(f'"endis5" must be set if "mode" is "sv".')
#
#        LOGGER_RPPLIST.info('Beginning initial fetch')
#        (relevant_qname_list, new_fetch_range) = self._initial_fetch(
#            bam, chrom, start0, end0, mode, endis5, 
#            fetch_padding_common, fetch_padding_sv,
#            new_fetch_padding, long_insert_threshold)
#
#        if relevant_qname_list is None:
#            self.qname_list = None
#            self.fetch_range = None
#            self.rpplist = None
#            self.start = None
#            self.end = None
#        else:
#            self.qname_list = relevant_qname_list
#            self.fetch_range = new_fetch_range
#    
#            LOGGER_RPPLIST.info('Beginning refined fetch')
#            fetchresult_dict = self._refined_fetch(bam, chrom, 
#                                                   new_fetch_range, 
#                                                   relevant_qname_list)
#    
#            LOGGER_RPPLIST.info('Beginning assembly into readpluspair')
#            for readlist in fetchresult_dict.values():
#                rpp = self._get_readpluspair(readlist, bam, fasta, chromdict)
#                if rpp is not None:
#                    self.append(rpp)
#
#            # sort
#            self.sort(key=(lambda rpp: min(rpp.rp1.read.reference_start,
#                                           rpp.rp2.read.reference_start)))
#
#    def update_alleleinfo(self, vcfspec):
#        for rpp in self:
#            rpp.update_alleleinfo(vcfspec)
#
#    def _readfilter_qname_nonsv(self, read, start0, end0):
#        """For use in "_initial_fetch" function.
#        The read is selected if True.
#        For decision of whether the read query_name will be included in the
#            "relevant_qname_list", for non-SV.
#        """
#
#        #softclip_range0 = readhandler.get_softclip_ends_range0(read)
#        #return (softclip_range0.stop > start0 and
#        #        softclip_range0.start < end0)
#
#        return True
#
#    def _readfilter_qname_sv(self, read, start0, end0, endis5):
#        """For use in "_initial_fetch" function.
#        The read is selected if True.
#        For decision of whether the read query_name will be included in the
#            "relevant_qname_list", for SV.
#        """
#
#        cond_overlapping = (read.reference_start < end0 and 
#                            read.reference_end > start0)
#        if endis5:
#            cond_distant = (read.is_reverse and read.reference_start >= end0)
#        else:
#            cond_distant = (read.is_forward and read.reference_end <= start0)
#
#        result = (cond_overlapping or cond_distant)
#
#        return result
#
#    def _readfilter_matesearch(self, read, long_insert_threshold):
#        """For use in "_initial_fetch" function.
#        The read is selected if True.
#        For decision of whether mate start position will be included in the 
#            "new_fetch_range".
#        """
#
#        return (read.reference_name == read.next_reference_name and
#                abs(read.template_length) < long_insert_threshold)
#
#    def _get_initial_fetcher(self, bam, chrom, start0, end0, mode, endis5,
#                             fetch_padding_common, fetch_padding_sv):
#        if mode == 'sv':
#            if endis5:
#                fetcher = readhandler.get_fetch(
#                    bam, chrom, 
#                    start=(start0 - fetch_padding_common),
#                    end=(end0 
#                         + fetch_padding_common
#                         + fetch_padding_sv),
#                    filter_fun=readhandler.check_bad_read)
#            else:
#                fetcher = readhandler.get_fetch(
#                    bam, chrom, 
#                    start=(start0 
#                           - fetch_padding_common
#                           - fetch_padding_sv),
#                    end=(end0 + fetch_padding_common),
#                    filter_fun=readhandler.check_bad_read)
#        elif mode in ('nonsv', 'view'):
#            fetcher = readhandler.get_fetch(
#                bam, chrom, 
#                start=(start0 - fetch_padding_common),
#                end=(end0 + fetch_padding_common),
#                filter_fun=readhandler.check_bad_read)
#
#        return fetcher
#
#    def _initial_fetch(self, bam, chrom, start0, end0, mode, endis5,
#                       fetch_padding_common, fetch_padding_sv,
#                       new_fetch_padding, long_insert_threshold):
#        """Read filtering done by readhandler.check_bad_read"""
#
#        relevant_qname_list = list()
#        start0_list = list()
#
#        # fetch start/end is padded with "initial_fetch_padding" to include 
#            # reads which overlap range(start0, end0) with soft-clipped 
#            # bases on IGV.
#
#        fetcher = self._get_initial_fetcher(
#            bam, chrom, start0, end0, mode, endis5,
#            fetch_padding_common, fetch_padding_sv)
#
#        for read in fetcher:
#            if mode == 'sv':
#                read_is_selected = self._readfilter_qname_sv(
#                    read, start0, end0, endis5)
#            elif mode == 'nonsv':
#                read_is_selected = self._readfilter_qname_nonsv(
#                    read, start0, end0)
#            elif mode == 'view':
#                read_is_selected = True
#
#            if read_is_selected:
#                relevant_qname_list.append(read.query_name)
#                start0_list.append(read.reference_start)
#                if self._readfilter_matesearch(read, long_insert_threshold):
#                    start0_list.append(read.next_reference_start)
#
#        if len(relevant_qname_list) == 0:
#            relevant_qname_list = None
#            new_fetch_range = None
#        else:
#            relevant_qname_list = list(set(relevant_qname_list))
#            new_fetch_range = range(min(start0_list) - new_fetch_padding, 
#                                    max(start0_list) + 1 + new_fetch_padding)
#
#        return relevant_qname_list, new_fetch_range
#
#    def _refined_fetch(self, bam, chrom, new_fetch_range, 
#                       relevant_qname_list):
#        """Read filtering is done by check_bad_read.
#        Store the fetched read only if its qname is included in the
#            "relevant_qname_list".
#        """
#
#        fetchresult_dict = dict()
#
#        for qname in relevant_qname_list:
#            fetchresult_dict[qname] = list()
#
#        for read in readhandler.get_fetch(
#                bam, chrom, start=new_fetch_range.start, 
#                end=new_fetch_range.stop, 
#                filter_fun=readhandler.check_bad_read):
#            if read.query_name in fetchresult_dict:
#                fetchresult_dict[read.query_name].append(read)
#
#        return fetchresult_dict
#
#    def _get_readpluspair(self, readlist, bam, fasta, chromdict):
#        # classify reads into primary and non-primary ones
#        readlist_primary = list()
#        readlist_nonprimary = list()
#        for read in readlist:
#            if readhandler.check_primary_alignment(read):
#                readlist_primary.append(read)
#            else:
#                readlist_nonprimary.append(read)
#
#        if len(readlist_primary) > 2:
#            e_msg_list = list()
#            e_msg_list.append('More than two primary reads found:')
#            for read in readlist_primary:
#                e_msg_list.append(read.to_string())
#            raise Exception('\n'.join(e_msg_list))
#        elif len(readlist_primary) == 0: 
#            # found reads are all non-primary
#            # supplementary reads overlapping pos0 can be missed
#            rpp = None
#        else:
#            rplist_nonprimary = [ReadPlus(x, fasta)
#                                 for x in readlist_nonprimary]
#            if len(readlist_primary) == 1: # mate read is somewhere far away
#                mate = readhandler.get_primary_mate(readlist_primary[0], bam)
#                if mate is None:
#                    raise Exception(
#                        f'Mate not found for this read:\n'
#                        f'{readlist_primary[0].to_string()}')
#                rplist_primary = [
#                    ReadPlus(readlist_primary[0], fasta), 
#                    ReadPlus(mate, fasta)]
#
#            else: # len(readlist_primary) == 2
#                rplist_primary = [ReadPlus(x, fasta) 
#                                  for x in readlist_primary]
#
#            rpp = ReadPlusPair(rplist_primary, rplist_nonprimary, chromdict)
#
#        return rpp


