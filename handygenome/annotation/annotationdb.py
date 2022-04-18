import re
import pprint
import sys
import collections
from operator import attrgetter


import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))

infoformat = importlib.import_module('.'.join([top_package_name, 'variantplus', 'infoformat']))
varianthandler = importlib.import_module('.'.join([top_package_name, 'variantplus', 'varianthandler']))

annotation = importlib.import_module('.'.join([top_package_name, 'annotation']))
ensembl_parser = importlib.import_module('.'.join([top_package_name, 'annotation', 'ensembl_parser']))
ensembl_rest = importlib.import_module('.'.join([top_package_name, 'annotation', 'ensembl_rest']))
veplib = importlib.import_module('.'.join([top_package_name, 'annotation', 'veplib']))
customfile = importlib.import_module('.'.join([top_package_name, 'annotation', 'customfile']))
rnalib = importlib.import_module('.'.join([top_package_name, 'annotation', 'rnalib']))

structvars = importlib.import_module('.'.join([top_package_name, 'svlib', 'structvars']))


# CONSTANTS

SEP_TOP = '$$$'
SEP_KEYVAL = ':::'
SEP_SEQ = '&&&'
SEP_SUBKEYVAL = '@@@'

ATOMIC_TYPES = (type(None), str, int, float, bool)


# POPFREQ CONSTANTS

def get_population_names(refver):
    pop_names = list()
    vcf = annotation.VCFS_DBSNP[refver]
    for key in vcf.header.info:
        if key.startswith('AF_'):
            pop_names.append(key[3:])

    return pop_names

POP_NAMES = set(get_population_names('hg19') + get_population_names('hg38'))


# COSMIC CONSTANTS

def get_cosmic_metadata():
    result = dict()
    for refver in ('hg19', 'hg38'):
        result[refver] = dict()
        for dbtype in ('coding', 'noncoding'):
            result[refver][dbtype] = dict()
            vcf = annotation.VCFS_COSMIC[refver][dbtype]
            for rec in vcf.header.records:
                if rec.key == 'cosmic_version':
                    result[refver][dbtype]['version'] = rec.value
                    continue
                elif rec.key.startswith('sample_count'):
                    result[refver][dbtype][rec.key] = int(rec.value)
                    continue

    return result

COSMIC_METADATAS = get_cosmic_metadata()


# MAIN

ANNOTDB_TYPES = ('plain', 'bnd1', 'bnd2')
ANNOTDB_ATTRNAMES = ('transcript', 'regulatory', 'motif', 'repeat', 
                     'popfreq', 'cosmic')
ANNOTDB_ATTRNAMES_FEATURES = ('transcript', 'regulatory', 'repeat')

INFOMETA_ID_PREFIX = {'plain': 'ANNOTS', 'bnd1': 'ANNOTS_BND1',
                      'bnd2': 'ANNOTS_BND2'}

INFOMETA_NUMS = {'transcript': '.', 'regulatory': '.', 'motif': '.',
                 'repeat': '.', 'popfreq': 1, 'cosmic': 1}
INFOMETA_DESCRIPTIONS = {
    'transcript': 'Transcript feature annotations encoded as a string',
    'regulatory': 'Regulatory element feature annotations encoded as a string',
    'motif': ('Transcription factor binding motif feature annotations '
              'encoded as a string'),
    'repeat': 'Repeat feature annotations encoded as a string',
    'popfreq': 'Population frequencies encoded as a string',
    'cosmic': 'COSMIC information encoded as a string',
    }

def get_infometas(ANNOTDB_TYPES, ANNOTDB_ATTRNAMES, INFOMETA_ID_PREFIX, 
                  INFOMETA_NUMS):
    """
    Returns: {
        'plain': { 
            'transcript': {'ID': ID, 'Type': Type, ...}, 
            'popfreq': {'ID': ID, 'Type': Type, ...}, 
            ...
            },
        'bnd1': { 
            'transcript': {'ID': ID, 'Type': Type, ...}, 
            'popfreq': {'ID': ID, 'Type': Type, ...}, 
            ... 
            },
        ...
        }
    """

    INFOMETAS = dict()
    for annotdb_type in ANNOTDB_TYPES:
        INFOMETAS[annotdb_type] = dict()
        for annotdb_attrname in ANNOTDB_ATTRNAMES:
            INFOMETAS[annotdb_type][annotdb_attrname] = {
                'ID': (INFOMETA_ID_PREFIX[annotdb_type] 
                       + '_' + annotdb_attrname),
                'Number': INFOMETA_NUMS[annotdb_attrname], 
                'Type': 'String',
                'Description': (f'{INFOMETA_DESCRIPTIONS[annotdb_attrname]}, '
                                f'for {annotdb_type}.'),
                }

    return INFOMETAS

INFOMETAS = get_infometas(ANNOTDB_TYPES, ANNOTDB_ATTRNAMES, 
                          INFOMETA_ID_PREFIX, INFOMETA_NUMS)


def add_infometas(pysamhdr):
    for annotdb_type in ANNOTDB_TYPES:
        for annotdb_attrname in ANNOTDB_ATTRNAMES:
            infometa = INFOMETAS[annotdb_type][annotdb_attrname]
            pysamhdr.add_meta(key = 'INFO', items = infometa.items())


class AnnotDB:
    """
    Attributes:
        featureDB
        transcript
        transcript_canonical
        regulatory
        motif
        repeat
        popfreq
        cosmic
    """

    def __init__(self, annotdb_type, refver, fasta, vr=None):
        assert annotdb_type in ANNOTDB_TYPES, (
            f'annotdb_type must be one of: {ANNOTDB_TYPES}')

        self.annotdb_type = annotdb_type
        self.refver = refver
        #assert self.refver in ensembl_rest.REFVERS_ALLOWED

        if fasta is None:
            self.fasta = pysam.FastaFile(
                common.DEFAULT_FASTA_PATHS[self.refver])
        else:
            self.fasta = fasta

        if vr is None:
            self.init_empty()
        else:
            self.load_vr(vr)

    def __repr__(self):
        return common.cpformat(self.get_showdict())

    def get_showdict(self):
        showdict = dict()
        showdict['popfreq'] = dict(self.popfreq)
        showdict['cosmic'] = dict(self.cosmic)

        showdict['transcript'] = self.transcript.get_showdict()
        showdict['transcript_canonical'] = \
            self.transcript_canonical.get_showdict()
        showdict['regulatory'] = self.regulatory.get_showdict()
        showdict['motif'] = self.motif.get_showdict()
        showdict['repeat'] = self.repeat.get_showdict()


        return showdict

    def init_empty(self):
        self.transcript = AnnotItemDict()
        self.regulatory = AnnotItemDict()
        self.motif = AnnotItemDict()
        self.repeat = AnnotItemDict()

        self.set_transcript_canonical()

        self.popfreq = Popfreq()
        self.cosmic = Cosmic()

    def load_vr(self, vr):
        self.transcript = AnnotItemDict(
            vr, infokey=INFOMETAS[self.annotdb_type]['transcript']['ID'])
        self.regulatory = AnnotItemDict(
            vr, infokey=INFOMETAS[self.annotdb_type]['regulatory']['ID'])
        self.motif = AnnotItemDict(
            vr, infokey=INFOMETAS[self.annotdb_type]['motif']['ID'])
        self.repeat = AnnotItemDict(
            vr, infokey=INFOMETAS[self.annotdb_type]['repeat']['ID'])

        self.set_transcript_canonical()

        self.popfreq = Popfreq(
            vr, infokey=INFOMETAS[self.annotdb_type]['popfreq']['ID'])
        self.cosmic = Cosmic(
            vr, infokey=INFOMETAS[self.annotdb_type]['cosmic']['ID'])

    def write(self, vr, addkey=False):
        if addkey:
            add_infometas(vr.header)

        self.transcript.write(
            vr, infokey=INFOMETAS[self.annotdb_type]['transcript']['ID'], 
            addkey=False)
        self.regulatory.write(
            vr, infokey=INFOMETAS[self.annotdb_type]['regulatory']['ID'], 
            addkey=False)
        self.motif.write(
            vr, infokey=INFOMETAS[self.annotdb_type]['motif']['ID'], 
            addkey=False)
        self.repeat.write(
            vr, infokey=INFOMETAS[self.annotdb_type]['repeat']['ID'], 
            addkey=False)

        self.popfreq.write(
            vr, infokey=INFOMETAS[self.annotdb_type]['popfreq']['ID'], 
            addkey=False)
        self.cosmic.write(
            vr, infokey=INFOMETAS[self.annotdb_type]['cosmic']['ID'], 
            addkey=False)

    # POPFREQ ###################################################

    def update_popfreq(self, vcfspec, dbsnp_vcf, search_equivs=True, 
                       overwrite=False):
        assert self.annotdb_type == 'plain', \
            f'This method is available only for annotdb_type "plain".'
        dbsnp_vr = fetch_dbsnp_vr(vcfspec, dbsnp_vcf, self.fasta, 
                                  search_equivs)
        self.popfreq.load_other(parse_dbsnp_vr(dbsnp_vr), overwrite=overwrite)

    # COSMIC #####################################################

    def update_cosmic(self, vcfspec, cosmic_coding_vcf, cosmic_noncoding_vcf,  
                      search_equivs=True, overwrite=False):
        assert self.annotdb_type == 'plain', \
            f'This method is available only for annotdb_type "plain".'
        cosmic_vr, cosmic_vr_noncoding = fetch_cosmic_vr(
            vcfspec, cosmic_coding_vcf, cosmic_noncoding_vcf, self.fasta, 
            search_equivs=search_equivs)
        self.cosmic.load_other(parse_cosmic_vr(cosmic_vr, cosmic_vr_noncoding), 
            overwrite = overwrite)

    # FEATURES ###################################################

    # post-cmdline-vep modifiers

    def update_features_postvep_plain(
            self, vcfspec, tabixfile_geneset, tabixfile_regulatory, 
            tabixfile_repeats, distance=veplib.DEFAULT_DISTANCE):
        """To be run with non-SV variants"""

        assert self.annotdb_type == 'plain'

        mttype = common.get_mttype(vcfspec.ref, vcfspec.alt)
        assert mttype != 'sv', f'Input variant must not be an SV.'

        # update using customfiles
        fetch_interval = common.Interval(chrom=vcfspec.chrom, 
                                         start1=(vcfspec.pos - distance), 
                                         end1=(vcfspec.pos + distance))
        self.update_customfile_transcript(fetch_interval, tabixfile_geneset)
        self.update_customfile_regulatory(fetch_interval, tabixfile_regulatory)
        self.update_customfile_repeat(fetch_interval, tabixfile_repeats)

        # modifier
        self.modifier_for_plain(vcfspec, mttype)
        self.set_transcript_canonical()

    def update_features_postvep_interval(
            self, interval, tabixfile_geneset, tabixfile_regulatory, 
            tabixfile_repeats, distance=veplib.DEFAULT_DISTANCE):
        assert self.annotdb_type == 'plain'

        # update using customfiles
        fetch_interval = common.Interval(chrom=interval.chrom, 
                                         start1=(interval.start1 - distance),
                                         end1=(interval.end1 + distance))
        self.update_customfile_transcript(fetch_interval, tabixfile_geneset)
        self.update_customfile_regulatory(fetch_interval, tabixfile_regulatory)
        self.update_customfile_repeat(fetch_interval, tabixfile_repeats)

        # modifier
        self.modifier_for_interval(interval)
        self.set_transcript_canonical()

    def update_features_postvep_bnd(
            self, chrom, pos, endis5, tabixfile_geneset, tabixfile_regulatory, 
            tabixfile_repeats, distance=veplib.DEFAULT_DISTANCE):
        assert self.annotdb_type in ('bnd1', 'bnd2')

        # update using customfiles
        if endis5:
            fetch_interval = common.Interval(chrom=chrom, start1=pos, 
                                             end1=(pos + distance))
        else:
            fetch_interval = common.Interval(chrom=chrom, 
                                             start1=(pos - distance), 
                                             end1=pos)

        self.update_customfile_transcript(fetch_interval, tabixfile_geneset)
        self.update_customfile_regulatory(fetch_interval, tabixfile_regulatory)
        self.update_customfile_repeat(fetch_interval, tabixfile_repeats)

        # modifier
        self.modifier_for_bnd(pos, endis5)
        self.set_transcript_canonical()

    # ensembl rest wrappers

    def update_features_ensembl_wrapper_plain(
            self, vcfspec, distance=veplib.DEFAULT_DISTANCE,
            vep=False, overlap=False, transcript=False, regulatory=False):
        """To be run with non-SV variants"""

        assert self.annotdb_type == 'plain'
        mttype = common.get_mttype(vcfspec.ref, vcfspec.alt)
        assert mttype != 'sv', f'Input variant must not be an SV.'

        if vep:
            self.update_ensembl_rest_vep_plain(vcfspec, distance)
        if overlap:
            self.update_ensembl_rest_overlap_plain(vcfspec, distance)
        if transcript:
            self.update_ensembl_rest_lookup_transcript()
        if regulatory:
            self.update_ensembl_rest_regulatory()

        self.modifier_for_plain(vcfspec, mttype)
        self.set_transcript_canonical()

    def update_features_ensembl_wrapper_interval(
            self, interval, distance=veplib.DEFAULT_DISTANCE,
            vep=False, overlap=False, transcript=False, regulatory=False):
        assert self.annotdb_type == 'plain'

        if vep:
            self.update_ensembl_rest_vep_interval(interval, distance)
        if overlap:
            self.update_ensembl_rest_overlap_interval(interval, distance)
        if transcript:
            self.update_ensembl_rest_lookup_transcript()
        if regulatory:
            self.update_ensembl_rest_regulatory()

        self.modifier_for_interval(interval)
        self.set_transcript_canonical()

    def update_features_bnd(
            self, chrom, pos, endis5, distance=veplib.DEFAULT_DISTANCE,
            vep=False, overlap=False, transcript=False, regulatory=False):
        assert self.annotdb_type in ('bnd1', 'bnd2')

        if vep:
            self.update_ensembl_rest_vep_bnd(chrom, pos, distance)
        if overlap:
            self.update_ensembl_rest_overlap_bnd(chrom, pos, endis5, distance)
        if transcript:
            self.update_ensembl_rest_lookup_transcript()
        if regulatory:
            self.update_ensembl_rest_regulatory()

        self.modifier_for_bnd(pos, endis5)
        self.set_transcript_canonical()


    # FEATURE UPDATE MODIFIERS ##########################################

    # for setting missing distances

    def set_missing_distances(self, distance_setter):
        """Sets distance for annotitems not derived from VEP (e.g. repeat 
        elements). Simulates the distance given by VEP.
        """

        for attrname in ANNOTDB_ATTRNAMES_FEATURES:
            annotitemdict = getattr(self, attrname)
            for key, annotitem in annotitemdict.items():
                if 'distance' not in annotitem: 
                    # annotitems derived from "ensembl rest overlap" 
                    #   do not have "distance" value
                    assert 'start0' in annotitem and 'end0' in annotitem, \
                        (f'"start0" and "end0" values are not set for this '
                         f'AnnotItem object:\n{annotitem}')
                    feature_range0 = range(annotitem['start0'], 
                                           annotitem['end0'])
                    annotitem['distance'] = distance_setter(feature_range0)

    def get_distance_setter_ins(self, vcfspec):
        pos0 = vcfspec.pos - 1
        def distance_setter(feature_range0):
            if pos0 <= feature_range0.start - 1:
                distance = feature_range0.start - pos0 - 1
            elif pos0 >= feature_range0.stop - 1:
                distance = feature_range0.stop - 1 - pos0
            else:
                distance = None

            return distance

        return distance_setter
                
    def get_distance_setter_snv(self, vcfspec):
        pos0 = vcfspec.pos - 1
        def distance_setter(feature_range0):
            if pos0 < feature_range0.start:
                distance = feature_range0.start - pos0
            elif pos0 >= feature_range0.stop:
                distance = feature_range0.stop - 1 - pos0
            else:
                distance = None

            return distance
                
        return distance_setter

    def get_distance_setter_interval(self, var_range0):
        def distance_setter(feature_range0):
            if var_range0.stop <= feature_range0.start:
                distance = feature_range0.start - var_range0.stop + 1
            elif var_range0.start >= feature_range0.stop:
                distance = feature_range0.stop - 1 - var_range0.start
            else:
                distance = None

            return distance
                
        return distance_setter

    def get_distance_setter_del(self, vcfspec):
        var_range0 = range(vcfspec.pos, vcfspec.pos + len(vcfspec.ref) - 1)
        return self.get_distance_setter_interval(var_range0)

    def get_distance_setter_bnd(self, pos):
        vcfspec = common.Vcfspec(None, pos, None, None)
        return self.get_distance_setter_snv(vcfspec)

    # modifier for plain
    def modifier_for_plain(self, vcfspec, mttype):
        if mttype == 'snv':
            distance_setter = self.get_distance_setter_snv(vcfspec)
        elif mttype == 'ins':
            distance_setter = self.get_distance_setter_ins(vcfspec)
        else: # del or cindel
            distance_setter = self.get_distance_setter_del(vcfspec)
        self.set_missing_distances(distance_setter)

    # modifier for interval
    def modifier_for_interval(self, interval):
        distance_setter = self.get_distance_setter_interval(interval.range0)
        self.set_missing_distances(distance_setter)
        self.set_is_enclosed(interval.range0)

    def set_is_enclosed(self, var_range0):
        """
        Sets "is_enclosed" attribute. Must be run after "start0" and "end0" 
        values have been set for all annotitems.
        """

        for attrname in ANNOTDB_ATTRNAMES_FEATURES:
            annotitemdict = getattr(self, attrname)
            for key, annotitem in annotitemdict.items():
                feature_range0 = range(annotitem['start0'], annotitem['end0'])
                annotitem['is_enclosed'] = (
                    feature_range0.start >= var_range0.start and
                    feature_range0.stop <= var_range0.stop
                )

    # modifier for breakend
    def modifier_for_bnd(self, pos, endis5):
        distance_setter = self.get_distance_setter_bnd(pos)
        self.set_missing_distances(distance_setter)

        self.settle_border_issues(pos, endis5)

    def settle_border_issues(self, pos, endis5):
        """
        Must be run after "set_missing_distances" and "lookup_transcript" or 
        "lookup_regulatory" has been run.

        1) Removes features entirely on the other side
        2) For features whose border coincide with the breakend border,
        distance is changed from None to 0.
        3) Sets "is_broken" attribute
        """

        on_the_other_side, on_the_border = self.get_feature_classifiers(endis5)

        for attrname in ANNOTDB_ATTRNAMES_FEATURES:
            annotitemdict = getattr(self, attrname)
            # step 1) and 2)
            keys_to_delete = list()
            for key, annotitem in annotitemdict.items():
                if on_the_other_side(annotitem):
                    keys_to_delete.append(key)
                else:
                    if annotitem['distance'] is None:
                        if on_the_border(annotitem, pos):
                            annotitem['distance'] = 0
            for key in keys_to_delete:
                del annotitemdict[key]
            # step 3)
            for key, annotitem in annotitemdict.items():
                annotitem['is_broken'] = (annotitem['distance'] is None)

    def get_feature_classifiers(self, endis5):
        """
        (endis5 == True)
                           THE OTHER SIDE    THE BREAKEND SIDE
                           - - - - - - - -  ----------------- 
FEATURE ON THE OTHER SIDE  = = = = = =     | = = = = = = = =  : FEATURE ON THE BORDER
                           - - - - - - - -  -----------------
                           A C G T A C G T   A C G T A C G T
                                             *(pos)
        ("= = =" denotes a feature)
        """

        if endis5:
            def on_the_other_side(annotitem):
                return (
                    (annotitem['distance'] is not None) and 
                    (annotitem['distance'] < 0))

            def on_the_border(annotitem, pos):
                return (annotitem['start1'] == pos)
        else:
            def on_the_other_side(annotitem):
                return (
                    (annotitem['distance'] is not None) and 
                    (annotitem['distance'] > 0))

            def on_the_border(annotitem, pos):
                return (annotitem['end1'] == pos)

        return on_the_other_side, on_the_border

    # CUSTOM FILE FETCH FUNCTIONS #####################################

    def update_customfile_repeat(self, fetch_interval, tabixfile_repeats):
        repeat_dict = customfile.fetch_repeat(
            fetch_interval.chrom, fetch_interval.start0, fetch_interval.end0, 
            tabixfile_repeats)
        self.repeat.load_other(repeat_dict, overwrite=False, create_new=True)

    def update_customfile_transcript(self, fetch_interval, tabixfile_geneset):
        transcript_dict = customfile.fetch_transcript(
            fetch_interval.chrom, fetch_interval.start0, fetch_interval.end0, 
            tabixfile_geneset)
        self.transcript.load_other(transcript_dict, overwrite=False, 
                                   create_new=True)

    def update_customfile_regulatory(self, fetch_interval, 
                                     tabixfile_regulatory):
        regulatory_dict = customfile.fetch_regulatory(
            fetch_interval.chrom, fetch_interval.start0, fetch_interval.end0, 
            tabixfile_regulatory)
        self.regulatory.load_other(regulatory_dict, overwrite=False, 
                                   create_new=True)

    # ENSEMBL UNIT FUNCTIONS ##########################################

    # vep
    def update_ensembl_rest_vep_plain(self, vcfspec,
                                      distance=veplib.DEFAULT_DISTANCE):
        self._update_ensembl_rest_vep_helper(
            vcfspec = vcfspec, hgvsg = None, distance = distance)

    def update_ensembl_rest_vep_interval(self, interval,
                                         distance=veplib.DEFAULT_DISTANCE):
        sv_del = structvars.Deletion(interval.chrom, interval.start1, 
                                     interval.end1, fasta=self.fasta)
        hgvsg = sv_del.get_hgvsg()
        self._update_ensembl_rest_vep_helper(
            vcfspec=None, hgvsg=hgvsg, distance=distance)

    def update_ensembl_rest_vep_bnd(self, chrom, pos,
                                    distance=veplib.DEFAULT_DISTANCE):
        ref = self.fasta.fetch(chrom, pos - 1, pos)
        vcfspec = common.Vcfspec(chrom, pos, ref, 'N')
        self._update_ensembl_rest_vep_helper(
            vcfspec=vcfspec, hgvsg=None, distance=distance)

    # overlap
    def update_ensembl_rest_overlap_plain(self, vcfspec,
                                          distance=veplib.DEFAULT_DISTANCE):
        chrom = vcfspec.chrom
        start1 = vcfspec.pos - distance
        end1 = vcfspec.pos + distance
        interval = common.Interval(chrom, start1, end1)
        self._update_ensembl_rest_overlap_helper(interval)

    def update_ensembl_rest_overlap_interval(self, interval, 
                                             distance=veplib.DEFAULT_DISTANCE):
        new_interval = common.Interval(chrom=interval.chrom,
                                       start1=(interval.start1 - distance),
                                       end1=(interval.end1 + distance))
        self._update_ensembl_rest_overlap_helper(new_interval)

    def update_ensembl_rest_overlap_bnd(self, chrom, pos, endis5, 
                                        distance=veplib.DEFAULT_DISTANCE):
        if endis5:
            start1 = pos
            end1 = pos + distance
        else:
            start1 = pos - distance
            end1 = pos
        interval = common.Interval(chrom, start1, end1)
        self._update_ensembl_rest_overlap_helper(interval)

    # lookup_transcript
    def update_ensembl_rest_lookup_transcript(self):
        parsed = ensembl_parser.parse_rest_lookup_transcript_post(
            ensembl_rest.lookup_id_post(tuple(self.transcript.keys()),
                                        refver=self.refver, expand=True),
            refver=self.refver, set_gene_name=True)
        self._load_ensembl_parsed(parsed, overwrite=False, create_new=False)

    # regulatory
    def update_ensembl_rest_regulatory(self):
        for ID in self.regulatory.keys():
            self._update_ensembl_rest_regulatory_helper(ID)

    # cmdline VEP
    def update_cmdline_vep(self, vr, overwrite=True, create_new=True):
        parsed = ensembl_parser.parse_cmdline_vep(vr)
        self._load_ensembl_parsed(parsed, overwrite=overwrite, 
                                  create_new=create_new)

    # helpers
    def _load_ensembl_parsed(self, parsed, overwrite=False, create_new=False):
        for annotdb_attrname, other_annotitemdict in parsed.items():
            annotitemdict = getattr(self, annotdb_attrname)
            annotitemdict.load_other(other_annotitemdict, overwrite=overwrite, 
                                     create_new=create_new)

    def _update_ensembl_rest_vep_helper(self, vcfspec, hgvsg, distance):
        parsed = ensembl_parser.parse_rest_vep(
            ensembl_rest.vep(
                refver=self.refver, vcfspec=vcfspec, hgvsg=hgvsg,
                distance=distance, with_CADD=True, with_Phenotypes=False,
                with_canonical=True, with_mane=True, with_miRNA=False,
                with_numbers=True, with_protein=True, with_ccds=True,
                with_hgvs=True))
        self._load_ensembl_parsed(parsed, overwrite=False, create_new=True)

    def _update_ensembl_rest_overlap_helper(self, interval):
        parsed = ensembl_parser.parse_rest_overlap( 
            ensembl_rest.overlap(
                chrom=interval.chrom, start1=interval.start1, 
                end1=interval.end1, refver=self.refver, transcript=False, 
                regulatory=False, motif=False, repeat=True),
            refver=self.refver, include_motif_without_evidence=False)
        self._load_ensembl_parsed(parsed, overwrite=False, create_new=True)

    def _update_ensembl_rest_lookup_transcript_helper(self, ID):
        parsed = ensembl_parser.parse_rest_lookup_transcript(
            ensembl_rest.lookup_id(ID, refver=self.refver, expand=True),
            refver=self.refver, set_gene_name=True)
        self._load_ensembl_parsed(parsed, overwrite=False, create_new=False)

    def _update_ensembl_rest_regulatory_helper(self, ID):
        parsed = ensembl_parser.parse_rest_regulatory(
            ensembl_rest.regulatory(ID, refver=self.refver))
        self._load_ensembl_parsed(parsed, overwrite=False, create_new=False)

    # FEATURE ADDITIONAL ATTRIBUTES GENERATION #####################

    def set_transcript_canonical(self):
        self.transcript_canonical = AnnotItemDict()
        for ID, annotitem in self.transcript.items():
            if annotitem['is_canonical']:
                self.transcript_canonical[ID] = annotitem

    def set_exon_intron_ranges(self):
        def make_exon_ranges(annotitem):
            exon_ranges = dict()
            for key, val in annotitem['exons'].items():
                keysp = key.split('_')
                valsp = [int(x) for x in val.split('_')]
                exon_number = int(keysp[0])
                exon_ranges[exon_number] = range(valsp[0], valsp[1])
            return exon_ranges

        def make_intron_ranges(exon_ranges, annotitem):
            if len(exon_ranges) == 1:
                intron_ranges = None
            else:
                intron_ranges = dict()
                ascending = sorted(exon_ranges.items(), key=lambda x: x[1][0])
                for idx in range(len(ascending) - 1):
                    intron_number = (ascending[idx][0] 
                                     if annotitem['is_forward'] else 
                                     ascending[idx+1][0])
                    range_item = range(ascending[idx][1][-1] + 1, 
                                       ascending[idx+1][1][0])
                    intron_ranges[intron_number] = range_item
                intron_ranges = dict(
                    sorted(intron_ranges.items(), key=lambda x: x[0]))

            return intron_ranges

        for ID, annotitem in self.transcript.items():
            if 'exons' in annotitem.keys():
                annotitem.exon_ranges = make_exon_ranges(annotitem)
                annotitem.intron_ranges = make_intron_ranges(
                    annotitem.exon_ranges, annotitem)


class AnnotItem(dict):
    infokey = None

    @staticmethod
    def decode(infostring, dic):
        """dic is updated with contents of infostring"""

        parsed = dict()
        for x in infostring.split(SEP_TOP):
            x_sp = x.split(SEP_KEYVAL)
            parsed[x_sp[0]] = x_sp[1]
        for key, val in parsed.items():
            if len(val.split(SEP_SUBKEYVAL)) == 1: # not a dict
                new_val = [common.str_to_nonstr(x) for x in val.split(SEP_SEQ)]
                if len(new_val) == 1:
                    new_val = new_val[0]
                else:
                    new_val = new_val[1:-1]
            else:
                tmp = list()
                for x in val.split(SEP_SEQ):
                    x_split = x.split(SEP_SUBKEYVAL)
                    tmp.append((x_split[0], common.str_to_nonstr(x_split[1])))
                new_val = dict(tmp)

            parsed[key] = new_val
        dic.update(parsed)

    def __init__(self, vr = None, infokey = None):
        super().__init__()
        if (vr is not None):
            self.load_vr(vr, infokey)

    def __repr__(self):
        return common.cpformat(dict(self))

    def __setitem__(self, key, val):
        self._setitem_sanitycheck(key, val)
        super().__setitem__(key, val)

    def _setitem_sanitycheck(self, key, val):
        try:
            assert isinstance(key, str), f'The key must be a str object.'
            if isinstance(val, dict):
                for subkey, subval in val.items():
                    assert isinstance(subkey, str), (
                        f'When the value is a dict, each subkey must be a '
                        f'str object.')
                    assert isinstance(subval, ATOMIC_TYPES), (
                        f'When the value is a dict, the type of each '
                        f'subvalue must be one of {str(ATOMIC_TYPES)}.')
            elif isinstance(val, list):
                for subval in val:
                    assert isinstance(subval, ATOMIC_TYPES), (
                        f'When the value is a list, the type of each item '
                        f'must be one of {str(ATOMIC_TYPES)}.')
            else:
                assert isinstance(val, ATOMIC_TYPES), (
                    f'When the value is neither list nor dict, its type '
                    f'must be one of {str(ATOMIC_TYPES)}.')
        except AssertionError as e:
            new_e = Exception(f'Invalid key/value for assignment to '
                               f'AnnotItem: key - {key} ; value - {value}')
            raise new_e from e

    def load_infostring(self, infostring):
        self.decode(infostring, self)

    def load_vr(self, vr, infokey = None):
        if infokey is None:
            infokey = self.__class__.infokey

        assert vr.header.info[infokey].number == 1
        assert vr.header.info[infokey].type == 'String'

        if not infoformat.check_NA_info(vr, infokey):
            infostring = infoformat.get_info(vr, infokey)
            self.load_infostring(infostring)

    def load_other(self, other, overwrite = False):
        """
        other: an AnnotItem object
        """

        if overwrite:
            for key, val in other.items():
                self[key] = val
        else:
            for key, val in other.items():
                try:
                    old_val = self[key]
                except KeyError:
                    self[key] = val
                else:
                    if old_val is None:
                        self[key] = val
                    else:
                        if val != old_val:
                            e_msg = (
                                f'Old value and new value are different; '
                                f'key: {key} ; old_val: {old_val} ; '
                                f'new_val: {val}')
                            raise Exception(e_msg)

    def encode(self):
        if len(self) == 0:
            infostring = None
        else:
            result = list()
            for key, val in self.items():
                if isinstance(val, (tuple, list)):
                    modified_val = (SEP_SEQ 
                                    + SEP_SEQ.join(map(str, val)) 
                                    + SEP_SEQ)
                elif isinstance(val, dict):
                    modified_val = SEP_SEQ.join(
                        SEP_SUBKEYVAL.join(map(str, x)) 
                        for x in val.items())
                else:
                    modified_val = val

                result.append(f'{key}{SEP_KEYVAL}{modified_val}')

            infostring = SEP_TOP.join(result)

        return infostring

    def write(self, vr, infokey=None, addkey=False):
        if infokey is None:
            infokey = self.__class__.infokey

        if addkey:
            add_infometas(vr.header)

        infoformat.set_info(vr, infokey, self.encode())


class AnnotItemDict(dict):
    def __init__(self, vr = None, infokey = None):
        super().__init__()
        if (vr is not None) and (infokey is not None):
            self.load_vr(vr, infokey)

    def __setitem__(self, key, val):
        assert isinstance(val, AnnotItem)
        super().__setitem__(key, val)

    def __repr__(self):
        return common.cpformat(self.get_showdict())

    def get_showdict(self):
        showdict = dict()
        for key, val in sorted(self.items(), key = lambda x: x[1]['start0']):
            showdict[key] = dict(val)

        return showdict

    def get_list(self):
        return sorted(self.values(), key = attrgetter('start0'))

    def load_vr(self, vr, infokey):
        assert vr.header.info[infokey].number == '.'
        assert vr.header.info[infokey].type == 'String'

        if not infoformat.check_NA_info(vr, infokey):
            infoval = infoformat.get_info(vr, infokey, collapse_tuple = False)
            for infostring in infoval:
                annotitem = AnnotItem()
                annotitem.load_infostring(infostring)
                self[annotitem['id']] = annotitem

    def load_other(self, other, overwrite = True, create_new = False):
        for ID, annotitem in other.items():
            if ID in self.keys():
                self[ID].load_other(annotitem, overwrite = overwrite)
            else:
                if create_new:
                    self[ID] = annotitem

    def write(self, vr, infokey, addkey = False):
        if addkey:
            add_infometas(vr.header)

        infoval = list()
        for annotitem in self.values():
            infoval.append(annotitem.encode())

        vr.info[infokey] = infoval


######################################


class Popfreq(AnnotItem):
    #infokey = INFOMETAS['nonSV']['popfreq']['ID']
    pass


class Cosmic(AnnotItem):
    #infokey = INFOMETAS['nonSV']['cosmic']['ID']
    pass


######################################


def fetch_dbsnp_vr(vcfspec, dbsnp_vcf, fasta, search_equivs=True):
    """
    Result may be None
    """

    dbsnp_vr = customfile.fetch_relevant_vr(vcfspec, dbsnp_vcf, fasta=fasta, 
                                            search_equivs=search_equivs)
    
    return dbsnp_vr


def fetch_cosmic_vr(vcfspec, cosmic_coding_vcf, cosmic_noncoding_vcf, fasta, 
                    search_equivs=True):
    """
    Result may be (None, None)
    """

    cosmic_vr = customfile.fetch_relevant_vr(
        vcfspec, cosmic_coding_vcf, fasta=fasta, 
        search_equivs=search_equivs)
    cosmic_vr_noncoding = customfile.fetch_relevant_vr(
        vcfspec, cosmic_noncoding_vcf, fasta=fasta, 
        search_equivs=search_equivs)

    return cosmic_vr, cosmic_vr_noncoding


######################################


def parse_dbsnp_vr(dbsnp_vr):
    popfreq = Popfreq()

    for pop in POP_NAMES:
        popfreq[pop] = None

    if dbsnp_vr is not None:
        popfreq['id'] = dbsnp_vr.info['dbSNP_ID']

        for key, val in dbsnp_vr.info.items():
            if key.startswith('AF_'):
                pop = key[3:]
                popfreq[pop] = infoformat.get_info(dbsnp_vr, key)

    return popfreq


def parse_cosmic_vr(cosmic_vr, cosmic_vr_noncoding):
    def handler(val):
        if val is None:
            result = None
        else:
            result = dict()
            #if not isinstance(val, (tuple, list)):
            #    val = [ val ]
            for x in val:
                sp = x.split(':::')
                result[sp[0]] = common.str_to_nonstr(sp[1])

        return result

    def set_vals_common(cosmic, vr):
        cosmic['id'] = vr.info['cosmic_ID']
        cosmic['occurrence'] = handler(vr.info['cosmic_occurrence'])
        cosmic['portion'] = handler(vr.info['cosmic_portion'])
        cosmic['total_occurrence'] = vr.info['cosmic_total_occurrence']
        cosmic['total_portion'] = vr.info['cosmic_total_portion']

        # below lines use 'infoformat.get_info' because each value may be missing
        cosmic['occurrence_somatic'] = handler(
            infoformat.get_info(vr, 'cosmic_occurrence_somatic', 
                                collapse_tuple=False))
        cosmic['portion_somatic'] = handler(
            infoformat.get_info(vr, 'cosmic_portion_somatic', 
                                collapse_tuple=False))
        cosmic['total_occurrence_somatic'] = infoformat.get_info(
            vr, 'cosmic_total_occurrence_somatic')
        cosmic['total_portion_somatic'] = infoformat.get_info(
            vr, 'cosmic_total_portion_somatic')

        cosmic['coding_score'] = infoformat.get_info(vr, 'cosmic_coding_score')

    def set_vals_coding(cosmic, vr):
        set_vals_common(cosmic, vr)
        cosmic['noncoding_score'] = None

    def set_vals_noncoding(cosmic, vr):
        set_vals_common(cosmic, vr)
        cosmic['noncoding_score'] = vr.info['cosmic_noncoding_score']

    def set_max_occur(cosmic):
        cosmic['max_occur_site'], cosmic['max_occur_count'] = max(
            cosmic['occurrence'].items(), key=(lambda x: x[1]))
        cosmic['max_occur_portion'] = (
            cosmic['portion'][cosmic['max_occur_site']])

        if cosmic['occurrence_somatic'] is None:
            cosmic['max_occur_site_somatic'] = None
            cosmic['max_occur_count_somatic'] = None
            cosmic['max_occur_portion_somatic'] = None
        else:
            cosmic['max_occur_site_somatic'], \
            cosmic['max_occur_count_somatic'] = max(
                cosmic['occurrence_somatic'].items(), key=(lambda x: x[1]))
            cosmic['max_occur_portion_somatic'] = (
                cosmic['portion_somatic'][cosmic['max_occur_site_somatic']])

    def set_empty(cosmic):
        cosmic['id'] = None
        cosmic['occurrence'] = None
        cosmic['portion'] = None
        cosmic['total_occurrence'] = None
        cosmic['total_portion'] = None

        # below lines use 'infoformat.get_info' because each value may be missing
        cosmic['occurrence_somatic'] = None
        cosmic['portion_somatic'] = None
        cosmic['total_occurrence_somatic'] = None
        cosmic['total_portion_somatic'] = None

        cosmic['coding_score'] = None
        cosmic['noncoding_score'] = None

        cosmic['max_occur_site'] = None
        cosmic['max_occur_count'] = None
        cosmic['max_occur_portion'] = None

        cosmic['max_occur_site_somatic'] = None
        cosmic['max_occur_count_somatic'] = None
        cosmic['max_occur_portion_somatic'] = None

    def main(cosmic_vr, cosmic_vr_noncoding):
        cosmic = Cosmic()

        if cosmic_vr is None and cosmic_vr_noncoding is None:
            set_empty(cosmic)
        else:
            # fetch coding database
            if cosmic_vr is not None:
                set_vals_common(cosmic, cosmic_vr)
                set_vals_coding(cosmic, cosmic_vr)
                set_max_occur(cosmic)

            # fetch noncoding database, only if coding data is unavailable
            # raises an exception if COSV IDs are different between coding 
            # and noncoding databases
            if cosmic_vr_noncoding is not None:
                if cosmic_vr is None:
                    set_vals_common(cosmic, cosmic_vr_noncoding)
                    set_vals_noncoding(cosmic, cosmic_vr_noncoding)
                    set_max_occur(cosmic)
                else:
                    if cosmic['id'] == cosmic_vr_noncoding.info['cosmic_ID']:
                        pass
                    else:
                        vcfspec = common.Vcfsepc(
                            cosmic_vr.chrom, cosmic_vr.pos, 
                            cosmic_vr.ref, cosmic_vr.alts[0])
                        raise Exception(
                            f'Different COSV IDs between coding and '
                            f'noncoding databases for the same '
                            f'mutation: {vcfspec}')

        return cosmic

    return main(cosmic_vr, cosmic_vr_noncoding)
