import re
import pprint

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
breakends_module = importlib.import_module('.'.join([top_package_name, 'variantplus', 'breakends']))
varianthandler = importlib.import_module('.'.join([top_package_name, 'variantplus', 'varianthandler']))
infoformat = importlib.import_module('.'.join([top_package_name, 'variantplus', 'infoformat']))
#veplib = importlib.import_module('.'.join([top_package_name, 'annotation', 'veplib']))
annotationdb = importlib.import_module('.'.join([top_package_name, 'annotation', 'annotationdb']))


class VariantPlus:
    """
    Attributes:
        vr
        fasta
        chromdict
        bnds
        is_bnd1
        VEPkeys
        is_VEPannot
        is_splitVEPannot
        VEPannot_Nfeatures
        VEPannot_list
        canonical_VEPannot_list
        VEPannot_dict
        canonical_VEPannot_dict
        is_sv
        is_bnd1
        #bnds
        bnds
        bnds_equivs
        bnds_maxspan_pos1
        bnds_maxspan_pos2
        homlen
        homseq
    """

    def __init__(self, vr, fasta=None, chromdict=None, refver=None,
                 add_annotdb_infometas=False, set_annotdb=True):
        """
        Args:
            vr: pysam.VariantRecord instance
            fasta: pysam.FastaFile instance
            chromdict: handygenome.common.Chromdict instance
        """

        assert len(vr.alts) == 1, (
            f'Multiallelic variant record is not allowed:\n{vr}')

        self.vr = vr
        self.refver = (common.infer_refver_vr(self.vr)
                       if refver is None else
                       refver)
        self.fasta = (
            pysam.FastaFile(common.DEFAULT_FASTA_PATH_DICT[self.refver])
            if fasta is None else
            fasta)
        self.chromdict = (common.ChromDict(fasta=self.fasta)
                          if chromdict is None else
                          chromdict)
        self.vcfspec = varianthandler.get_vcfspec(self.vr)

        # SV attrs
        self.is_sv = varianthandler.check_SV(self.vr)
        self.set_bnds_attributes() # is_bnd1, bnds

        # annotdb
        if add_annotdb_infometas:
            annotationdb.add_infometas(self.vr.header)
        if set_annotdb:
            self.set_annotdb()

    def __repr__(self):
        vr_string = '\t'.join(str(self.vr).split('\t')[:5])
        return f'<VariantPlus object ({vr_string})>'

    def set_annotdb(self):
        if self.is_sv:
            if self.bnds.chrom_bnd1 != self.bnds.chrom_bnd2:
                self.annotdb_bnd1 = annotationdb.AnnotDB(
                    'bnd1', self.refver, self.fasta, self.chromdict,
                    vr=self.vr)
                self.annotdb_bnd2 = annotationdb.AnnotDB(
                    'bnd2', self.refver, self.fasta, self.chromdict,
                    vr=self.vr)
            else:
                self.annotdb = annotationdb.AnnotDB(
                    'plain', self.refver, self.fasta, self.chromdict,
                    vr=self.vr)
        else:
            self.annotdb = annotationdb.AnnotDB(
                'plain', self.refver, self.fasta, self.chromdict, vr=self.vr)

    ##################################################################

    def get_gene_names(self, canonical_only=True):
        if canonical_only:
            return [feature['gene_name'] 
                    for feature in self.annotdb.transcript_canon_ovlp.values()]
        else:
            return [feature['gene_name'] 
                    for feature in self.annotdb.transcript_ovlp.values()]

    def check_intergenic(self):
        return len(self.annotdb.transcript) == 0

    def get_info(self, key):
        return infoformat.get_value_info(self.vr, key)

    def get_format(self, sampleid, key):
        return infoformat.get_value_format(self.vr, sampleid, key)

    def check_NA_info(self, key):
        return infoformat.check_NA_info(self.vr, key)

    def check_NA_format(self, sampleid, key):
        return infoformat.check_NA_format(self.vr, sampleid, key)

    def show_info(self):
        infoformat.show_info(self.vr)

    def show_format(self):
        infoformat.show_format(self.vr)

    ##################################################################

#    def update_popfreq(self):
#        self.annotdb.update_popfreq(self.vr, search_equivs = True)
#        self.annotdb.write_popfreq(self.vr, addkey = False)
#    
#    def update_cosmic(self):
#        self.annotdb.update_cosmic(self.vr, search_equivs = True)
#        self.annotdb.write_cosmic(self.vr, addkey = False)

#    def update_features(self, distance=5000, vep=True, overlap=True, 
#                        transcript=True, regulatory=True):
#        self.annotdb.update_features(
#            self.vr, hg19=None, distance=distance,
#            vep=vep, overlap=overlap, 
#            transcript=transcript, regulatory=regulatory)
#        self.annotdb.write_features(self.vr, addkey = False)

    ##################################################################

    def get_vr_bnd2(self):
        assert self.bnds is not None, f'"bnds" attribute must be set.'

        vr_bnd2 = self.vr.header.new_record()
        vr_bnd2.id = self.bnds.get_id_bnd2()
        vcfspec_bnd2 = self.bnds.get_vcfspec_bnd2()
        varianthandler.apply_vcfspec(vr_bnd2, vcfspec_bnd2)
        vr_bnd2.info['MATEID'] = self.vr.id

        return vr_bnd2

    ##################################################################

    def set_bnds_attributes(self):
        if self.is_sv:
            vr_svinfo = breakends_module.get_vr_svinfo_standard_vr(
                    self.vr, self.fasta, self.chromdict,
                    )
            self.is_bnd1 = vr_svinfo['is_bnd1']
            self.bnds = breakends_module.get_bnds_from_vr_svinfo(
                self.vr, vr_svinfo, self.fasta, self.chromdict)
        else:
            self.is_bnd1 = None
            self.bnds = None

    def _refine_vr_InfoFormatValues(self):
        infoformat.refine_vr_InfoFormatValue(self.vr)

    ##################################################################

    def _set_VEPkeys(self, VEPkeys):
        if VEPkeys is None:
            self.VEPkeys = veplib.get_VEPkeys(self.vr.header)
        else:
            self.VEPkeys = VEPkeys

    def _set_vep_attributes(self):
        self._set_vep_attributes_annotformat()
        self._set_VEPannot_Nfeatures()
        self._set_VEPannot_spec()

    def _set_vep_attributes_annotformat(self):
        self.is_VEPannot, self.is_splitVEPannot = veplib.get_VEPannot_format(self.vr, self.VEPkeys)

    def _set_VEPannot_Nfeatures(self):
        if self.is_VEPannot:
            if self.is_splitVEPannot:
                val = self.vr.info['CSQ_' + self.VEPkeys[0]]
            else:
                val = self.vr.info['CSQ']

            if isinstance(val, tuple):
                self.VEPannot_Nfeatures = len(val)
            else:
                self.VEPannot_Nfeatures = 1
        else:
            self.VEPannot_Nfeatures = None

    def _set_VEPannot_spec(self):
        if self.is_VEPannot:
            VEPannot_list = self._make_VEPannot_list()
    
            # make canonical list
            canonical_VEPannot_list = [ x for x in VEPannot_list if x['CANONICAL'] == 'YES' ]
    
            # make VEPannot_dict from VEPannot_list
            VEPannot_dict = dict()
            canonical_VEPannot_dict = dict()
            for key in self.VEPkeys:
                VEPannot_dict[key] = tuple(( x[key] for x in VEPannot_list ))
                canonical_VEPannot_dict[key] = tuple(( x[key] for x in canonical_VEPannot_list ))
    
            self.VEPannot_list = VEPannot_list
            self.canonical_VEPannot_list = canonical_VEPannot_list
            self.VEPannot_dict = VEPannot_dict
            self.canonical_VEPannot_dict = canonical_VEPannot_dict
        else:
            self.VEPannot_list = None
            self.canonical_VEPannot_list = None
            self.VEPannot_dict = None
            self.canonical_VEPannot_dict = None

    def _make_VEPannot_list(self):
        VEPannot_list = list()
        for idx in range(self.VEPannot_Nfeatures):
            VEPannot_list.append(dict())

        # fill in values
        if self.is_splitVEPannot:
            for VEPkey in self.VEPkeys:
                val = infoformat.get_value_info(self.vr, 'CSQ_' + VEPkey)

                if not isinstance(val, tuple):
                    val = tuple([val])
                # val: ('GENE1', 'GENE2') (for CSQ_SYMBOL)
                for idx, subval in enumerate(val):
                    if subval is None:
                        VEPannot_list[idx][VEPkey] = None
                    else:
                        VEPannot_list[idx][VEPkey] = common.str_to_nonstr(subval)
        else:
            val = self.vr.info['CSQ']
            if not isinstance(val, tuple):
                val = tuple([val])
            # val: ('AA|AA|AA', 'BB|BB|BB')
            for idx, subval in enumerate(val):
                # subval: 'AA|AA|AA'
                for VEPkey, subsubval in zip(self.VEPkeys, subval.split('|')):
                    subsubval = re.sub('^[^:]+:::', '', subsubval) # for 'withkey' form of merged VEPannot format
                    # subsubval: 'AA'
                    if subsubval in infoformat.NA_VALUES:
                        VEPannot_list[idx][VEPkey] = None
                    else:
                        VEPannot_list[idx][VEPkey] = common.str_to_nonstr(subsubval)

        return VEPannot_list


def get_vp_sortkey(chromdict):
    vr_sortkey = common.get_vr_sortkey(chromdict)
    def vp_sortkey(vp):
        return vr_sortkey(vp.vr)

    return vp_sortkey
