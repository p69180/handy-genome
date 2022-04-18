import warnings

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))

variantplus = importlib.import_module('.'.join([top_package_name, 'variantplus', 'variantplus']))
variantpluspair = importlib.import_module('.'.join([top_package_name, 'variantplus', 'variantpluspair']))
varianthandler = importlib.import_module('.'.join([top_package_name, 'variantplus', 'varianthandler']))
headerhandler = importlib.import_module('.'.join([top_package_name, 'variantplus', 'headerhandler']))
initvcf = importlib.import_module('.'.join([top_package_name, 'vcfeditor', 'initvcf']))

breakends_module = importlib.import_module('.'.join([top_package_name, 'variantplus', 'breakends']))

#veplib = importlib.import_module('.'.join([top_package_name, 'annotation', 'veplib']))
annotationdb = importlib.import_module('.'.join([top_package_name, 'annotation', 'annotationdb']))


GREMLIN_HEADER = [
        'CHR1',
        'POS1',
        'CHR2',
        'POS2',
        'CT',
        'SVTYPE',
        'tumor_var_read',
        'tumor_var_split_read',
        'tumor_var_sa_read',
        'normal_var_read',
        'normal_var_split_read',
        'normal_same_clip',
        'tumor_ref_read_bpt1',
        'tumor_ref_read_bpt2',
        'p_value_ref_var_read',
        'tumor_vaf_bpt1',
        'tumor_vaf_bpt2',
        'tumor_var_mapq_bpt1',
        'tumor_var_mapq_bpt2',
        'tumor_depth_change_bpt1',
        'tumor_depth_change_bpt2',
        'normal_other_var_cluster_bpt1',
        'normal_other_var_cluster_bpt2',
        'normal_depth_bpt1',
        'normal_depth_bpt2',
        'normal_panel',
        'score',
        'mate_chr_diff_bp1',
        'mate_chr_same_bp1',
        'mate_chr_diff_prop_bp1',
        'mate_chr_diff_bp2',
        'mate_chr_same_bp2',
        'mate_chr_diff_prop_bp2',
        ]


class VcfPlus:
    """
    Attributes:
        vcf
        header
        has_good_contigs
        VEPkeys
        vplist
        vplist_sorted
        has_id
        has_duplicate_id
        has_mateid
        is_VEPannot
        vp_dict
        mateid_dict
        mateid_pairs
        vpp_list
        vpp_list_sorted
    """

    def __init__(self, vcf_path, fasta=None, chromdict=None, init_vp=True):
        """
        Args:
            fasta : pysam.FastaFile object
            vcf_path : Path to a vcf file
        """

        self.vcf_path = vcf_path
        self.set_fasta_chromdict(vcf_path, fasta, chromdict)
        self.init_vcf(vcf_path)

        if init_vp:
            self.init_vp_containers()

    ###########################################################

    def set_fasta_chromdict(self, vcf_path, fasta, chromdict):
        if chromdict is None:
            with pysam.VariantFile(vcf_path, 'r') as in_vcf:
                self.chromdict = common.ChromDict(vcfheader=in_vcf.header)
        else:
            self.chromdict = chromdict

        if fasta is None:
            refver = common.infer_refver(self.chromdict)
            self.fasta = pysam.FastaFile(
                common.DEFAULT_FASTA_PATH_DICT[refver])
        else:
            self.fasta = fasta

    def init_vcf(self, vcf_path):
        self.vcf = pysam.VariantFile(self.vcf_path, 'r')
        #self.header = self.vcf.header
        annotationdb.add_infometas(self.vcf.header)

    def init_vp_containers(self, vplist = None):
        self.set_vplist(vplist)
        #if check_VEPannot:
        #    self.set_is_VEPannot()
        #if set_details:
            #self.set_ID_attributes()
            #self.set_more_vp_containers()

    ###########################################################

    def set_header(self, header = None):
        if header is None:
            if self.vcf is None:
                self.header = None
            else:
                self.header = self.vcf.header
        else:
            self.header = header

    def _update_contig(self):
        if not check_having_good_contigs(self.vcf):
            self.vcf = get_updated_pysamvcf(self.vcf, chromdict = self.chromdict)

    def _set_has_good_contigs(self):
        self.has_good_contigs = check_having_good_contigs(self.vcf)
    
    def set_VEPkeys(self):
        if self.header is None:
            self.VEPkeys = None
        else:
            self.VEPkeys = veplib.get_VEPkeys(self.header)

    ###########################################################

    def set_vplist(self, vplist=None):
        """SV variant records for bnd2 are not loaded"""

        if vplist is None:
            self.vplist = list()
            for vr in self.vcf.fetch():
                vp = variantplus.VariantPlus(vr, self.fasta, self.chromdict)
                if vp.is_sv:
                    if vp.is_bnd1:
                        self.vplist.append(vp)
                else:
                    self.vplist.append(vp)
        else:
            self.vplist = vplist

        #self.vplist_sorted = self.check_vplist_is_sorted()

    def write_vcf(
            self, 
            outfile_path, 
            mode_bcftools = 'z', 
            mode_pysam = None, 
            nosort = False,
            ):
        mode_pysam = common.write_mode_arghandler(mode_bcftools, mode_pysam)

        out_vrlist = list()
        for vp in self.vplist:
            if vp.is_sv:
                vr_pos1, vr_pos2 = vp.get_vrs_for_write()
                out_vrlist.append(vr_pos1)
                out_vrlist.append(vr_pos2)
            else:
                out_vrlist.append(vp.vr)

        out_vrlist.sort(key = lambda vr: varianthandler.vr_sortkey(vr, self.chromdict))

        with pysam.VariantFile(outfile_path, mode = mode_pysam, header = self.vcf.header) as out_vcf:
            for vr in out_vrlist:
                out_vcf.write(vr)

    ###########################################################

    def set_ID_attributes(self):
        self.set_has_duplicate_pos_alleles()
        self.set_has_id()
        self.set_has_duplicate_id()
        self.set_has_mateid()
        self._sanity_check_dupID()

    def set_has_duplicate_pos_alleles(self):
        self.has_duplicate_pos_alleles = \
                varianthandler.check_has_duplicate_pos_alleles(vp.pysamvr for vp in self.vplist)

    def set_has_id(self):
        self.has_id = varianthandler.check_has_id(vp.pysamvr for vp in self.vplist)

    def set_has_duplicate_id(self):
        if self.has_id:
            self.has_duplicate_id = varianthandler.check_has_duplicate_id(vp.pysamvr for vp in self.vplist)
        else:
            self.has_duplicate_id = False

    def set_has_mateid(self):
        if self.has_id:
            self.has_mateid = varianthandler.check_has_mateid(vp.pysamvr for vp in self.vplist)
        else:
            self.has_mateid = False

    def _sanity_check_dupID(self):
        if self.has_duplicate_id:
            raise Exception(f'VCF file {self.vcf_path} has duplicate IDs.')

    ###########################################################

    def set_more_vp_containers(self):
        self.set_id_dict()
        self.set_mateid_dict()
        self._sanity_check_mateid()
        self.set_vpp_list()

    def set_id_dict(self):
        self.id_dict = dict()
        for vp in self.vplist:
            if vp.vr.id is not None:
                self.id_dict[vp.vr.id] = vp

    def set_mateid_dict(self):
        self.mateid_dict = dict()
        for vp in self.vplist:
            if vp.vr.id is not None:
                if vp.check_NA_info('MATEID'):
                    self.mateid_dict[vp.vr.id] = None
                else:
                    self.mateid_dict[vp.vr.id] = vp.get_value_info('MATEID')

    def _sanity_check_mateid(self):
        all_IDs = set(self.id_dict.keys())
        for ID, MATEID in self.mateid_dict.items():
            if MATEID is not None:
                if MATEID in all_IDs:
                    mate_of_mate = self.mateid_dict[MATEID]
                    if mate_of_mate != ID:
                        raise Exception(f'MATEID validity error. {ID} -> {MATEID} but {MATEID} -> {mate_of_mate}')
                else:
                    raise Exception(f'MATEID validity error. "{MATEID}", which is the mate of "{ID}", is not in the vcf.')

    def set_vpp_list(self):
        def get_vpp_from_one_SVvp(vp):
            new_vp1 = vp.make_maxspan_vp(pos1 = True)
            new_vp2 = vp.make_maxspan_vp(pos2 = True)
            return variantpluspair.VariantPlusPair([new_vp1, new_vp2])

        self.vpp_list = list()
        checked_IDs = set()

        for vp in self.vplist:
            if vp.pysamvr.id in checked_IDs:
                continue
            else:
                if vp.check_NA_info('MATEID'):
                    if varianthandler.check_SV(vp.pysamvr):
                        vpp = get_vpp_from_one_SVvp(vp)
                    else:
                        vpp = variantpluspair.VariantPlusPair([vp])

                    self.vpp_list.append(vpp)
                    checked_IDs.add(vp.pysamvr.id)
                else:
                    mateid = self.mateid_dict[vp.pysamvr.id]
                    mate_vp = self.id_dict[mateid]

                    if vp.bnds.sameseq(mate_vp.bnds):
                        vpp = variantpluspair.VariantPlusPair([vp, mate_vp])
                        self.vpp_list.append(vpp)
                    else:
                        vpp1 = get_vpp_from_one_SVvp(vp)
                        self.vpp_list.append(vpp1)
                        vpp2 = get_vpp_from_one_SVvp(mate_vp)
                        self.vpp_list.append(vpp2)

                    checked_IDs.add(vp.pysamvr.id)
                    checked_IDs.add(mateid)

        #self.vp_pair_list_sorted = self.check_vp_pair_list_is_sorted()

    ###########################################################

    def check_vplist_is_sorted(self):
        for idx in range(len(self.vplist) - 1):
            vp_pre = self.vplist[idx]
            vp_post = self.vplist[idx + 1]
            if common.get_order(
                    vp_pre.vcfspec.chrom,
                    vp_pre.vcfspec.pos,
                    vp_post.vcfspec.chrom,
                    vp_post.vcfspec.pos,
                    self.chromdict,
                    ) > 0:
                return False

        return True

    def sort_vplist(self):
        self.vplist.sort(
                key = lambda vp: common.coord_sortkey(
                    vp.pysamvr.contig,
                    vp.pysamvr.pos,
                    self.chromdict,
                    )
                )

    def check_vp_pair_list_is_sorted(self):
        for idx in range(len(self.vp_pair_list) - 1):
            vp_pair_pre = self.vp_pair_list[idx]
            vp_pair_post = self.vp_pair_list[idx + 1]
            if common.get_order(
                    vp_pair_pre['vp1'].pysamvr.contig,
                    vp_pair_pre['vp1'].pysamvr.pos,
                    vp_pair_post['vp1'].pysamvr.contig,
                    vp_pair_post['vp1'].pysamvr.pos,
                    self.chromdict,
                    ) > 0:
                return False
        return True

    def sort_vp_pair_list(self):
        """
        Sort by coordinate of breakend 1
        """
        self.vp_pair_list.sort(
                key = lambda x: varianthandler.pysamvr_sortkey(x[0].pysamvr, self.chromdict)
                )

    ###########################################################

                

    def write_gremlin(self, outfile_path):
        NAchar = '.'
        with open(outfile_path, 'w') as outfile:
            outfile.write('#' + '\t'.join(GREMLIN_HEADER) + '\n')
            for pair in self.vp_pair_list:
                vp = pair['vp1']
                line = list()

                line.append(vp.bnds.chrom1) # CHR1
                line.append(str(vp.bnds.pos1)) # POS1
                line.append(vp.bnds.chrom2) # CHR1
                line.append(str(vp.bnds.pos2)) # POS1

                CT1 = '5' if vp.bnds.endtype1_is5 else '3'
                CT2 = '5' if vp.bnds.endtype2_is5 else '3'
                CT = CT1 + 'to' + CT2
                line.append(CT) # CT

                line.append(vp.bnds.SVTYPE) # SVTYPE

                for field in GREMLIN_HEADER[6:]:
                    if field in vp.pysamvr.info.keys():
                        val = str( vp.pysamvr.info[field] )
                    else:
                        val = NAchar
                    line.append(val)

                outfile.write('\t'.join(line) + '\n')

    ###########################################################



def check_having_good_contigs(vcf):
    """
    Checks if CHROM names of all variant records are included in metadata contig names.
    Returns True if so, or returns False.

    Args:
        vcf: pysam.VariantFile instance
    """

    good = True

    chromdict = common.ChromDict(pysamhdr = vcf.header)
    for vr in vcf.fetch():
        if vr.contig in chromdict.contigs:
            if vr.pos in range(1, chromdict[vr.contig]+1):
                continue
            else:
                good = False
                break
        else:
            good = False
            break

    return good



######################################


def get_updated_pysamvcf(
        vcf, chromdict = None, samples = None, pysamhdr = None, vcfver = common.DEFAULT_VCFVER,
        ):
    """
    Writes a new temporary vcf file (removed at the end) with the same contents as input pysam.VariantFile except contig headers.
    Returns a pysam.VariantFile object generated with the newly written vcf file.
    """

    tf_path = common.get_tempfile_path(delete = False)
    write_updated_vcf(tf_path, vcf, chromdict, samples, pysamhdr, vcfver)
    new_vcf = pysam.VariantFile(tf_path, 'r')
    os.remove(tf_path)

    return new_vcf


def write_updated_vcf(
        outfile_path, vcf, chromdict = None, samples = None, pysamhdr = None, vcfver = common.DEFAULT_VCFVER, 
        mode_bcftools = 'z', mode_pysam = None,
        ):
    mode_pysam = common.write_mode_arghandler(mode_bcftools, mode_pysam)
    if pysamhdr is None:
        header = initvcf.create_header(chromdict, samples, vcf.header, vcfver)
    else:
        header = initvcf.create_header(chromdict, samples, pysamhdr, vcfver)

    with pysam.VariantFile(outfile_path, mode = mode_pysam, header = header) as out_vcf:
        for vr in vcf.fetch():
            out_vcf.write(varianthandler.reform_samples(vr, pysamhdr = header))
