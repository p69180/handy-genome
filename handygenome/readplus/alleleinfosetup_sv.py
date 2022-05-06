import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
breakends = importlib.import_module('.'.join([top_package_name, 'variantplus', 'breakends']))


#def check_rp_spans_bnds_unused(rp, bnds, pos_range0_bnd1, pos_range0_bnd2):
#    def subfun(rp, chrom_bnd, endis5, pos_range0):
#        if rp.read.reference_name == chrom_bnd:
#            if endis5:
#                if (((pos_range0.start - 1) in rp.pairs_dict['refpos0']) and
#                    ((pos_range0.stop - 1) in rp.pairs_dict['refpos0'])):
#                    traverses = True
#                    clippedat = False
#                else:
#                    traverses = False
#                    if rp.read.reference_start in pos_range0:
#                        clippedat = not (rp.pairs_dict['refpos0'][0] == 
#                            rp.read.reference_start)
#
#                        # this checks all preceding cigars are softclip
#    #                    refstart_index = rp.pairs_dict['refpos0'].index(
#    #                        rp.read.reference_start)
#    #                    cigarops_before = rp.pairs_dict['cigarop'][:refstart_index]
#    #                    clippedat = set(cigarops_before).issubset({4})
#                    else:
#                        clippedat = False
#            else:
#                if ((pos_range0.start in rp.pairs_dict['refpos0']) and
#                    (pos_range0.stop in rp.pairs_dict['refpos0'])):
#                    traverses = True
#                    clippedat = False
#                else:
#                    traverses = False
#                    if (rp.read.reference_end - 1) in pos_range0:
#                        clippedat = not (rp.pairs_dict['refpos0'][-1] ==
#                            rp.read.reference_end - 1)
#                    else:
#                        clippedat = False
#        else:
#            traverses = False
#            clippedat = False
#
#        return traverses, clippedat
#
#    traverses_bnd1, clippedat_bnd1 = subfun(rp, bnds.chrom_bnd1, 
#                                            bnds.endis5_bnd1, pos_range0_bnd1)
#    traverses_bnd2, clippedat_bnd2 = subfun(rp, bnds.chrom_bnd2, 
#                                            bnds.endis5_bnd2, pos_range0_bnd2)
#
#    return traverses_bnd1, clippedat_bnd1, traverses_bnd2, clippedat_bnd2


def check_read_spans_bnd(read, chrom_bnd, pos_range0_bnd, endis5_bnd):
    if read.reference_name == chrom_bnd:
        if endis5_bnd:
            if (
                    read.reference_end <= pos_range0_bnd.start or 
                    read.reference_start >= pos_range0_bnd.stop):
                clippedat = False
                traverses = False

        else:

            if endis5_bnd:
                if read.reference_start < pos_range0_bnd.start:
                    clippedat = False
                    traverses = True
                else:
                    if read.cigartuples[0][0] == 4:
                        clippedat = True
                        traverses = False
                    else:
                        clippedat = False
                        traverses = False
            else:
                if read.reference_end > pos_range0_bnd.stop:
                    clippedat = False
                    traverses = True
                else:
                    if read.cigartuples[-1][0] == 4:
                        clippedat = True
                        traverses = False
                    else:
                        clippedat = False
                        traverses = False
    else:
        clippedat = False
        traverses = False

    return clippedat, traverses
    

def check_rp_spans_bnd(rp, chrom_bnd, pos_bnd_range0, endis5):
    if rp.read.reference_name == chrom_bnd:
        if endis5:
            if rp.read.reference_end >= pos_bnd_range0.stop:
                if rp.read.reference_start <= (pos_bnd_range0.start - 1):
                    spans = True
                else:
                    if rp.pairs_dict['refpos0'][0] != rp.read.reference_start:
                        spans = True
                    else:
                        spans = False
            else:
                spans = False
        else:
            if rp.read.reference_start <= pos_bnd_range0.start:
                if rp.read.reference_end > pos_bnd_range0.stop:
                    spans = True
                else:
                    if rp.pairs_dict['refpos0'][-1] != (
                            rp.read.reference_end - 1):
                        spans = True
                    else:
                        spans = False
            else:
                spans = False
    else:
        spans = False

    return spans


#def get_supporting_bnds(rp):
#
#    def get_suppl_info(rp, SAitem):
#        """suppl_is5prime: 
#            Indicates whether the aligned portion in the supplementary 
#            alignment belongs to the 5-prime side (toward the first base) 
#            of the read.
#        suppl_isleft:
#            Indicates whether the aligned portion in the supplementary 
#            alignment belongs to the leftmost side, with regard to the
#            plus strand of the reference sequence.
#        """
#
#        if SAitem['is_forward']:
#            if (
#                    SAitem['cigartuples'][0][0] == 4 and 
#                    SAitem['cigartuples'][-1][0] == 0):
#                suppl_is5prime = False
#                suppl_isleft = False
#            elif (
#                    SAitem['cigartuples'][0][0] == 0 and 
#                    SAitem['cigartuples'][-1][0] == 4):
#                suppl_is5prime = True
#                suppl_isleft = True
#        else:
#            if (
#                    SAitem['cigartuples'][0][0] == 4 and 
#                    SAitem['cigartuples'][-1][0] == 0):
#                suppl_is5prime = True
#                suppl_isleft = False
#            elif (
#                    SAitem['cigartuples'][0][0] == 0 and 
#                    SAitem['cigartuples'][-1][0] == 4):
#                suppl_is5prime = False
#                suppl_isleft = True
#
#        return suppl_is5prime, suppl_isleft
#
#    def get_primary_clip_info(rp, suppl_is5prime):
#        """The softclip length in the primary alignment, 
#            corresponding to the supplementary alignment.
#        """
#
#        if suppl_is5prime:
#            clip_cigartuple = (rp.read.cigartuples[0]
#                               if rp.read.is_forward else
#                               rp.read.cigartuples[-1])
#        else:
#            clip_cigartuple = (rp.read.cigartuples[-1]
#                               if rp.read.is_forward else
#                               rp.read.cigartuples[0])
#
#        if clip_cigartuple[0] != 4:
#            raise Exception(
#                f'The cigarop of the primary read on the SA side is not '
#                f'a softclip.')
#
#        primary_cliplen = clip_cigartuple[1]
#
#        return primary_cliplen
#
#    # main
#    if rp.SAlist is None:
#        return None
#    
#    for SAitem in rp.SAlist:
#        suppl_is5prime, suppl_isleft = get_suppl_info(rp, SAitem)
#        primary_cliplen = get_primary_cliplen(rp, suppl_is5prime)
                
                
        
        
        
