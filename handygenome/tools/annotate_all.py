import re
import os
import stat
import argparse
import subprocess
import shutil
import textwrap

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
workflow = importlib.import_module('.'.join([top_package_name, 'workflow']))
toolsetup = importlib.import_module('.'.join([top_package_name, 'workflow', 'toolsetup']))
split = importlib.import_module('.'.join([top_package_name, 'vcfeditor', 'split']))
concat = importlib.import_module('.'.join([top_package_name, 'vcfeditor', 'concat']))
veplib = importlib.import_module('.'.join([top_package_name, 'annotation', 'veplib']))
varianthandler = importlib.import_module('.'.join([top_package_name, 'variantplus', 'varianthandler']))
variantplus = importlib.import_module('.'.join([top_package_name, 'variantplus', 'variantplus']))
breakends = importlib.import_module('.'.join([top_package_name, 'variantplus', 'breakends']))
annotation = importlib.import_module('.'.join([top_package_name, 'annotation']))
annotationdb = importlib.import_module('.'.join([top_package_name, 'annotation', 'annotationdb']))
ensembl_parser = importlib.import_module('.'.join([top_package_name, 'annotation', 'ensembl_parser']))
indexing = importlib.import_module('.'.join([top_package_name, 'vcfeditor', 'indexing']))


def unit_job(split_infile_path, split_outfile_path, fasta_path, species,
             assembly, distance, refver, do_features, do_cosmic, do_popfreq):
    # basic setup
    fasta = pysam.FastaFile(fasta_path)
    chromdict = common.ChromDict(fasta=fasta)
    # setup for features 
    tabixfile_geneset = annotation.TABIXFILES_GENESET[refver]
    tabixfile_regulatory = annotation.TABIXFILES_REGULATORY[refver]
    tabixfile_repeats = annotation.TABIXFILES_REPEATS[refver]
    # run VEP
    if do_features:
        vep_vr_dict = setup_vep(split_infile_path, fasta_path, fasta, 
                                chromdict, species, assembly, distance)
    else:
        vep_vr_dict = None
    # setup for popfreq
    dbsnp_vcf = annotation.VCFS_DBSNP[refver]
    # setup for COSMIC
    cosmic_coding_vcf = annotation.VCFS_COSMIC[refver]['coding']
    cosmic_noncoding_vcf = annotation.VCFS_COSMIC[refver]['noncoding']

    add_annotations(
        split_infile_path, split_outfile_path, vep_vr_dict, 
        fasta, distance, refver, chromdict,
        dbsnp_vcf, cosmic_coding_vcf, cosmic_noncoding_vcf, 
        tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats,
        do_features, do_cosmic, do_popfreq)


def setup_vep(split_infile_path, fasta_path, fasta, chromdict, species,
              assembly, distance):
    vepinput_path = re.sub('\.vcf\.gz$', '.vepinput.vcf.gz', 
                           split_infile_path)
    vepoutput_path = re.sub('\.vcf\.gz$', '.vepoutput.vcf', 
                            split_infile_path)
        # VEP output is not compressed
    make_vepinput(split_infile_path, vepinput_path, fasta, chromdict)
    run_vep(vepinput_path, vepoutput_path, fasta_path, species, 
            assembly, distance)

    # get vep_vr_dict
    vep_vr_dict = dict()
    with pysam.VariantFile(vepoutput_path) as in_vcf:
        for vr in in_vcf.fetch():
            vep_vr_dict[vr.id] = vr

    return vep_vr_dict


def make_vepinput(split_infile_path, vepinput_path, fasta, chromdict):
    def handle_nonsv(vr, in_vcf, out_vr_list):
        new_vr = in_vcf.header.new_record()
        new_vr.contig = vr.contig
        new_vr.pos = vr.pos
        new_vr.alleles = vr.alleles
        new_vr.id = varianthandler.get_vcfspec(vr).get_id()
        out_vr_list.append(new_vr)

    def handle_sv(vr, in_vcf, out_vr_list, bnds, fasta):
        vr_bnd1 = in_vcf.header.new_record()
        vr_bnd1.contig = bnds.chrom_bnd1
        vr_bnd1.pos = bnds.pos_bnd1
        vr_bnd1.ref = fasta.fetch(bnds.chrom_bnd1, bnds.pos_bnd1 - 1, 
                                  bnds.pos_bnd1)
        vr_bnd1.alts = ('N',)
        vr_bnd1.id = bnds.get_id_bnd1()
        out_vr_list.append(vr_bnd1)

        vr_bnd2 = in_vcf.header.new_record()
        vr_bnd2.contig = bnds.chrom_bnd2
        vr_bnd2.pos = bnds.pos_bnd2
        vr_bnd2.ref = fasta.fetch(bnds.chrom_bnd2, bnds.pos_bnd2 - 1, 
                                  bnds.pos_bnd2)
        vr_bnd2.alts = ('N',)
        vr_bnd2.id = bnds.get_id_bnd2()
        out_vr_list.append(vr_bnd2)

    def handle_cpgmet(vr, in_vcf, out_vr_list):
        new_vr = in_vcf.header.new_record()
        new_vr.contig = vr.contig
        new_vr.pos = vr.pos
        new_vr.ref = vr.ref
        new_vr.alts = 'N'
        new_vr.id = varianthandler.get_vcfspec(vr).get_id()
        out_vr_list.append(new_vr)

    with pysam.VariantFile(split_infile_path) as in_vcf:
        with pysam.VariantFile(vepinput_path, mode='wz', 
                               header=in_vcf.header) as out_vcf:
            out_vr_list = list()
            for vr in in_vcf.fetch():
                if varianthandler.check_SV(vr):
                    vr_svinfo = breakends.get_vr_svinfo_standard_vr(
                        vr, fasta, chromdict)
                    if vr_svinfo['is_bnd1']:
                        bnds = breakends.get_bnds_from_vr_svinfo(
                            vr, vr_svinfo, fasta, chromdict)
                        handle_sv(vr, in_vcf, out_vr_list, bnds, fasta)
                elif varianthandler.check_cpgmet(vr):
                    handle_cpgmet(vr, in_vcf, out_vr_list)
                else:
                    handle_nonsv(vr, in_vcf, out_vr_list)

            out_vr_list.sort(key=common.get_vr_sortkey(chromdict))

            for vr in out_vr_list:
                out_vcf.write(vr)


def run_vep(vepinput_path, vepoutput_path, fasta_path, species, assembly,
            distance):
    vep_args = veplib.get_vep_args(vepinput_path, vepoutput_path,
                                   fasta_path, species, assembly,
                                   distance=distance,
                                   force_overwrite=False)
    p = subprocess.run(vep_args, text=True, capture_output=True)
    if p.returncode != 0:
        raise Exception(f'VEP run failed. Error message:\n{p.stderr}')


def add_annotations(
        split_infile_path, split_outfile_path, vep_vr_dict, 
        fasta, distance, refver, chromdict,
        dbsnp_vcf, cosmic_coding_vcf, cosmic_noncoding_vcf,
        tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats,
        do_features, do_cosmic, do_popfreq):
    def add_vep_sv(vep_vr_bnd1, vep_vr_bnd2, bnds, annotdb_bnd1, annotdb_bnd2,
                   tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats, 
                   distance):
        annotdb_bnd1.update_cmdline_vep(
            vep_vr_bnd1, overwrite=True, create_new=True)
        annotdb_bnd1.update_features_postvep_bnd(
            bnds.chrom_bnd1, bnds.pos_bnd1, bnds.endis5_bnd1, 
            tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats,
            distance=distance)

        annotdb_bnd2.update_cmdline_vep(
            vep_vr_bnd2, overwrite=True, create_new=True)
        annotdb_bnd2.update_features_postvep_bnd(
            bnds.chrom_bnd2, bnds.pos_bnd2, bnds.endis5_bnd2, 
            tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats,
            distance=distance)

    def add_vep_nonsv(vep_vr, vcfspec, annotdb, tabixfile_geneset, 
                      tabixfile_regulatory, tabixfile_repeats, distance):
        annotdb.update_cmdline_vep(vep_vr, overwrite=True, create_new=True)
        annotdb.update_features_postvep_plain(
            vcfspec, tabixfile_geneset, tabixfile_regulatory, 
            tabixfile_repeats, distance)

    ###################

    def add_to_vr_cpgmet(
            vr, vep_vr_dict, distance, refver, fasta, chromdict,
            tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats,
            do_features, do_cosmic, do_popfreq):
        vp = variantplus.VariantPlus(vr=vr, refver=refver, fasta=fasta,
                                     chromdict=chromdict, set_annotdb=True,
                                     set_readstats=False)
        annotdb = vp.annotdb
        vcfspec = vp.vcfspec

        if do_features:
            vep_vr = vep_vr_dict[vcfspec.get_id()] 
            add_vep_nonsv(vep_vr, vcfspec, annotdb, tabixfile_geneset, 
                          tabixfile_regulatory, tabixfile_repeats, distance)

        # write
        annotdb.write(vr, add_meta=False)

    def add_to_vr_sv(
            vr, vep_vr_dict, bnds, distance, refver, fasta, chromdict,
            tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats,
            do_features, do_cosmic, do_popfreq):
        vp = variantplus.VariantPlus(vr=vr, refver=refver, fasta=fasta,
                                     chromdict=chromdict, set_annotdb=True,
                                     set_readstats=False)
        annotdb_bnd1 = vp.annotdb_bnd1
        annotdb_bnd2 = vp.annotdb_bnd2

        if do_features:
            vep_vr_bnd1 = vep_vr_dict[bnds.get_id_bnd1()]
            vep_vr_bnd2 = vep_vr_dict[bnds.get_id_bnd2()]
            add_vep_sv(vep_vr_bnd1, vep_vr_bnd2, bnds, annotdb_bnd1, 
                       annotdb_bnd2, tabixfile_geneset, tabixfile_regulatory, 
                       tabixfile_repeats, distance)

        # write
        annotdb_bnd1.write(vr, add_meta=False)
        annotdb_bnd2.write(vr, add_meta=False)

    def add_to_vr_nonsv(
            vr, vep_vr_dict, distance, refver, fasta, chromdict,
            dbsnp_vcf, cosmic_coding_vcf, cosmic_noncoding_vcf, 
            tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats,
            do_features, do_cosmic, do_popfreq):
        vp = variantplus.VariantPlus(vr=vr, refver=refver, fasta=fasta,
                                     chromdict=chromdict, set_annotdb=True,
                                     set_readstats=False)
        annotdb = vp.annotdb
        vcfspec = vp.vcfspec

        if do_features:
            vep_vr = vep_vr_dict[vcfspec.get_id()] 
            add_vep_nonsv(vep_vr, vcfspec, annotdb, tabixfile_geneset, 
                          tabixfile_regulatory, tabixfile_repeats, distance)
        if do_popfreq:
            annotdb.update_popfreq(vcfspec, dbsnp_vcf, search_equivs=True, 
                                   overwrite=True)
        if do_cosmic:
            annotdb.update_cosmic(vcfspec, cosmic_coding_vcf, 
                                  cosmic_noncoding_vcf,
                                  search_equivs=True, overwrite=True)
        # write
        annotdb.write(vr, add_meta=False)

    def add_to_vr(
            vr, vep_vr_dict, distance, fasta, chromdict,
            dbsnp_vcf, cosmic_coding_vcf, cosmic_noncoding_vcf, 
            tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats,
            do_features, do_cosmic, do_popfreq):
        if varianthandler.check_SV(vr):
            vr_svinfo = breakends.get_vr_svinfo_standard_vr(vr, fasta, 
                                                            chromdict)
            if vr_svinfo['is_bnd1']:
                bnds = breakends.get_bnds_from_vr_svinfo(vr, vr_svinfo, fasta, 
                                                         chromdict)
                add_to_vr_sv(
                    vr, vep_vr_dict, bnds, distance, refver, fasta, chromdict,
                    tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats,
                    do_features, do_cosmic, do_popfreq)
        elif varianthandler.check_cpgmet(vr):
            add_to_vr_cpgmet(
                vr, vep_vr_dict, distance, refver, fasta, chromdict,
                tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats,
                do_features, do_cosmic, do_popfreq)
        else:
            add_to_vr_nonsv(
                vr, vep_vr_dict, distance, refver, fasta, chromdict,
                dbsnp_vcf, cosmic_coding_vcf, cosmic_noncoding_vcf, 
                tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats,
                do_features, do_cosmic, do_popfreq)

    # main
    with pysam.VariantFile(split_infile_path) as in_vcf:
        annotationdb.add_infometas(in_vcf.header)
        with pysam.VariantFile(split_outfile_path, mode='wz', 
                               header=in_vcf.header) as out_vcf:
            for vr in in_vcf.fetch():
                add_to_vr(
                    vr, vep_vr_dict, distance, fasta, chromdict,
                    dbsnp_vcf, cosmic_coding_vcf, cosmic_noncoding_vcf, 
                    tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats,
                    do_features, do_cosmic, do_popfreq)
                out_vcf.write(vr)


#############################################################

def argument_parser(cmdargs):
    def sanity_check(args):
        pass

    def postprocess(args):
        args.species = veplib.REFVER_TO_VEPARGS[args.refver]['species']
        args.assembly = veplib.REFVER_TO_VEPARGS[args.refver]['assembly']

        if not any([args.do_features, args.do_cosmic, args.do_popfreq]):
            args.do_features = True
            args.do_cosmic = True
            args.do_popfreq = True

    parser_dict = workflow.init_parser(
        description=(
            f'Adds transcript, regulatory elements, repeat sequences, COSMIC '
            f'information, and population frequency to the input vcf file. '
            f'If one or more of "--do-features", "--do-cosmic", '
            f'or "--do-popfreq" is set, only those items are added; '
            f'if none of them is set, all is done.'))

    # required
    workflow.add_infile_arg(parser_dict['required'], required=True)
    workflow.add_outfile_arg(parser_dict['required'], required=True, 
                             must_not_exist='ask')
    workflow.add_refver_arg(parser_dict['required'], required=True, 
                            choices=('hg19', 'hg38'))

    # optional
    workflow.add_outfmt_arg(parser_dict['optional'], required=False)
    parser_dict['optional'].add_argument(
        '--distance', required=False, default=5000,
        help=(f'Distance (in bp) from the mutation location to fetch feature '
              f'information. Used for VEP and custom file annotation.'))
    workflow.add_logging_args(parser_dict)
    workflow.add_scheduler_args(parser_dict, default_parallel=1, 
                                default_sched='slurm')

    # flag
    workflow.add_rmtmp_arg(parser_dict)
    workflow.add_index_arg(parser_dict)
    parser_dict['flag'].add_argument(
        '--do-features', dest='do_features', action='store_true',
        help=(f'If set, runs vep and annotates geneset, regulatory '
              f'element, and repeat element information.'))
    parser_dict['flag'].add_argument(
        '--do-cosmic', dest='do_cosmic', action='store_true',
        help=f'If set, annotates COSMIC information.')
    parser_dict['flag'].add_argument(
        '--do-popfreq', dest='do_popfreq', action='store_true',
        help=(f'If set, annotates population frequency information from '
             f'dbSNP database.'))

    args = parser_dict['main'].parse_args(cmdargs)
    sanity_check(args)
    postprocess(args)

    return args


def make_tmpdir(infile_path):
    tmpdir_paths = workflow.get_tmpdir_paths(
        ['scripts', 'logs', 'split_infiles', 'split_outfiles'],
        prefix=(os.path.basename(infile_path) + '_vepannot'),
        where=os.path.dirname(infile_path))

    return tmpdir_paths


def split_infile(infile_path, tmpdir_paths, parallel):
    split_filenames = split.main(vcf_path=infile_path, 
                                 outdir=tmpdir_paths['split_infiles'], 
                                 n_file=parallel, mode_bcftools='z', 
                                 prefix='', suffix='.vcf.gz')

    return split_filenames


def write_jobscripts(tmpdir_paths, split_infile_path_list, fasta_path, 
                     species, assembly, distance, refver,
                     do_features, do_cosmic, do_popfreq,
                     jobname_prefix='annotate_all', ncore_perjob=2):
    jobscript_path_list = list()
    split_outfile_path_list = list()

    for zidx, split_infile_path in common.zenumerate(split_infile_path_list):
        basename = os.path.basename(split_infile_path)

        split_outfile_path = os.path.join(tmpdir_paths['split_outfiles'], 
                                          basename)
        split_outfile_path_list.append(split_outfile_path)

        jobscript_path = os.path.join(tmpdir_paths['scripts'], 
                                      f'{zidx}.sbatch')
        jobscript_path_list.append(jobscript_path)
        logpath = os.path.join(tmpdir_paths['logs'], 
                               os.path.basename(jobscript_path) + '.log')
        success_logpath = re.sub('\.log$', '.success', logpath)
        failure_logpath = re.sub('\.log$', '.failure', logpath)

        script_contents = textwrap.dedent(f"""\
            #!{common.PYTHON}

            #SBATCH -N 1
            #SBATCH -n 1

            #SBATCH -c {ncore_perjob}
            #SBATCH -o {os.devnull}
            #SBATCH -e {os.devnull}
            #SBATCH -J {jobname_prefix}_{zidx}

            import os
            import contextlib
            import traceback
            import sys
            sys.path.append('{common.PACKAGE_LOCATION}')
            from {__name__} import unit_job

            log = open('{logpath}', 'w')
            with \\
                    contextlib.redirect_stdout(log), \\
                    contextlib.redirect_stderr(log):
                try:
                    unit_job(
                        split_infile_path='{split_infile_path}',
                        split_outfile_path='{split_outfile_path}',
                        fasta_path='{fasta_path}', 
                        species='{species}', 
                        assembly='{assembly}', 
                        distance={distance},
                        refver='{refver}',
                        do_features={do_features},
                        do_cosmic={do_cosmic},
                        do_popfreq={do_popfreq},
                        )
                except:
                    print(traceback.format_exc())
                    success = False
                else:
                    success = True
            log.close()

            if success:
                os.rename('{logpath}', '{success_logpath}')
            else:
                os.rename('{logpath}', '{failure_logpath}')
                raise SystemExit(1)
            """)

        with open(jobscript_path, 'w') as outfile:
            outfile.write(script_contents)

        os.chmod(jobscript_path, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)

    return jobscript_path_list, split_outfile_path_list


def concat_results(split_outfile_path_list, outfile_path, mode_pysam):
    concat.main(infile_path_list = split_outfile_path_list,
                outfile_path = outfile_path, mode_pysam = mode_pysam,
                outfile_must_not_exist = 'no')


def main(cmdargs):
    args = argument_parser(cmdargs)

    # make temporary directory tree
    tmpdir_paths = make_tmpdir(args.infile_path)

    # setup logger
    logger = toolsetup.setup_logger(args=args, 
                                    tmpdir_root=tmpdir_paths['root'],
                                    with_genlog=True)
    logger.info('Beginning')

    # setup other parameters
    fasta_path = common.DEFAULT_FASTA_PATHS[args.refver]

    # split the input file and run jobs
    logger.info('Splitting the input file')
    split_infile_path_list = split_infile(args.infile_path, tmpdir_paths, 
                                          args.parallel)
    jobscript_path_list, split_outfile_path_list = write_jobscripts(
        tmpdir_paths, split_infile_path_list, fasta_path, args.species, 
        args.assembly, args.distance, args.refver,
        args.do_features, args.do_cosmic, args.do_popfreq)
    logger.info('Running annotation jobs for each split file')
    workflow.run_jobs(jobscript_path_list, sched=args.sched, 
                      intv_check=args.intv_check, 
                      intv_submit=args.intv_submit, 
                      logger=logger, log_dir=tmpdir_paths['logs'],
                      raise_on_failure=True)

    # concatenate split files
    logger.info('Merging split files')
    concat_results(split_outfile_path_list, args.outfile_path, 
                   args.mode_pysam)

    # rmtmp, index
    if not args.donot_rm_tmp:
        shutil.rmtree(tmpdir_paths['root'])
    if not args.donot_index:
        indexing.index_vcf(args.outfile_path)

    logger.info('All successfully finished')

