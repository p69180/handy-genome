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
split = importlib.import_module('.'.join([top_package_name, 'vcfeditor', 'split']))
concat = importlib.import_module('.'.join([top_package_name, 'vcfeditor', 'concat']))
veplib = importlib.import_module('.'.join([top_package_name, 'annotation', 'veplib']))
varianthandler = importlib.import_module('.'.join([top_package_name, 'variantplus', 'varianthandler']))
breakends = importlib.import_module('.'.join([top_package_name, 'variantplus', 'breakends']))
annotation = importlib.import_module('.'.join([top_package_name, 'annotation']))
annotationdb = importlib.import_module('.'.join([top_package_name, 'annotation', 'annotationdb']))
ensembl_parser = importlib.import_module('.'.join([top_package_name, 'annotation', 'ensembl_parser']))


# unit job

def unit_job(split_infile_path, split_outfile_path, fasta_path, species,
             assembly, distance, refver):
    vepinput_path = re.sub('\.vcf\.gz$', '.vepinput.vcf.gz', split_infile_path)
    vepoutput_path = re.sub('\.vcf\.gz$', '.vepoutput.vcf.gz', 
                            split_infile_path)

    fasta = pysam.FastaFile(fasta_path)
    chromdict = common.ChromDict(fasta=fasta)
    dbsnp_vcf = annotation.VCFS_DBSNP[refver]
    cosmic_coding_vcf = annotation.VCFS_COSMIC[refver]['coding']
    cosmic_noncoding_vcf = annotation.VCFS_COSMIC[refver]['noncoding']
    tabixfile_geneset = annotation.TABIXFILES_GENESET[refver]
    tabixfile_regulatory = annotation.TABIXFILES_REGULATORY[refver]
    tabixfile_repeats = annotation.TABIXFILES_REPEATS[refver]

    make_vepinput(split_infile_path, vepinput_path, fasta, chromdict)
    run_vep(vepinput_path, vepoutput_path, fasta_path, species, assembly,
            distance)
    add_annotations(
        split_infile_path, split_outfile_path, vepoutput_path, 
        fasta, distance, refver, chromdict,
        dbsnp_vcf, cosmic_coding_vcf, 
        cosmic_noncoding_vcf, tabixfile_geneset, 
        tabixfile_regulatory, tabixfile_repeats)


def make_vepinput(split_infile_path, vepinput_path, fasta, chromdict):
    def handle_nonsv(vr, in_vcf, out_vr_list):
        new_vr = in_vcf.header.new_record()
        new_vr.contig = vr.contig
        new_vr.pos = vr.pos
        new_vr.alleles = vr.alleles
        new_vr.id = varianthandler.get_vcfspec(vr).get_id()
        out_vr_list.append(new_vr)

    def handle_sv(vr, in_vcf, out_vr_list, bnds, fasta):
        vr_pos1 = in_vcf.header.new_record()
        vr_pos1.contig = bnds.chrom1
        vr_pos1.pos = bnds.pos1
        vr_pos1.ref = fasta.fetch(bnds.chrom1, bnds.pos1 - 1, bnds.pos1)
        vr_pos1.alts = ('N',)
        vr_pos1.id = bnds.get_id_pos1()
        out_vr_list.append(vr_pos1)

        vr_pos2 = in_vcf.header.new_record()
        vr_pos2.contig = bnds.chrom2
        vr_pos2.pos = bnds.pos2
        vr_pos2.ref = fasta.fetch(bnds.chrom2, bnds.pos2 - 1, bnds.pos2)
        vr_pos2.alts = ('N',)
        vr_pos2.id = bnds.get_id_pos2()
        out_vr_list.append(vr_pos2)

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
        split_infile_path, split_outfile_path, vepoutput_path, 
        fasta, distance, refver, chromdict,
        dbsnp_vcf, cosmic_coding_vcf, cosmic_noncoding_vcf,
        tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats):
    def get_vep_vr_dict(vepoutput_path):
        vep_vr_dict = dict()
        with pysam.VariantFile(vepoutput_path) as in_vcf:
            for vr in in_vcf.fetch():
                vep_vr_dict[vr.id] = vr

        return vep_vr_dict

    def add_vep_sv(vep_vr_dict, bnds, annotdb_bnd1, annotdb_bnd2,
                   tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats, 
                   distance):
        vep_vr_bnd1 = vep_vr_dict[bnds.get_id_pos1()]
        annotdb_bnd1.update_cmdline_vep(
            vep_vr_bnd1, overwrite=False, create_new=True)
        annotdb_bnd1.update_features_postvep_bnd(
            bnds.chrom1, bnds.pos1, bnds.pos1_endis5, 
            tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats,
            distance=distance)

        vep_vr_bnd2 = vep_vr_dict[bnds.get_id_pos2()]
        annotdb_bnd2.update_cmdline_vep(
            vep_vr_bnd2, overwrite=True, create_new=True)
        annotdb_bnd2.update_features_postvep_bnd(
            bnds.chrom2, bnds.pos2, bnds.pos2_endis5, 
            tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats,
            distance=distance)

    def add_vep_nonsv(vep_vr_dict, vcfspec, annotdb, tabixfile_geneset, 
                      tabixfile_regulatory, tabixfile_repeats, distance):
        vep_vr = vep_vr_dict[vcfspec.get_id()] 
        annotdb.update_cmdline_vep(vep_vr, overwrite=False, create_new=True)
        annotdb.update_features_postvep_plain(
            vcfspec, tabixfile_geneset, tabixfile_regulatory, 
            tabixfile_repeats, distance)

    ###################

    def add_annotations_to_vr_sv(
            vr, vep_vr_dict, bnds, distance, refver, fasta,
            tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats):
        annotdb_bnd1 = annotationdb.AnnotDB('bnd1', refver, fasta)
        annotdb_bnd2 = annotationdb.AnnotDB('bnd2', refver, fasta)
        add_vep_sv(vep_vr_dict, bnds, annotdb_bnd1, annotdb_bnd2,
                   tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats, 
                   distance)
        annotdb_bnd1.write(vr, addkey=False)
        annotdb_bnd2.write(vr, addkey=False)

    def add_annotations_to_vr_nonsv(
            vr, vep_vr_dict, distance, refver, fasta,
            dbsnp_vcf, cosmic_coding_vcf, cosmic_noncoding_vcf, 
            tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats):
        annotdb = annotationdb.AnnotDB('plain', refver, fasta)
        vcfspec = varianthandler.get_vcfspec(vr)
        # vep
        add_vep_nonsv(vep_vr_dict, vcfspec, annotdb, tabixfile_geneset, 
                      tabixfile_regulatory, tabixfile_repeats, distance)
        # popfreq
        annotdb.update_popfreq(vcfspec, dbsnp_vcf, search_equivs=True, 
                               overwrite=False)
        # cosmic
        annotdb.update_cosmic(vcfspec, cosmic_coding_vcf, cosmic_noncoding_vcf,  
                              search_equivs=True, overwrite=False)

        annotdb.write(vr, addkey=False)

    def add_annotations_to_vr(
            vr, vep_vr_dict, distance, fasta, chromdict,
            dbsnp_vcf, cosmic_coding_vcf, cosmic_noncoding_vcf, 
            tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats):
        if varianthandler.check_SV(vr):
            vr_svinfo = breakends.get_vr_svinfo_standard_vr(vr, fasta, 
                                                            chromdict)
            if vr_svinfo['is_bnd1']:
                bnds = breakends.get_bnds_from_vr_svinfo(vr, vr_svinfo, fasta, 
                                                         chromdict)
                add_annotations_to_vr_sv(
                    vr, vep_vr_dict, bnds, distance, refver, fasta,
                    tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats)
        else:
            add_annotations_to_vr_nonsv(
                vr, vep_vr_dict, distance, refver, fasta,
                dbsnp_vcf, cosmic_coding_vcf, cosmic_noncoding_vcf, 
                tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats)

    # main
    vep_vr_dict = get_vep_vr_dict(vepoutput_path)
    with pysam.VariantFile(split_infile_path) as in_vcf:
        annotationdb.add_infometas(in_vcf.header)
        with pysam.VariantFile(split_outfile_path, mode='wz', 
                               header=in_vcf.header) as out_vcf:
            for vr in in_vcf.fetch():
                add_annotations_to_vr(
                    vr, vep_vr_dict, distance, fasta, chromdict,
                    dbsnp_vcf, cosmic_coding_vcf, cosmic_noncoding_vcf, 
                    tabixfile_geneset, tabixfile_regulatory, tabixfile_repeats)
                out_vcf.write(vr)


#############################################################

def argument_parser(cmdargs):
    def sanity_check(args):
        pass

    def modifier(args):
        args.species = veplib.REFVER_TO_VEPARGS[args.refver]['species']
        args.assembly = veplib.REFVER_TO_VEPARGS[args.refver]['assembly']

    parser_dict = workflow.init_parser()

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
    workflow.add_rmtmp_arg(parser_dict)

    args = parser_dict['main'].parse_args(cmdargs)
    sanity_check(args)
    modifier(args)

    return args


def make_tmpdir(infile_path):
    tmpdir_paths = common.get_tmpdir_paths(
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
                     jobname_prefix='vep_annot', ncore_perjob=2):
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

        with open(jobscript_path, 'w') as outfile:
            outfile.write(textwrap.dedent(f"""\
                #!{common.PYTHON}

                #SBATCH -N 1
                #SBATCH -n 1

                #SBATCH -c {ncore_perjob}
                #SBATCH -o {logpath}
                #SBATCH -e {logpath}
                #SBATCH -J {jobname_prefix}_{zidx}

                import sys
                sys.path.append('{common.PACKAGE_LOCATION}')
                from {__name__} import unit_job

                unit_job(
                    split_infile_path='{split_infile_path}',
                    split_outfile_path='{split_outfile_path}',
                    fasta_path='{fasta_path}', 
                    species='{species}', 
                    assembly='{assembly}', 
                    distance={distance},
                    refver='{refver}',
                    )"""))

        os.chmod(jobscript_path, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)

    return jobscript_path_list, split_outfile_path_list


def concat_results(split_outfile_path_list, outfile_path, mode_pysam):
    concat.main(infile_path_list = split_outfile_path_list,
                outfile_path = outfile_path, mode_pysam = mode_pysam,
                outfile_must_not_exist = 'no')


def main(cmdargs):
    args = argument_parser(cmdargs)

    # setups logger
    logger = workflow.get_logger(name='vepannot', stderr=(not args.silent),
                                 filename=args.log, append=False)
    logger.info('Beginning')

    # setups other parameters
    fasta_path = common.DEFAULT_FASTA_PATH_DICT[args.refver]

    # splits the input file and run jobs
    tmpdir_paths = make_tmpdir(args.infile_path)
    logger.info('Splitting the input file')
    split_infile_path_list = split_infile(args.infile_path, tmpdir_paths, 
                                          args.parallel)
    jobscript_path_list, split_outfile_path_list = write_jobscripts(
        tmpdir_paths, split_infile_path_list, fasta_path, args.species, 
        args.assembly, args.distance, args.refver)
    logger.info('Running annotation jobs for each split file')
    workflow.run_jobs(jobscript_path_list, sched=args.sched, 
                      intv_check=args.intv_check, intv_submit=args.intv_submit, 
                      logger=logger, log_dir=tmpdir_paths['logs'],
                      raise_on_failure=True)

    # concatenates split files
    logger.info('Merging split files')
    concat_results(split_outfile_path_list, args.outfile_path, args.mode_pysam)

    if not args.dont_rm_tmp:
        shutil.rmtree(tmpdir_paths['top'])

    logger.info('All successfully finished')

