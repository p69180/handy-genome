import re
import os
import stat
import shutil

import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
workflow = importlib.import_module('.'.join([top_package_name, 'workflow']))
annotationdb = importlib.import_module('.'.join([top_package_name, 'annotation', 'annotationdb']))
split = importlib.import_module('.'.join([top_package_name, 'vcfeditor', 'split']))
concat = importlib.import_module('.'.join([top_package_name, 'vcfeditor', 'concat']))


def unit_job(
        infile_path,
        outfile_path,
        refver,
        mode_bcftools = common.DEFAULT_MODE_BCFTOOLS, 
        mode_pysam = None, 
        ):
    mode_pysam = common.write_mode_arghandler(mode_bcftools, mode_pysam)

    with \
        pysam.FastaFile(common.DEFAULT_FASTA_PATH_DICT[refver]) as fasta, \
        pysam.VariantFile(infile_path) as in_vcf:

        annotationdb.add_infometas(in_vcf.header)
        with pysam.VariantFile(outfile_path, mode = mode_pysam, header = in_vcf.header) as out_vcf:
            for vr in in_vcf.fetch():
                cosmic_vr, cosmic_vr_noncoding = \
                    annotationdb.fetch_cosmic_vr(vr, refver = refver, fasta = fasta, search_equivs = True)
                cosmic = annotationdb.parse_cosmic_vr(cosmic_vr, cosmic_vr_noncoding)
                cosmic.write(vr, addkey = False)

                out_vcf.write(vr)


def argument_parser(cmdargs):
    def sanity_check(args):
        pass

    parser_dict = workflow.init_parser()

    workflow.add_infile_arg(parser_dict['required'], required = True)
    workflow.add_outfile_arg(parser_dict['required'], required = True, must_not_exist = 'ask')
    workflow.add_refver_arg(parser_dict['required'], required = True, choices = 'human')
    workflow.add_outfmt_arg(parser_dict['optional'], required = False)

    workflow.add_logging_args(parser_dict)
    workflow.add_scheduler_args(
            parser_dict,
            default_parallel = 1,
            default_sched = 'slurm',
            )
    workflow.add_rmtmp_arg(parser_dict)

    args = parser_dict['main'].parse_args(cmdargs)
    sanity_check(args)

    return args


def make_tmpdir(infile_path):
    tmpdir_paths = common.get_tmpdir_paths(
            ['scripts', 'logs', 'split_infiles', 'split_outfiles'], 
            prefix = os.path.basename(infile_path) + '_annotate_cosmic',
            where = os.path.dirname(infile_path),
            )

    return tmpdir_paths


def split_infile(infile_path, tmpdir_paths, parallel):
    split_filenames = split.main(
            vcf_path = infile_path,
            outdir = tmpdir_paths['split_infiles'],
            n_file = parallel,
            mode_bcftools = 'z',
            prefix = '',
            suffix = '.vcf.gz',
            )

    return split_filenames


def write_jobscripts(
        tmpdir_paths, 
        split_filenames, 
        refver,
        mode_pysam,
        jobname_prefix = 'annotate_cosmic',
        ncore_perjob = 1, 
        ):
    jobscript_paths = list()
    split_outfile_names = list()

    for zidx, infile_path in common.zenumerate(split_filenames):
        outfile_path = os.path.join(
                tmpdir_paths['split_outfiles'],
                re.sub('.vcf.gz$', '.cosmic.vcf.gz', os.path.basename(infile_path)),
                )
        split_outfile_names.append(outfile_path)

        script_path = os.path.join(tmpdir_paths['scripts'], f'{zidx}.sbatch')
        jobscript_paths.append(script_path)

        logpath = os.path.join(tmpdir_paths['logs'], os.path.basename(script_path) + '.log')

        with open(script_path, 'w') as outfile:
            outfile.write(f'''\
#!{common.PYTHON}

#SBATCH -N 1
#SBATCH -n 1

#SBATCH -c {ncore_perjob}
#SBATCH -o {logpath}
#SBATCH -e {logpath}
#SBATCH -J {jobname_prefix}_{zidx}

import sys
sys.path.append('{common.PACKAGE_LOCATION}')
from {top_package_name}.tools.annotate_cosmic import unit_job

unit_job(
    infile_path = '{infile_path}',
    outfile_path = '{outfile_path}',
    refver = '{refver}',
    mode_pysam = '{mode_pysam}', 
    )
''')
        os.chmod(script_path, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)

    return jobscript_paths, split_outfile_names


def concat_results(split_outfile_names, outfile_path, mode_pysam):
    concat.main(
        infile_path_list = split_outfile_names,
        outfile_path = outfile_path,
        mode_pysam = mode_pysam,
        outfile_must_not_exist = 'no',
        )


def main(cmdargs):
    args = argument_parser(cmdargs)
    logger = workflow.get_logger(
            name = 'annotate_cosmic',
            stderr = (not args.silent),
            filename = args.log,
            append = False,
            )

    tmpdir_paths = make_tmpdir(args.infile_path)
    split_filenames = split_infile(args.infile_path, tmpdir_paths, args.parallel)

    jobscript_paths, split_outfile_names = write_jobscripts(
            tmpdir_paths, 
            split_filenames, 
            refver = args.refver,
            mode_pysam = args.mode_pysam,
            )
    success, exitcode_list = workflow.run_jobs(
            jobscript_paths, 
            args.sched, 
            args.intv_check, 
            args.intv_submit, 
            logger, 
            tmpdir_paths['logs'],
            raise_on_failure = True,
            )

    concat_results(split_outfile_names, args.outfile_path, args.mode_pysam)

    if not args.dont_rm_tmp:
        shutil.rmtree(tmpdir_paths['top'])

    logger.info('ALL SUCCESSFULLY FINISHED')

