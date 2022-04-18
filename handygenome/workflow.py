import sys
import os
import re
import argparse
import subprocess
import time
import logging
import uuid
import textwrap

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))


SLURMBIN = '/usr/local/slurm/bin'
SQUEUE = os.path.join(SLURMBIN, 'squeue')
SBATCH = os.path.join(SLURMBIN, 'sbatch')
SCONTROL = os.path.join(SLURMBIN, 'scontrol')
SINFO = os.path.join(SLURMBIN, 'sinfo')

PBSBIN = '/usr/local/pbs/bin/'
PBSNODES = os.path.join(PBSBIN, 'pbsnodes')
QSTAT = os.path.join(PBSBIN, 'qstat')

RE_PATS = dict()
RE_PATS['scontrol_split_job'] = re.compile('^([^=]+)(=(.*))?$')
RE_PATS['scontrol_split_nodes'] = re.compile(r'(?<=\S)\s+(?=\S+=)')

DEFAULT_INTV_CHECK = 60 # seconds
DEFAULT_INTV_SUBMIT = 1 # seconds


# logging #

DEFAULT_DATEFMT = '%Z %Y-%m-%d %H:%M:%S' # KST 2022-03-23 22:12:34
DEFAULT_LOG_FORMATTERS = {
    'without_name': logging.Formatter(
        fmt='[%(asctime)s %(levelname)s] %(message)s', 
        datefmt=DEFAULT_DATEFMT),
    'with_name': logging.Formatter(
        fmt='[%(asctime)s %(levelname)s] %(name)s: %(message)s', 
        datefmt=DEFAULT_DATEFMT),
    }


def get_logger(name=None, formatter=None, level='info', stderr=True, 
               filename=None, append=False):
    if name is None:
        name = str(uuid.uuid4())

    if formatter is None:
        formatter = DEFAULT_LOG_FORMATTERS['with_name']

    loglevel = getattr(logging, level.upper())

    logger = logging.getLogger(name)
    logger.setLevel(loglevel)
    logger.propagate = False

    if stderr:
        sh = logging.StreamHandler()
        sh.setLevel(loglevel)
        sh.setFormatter(formatter)
        logger.addHandler(sh)

    if filename is not None:
        fh = logging.FileHandler(filename, mode=('a' if append else 'w'))
        fh.setLevel(loglevel)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


################################################################

def get_split_filenames(n_file, outdir, prefix, suffix):
    """
    Returns:
        A list of file names each of which looks like 
            <outdir>/<prefix>000<suffix>.
        File creation is not done.
    """

    result = list()
    padded_indices = common.get_padded_indices(n_file)
    result = [os.path.join(outdir, prefix + idx_pad + suffix)
              for idx_pad in padded_indices]
    
    return result

get_indexed_filenames = get_split_filenames  # alias


# ARGUMENT PARSER SETUP FUNCTIONS

class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter):
    pass


def get_basic_parser():
    parser_dict = init_parser()

    add_infile_arg(parser_dict['required'], required=True)
    add_outfile_arg(parser_dict['required'], required=True, 
                    must_not_exist=True)
    add_fasta_arg(parser_dict['required'], required=True)

    add_refver_arg(parser_dict['optional'], required=False, choices='all')
    add_outfmt_arg(parser_dict['optional'], required=False)
    add_scheduler_args(parser_dict['optional'], default_parallel=1,
                       default_sched='slurm',
                       default_check_interval=DEFAULT_INTV_CHECK,
                       default_submit_interval=DEFAULT_INTV_SUBMIT)

    return parser_dict['main']


def init_parser(description=None):
    parser = argparse.ArgumentParser(description=description, 
                                     formatter_class=CustomFormatter,
                                     add_help=False)
    required = parser.add_argument_group(
        title='REQUIRED', description='Required ones. Accepts 1 argument.')
    optional = parser.add_argument_group(
        title='OPTIONAL', description='Optional ones. Accepts 1 argument.')
    flag = parser.add_argument_group(
        title='FLAG', description='Optional ones which accept 0 argument.')

    add_help_arg(flag)

    parser_dict = {'main': parser, 'required': required, 
                   'optional': optional, 'flag': flag}

    return parser_dict


def add_infile_arg(parser, required=True, help=f'Input vcf file path.'):
    parser.add_argument(
        '-i', '--infile', dest='infile_path', required=required, 
        type=arghandler_infile, metavar='<input file path>', help=help)


def add_infilelist_arg(
        parser, required=True, 
        help=f'One or more input file paths, separated with whitespaces.'):
    parser.add_argument(
        '--infilelist', dest='infile_path_list', required=required, 
        nargs='+', type=arghandler_infile, metavar='<input file path>', 
        help = help)


def add_indir_arg(parser, required=True, help=f'Input directory path.'):
    parser.add_argument(
        '--dir', dest='indir_path', required=required, type=arghandler_indir,
        metavar='<input directory path>', help=help)


def add_outfile_arg(
        parser, required=True, must_not_exist='ask',
        help='Output vcf file path. Must not exist in advance.'):
    handler = get_arghandler_outfile(must_not_exist)
    parser.add_argument(
        '-o', '--outfile', dest='outfile_path', required=required, 
        type=handler, metavar='<output file path>', help = help)


def add_outdir_arg(parser, required=True):
    parser.add_argument(
        '--outdir', dest='outdir_path', required=True,
        type=arghandler_outdir, metavar='<output directory>',
        help=('Output directory path. It will be created if it does not '
              'exist; otherwise it must be an empty directory.'))


def add_fasta_arg(parser, required=True):
    available_versions = ", ".join(common.DEFAULT_FASTA_PATH_DICT.keys())
    parser.add_argument(
        '-f', '--fasta', dest='fasta_path', required=required,
        type=arghandler_fasta, metavar='<fasta path>', 
        help=('A fasta file path or, alternatively, a reference genome '
              f'version(one of {available_versions}), in which case a preset '
              'fasta file for the reference version is used.'))


def add_refver_arg(parser, required=True, choices='all'):
    if choices == 'all':
        allowed_vals = tuple( common.CHR1_LENGTH_DICT.keys() )
    elif choices == 'mouse':
        allowed_vals = ('mm9', 'mm10', 'mm39')
    elif choices == 'human':
        allowed_vals = ('hg18', 'hg19', 'hg38')
    else:
        allowed_vals = choices

    helpmsg = f'Reference genome version. Must be one of {allowed_vals}.'
    parser.add_argument(
        '--refver', dest='refver', required=required, default=None, 
        choices=allowed_vals, metavar='<reference genome version>', 
        help = helpmsg)


def add_outfmt_arg(parser, required=False, default='z'):
    assert default in ('v', 'z', 'u', 'b')
    parser.add_argument(
        '-O', '--output-format', dest='mode_pysam', required=required, 
        default=default, choices=('v', 'z', 'u', 'b'), 
        type=(lambda x: common.PYSAM_FORMAT_DICT[x]),
        metavar='<output format>',
        help = ('Output vcf file format (bcftools style). Must be one of: '
                'v (uncompressed VCF), z (compressed VCF), '
                'u (uncompressed BCF), b (compressed BCF).'))


def add_parallel_arg(parser, required=False, default=1):
    parser.add_argument(
        '-p', '--parallel', dest='parallel', required=required,
        default=default, type=int,
        metavar='<number of parallelization>', 
        help=f'Number of parallelization.')


def add_sched_arg(parser, required=False, default='slurm'):
    assert default in ('slurm', 'local')
    parser.add_argument(
        '--sched', dest='sched', required=required, default=default, 
        choices=('slurm', 'local'), metavar='<"slurm" or "local">', 
        help = f'Parallelization method. Must be "slurm" or "local".')


def add_checkint_arg(parser, required=False, default=DEFAULT_INTV_CHECK):
    parser.add_argument(
        '--slurm-check-interval', dest='intv_check', required=required,
        default=default, type=int, metavar='<slurm check interval>', 
        help=f'Slurm job status check interval in seconds.')


def add_submitint_arg(parser, required=False, default=DEFAULT_INTV_SUBMIT):
    parser.add_argument(
        '--slurm-submit-interval', dest='intv_submit', required=required,
        default=default, type=int, metavar='<slurm submit interval>', 
        help=f'Slurm job submission interval in seconds.')


def add_help_arg(parser):
    parser.add_argument('-h', '--help', action='help',
                        help='show this help message and exit')


# functions which receive 'parser_dict' argument


def add_logging_args(parser_dict):
    parser_dict['optional'].add_argument(
        '--log', required=False, 
        help=('If used, progress messages will be written to this file. '
              'Existing file will be truncated.'))
    parser_dict['flag'].add_argument(
        '--silent', action='store_true',
        help='If set, progress messages will not be printed to the terminal.')


def add_scheduler_args(parser_dict, default_parallel=1, default_sched='slurm',
                       default_check_interval=DEFAULT_INTV_CHECK,
                       default_submit_interval=DEFAULT_INTV_SUBMIT):
    add_parallel_arg(parser_dict['optional'], required=False, 
                     default=default_parallel)
    add_sched_arg(parser_dict['optional'], required=False, 
                  default=default_sched)
    add_checkint_arg(parser_dict['optional'], required=False, 
                     default=default_check_interval)
    add_submitint_arg(parser_dict['optional'], required=False, 
                      default=default_submit_interval)


def add_rmtmp_arg(parser_dict):
    parser_dict['flag'].add_argument(
        '--dont-remove-tmp', dest='dont_rm_tmp', action='store_true',
        help=('If set, temporary files are not removed (which are removed '
              'by default).'))


############################################


def arghandler_infile(arg):
    arg = os.path.abspath(arg)

    try:
        common.check_infile_validity(arg)
    except Exception as e:
        raise argparse.ArgumentTypeError(str(e))
    else:
        return arg


def arghandler_infile_list(arg):
    """Args:
        arg: A list of input file paths 
    """

    return [arghandler_infile(x) for x in arg]


def arghandler_indir(arg):
    arg = os.path.abspath(arg)

    try:
        common.check_indir_validity(arg)
    except Exception as e:
        raise argparse.ArgumentTypeError(str(e))
    else:
        return arg


def arghandler_outfile_ask(arg):
    arg = os.path.abspath(arg)

    try:
        exists = common.check_outfile_validity(arg, must_not_exist=False)
    except Exception as e:
        raise argparse.ArgumentTypeError(str(e))
    else:
        if exists:
            msg = (f'Specified output file "{arg}" already exists. Would '
                       f'you like to proceed with overwriting? (y/n) ')
            ans = input(msg)
            while True:
                if ans == 'y':
                    return arg
                elif ans == 'n':
                    sys.exit(1)
                else:
                    ans = input('Please enter "y" or "n". ')
                    continue
        else:
            return arg


def arghandler_outfile_mayexist(arg):
    arg = os.path.abspath(arg)

    try:
        common.check_outfile_validity(arg, must_not_exist=False)
    except Exception as e:
        raise argparse.ArgumentTypeError(str(e))
    else:
        return arg


def arghandler_outfile_mustbeabsent(arg):
    arg = os.path.abspath(arg)

    try:
        common.check_outfile_validity(arg, must_not_exist=True)
    except Exception as e:
        raise argparse.ArgumentTypeError(str(e))
    else:
        return arg


def get_arghandler_outfile(must_not_exist):
    e_msg = '"must_not_exist" argument must be on of "yes", "no", or "ask".'
    assert must_not_exist in ('yes', 'no', 'ask'), e_msg
    if must_not_exist == 'yes':
        return arghandler_outfile_mustbeabsent
    elif must_not_exist == 'no':
        return arghandler_outfile_mayexist
    elif must_not_exist == 'ask':
        return arghandler_outfile_ask


def arghandler_outdir(arg):
    arg = os.path.abspath(arg)

    try:
        exists = common.check_outdir_validity(arg)
    except Exception as e:
        raise argparse.ArgumentTypeError(str(e))
    else:
        if not exists:
            os.mkdir(arg)
        return arg


def arghandler_fasta(arg):
    if arg in common.DEFAULT_FASTA_PATH_DICT.keys():
        return common.DEFAULT_FASTA_PATH_DICT[arg]
    else:
        return arghandler_infile(arg)


##########################################################




"""
SLURM JOB HANDLING
"""


class Job:
    """
    Designed for slurm

    Attributes:
        jobscript_path: Path of job script file to run.
        jobscript_string: Multi-line string to be written to stdin of sbatch.
        jobid: JobId as integer. None before job submission.
        submitted: True if job has been submitted, False otherwise.
        exitcode
        success: 
            None: 
                1) before job has been finished
                2) job has been finished but its state cannot be checked 
                    because job id has disappeared from slurm
            True if the job has finished successfully
            False if the job has finished not successfully
        no_jobid: 
            None: before job has been submitted
            True: job id is no longer remaining on slurm
            False: job id exists on slurm
        scontrol_result: Dictionary which contains result of 
            "scontrol show job <jobid>" command.
        sbatch_err_info: An informative string created when sbatch returns 
            a nonzero exit code.
        stdout_path: job stdout path
        stderr_path: job stderr path
        status
        logger
    """

    jobstates_pending = ('CONFIGURING', 'PENDING')
    jobstates_running = ('COMPLETING', 'RUNNING', 'RESV_DEL_HOLD', 
                         'REQUEUE_FED', 'REQUEUE_HOLD', 'RESIZING', 
                         'SIGNALING', 'SPECIAL_EXIT', 'STAGE_OUT', 'STOPPED', 
                         'SUSPENDED')
    jobstates_finished = ('CANCELLED', 'COMPLETED', 'BOOT_FAIL', 'DEADLINE', 
                          'FAILED', 'NODE_FAIL', 'OUT_OF_MEMORY', 'PREEMPTED', 
                          'TIMEOUT')
    jobstates_unknown = ('REVOKED',)

    def __init__(self, jobscript_path=None, jobscript_string=None, 
                 verbose=True, logger=None):
        common.check_num_None(1, (jobscript_path, jobscript_string), 
                              ('jobscript_path', 'jobscript_string'))

        self.jobscript_path = jobscript_path
        self.jobscript_string = jobscript_string
        self.jobid = None
        self.submitted = False
        self.exitcode = None
        self.success = None
        self.no_jobid = None
        self.scontrol_result = None
        self.sbatch_err_info = None
        self.stdout_path = None
        self.stderr_path = None

        self.status = {'pending': False, 'running': False, 'finished': False}

        if logger is None:
            self.logger = get_logger(
                formatter=DEFAULT_LOG_FORMATTERS['without_name'], 
                stderr=verbose)
        else:
            self.logger = logger

    def submit(self):
        if self.jobscript_path is None:
            self._submit_string()
        else:
            self._submit_path()

    def check_pending(self):
        if self.submitted:
            self.update()
            return self.status['pending']
        else:
            return False

    def check_running(self):
        if self.submitted:
            self.update()
            return self.status['running']
        else:
            return False

    def check_finished(self):
        if self.submitted:
            self.update()
            return self.status['finished']
        else:
            return False

    def update(self):
        assert self.submitted

        scontrol_result, no_jobid = get_scontrol_job_result(self.jobid)[:2]
        self.no_jobid = no_jobid

        if self.no_jobid:
            # self.success and self.scontrol_result are not touched in this case
            self.status = {'pending': False, 'running': False, 
                           'finished': True}
        else:
            self.scontrol_result = scontrol_result

            # update self.stderr/stdout path
            self.stdout_path = scontrol_result['StdOut']
            self.stderr_path = scontrol_result['StdErr']

            # update self.status and self.success
            jobstate = self.scontrol_result['JobState']
            if jobstate in Job.jobstates_pending:
                self.status = {
                    'pending': True, 'running': False, 'finished': False}
            elif jobstate in Job.jobstates_running:
                self.status = {
                    'pending': False, 'running': True, 'finished': False}
            elif jobstate in Job.jobstates_finished:
                self.status = {
                    'pending': False, 'running': False, 'finished': True}
            else:
                raise Exception(textwrap.dedent(f"""\
                    Unable to interpret JobState.
                    JobState: {jobstate}
                    scontrol result: {self.scontrol_result}"""))

            # update self.exitcode, self.success
            if self.status['finished']:
                self.success = (jobstate == 'COMPLETED')
                self.exitcode = int(scontrol_result['ExitCode'].split(':')[0])

    ##########################

    def _submit_string(self):
        p = subprocess.run([SBATCH], input=self.jobscript_string,
                           capture_output=True, text=True)

        if p.returncode != 0:
            self.sbatch_err_info = textwrap.dedent(f"""\
                Job submission using sbatch failed.

                string written to the stdin of sbatch: {self.jobscript_string}
                sbatch stdout: {p.stdout}
                sbatch stderr: {p.stderr}
                sbatch exit code: {p.returncode}""")
            raise Exception(self.sbatch_err_info)
        else:
            self.jobid = int(p.stdout.split()[-1])
            self.submitted = True
            self.logger.info(f'Submitted a job: JobID - {self.jobid}')

    def _submit_path(self):
        p = subprocess.run([SBATCH, self.jobscript_path],
                           capture_output=True, text=True)

        if p.returncode != 0:
            self.sbatch_err_info = textwrap.dedent(f"""\
                Job submission using sbatch failed.

                job script path: {self.jobscript_path}
                sbatch stdout: {p.stdout}
                sbatch stderr: {p.stderr}
                sbatch exit code: {p.returncode}""")
            raise Exception(self.sbatch_err_info)
        else:
            self.jobid = int(p.stdout.split()[-1])
            self.submitted = True
            self.logger.info(f'Submitted a job: JobID - {self.jobid}')


class JobList(list):
    """
    Attributes:
        logger
        success
        jobs_notsubmit
        jobs_pending
        jobs_running
        jobs_finished
        jobs_success
        jobs_failure
        jobs_unknown
        
    """

    def __init__(self, jobscript_path_list, verbose=True, logpath=None, 
                 logger=None):
        for jobscript_path in jobscript_path_list:
            self.append(Job(jobscript_path=jobscript_path, verbose=verbose, 
                            logger=logger))

        if logger is None:
            self.logger = get_logger(
                formatter=DEFAULT_LOG_FORMATTERS['without_name'], 
                stderr=verbose, filename=logpath)
        else:
            self.logger = logger

        self.success = None

    def submit(self, intv_submit=DEFAULT_INTV_SUBMIT):
        msg = 'One or more jobs have been already submitted.'
        assert not any(job.submitted for job in self), msg

        for job in self:
            job.submit()
            time.sleep(intv_submit)

    def wait(self, intv_check=DEFAULT_INTV_CHECK, edit_log_suffix=True):
        def log_status():
            pending = ", ".join(str(job.jobid) for job in self.jobs_pending)
            running = ", ".join(str(job.jobid) for job in self.jobs_running)
            finished = ", ".join(str(job.jobid) for job in self.jobs_finished)
            msg = textwrap.dedent(f"""\
                Current job status:
                  Not submitted yet: {len(self.jobs_notsubmit)}
                  Pending: {len(self.jobs_pending)} (JobIDs: {pending})
                  Running: {len(self.jobs_running)} (JobIDs: {running})
                  Finished: {len(self.jobs_finished)} (JobIDs: {finished})""")
            self.logger.info(msg)

        def log_epilogue():
            msg = f'''
    All finished.
    Successful jobs: {len(self.jobs_success)} (JobIDs: {", ".join(str(job.jobid) for job in self.jobs_success)})
    Failed jobs: {len(self.jobs_failure)} (JobIDs: {", ".join(str(job.jobid) for job in self.jobs_failure)})
    Jobs with unknown exit statuses: {len(self.jobs_unknown)} (JobIDs: {", ".join(str(job.jobid) for job in self.jobs_unknown)})
'''
            self.logger.info(msg)

        def main():
            while True:
                self.update()
                self.set_sublists_statuses()
                log_status()
                if all( job.status['finished'] for job in self ):
                    break
                else:
                    time.sleep(intv_check)
                    continue

            self.set_sublists_exitcodes()
            log_epilogue()
            self.success = all(job.success for job in self)
            if edit_log_suffix:
                self._edit_log_suffix()

        main()

    def get_exitcodes(self):
        return [ job.exitcode for job in self ]

    def update(self):
        for job in self:
            job.update()

    def set_sublists_statuses(self):
        self.jobs_notsubmit = [ job for job in self if not job.submitted ]
        self.jobs_pending = [ job for job in self if job.status['pending'] ]
        self.jobs_running = [ job for job in self if job.status['running'] ]
        self.jobs_finished = [ job for job in self if job.status['finished'] ]

    def set_sublists_exitcodes(self):
        self.jobs_success = [ job for job in self if job.success is True ]
        self.jobs_failure = [ job for job in self if job.success is False ]
        self.jobs_unknown = [ job for job in self if job.success is None ]

    def _edit_log_suffix(self):
        for job in self:
            if job.stderr_path is None:
                continue

            if job.exitcode == 0:
                new_logpath = re.sub('\.log$', '.success', job.stderr_path)
            else:
                new_logpath = re.sub('\.log$', '.failure', job.stderr_path)

            os.rename(job.stderr_path, new_logpath)


def get_scontrol_job_result(jobid):
    """
    Returns:
        A tuple (scontrol_result, no_jobid, returncode, stderr)

        scontrol_result: A dict with keys: 'JobId', 'JobName', 'UserId', 
            'GroupId', 'MCS_label', 'Priority', 'Nice', 'Account', 'QOS', 
            'JobState', 'Reason', 'Dependency', 'Requeue', 'Restarts', 
            'BatchFlag', 'Reboot', 'ExitCode', 'RunTime', 'TimeLimit', 
            'TimeMin', 'SubmitTime', 'EligibleTime', 'AccrueTime', 
            'StartTime', 'EndTime', 'Deadline', 'SuspendTime', 
            'SecsPreSuspend', 'LastSchedEval', 'Partition', 'AllocNode:Sid', 
            'ReqNodeList', 'ExcNodeList', 'NodeList', 'BatchHost', 'NumNodes', 
            'NumCPUs', 'NumTasks', 'CPUs/Task', 'ReqB:S:C:T', 'TRES', 
            'Socks/Node', 'NtasksPerN:B:S:C', 'CoreSpec', 'MinCPUsNode', 
            'MinMemoryNode', 'MinTmpDiskNode', 'Features', 'DelayBoot', 
            'OverSubscribe', 'Contiguous', 'Licenses', 'Network', 'Command', 
            'WorkDir', 'StdErr', 'StdIn', 'StdOut', 'Power', 'NtasksPerTRES:0'
        no_jobid: True or False 
        returncode
        stderr
    """

    cmdargs = [ SCONTROL, '-o', 'show', 'job', str(jobid) ]
    p = subprocess.run(args = cmdargs, text = True, capture_output = True)

    returncode = p.returncode
    stderr = p.stderr

    if p.returncode != 0:
        if p.stderr.strip() == 'slurm_load_jobs error: Invalid job id specified':
            scontrol_result = None
            no_jobid = True
        else:
            e_msg = f'''\
scontrol show job finished with nonzero exit code.
commands: {cmdargs}
stdout: {p.stdout}
stderr: {p.stderr}
exit code: {p.returncode}'''
            raise Exception(e_msg)

    else:
        no_jobid = False

        scontrol_result = dict()
        for word in p.stdout.split():
            wordsp = word.split('=', maxsplit = 1)
            if len(wordsp) == 1:
                key = wordsp[0]
                val = None
            elif len(wordsp) == 2:
                key, val = wordsp
                """
            else:
                e_msg = f'''\
Field with more than 2 "=" character found from "scontrol show job" output.
Field: {word}
scontrol command: {cmdargs}'''
                raise Exception(e_msg)
                """
    
            scontrol_result[key] = val
    
    return scontrol_result, no_jobid, returncode, stderr


###############################


def run_jobs(jobscript_paths, sched, intv_check, intv_submit, logger, 
             log_dir, raise_on_failure=True):
    assert sched in ('local', 'slurm')

    if sched == 'local':
        success, exitcode_list = run_jobs_local(jobscript_paths)
    elif sched == 'slurm':
        joblist = JobList(jobscript_paths, logger=logger)
        joblist.submit(intv_submit=intv_submit)
        joblist.wait(intv_check=intv_check)

        success = joblist.success
        exitcode_list = joblist.get_exitcodes()

    if raise_on_failure:
        if not success:
            raise SystemExit(
                (f'One or more jobs have finished unsuccessfully. Refer to '
                 f'log files in {log_dir}'))
        else:
            return None
    else:
        return success, exitcode_list


def run_jobs_local(jobscript_path_list):
    plist = list()
    for jobscript_path in jobscript_path_list:
        p = subprocess.Popen([jobscript_path], stdout=subprocess.PIPE, 
                             stderr=subprocess.PIPE, text=True)
        plist.append(p)

    for p in plist:
        p.wait()

    exitcode_list = [p.returncode for p in plist]
    success = all(x == 0 for x in exitcode_list)

    return success, exitcode_list


"""
JOB SCRIPT GENERATION
"""

def make_multiline_command(lines, leading_taps = 0):
    new_lines = list()
    new_lines.append( '\t'*leading_taps + lines[0] + ' \\' )
    new_lines.extend( [ '\t'*(leading_taps+1) + x + ' \\' for x in lines[1:-1] ] )
    new_lines.append( '\t'*(leading_taps+1) + lines[-1] )

    return '\n'.join(new_lines)


def make_jobscript_string(lines, shell = False, python = False, **kwargs):
    string_list = list()


    string_list.append('#!' + common.BASH)

    if 'N' not in kwargs:
        kwargs['N'] = 1
    if 'n' not in kwargs:
        kwargs['n'] = 1
    if 'o' not in kwargs:
        kwargs['o'] = '/dev/null'
    if 'c' not in kwargs:
        kwargs['c'] = 1

    for k,v in kwargs.items():
        if v is not None:
            string_list.append(f'#SBATCH -{k} {v}')

    string_list.append('')

    for line in lines:
        string_list.append(line)

    return '\n'.join(string_list)


def make_jobscript(jobscript_path, lines, **kwargs):
    common.check_outfile_validity(jobscript_path)
    jobscript_string = make_jobscript_string(lines, **kwargs)
    with open(jobscript_path, 'w') as outfile:
        outfile.write(jobscript_string)

