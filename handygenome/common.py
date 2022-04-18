import sys
import os
import re
import datetime
import tempfile
import argparse
import collections
import json
import pprint
import shutil
import urllib.request
import urllib.parse
import urllib.error

import pysam


'''
Reference genome version aliases

MGSCv37 == mm9
GRCm38 == mm10
GRCm39 == mm39
NCBI36 == hg18
GRCh37 == hg19
GRCh38 == hg38
'''

'''
pysam file mode string

- Valid mode string pattern : ^[rwa]([bzu]?[0-9]?|[0-9]?[bzu]?)$
- Conflicting characters: (b,z,u)

- 'wb' : Compressed BCF (level 6), regardless of file name. 
- 'wb[0-9]' : BCF with level of compression indicated by the number, regardless of file name. 

- 'wz' : Compressed VCF (level 6), regardless of file name.
- 'wz[0-9]' : VCF with level of compression indicated by the number, regardless of file name. 

- 'wu' : Uncompressed VCF, regardless of file name.
- 'wu[0-9]' : Uncompressed VCF, regardless of file name.

- 'w[0-9]' : Uncompressed VCF, regardless of file name. 

- 'w' :
    *.vcf : uncompressed VCF
    *.vcf.gz : compressed VCF (level 6)
    *.bcf : compressed BCF (level 6)
    *.bcf.gz : compressed VCF (level 6)
'''

PACKAGE_LOCATION = '/home/users/pjh/scripts/python_genome_packages'

THRESHOLD_TEMPLATE_LENGTH = 1000
DEFAULT_VCFVER = '4.3'
#DEFAULT_INTV_CHECK = 60
#DEFAULT_INTV_SUBMIT = 1

# re patterns
RE_PATS = {
    'int' : re.compile('-?[0-9]+'),
    'float' : re.compile('(-?[0-9]+\.[0-9]+)|(-?[0-9]+(\.[0-9]+)?e-?[0-9]+)'),
    'nucleobases' : re.compile('[ACGTNacgtn]+'),
    'numbered_chromosome' : re.compile('(chr)?[0-9]+'),
    'assembled_chromosome' : re.compile('(chr)?([0-9]+|X|Y)'),

    'alt_bndstring_1' : re.compile('^(?P<t>[^\[\]]+)(?P<bracket1>\[|\])(?P<matechrom>[^:]+):(?P<matepos>[0-9]+)(?P<bracket2>\[|\])$'),
    'alt_bndstring_2' : re.compile('^(?P<bracket1>\[|\])(?P<matechrom>[^:]+):(?P<matepos>[0-9]+)(?P<bracket2>\[|\])(?P<t>[^\[\]]+)$'),
        # these bndstring patterns assume that "t" portion must not be blank
}

# executable paths
BASH = '/usr/bin/bash'
BWA = '/usr/local/bin/bwa'
PERL = '/home/users/pjh/scripts/conda_wrapper/perl'
PYTHON = '/home/users/pjh/tools/miniconda/210821/miniconda3/envs/genome_v5/bin/python'

VEP_V102 = '/home/users/pjh/scripts/conda_wrapper/vep_v102' # for mm10
#VEP_V104 = '/home/users/pjh/scripts/conda_wrapper/vep_v104'
VEP_V105 = '/home/users/pjh/scripts/conda_wrapper/vep_v105'
VEP = VEP_V105
VEP_MM10 = VEP_V102

CONDABIN_PATH = '/home/users/pjh/conda_bin'
SAMTOOLS = os.path.join(CONDABIN_PATH, 'samtools')
BCFTOOLS = os.path.join(CONDABIN_PATH, 'bcftools')
BEDTOOLS = os.path.join(CONDABIN_PATH, 'bedtools')

# vep cache directory
VEP_CACHE_DIR = '/home/users/pjh/.vep'

# chr1 lengths
CHR1_LENGTH_DICT = {
    'mm9': 197_195_432,
    'mm10': 195_471_971,
    'mm39': 195_154_279,
    'hg18': 247_249_719,
    'hg19': 249_250_621,
    'hg38': 248_956_422,
    }
CHR1_LENGTH_DICT_REV = { val : key for key, val in CHR1_LENGTH_DICT.items() }

# default fasta paths
DEFAULT_FASTA_PATHS = {
    'hg18': '/home/users/pjh/References/reference_genome/NCBI36/ucsc/custom_files/hg18.fa',
    'hg19': '/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta',
    'hg19_hs37d5': '/home/users/sypark/02_Reference/15_pcawg/genome.fa',
    'hg38': '/home/users/data/01_reference/human_g1k_v38/Homo_sapiens_assembly38.fasta',

    'mm9': '/home/users/pjh/References/reference_genome/mm9/ucsc/custom_files/mm9.fa',
    'mm10': '/home/users/pjh/References/reference_genome/GRCm38/ucsc/custom_files/mm10.fa',
    'mm39': '/home/users/pjh/References/reference_genome/GRCm39/ucsc/custom_files/mm39.fa',
    }
DEFAULT_FASTA_PATH_DICT = DEFAULT_FASTA_PATHS

BCFTOOLS_FORMAT_DICT = { 'v':'w', 'z':'wz', 'u':'wbu', 'b':'wb' }
CYVCF2_FORMAT_DICT = BCFTOOLS_FORMAT_DICT
PYSAM_FORMAT_DICT = { 'v':'wz0', 'z':'wz', 'u':'wb0', 'b':'wb' }
PYSAM_MODE_DICT = PYSAM_FORMAT_DICT
DEFAULT_MODE_BCFTOOLS = 'z'

# http function constants
HTTP_HEADER_POST = {"Content-Type": "application/json", 
                    "Accept": "application/json"}
HTTP_HEADER_GET = {'Content-Type': 'application/json'}

# color escape sequences
COLORS = {
    'red':     '\033[38;5;196m',
    'magenta': '\033[38;5;201m',
    'pink':    '\033[38;5;213m',
    'orange':  '\033[38;5;9m',
    'yellow':  '\033[38;5;214m',
    'gold':    '\033[38;5;11m',
    'green':   '\033[38;5;40m',
    'blue':    '\033[38;5;33m',
    'cyan':    '\033[38;5;14m',
    'purple':  '\033[38;5;93m',
    'gray':    '\033[38;5;8m',
    'white':   '\033[38;5;15m',
    'end':     '\033[0m',
    }


class ChromDict(collections.OrderedDict):
    def __init__(self, fasta_path=None, fasta=None, bam_path=None, bam=None, 
                 vcfheader=None, custom=None, refver=None):
        """
        Args:
            fasta: pysam.FastaFile object
            bam: pysam.AlignmentFile object
            vcfheader: pysam.VariantHeader object
            custom: {'contigs': ['contig1', 'contig2', ...], 
                     'lengths': [length1, length2, ...] }
        """

        # sanity check
        check_num_notNone(
            1, 
            (fasta_path, fasta, bam_path, bam, vcfheader, custom, refver), 
            ('fasta_path', 'fasta', 'bam_path', 'bam', 'vcfheader', 
             'custom', 'refver'),
            )
        if refver is not None:
            check_arg_choices(refver, 'refver', DEFAULT_FASTA_PATH_DICT.keys())

        # set self dict
        if vcfheader is not None:
            for contig in vcfheader.contigs.values():
                self[contig.name] = contig.length
        elif custom is not None:
            for chrom, length in zip(custom['contigs'], custom['lengths']):
                self[chrom] = length
        else:
            if fasta_path is not None:
                wrapper = pysam.FastaFile(fasta_path)
            elif fasta is not None:
                wrapper = fasta
            elif bam_path is not None:
                wrapper = pysam.AlignmentFile(bam_path)
            elif bam is not None:
                wrapper = bam
            elif refver is not None:
                 wrapper = pysam.FastaFile(DEFAULT_FASTA_PATH_DICT[refver])
    
            for chrom, length in zip(wrapper.references, wrapper.lengths):
                self[chrom] = length
    
            if ((fasta_path is not None) or 
                (bam_path is not None) or 
                (refver is not None)):
                wrapper.close()

        # set contigs, lengths
        self.contigs = list(self.keys())
        self.lengths = list(self.values())

class Vcfspec(
    collections.namedtuple('Vcfspec', 
                           field_names = ('chrom', 'pos', 'ref', 'alt'), 
                           defaults = (None, None, None, None))):

    def get_id(self):
        return '_'.join([self.chrom, str(self.pos), self.ref, self.alt])


class Interval:
    """
    Attributes:
        chrom
        start0
        end0
        start1
        end1
        range0
        length
    """

    def __init__(self, chrom, start1 = None, end1 = None, start0 = None, end0 = None):
        """
        'chrom' is mandatory
        ('start1' and 'end1') or ('start0' and 'end0') must be given(for coordinate setting).

        'start0' and 'end0' are 0-based half-open system.
        'start1' and 'end1' are 1-based closed system.
        """

        self.chrom = chrom

        # start0, end0, start1, end1
        if (start1 is not None) and (end1 is not None):
            self.start0 = start1 - 1
            self.end0 = end1
        else:
            self.start0 = start0
            self.end0 = end0

        self.start1 = self.start0 + 1
        self.end1 = self.end0

        # range
        self.range0 = range(self.start0, self.end0)

        # length
        self.length = self.end0 - self.start0

    def __repr__(self):
        return f'<Interval> {self.chrom}:{self.start1:,}-{self.end1:,}'


class ColorsQQ:
    mine="\033[48;5;6m"
    busy="\033[48;5;244m"
    free="\033[48;5;238m"
    end="\033[0m"
    nor="\033[48;5;160m"
    nor2="\033[48;5;52m"
    danger1="\033[38;5;208m"
    danger2="\033[38;5;196m"
    darkgray="\033[38;5;240m"


def visualize_colors():
    for i in range(256):
        print(f'{i:<3d} \033[38;5;{i}m\\033[38;5;{i}m\033[0m')


###################################################


def cpformat(obj):
    result = pprint.pformat(obj)
    result = re.sub('(True)', COLORS['green'] + '\\1' + COLORS['end'], result)
    result = re.sub('(False)', COLORS['red'] + '\\1' + COLORS['end'], result)
    result = re.sub('(None)', COLORS['purple'] + '\\1' + COLORS['end'], result)
    return result


def cpprint(obj):
    print(cpformat(obj))


###################################################


def timer(func):
    """Print the runtime of the decorated function"""

    import functools
    import time

    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        start_time = time.perf_counter()    # 1
        value = func(*args, **kwargs)
        end_time = time.perf_counter()      # 2
        run_time = end_time - start_time    # 3
        print(f"Finished {func.__name__!r} in {run_time:.4f} secs")
        return value
    return wrapper_timer


###################################################

# argument sanity checks

def check_num_None(n, values, names):
    None_count = sum(x is None for x in values)
    if None_count != n:
        names = [f"'{x}'" for x in names]
        raise Exception(f'There must be exactly {n} None among {", ".join(names)}.')


def check_num_notNone(n, values, names):
    notNone_count = sum(x is not None for x in values)
    if notNone_count != n:
        names = [f"'{x}'" for x in names]
        raise Exception(f'Exactly {n} of {", ".join(names)} must be set as a non-None value.')


def check_arg_choices(arg, argname, choices):
    if arg not in choices:
        raise ValueError(f"'{argname}' must be one of {', '.join(choices)}")


###################################################


def zenumerate(iterable):
    length = len(tuple(iterable))
    width = len(str(length-1))
    for idx, item in enumerate(iterable):
        yield str(idx).zfill(width), item


def round_sig(num, n):
    """
    round 'num' with 'n' significant digits
    """

    assert n > 0
    from math import log10, floor
    return round(num, -floor(log10(abs(num))) + (n-1))


def str_to_nonstr(val):
    assert isinstance(val, str)

    if val == 'None':
        return None
    elif val == 'True':
        return True
    elif val == 'False':
        return False
    elif RE_PATS['int'].fullmatch(val) is not None:
        return int(val)
    elif RE_PATS['float'].fullmatch(val) is not None:
        return float(val)
    else:
        return val


def nonstr_to_str(val):
    assert not isinstance(val, str)

    if val is None:
        return None
    else:
        return str(val)


def get_datestring():
    """
    Returns a string like '2021-12-06 11:49:55'
    """

    return str( datetime.datetime.now() ).split('.')[0]


def get_timestamp():
    """
    Returns a string like 'KST 2021-12-06 11:51:36'
    """

    dt = datetime.datetime.now().astimezone()
    timestamp = f'{str(dt.tzinfo)} {str(dt).split(".")[0]}'

    return timestamp


def print_err(*args, stderr = True, files = None, **kwargs):
    """
    Args:
        stderr: (Bool) Whether to write to stderr
        files: A list of file paths to which message is written
    """

    if stderr:
        print(*args, file = sys.stderr, flush = True, **kwargs)

    if files is not None:
        for fname in files:
            with open(fname, 'a') as f:
                print(*args, file = f, flush = True, **kwargs)


def printerr(*args, **kwargs):
    print_err(*args, **kwargs)


def printerrdate(*args, **kwargs):
    datestring = get_datestring()
    print_err(f'[{datestring}]', *args, **kwargs)


def print_timestamp(*args, **kwargs):
    timestamp = get_timestamp()
    print_err(f'[{timestamp}]', *args, **kwargs)


def coord_sortkey(chrom, pos, chromdict): 
    """
    Args:
        pos: 1-based
        chromdict: ChromDict class instance
    """

    #return sum(chromdict.lengths[ : chromdict.contigs.index(chrom) ]) + pos
    return (chromdict.contigs.index(chrom), pos)


def get_vr_sortkey(chromdict):
    def sortkey(vr):
        return coord_sortkey(vr.chrom, vr.pos, chromdict)

    return sortkey


def get_vcfspec_sortkey(chromdict):
    def sortkey(vcfspec):
        return coord_sortkey(vcfspec.chrom, vcfspec.pos, chromdict)

    return sortkey


def get_vcfspec_order(vcfspec1, vcfspec2, chromdict):
    """
    Returns:
        0: equal
        negative integer: vcfspec1 comes first
        positive integer: vcfspec2 comes first
    """

    return (coord_sortkey(vcfspec1.chrom, vcfspec1.pos, chromdict) 
            - coord_sortkey(vcfspec2.chrom, vcfspec2.pos, chromdict))


def get_order(chrom1, pos1, chrom2, pos2, chromdict):
    """
    Returns:
        0: equal; -1: chrom1/pos1 comes first; 1: chrom2/pos2 comes first
    """

    if chrom1 == chrom2 and pos1 == pos2:
        return 0
    else:
        if (coord_sortkey(chrom1, pos1, chromdict) 
            < coord_sortkey(chrom2, pos2, chromdict)):
            return -1
        else:
            return 1

###################################################


def rm_newline(line):
    return re.sub('(\r)?\n$', '', line)


def get_linesp(line, sep = '\t'):
    return rm_newline(line).split(sep)


###################################################


def get_mttype(ref, alt):
    if RE_PATS['nucleobases'].match(alt) is not None:
        if len(ref) == len(alt):
            if len(ref) == 1:
                mttype = 'snv'
            else:
                mttype = 'mnv'
        else:
            if len(ref) == 1:
                mttype = 'ins'
            elif len(alt) == 1:
                mttype = 'del'
            else:
                mttype = 'cindel'
    else:
        mttype = 'sv'

    return mttype


def get_indelseq(ref, alt):
    mttype = get_mttype(ref, alt)
    if mttype == 'ins':
        indelseq = alt[1:]
    elif mttype == 'del':
        indelseq = ref[1:]
    else:
        indelseq = None
    
    return indelseq


###################################################

# filename-related functions


def listdir(path):
    return sorted(os.path.join(path, x) for x in os.listdir(path))


def get_tempfile_path(prefix = None, suffix = None, where = os.getcwd(), delete = False, is_dir = False):
    if delete:
        fd, path = tempfile.mkstemp(prefix = prefix, suffix = suffix, dir = where)
        os.close(fd)
        os.remove(path)
    else:
        if is_dir:
            path = tempfile.mkdtemp(prefix = prefix, suffix = suffix, dir = where)
        else:
            fd, path = tempfile.mkstemp(prefix = prefix, suffix = suffix, dir = where)
            os.close(fd)

    return path

get_tmpfile_path = get_tempfile_path


def get_tempdir_paths(subdirs, prefix = None, suffix = None, where = os.getcwd(), top_path = None):
    """
    Args:
        subdirs: An iterable which contains subdirectory basenames
    """

    assert 'top' not in subdirs

    if top_path is None:
        top_path = os.path.abspath(
                get_tempfile_path(prefix = prefix, suffix = suffix, where = where, delete = False, is_dir = True)
                )
    else:
        exists = check_outdir_validity(top_path)
        if not exists:
            os.mkdir(top_path)

    tmpdir_paths = dict()
    tmpdir_paths['top'] = top_path
    for basename in subdirs:
        tmpdir_paths[basename] = os.path.join(tmpdir_paths['top'], basename)
        os.mkdir(tmpdir_paths[basename])

    return tmpdir_paths

get_tmpdir_paths = get_tempdir_paths


def get_padded_indices(n):
    """Begins with 0"""

    width = len(str(n-1))
    result = [ str(idx).zfill(width) for idx in range(n) ]

    return result


def check_outdir_validity(outdir):
    """
    Return: 
        True if outdir already exists, otherwise False.

    Raise:
        If 1) outdir is an existing non-directory file, or 2) outdir is an existing non-empty directory
    """

    if not os.path.exists(outdir):
        if os.access(os.path.dirname(outdir), os.W_OK | os.X_OK):
            return False
        else:
            raise Exception(f"You do not have w/x permission on 'dirname' of '{outdir}'.")

    elif os.path.isfile(outdir):
        raise Exception(f"'{outdir}' is an existing regular file.")

    elif os.path.isdir(outdir):
        if os.access(outdir, os.W_OK | os.X_OK):
            if len(os.listdir(outdir)) != 0:
                raise Exception(f"'{outdir}' is an existing directory and is not empty.")
            else:
                return True
        else:
            raise Exception(f"'{outdir}' is an existing directory which you do not have w/x permission on.")

    else:
        raise Exception(f"'{outdir}' is an existing file which is neither a regular file nor a directory.")


def check_outfile_validity(outfile, must_not_exist = False):
    """
    Return: 
        True if outfile already exists, otherwise False.
    Raise:
        permission problem
    """

    if not os.path.exists(outfile):
        if os.access(os.path.dirname(outfile), os.W_OK | os.X_OK):
            return False
        else:
            raise Exception(f"You do not have w/x permission on dirname of '{outfile}'.")

    elif os.path.isfile(outfile):
        if must_not_exist:
            raise Exception(f"Specified file '{outfile}' must not exist in advance.")
        else:
            return True

    elif os.path.isdir(outfile):
        raise Exception(f"'{outfile}' is an existing directory.")

    else:
        raise Exception(f"'{outfile}' is an existing file which is neither a regular file nor a directory.")


def check_infile_validity(infile):
    """
    Raise:
        If infile does not exist or user does not have read permission
    """

    if os.path.exists(infile):
        if os.path.isfile(infile):
            if not os.access(infile, os.R_OK):
                raise Exception(f"You do not have read permission on '{infile}'.")
        else:
            raise Exception(f"'{infile}' is not a regular file.")
    else:
        raise Exception(f"'{infile}' does not exist.")


def check_indir_validity(indir):
    if os.path.exists(indir):
        if os.path.isdir(indir):
            if not os.access(indir, os.R_OK | os.X_OK):
                raise Exception(f"You do not have read and execution permission on '{indir}'.")
        else:
            raise Exception(f"'{indir}' is not a directory.")
    else:
        raise Exception(f"'{indir}' does not exist.")


###################################################


def printwidth_get_width_list(df):
    '''
    df: [
    [line1_field1, line1_field2, ... ],
    [line2_field1, line2_field2, ... ],
    ...,
    ]
    '''
    width_list = list()
    for i in range(len(df[0])):
        width_list.append(list())

    for line in df:
        for idx, field in enumerate(line):
            width_list[idx].append(len(field))

    for idx, e in enumerate(width_list):
        width_list[idx] = max(e)

    return width_list


def printwidth_print_line(line, width_list, margin, target):
    printresult = ''
    for idx, e in enumerate(line):
        printresult += f'{e:>{width_list[idx] + margin}s}'
    if target == 'out':
        print(printresult, flush = True)
    elif target == 'err':
        print(printresult, flush = True, file = sys.stderr)


def printwidth(df, margin = 2, target = 'out'):
    for line in df:
        for idx, e in enumerate(line):
            line[idx] = str(e)

    width_list = printwidth_get_width_list(df)
    for line in df:
        printwidth_print_line(line, width_list, margin, target)


###################################################


def infer_refver(chromdict):
    if '1' in chromdict.contigs:
        chr1_length = chromdict.lengths[chromdict.contigs.index('1')]
    elif 'chr1' in chromdict.contigs:
        chr1_length = chromdict.lengths[chromdict.contigs.index('chr1')]
    else:
        raise Exception('"1" and "chr1" both not included in the chromosome name list.')
    
    if chr1_length in CHR1_LENGTH_DICT_REV:
        refver = CHR1_LENGTH_DICT_REV[chr1_length]
    else:
        refver = None # unknown reference genome
    
    return refver


def infer_refver_vr(vr):
    return infer_refver(ChromDict(vcfheader=vr.header))


###################################################


def http_get(url, params = None, headers = dict(), text = False):
    if params is not None:
        url = url + '?' + urllib.parse.urlencode(params)

    req = urllib.request.Request(url, headers = headers, method = 'GET')
    try:
        with urllib.request.urlopen(req) as response:
            if text:
                result = str(response.read(), 'utf-8')
            else:
                result = json.loads(response.read())
    except urllib.error.HTTPError as e:
        print(str(e.read(), 'utf-8'))
        raise

    return result


def http_post(url, data, params = None, headers = dict(), text = False):
    data = json.dumps(data).encode('ascii')

    if params is not None:
        url = url + '?' + urllib.parse.urlencode(params)

    req = urllib.request.Request(url, data = data, headers = headers, method = 'POST')
    try:
        with urllib.request.urlopen(req) as response:
            if text:
                result = str(response.read(), 'utf-8')
            else:
                result = json.loads(response.read())
    except urllib.error.HTTPError as e:
        print(str(e.read(), 'utf-8'))
        raise

    return result


def download(url, path):
    with urllib.request.urlopen(url) as response:
        with open(path, 'wb') as outfile:
            shutil.copyfileobj(response, outfile)


###################################################


def write_mode_arghandler(mode_bcftools, mode_pysam):
    """
    Args:
        bcftools_mode: bcftools style mode string
        mode_pysam: pysam style mode string

    If 'mode_pysam' is set as a non-None value,
    it overrides 'mode_bcftools'.
    """

    if mode_pysam is None:
        return PYSAM_MODE_DICT[mode_bcftools]
    else:
        return mode_pysam


###################################################


def get_different_base(base):
    assert base in 'ACTG'
    if base == 'A':
        return 'C'
    else:
        return 'A'


