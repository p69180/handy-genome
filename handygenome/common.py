import sys
import os
import re
import time
import datetime
import tempfile
import collections
import json
import pprint
import shutil
import urllib.request
import urllib.parse
import urllib.error
import io
import contextlib
import inspect
import functools

import pysam


"""
Reference genome version aliases

MGSCv37 == mm9
GRCm38 == mm10
GRCm39 == mm39
NCBI36 == hg18
GRCh37 == hg19
GRCh38 == hg38
"""

"""
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
"""

PACKAGE_LOCATION = '/home/users/pjh/scripts/python_genome_packages'
DEFAULT_VCFVER = '4.3'

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

# SV symbolic allele strings
SV_ALTS = ('DEL', 'INS', 'DUP', 'INV', 'CNV', 'BND')
CPGMET_ALT = 'CPGMET'

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


# colors
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


def cpformat(obj):
    result = pprint.pformat(obj)
    result = re.sub('(True)', COLORS['green'] + '\\1' + COLORS['end'], result)
    result = re.sub('(False)', COLORS['red'] + '\\1' + COLORS['end'], result)
    result = re.sub('(None)', COLORS['purple'] + '\\1' + COLORS['end'], result)
    return result


def cpprint(obj):
    print(cpformat(obj))


# timer decorator

def deco_timer(func):
    """Print the runtime of the decorated function"""

    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        start_time = time.perf_counter()    # 1
        value = func(*args, **kwargs)
        end_time = time.perf_counter()      # 2
        run_time = end_time - start_time
        print(run_time)

        return value

    return wrapper_timer


###################################################

# argument sanity check decorators

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


def get_deco_num_set(names, n):
    def decorator(func):
        sig = inspect.signature(func)
        if not set(names).issubset(sig.parameters.keys()):
            raise Exception(
                f'The names of parameters given to '
                f'"get_deco_num_set" function '
                f'is not included in the parameter names of '
                f'the function "{func.__name__}".')

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            n_set = sum((name in ba.arguments) for name in names)
            if n_set != n:
                raise ValueError(
                    f'For the function "{func.__name__}", the '
                    #f'number of parameters set from arguments '
                    f'number of parameters being set, '
                    f'among {tuple(names)}, must be {n}.')

            return func(*args, **kwargs)

        return wrapper

    return decorator


def get_deco_num_set_differently(names, n):
    def decorator(func):
        sig = inspect.signature(func)
        if not set(names).issubset(sig.parameters.keys()):
            raise Exception(
                f'The names of parameters given to '
                f'"get_deco_num_set_differently" function '
                f'is not included in the parameter names of '
                f'the function "{func.__name__}".')

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            ba.apply_defaults()
            n_diff = sum((sig.parameters[name].default != ba.arguments[name])
                        for name in names)
            if n_diff != n:
                raise ValueError(
                    f'For the function "{func.__name__}", the '
                    f'number of parameters, among {tuple(names)}, '
                    f'being set as a value different from the default, '
                    f'must be {n}.')

            return func(*args, **kwargs)

        return wrapper

    return decorator


def get_deco_arg_choices(mapping):
    """Args:
        mapping: {'argname': (valid_value1, valid_value2, ...), ...}
    """

    def decorator(func):
        sig = inspect.signature(func)
        if not set(mapping.keys()).issubset(sig.parameters.keys()):
            raise Exception(
                f'The names of parameters given to '
                f'"get_deco_check_arg_choices" function '
                f'is not included in the parameter names of '
                f'the function "{func.__name__}".')

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            ba.apply_defaults()
            for key, val in mapping.items():
                if ba.arguments[key] not in val:
                    raise ValueError(
                        f'For the function "{func.__name__}", '
                        f'the parameter "{key}" must be one of these values: '
                        f'{tuple(val)}.')

            return func(*args, **kwargs)

        return wrapper

    return decorator


###################################################

class ChromDict(collections.OrderedDict):
    @get_deco_num_set(('fasta_path', 'fasta', 'bam_path', 'bam', 
                       'vcfheader', 'bamheader', 'custom', 'refver'), 1)
    def __init__(self, fasta_path=None, fasta=None, bam_path=None, bam=None, 
                 vcfheader=None, bamheader=None, custom=None, refver=None):
        """
        Args:
            fasta: pysam.FastaFile object
            bam: pysam.AlignmentFile object
            vcfheader: pysam.VariantHeader object
            bamheader: pysam.AlignmentHeader object
            custom: {'contigs': ['contig1', 'contig2', ...], 
                     'lengths': [length1, length2, ...] }
        """

        # sanity check
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
            elif bamheader is not None:
                wrapper = bamheader
            elif refver is not None:
                wrapper = pysam.FastaFile(DEFAULT_FASTA_PATH_DICT[refver])
    
            for chrom, length in zip(wrapper.references, wrapper.lengths):
                self[chrom] = length
    
            if (
                    (fasta_path is not None) or 
                    (bam_path is not None) or 
                    (refver is not None)):
                wrapper.close()

        # set contigs, lengths
        self.contigs = list(self.keys())
        self.lengths = list(self.values())


class Vcfspec:
    def __init__(self, chrom=None, pos=None, ref=None, alts=None):
        if alts is not None:
            if not isinstance(alts, (tuple, list)):
                raise Exception(f'"alts" argument must be a tuple or a list.')

        self.chrom = chrom
        self.pos = pos
        self.pos0 = pos - 1
        self.ref = ref
        self.alts = tuple(alts)

    def __repr__(self):
        if len(self.alts) == 1:
            altstring = str(self.alts[0])
        else:
            altstring = str(list(self.alts))
        return (f'<Vcfspec ({self.chrom}:{self.pos} '
                f'{self.ref}>{altstring})>')

    def get_id(self):
        return '_'.join([self.chrom, 
                         str(self.pos), 
                         self.ref, 
                         '|'.join(self.alt)])

    def get_mutation_type(self, idx=0):
        return get_mttype(self.ref, self.alts[idx])

    def get_mttype_firstalt(self):
        return self.get_mutation_type(0)

    def get_tuple(self):
        return (self.chrom, self.pos, self.ref, self.alts)

    def __hash__(self):
        return hash(self.get_tuple())

    ### ranges
    def get_REF_range0(self):
        return range(self.pos0, self.pos0 + len(self.ref))

    @functools.cache
    def get_preflank_range0(self, idx=0, flanklen=1):
        assert flanklen >= 1, f'"flanklen" argument must be at least 1.'

        if self.alts[idx][0] == self.ref[0]:
            flanklen = flanklen - 1
        return range(self.pos0 - flanklen, self.pos0)

    @functools.cache
    def get_postflank_range0(self, flanklen=1):
        assert flanklen >= 1, f'"flanklen" argument must be at least 1.'

        return range(self.pos0 + len(self.ref),
                     self.pos0 + len(self.ref) + flanklen)


def check_vcfspec_monoallele(vcfspec):
    if len(vcfspec.alts) != 1:
        raise Exception('The input vcfspec must be monoalleleic.')


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

    def __init__(self, chrom, start1=None, end1=None, start0=None, end0=None):
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

    return (chromdict.contigs.index(chrom), pos)


def get_read_sortkey(chromdict):
    def sortkey(read):
        return coord_sortkey(read.reference_name, read.reference_start, 
                             chromdict)

    return sortkey


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


def compare_coords(chrom1, pos1, chrom2, pos2, chromdict):
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


def get_linesp(line, sep='\t'):
    return rm_newline(line).split(sep)


###################################################


def get_mttype(ref, alt):
    if RE_PATS['nucleobases'].fullmatch(alt) is None:
        if any(
                (re.fullmatch(f'<{x}(:.+)?>', alt) is not None)
                for x in SV_ALTS):
            mttype = 'sv'
        elif (
                (RE_PATS['alt_bndstring_1'].fullmatch(alt) is not None) or 
                (RE_PATS['alt_bndstring_2'].fullmatch(alt) is not None)):
            mttype = 'sv'
        elif alt == f'<{CPGMET_ALT}>':
            mttype = 'cpgmet'
        else:
            raise Exception(f'Unexpected symbolic ALT allele: {alt}')
    else:
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
                mttype = 'delins'

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

def listdir(path):
    return sorted(os.path.join(path, x) for x in os.listdir(path))


def get_padded_indices(n):
    """Begins with 0"""

    width = len(str(n-1))
    result = [str(idx).zfill(width) for idx in range(n)]

    return result


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

@get_deco_num_set(('chromdict', 'vcfheader', 'bamheader'), 1)
def infer_refver(chromdict=None, vcfheader=None, bamheader=None):
    if chromdict is not None:
        return infer_refver_chromdict(chromdict)
    elif vcfheader is not None:
        return infer_refver_vcfheader(vcfheader)
    elif bamheader is not None:
        return infer_refver_bamheader(bamheader)


def infer_refver_chromdict(chromdict):
    if '1' in chromdict.contigs:
        chr1_length = chromdict.lengths[chromdict.contigs.index('1')]
    elif 'chr1' in chromdict.contigs:
        chr1_length = chromdict.lengths[chromdict.contigs.index('chr1')]
    else:
        raise Exception('"1" and "chr1" both not included in the chromosome '
                        'name list.')
    
    if chr1_length in CHR1_LENGTH_DICT_REV:
        refver = CHR1_LENGTH_DICT_REV[chr1_length]
    else:
        raise Exception(f'Cannot infer refver: unknown chr1 length')
        #refver = None # unknown reference genome
    
    return refver


def infer_refver_vcfheader(vcfheader):
    return infer_refver_chromdict(ChromDict(vcfheader=vcfheader))


def infer_refver_bamheader(bamheader):
    return infer_refver_chromdict(ChromDict(bamheader=bamheader))


def infer_refver_vr(vr):
    return infer_refver_chromdict(ChromDict(vcfheader=vr.header))


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


###################################################

# deprecated; this does not work
def get_vcf_noerr(*args, **kwargs):
    with contextlib.redirect_stderr(io.StringIO()) as err, \
            contextlib.redirect_stdout(io.StringIO()) as out:
        vcf = pysam.VariantFile(*args, **kwargs)

    for buf in (err, out):
        msg = buf.getvalue()
        if not msg.startswith('[E::idx_find_and_load] '
                              'Could not retrieve index file for'):
            print(msg, end='', file=sys.stderr)

    return vcf
        
