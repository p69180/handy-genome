import os
import tempfile
import subprocess
import gzip
import collections
import re

import Bio.SeqIO
import Bio.SeqRecord
import Bio.Seq
import pysam

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))


def run_bwa(seq, refver):
    ref_path = common.DEFAULT_FASTA_PATHS[refver]
    seqrec = Bio.SeqRecord.SeqRecord(seq=Bio.Seq.Seq(seq), id='query')
    with tempfile.TemporaryDirectory(dir=os.getcwd()) as tmpdir:
        query_path = os.path.join(tmpdir, 'input.fasta')
        sam_path = os.path.join(tmpdir, 'output.sam')
        Bio.SeqIO.write([seqrec], query_path, 'fasta')
        p = subprocess.run([common.BWA, 'mem', 
                            '-Y', 
                            '-M', 
                            '-t', '1',
                            '-o', sam_path,
                            ref_path,
                            query_path])
        readlist = list(pysam.AlignmentFile(sam_path).fetch())

    return readlist

