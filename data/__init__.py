import biox

"""
biox.data module for managing FASTQ and FASTA data files.
"""

# modules
from TabReader import *
from Fastq import *
from Fasta import *
from Bedgraph import *
from Gff3 import *
from Gtf import *
from Gene import *
from GeneFeature import *
from Bam import *

def fastq_qminmax(filename):
    """
    Returns min and max quality ord value from fastq file.
    """
    qmin = ()
    qmax = 0
    f = biox.data.Fastq(filename)
    while f.read():
        qmin = min(qmin, ord(f.quality[-1]))
        qmax = max(qmax, ord(f.quality[0]))        
    return qmin, qmax
    
def fasta_check(filename, allowed_chars=["A", "C", "T", "G", "N"]):
    """
    Checks if FASTA file format is valid.
    """
    f = biox.data.Fasta(filename)
    valid_fasta = True
    while f.read():
        seq_len = 0
        if len(f.sequence)==0:
            return False, "Sequence with ID %s has length 0" % (f.id)
        for allowed in allowed_chars:
            seq_len += f.sequence.upper().count(allowed)
        if seq_len!=len(f.sequence):
            return False, "Sequence with characters other than %s" % str(allowed_chars)
    return True, "File in FASTA format"
