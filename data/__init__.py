import biox

"""
biox.data module for managing FASTQ and FASTA data files.
"""

# modules
from TabReader import *
from Fastq import *
from Fasta import *

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