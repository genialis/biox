import biox

# modules
from TabReader import *
from Fastq import *
from Fasta import *

"""
Returns min and max quality ord value from fastq file.
"""
def fastq_qminmax(filename):
    qmin = ()
    qmax = 0
    f = biox.data.Fastq(filename)
    while f.read():
        qmin = min(qmin, ord(f.quality[-1]))
        qmax = max(qmax, ord(f.quality[0]))        
    return qmin, qmax