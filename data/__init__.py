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
    valid_format = True
    reason = "File is in FASTA format"
    part = "id"
    f = open(filename)
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "")
        if part=="sequence" and r.startswith(">"):
            part = "id"
        if part=="id" and not r.startswith(">"):
            valid_format = False
            reason = "No ID"
            break
        if part=="sequence" and r.startswith(">"):
            valid_format = False
            reason = "ID without sequence"
            break
        if part=="sequence":
            r = r.upper()
            seq_len = 0
            for allowed in allowed_chars:
                seq_len += r.count(allowed)
            if seq_len != len(r):
                valid_format = False
                reason = "Sequence contains other characters than %s" % str(allowed_chars)
                break
        if part=="id" and r.startswith(">"):
            part = "sequence"
        if part=="sequence" and r=="":
            valid_format = False
        r = f.readline()
    if part=="sequence" and r=="":
        valid_format = False
        reason = "ID without sequence"
    return valid_format, reason
