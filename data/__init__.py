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

def prepare_fasta_gbrowse(filename):
    f = biox.data.Fasta(filename)
    f_gff3 = open("chromosomes.gff3", "wt")
    while f.read():
        fout = open("%s.fa" % f.id, "wt")
        fout.write(">%s\n" % f.id)
        sequences = [f.sequence[i:i+50] for i in range(0, len(f.sequence), 50)]
        fout.write("\n".join(sequences) + "\n")
        fout.close()
        row = [f.id, "biox", "chromosome", 1, len(f.sequence), ".", ".", ".", "ID=%s;Name=%s" % (f.id, f.id)]
        f_gff3.write("\t".join(str(x) for x in row) + "\n")
    f_gff3.close()
    
def prepare_fasta_mapability(input, output, L = 35):
    f_output = open(output, "wt")
    id_seq = 1
    f_input = biox.data.Fasta(input)
    while f_input.read():
      for i in xrange(0, len(f_input.sequence)-L):
        f_output.write(">%s\n%s\n" % (id_seq, f_input.sequence[i:i+L]))
        id_seq += 1
    f_output.close()
