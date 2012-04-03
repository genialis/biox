"""
Biox
====

Biox is a set of python modules, designed to ease next-generation sequencing data analysis.
"""

import biox
import data
import map
import utils
from config import *

def gff3_from_fasta(fasta_file):
    f = biox.data.Fasta(fasta_file)
    while f.read():
        row = [f.id, "biox", "chromosome", "1", len(f.sequence), ".", ".", ".", "ID=%s;Name=%s" % (f.id, f.id)]
        print "\t".join(str(x) for x in row)
