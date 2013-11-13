import math

import biox
import data
import map
import utils
import expression


def gff3_from_fasta(fasta_file):
    f = biox.data.Fasta(fasta_file)
    while f.read():
        row = [f.id, "biox", "chromosome", "1", len(f.sequence), ".", ".", ".", "ID=%s;Name=%s" % (f.id, f.id)]
        print "\t".join(str(x) for x in row)


bowtie_folder = ""
bowtie_index_folder = ""
samtools_folder = ""
os_shell = ""
bedtools_folder = ""
kentutils_folder = ""

try:
    from config import *
except ImportError:
    log.error("Missing settings file config.py.")
