import biox

f = biox.data.fasta("data.fastq")
while f.read():
    print f.id
    print f.sequence
