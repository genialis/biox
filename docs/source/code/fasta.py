import biox

f = biox.data.fasta("data.fasta")
while f.read():
    print f.id
    print f.sequence
