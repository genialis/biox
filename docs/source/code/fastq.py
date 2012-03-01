import biox

f = biox.data.fastq("data.fastq")
while f.read():
    print f.id
    print f.sequence
