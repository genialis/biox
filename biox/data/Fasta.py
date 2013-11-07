import biox
import gzip
import os

class Fasta:

    """
    This class handles FASTA files.
    """

    def __init__(self, filename):
        self.filename = filename
        if os.path.exists(self.filename+".index"):
            self.read_index()
        else:
            self.create_index()
            self.read_index()
        self.open_file()
        
    def open_file(self):
        self.id = None
        self.next_id = None
        if self.filename.endswith(".gz"):
            self.f = gzip.open(self.filename, "rt")
        else:
            self.f = open(self.filename, "rt")
        
    def read(self):
        self.sequence = []
        if self.next_id!=None:
            self.id = self.next_id
            self.next_id = None
        r = self.f.readline()
        if not r:
            return False
        while r:
            r = r.replace("\r", "").replace("\n", "")
            if r=="":
                r = self.f.readline()
                continue
            if r[0]==">" and self.id==None:
                self.id = r[1:]
                r = self.f.readline()
                continue
            elif r[0]==">":
                self.next_id = r[1:]
                self.sequence = "".join(self.sequence)
                return True
            self.sequence.append(r)
            r = self.f.readline()
        self.sequence = "".join(self.sequence)
        return True
        
    def create_index(self):
        self.open_file()
        f_linear = open(self.filename + ".linear", "wt")
        f_index = open(self.filename + ".index", "wt")
        while self.read():
            f_linear.write(self.sequence)
            f_index.write("%s\t%s\n" % (self.id, len(self.sequence)))
        f_linear.close()
        f_index.close()
        
    def read_index(self):
        self.index = {}
        f = open(self.filename+".index", "rt")
        r = f.readline()
        start = 0
        while r:
            r = r.replace("\r", "").replace("\n", "").split("\t")
            self.index[r[0]] = start
            start += int(r[1])
            r = f.readline()

    def get_seq(self, chr, start, stop, reverse_complement=False):
        f = open(self.filename+".linear", "rb")
        num_bytes = stop - start + 1
        pos_from = self.index[chr] + start - 1
        pos_to = pos_from + num_bytes
        f.seek(pos_from)
        seq = f.read(num_bytes)
        f.close()
        if reverse_complement:
            seq = biox.data.Sequence.rev_comp(seq)
        return seq
