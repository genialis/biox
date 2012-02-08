"""
Reads FASTA files

Example:

f = FASTAreader(file_name)
while f.read():
    print f.id
"""

import gzip

class Fasta:

    def __init__(self, file_name):
        if file_name.endswith(".gz"):
            self.f = gzip.open(file_name, "rt")
        else:
            self.f = open(file_name, "rt")
        self.id = ""
        self.next_id = None
        
    def read(self):
        self.sequence = []
        if self.next_id!=None:
            self.id = self.next_id
            self.next_id = None
        r = self.f.readline()
        if not r:
            return False
        while r:
            r = r.rstrip("\r").rstrip("\n")
            if r[0]==">" and self.id=="":
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
