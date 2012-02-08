"""
Reads FASTQ files

Example:

f = FASTQreader(file_name)
while f.read():
    print f.identifier
"""

import gzip

class Fastq:
    def read_line(self):
        return self.f.readline().rstrip("\r").rstrip("\n")
        
    def read(self):
        self.identifier = self.read_line()
        self.sequence = self.read_line()
        self.plus = self.read_line() # +
        self.quality = self.read_line()
        self.uncut_sequence = self.sequence
        self.uncut_quality = self.quality
        self.records_read += 1
        if self.identifier == "":
                return False
        if self.cut_bad:
            qual = self.quality.rstrip("B")
            if len(qual) <> len(self.quality):
                self.sequence = self.sequence[:len(qual)]
                self.quality = qual
        return True

    def __init__(self, file_name, cut_bad=False):
        self.cut_bad = cut_bad
        if file_name.endswith(".gz"):
            self.f = gzip.open(file_name, "rt")
        else:
            self.f = open(file_name, "rt")
        self.identifier = ""
        self.sequence = ""
        self.plus = ""
        self.quality = ""
        self.records_read = 0