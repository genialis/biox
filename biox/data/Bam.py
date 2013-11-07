import biox
from os.path import join as pjoin

class Bam():
    
    def __init__(self, filename):
        self.samtools_exec = pjoin(biox.samtools_folder, "samtools")
        self.filename = filename
        
    def get_coverage(self, chr, strand, start, stop):
        raw_count = 0
        command = "{samtools_exec} view -F 4 -q 30 {filename} {chr}:{start}-{stop}".format(samtools_exec = self.samtools_exec, filename = self.filename, chr = chr, start = start, stop = stop)
        output, err = biox.utils.cmd(command)
        output = output.split("\n")
        for line in output:
            line = line.split("\t")
            if len(line)>3:
                if int(line[3])<=start and int(line[7])<=stop:
                    raw_count += 1
        return raw_count
