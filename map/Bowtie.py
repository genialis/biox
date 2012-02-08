import biox
from os.path import join as pjoin

class Bowtie():

    def __init__(self):
        self.bowtie_exec = pjoin(biox.map.bowtie_folder, "bowtie")
        self.bowtie_build_exec = pjoin(biox.map.bowtie_folder, "bowtie-build")
        self.trim3 = False
        self.sam = True
        self.mode_n = True
        self.mode_v = False
        self.n = self.v = 2
        
    def set_mode_n(self):
        self.mode_n = True
        self.mode_v = False
    
    def set_mode_v(self):
        self.mode_v = True
        self.mode_n = False
        
    def set_m(self, m):
        self.m = m
        
    def set_v(self, v):
        self.v = v
        self.set_mode_v()
    
    def set_n(self, n):
        self.n = n
        self.set_mode_n()
    
    def set_a(self, a):
        self.a = a
        
    def set_sam(self, sam):
        self.sam = sam
        
    def set_mode_fasta(self):
        self.mode_fasta = True
        self.set_mode_v()
        
    def trim_3(self, trim3, trim3_iter=None, trim3_step=None):
        self.trim3 = trim3
        self.trim3_iter = trim3_iter
        self.trim3_step - trim3_step
        
    def map(self, index, input, output):
        sam_par = "--sam" if self.sam else ""
        n_par = "-n %s" % self.n if self.mode_n else ""
        v_par = "-v %s" % self.v if self.mode_v else ""
        fasta_par = "-f" if self.mode_fasta==True else ""
        index_par = pjoin(biox.map.bowtie_index_folder, index)
        if self.trim3==False:
            command = "{bowtie_exec} {n_par} {v_par} {sam_par} {index_par} {fasta_par} {input} 1>{output}".format \
            (bowtie_exec = self.bowtie_exec, fasta_par = fasta_par, index_par = index_par, input = input, output = output, sam_par = sam_par, n_par = n_par, v_par = v_par)
            str, err = biox.utils.cmd(command)
        return True
        
    def make_index(self, fasta, index_name):
        output = pjoin(biox.map.bowtie_index_folder, index_name)
        command = "{bowtie_build_exec} {fasta} {output}".format \
        (bowtie_build_exec = self.bowtie_build_exec, fasta = fasta, output = output)
        std, err = biox.utils.cmd(command)
