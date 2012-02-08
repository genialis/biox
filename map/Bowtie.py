import biox
from os.path import join as pjoin
import os
import sys

class Bowtie():

    def __init__(self):
        self.bowtie_exec = pjoin(biox.bowtie_folder, "bowtie")
        self.bowtie_build_exec = pjoin(biox.bowtie_folder, "bowtie-build")
        self.samtools_exec = pjoin(biox.map.samtools_folder, "samtools")
        self.trim3 = False
        self.sam = True
        self.mode_n = True
        self.mode_v = False
        self.un_enabled = None
        self.max_enabled = None
        self.n = self.v = 2
        self.m = None
        
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
        
    def enable_trim3(self, trim3_iter=None, trim3_step=None):
        self.trim3 = True
        self.trim3_iter = trim3_iter if trim3_iter!=None else 5
        self.trim3_step = trim3_step if trim3_step!=None else 2
        
    def disable_trim3(self):
        self.trim3 = False
        
    def enable_un(self):
        self.un_enabled = True
        
    def disable_un(self):
        self.un_enabled = False
        
    def enable_max(self):
        self.max_enabled = True
        
    def disable_max(self):
        self.max_enabled = False
        
    def read_statfile(self, file):
        reads = mapped = 0
        f = open(file, "rt")
        r = f.readline()
        while r:
            r = r.rstrip("\r").rstrip("\n")
            if r.startswith("# reads processed: "):
                reads = int(r.split("# reads processed: ")[1])
            if r.startswith("# reads with at least one reported alignment: "):
                mapped = int(r.split("# reads with at least one reported alignment: ")[1].split(" ")[0])
            r = f.readline()
        return (reads, mapped)
        
    def map(self, index, input, output, stats=None, bam=True, delete_temp=True):
        sam_par = "--sam" if self.sam else ""
        n_par = "-n %s" % self.n if self.mode_n else ""
        v_par = "-v %s" % self.v if self.mode_v else ""
        m_par = "-m %s" % self.m if self.m != None else ""
        fasta_par = "-f" if self.mode_fasta==True else ""
        index_par = pjoin(biox.map.bowtie_index_folder, index)
        sam_files = []
        stat_files = []
        un_files = []
        bam_files = []
        if self.trim3==False:
            stats_par = "%s.stats" % (output) if stats==True else "/dev/null" 
            if type(input)==list:
                if len(input)==2:
                    input_par = "-1 %s -2 %s" % (input[0], input[1])
                else:
                    input_par = input[0]
            else:
                input_par = input
            output_par = "%s.sam" % output
            un_par = "--un "+output+".unmapped" if self.un_enabled else ""
            max_par = "--max "+output+".maxmulti" if self.max_enabled else ""
            command = "{bowtie_exec} {un_par} {max_par} {n_par} {v_par} {sam_par} {m_par} {index_par} {fasta_par} {input_par} 1>{output_par} 2>{output_stats}".format \
            (bowtie_exec = self.bowtie_exec, fasta_par = fasta_par, index_par = index_par, input_par = input_par, output_par = output_par, \
            sam_par = sam_par, n_par = n_par, v_par = v_par, output_stats = stats_par, m_par = m_par, \
            un_par = un_par, max_par = max_par)
            str, err = biox.utils.cmd(command)
            sam_files.append(output_par)
            if self.un_enabled:
                un_files.append(output+".unmapped")
            if bam: # create bam file
                bam_output = "%s.bam" % output
                command = "{samtools_exec} view -F 4 -bS {sam} > {bam}".format(samtools_exec = self.samtools_exec, sam = output_par, bam = bam_output)
                str, err = biox.utils.cmd(command)
            if stats:
                stat_files.append(stats_par)
                reads, mapped = self.read_statfile(stats_par)
                stats_par_tab = "%s.stats.tab" % (output)
                f = open(stats_par_tab, "wt")
                f.write("reads processed\treads mapped\n")
                f.write("%s\t%s\n" % (reads, mapped))
                f.close()
        else:
            reads = {}
            mapped = {}
            for t in range(0, self.trim3_iter+1):
                trim3 = t * self.trim3_step
                stats_par = "%s.trim%s.stats" % (output, trim3) if stats==True else "/dev/null" 
                
                if t==0:
                    if type(input)==list:
                        if len(input)==2:
                            input_par = "-1 %s -2 %s" % (input[0], input[1])
                        elif len(input)==1:
                            input_par = input[0]
                    else:
                        input_par = input
                else:
                    if type(input)==list:
                        if len(input)==2:
                            input_par = "-1 %s_1 -2 %s_2" % (un_par, un_par)
                    else:
                        input_par = un_par
                        
                # use _ instead of . here because pair-end mapping naming of
                # bowtie is not working correctly when using dot
                # eg.: output.trim0.unmapped -> output.trim0_1.unmapped, output.trim0_2.unmapped
                # instead of output.trim0.unmapped_1, output.trim0.unmapped_2
                un_par = output+"_trim%s_unmapped" % trim3
                    
                max_par = "--max "+output+".trim%s.maxmulti" % trim3 if self.max_enabled else ""
                trim3_par = "--trim3 %s" % trim3
                output_par = output+".trim%s.sam" % trim3
                command = "{bowtie_exec} --un {un_par} {max_par} {trim3_par} {n_par} {v_par} {sam_par} {m_par} {index_par} {fasta_par} {input_par} 1>{output_par} 2>{output_stats}".format \
                (bowtie_exec = self.bowtie_exec, fasta_par = fasta_par, index_par = index_par, input_par = input_par, output_par = output_par, \
                sam_par = sam_par, n_par = n_par, v_par = v_par, output_stats = stats_par, m_par = m_par, trim3_par = trim3_par, \
                un_par = un_par, max_par = max_par)
                str, err = biox.utils.cmd(command)
                sam_files.append(output_par)
                if type(input)==list:
                    if len(input)==2:
                        un_files.append("%s_1" % un_par)
                        un_files.append("%s_2" % un_par)
                    elif len(input)==1:
                        un_files.append(un_par)
                else:
                        un_files.append(un_par)
                if bam:
                    bam_file = output_par[:-3]+"bam"
                    command = "{samtools_exec} view -F 4 -bS {sam} > {bam}".format(samtools_exec = self.samtools_exec, sam = output_par, bam = bam_file)
                    bam_files.append(output_par[:-3]+"bam")
                    str, err = biox.utils.cmd(command)
                if stats:
                    stat_files.append(stats_par)
                    reads[trim3], mapped[trim3] = self.read_statfile(stats_par)
            # merge bam files
            if bam:
                bam_output = "%s.bam" % output
                command = "{samtools_exec} merge -h {header_sam} {bam_output} {bam_files}".format(samtools_exec = self.samtools_exec, \
                header_sam = output_par, bam_output = bam_output, bam_files = " ".join(bam_files))
                str, err = biox.utils.cmd(command)
            # create trim statistics
            trims = reads.keys()
            trims.sort()
            stats_par_tab = "%s.stats.tab" % (output)
            f = open(stats_par_tab, "wt")
            f.write("trim3 size\treads processed\treads mapped\n")
            for trim in trims:
                f.write("%s\t%s\t%s\n" % (trim, reads.get(trim, 0), mapped.get(trim, 0)))
            f.close()
                
        if bam: # sort and index bam file
            bam_sorted = bam_output[:-4]+"_sorted"
            command = "{samtools_exec} sort {bam} {bam_sorted}".format(samtools_exec = self.samtools_exec, bam = bam_output, bam_sorted = bam_sorted)
            str, err = biox.utils.cmd(command)
            command = "{samtools_exec} index {bam_sorted}".format(samtools_exec = self.samtools_exec, bam_sorted = bam_sorted+".bam")
            str, err = biox.utils.cmd(command)
            os.rename(bam_sorted+".bam", output+".bam")
            os.rename(bam_sorted+".bam.bai", output+".bam.bai")

        if delete_temp:
            for stat_file in stat_files:
                os.remove(stat_file)
        if delete_temp and bam: # remove sam
            for sam_file in sam_files:
                os.remove(sam_file)
            for bam_file in bam_files:
                os.remove(bam_file)
        if delete_temp and self.un_enabled in [False, None]:
            for un_file in un_files:
                os.remove(un_file)
        return True
        
    def make_index(self, fasta, index_name):
        output = pjoin(biox.map.bowtie_index_folder, index_name)
        command = "{bowtie_build_exec} {fasta} {output}".format \
        (bowtie_build_exec = self.bowtie_build_exec, fasta = fasta, output = output)
        std, err = biox.utils.cmd(command)
