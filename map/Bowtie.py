import biox
from os.path import join as pjoin
import os
import sys
import locale

locale.setlocale(locale.LC_ALL, '')

class Bowtie():

    """
    Mapping of reads to a reference genome
    """

    def __init__(self):
        self.bowtie_exec = pjoin(biox.bowtie_folder, "bowtie")
        self.bowtie_build_exec = pjoin(biox.bowtie_folder, "bowtie-build")
        self.samtools_exec = pjoin(biox.samtools_folder, "samtools")
        self.trim3 = False
        self.sam = True
        self.processors = 2
        self.un_enabled = False
        self.max_enabled = False
        self.mode_fasta = False
        self.m = 1
        self.u = False
        self.l = False
        self.strata = False
        self.enable_n()
        
    def enable_n(self, n=2):
        """
        Enable mode n of bowtie mapper.
        """
        self.n = n
        self.v = False
        
    def set_l(self, l):
        """
        Set seed length (l) and enable n mode.
        """
        self.l = l
        self.enable_n()
        
    def enable_strata(self, strata):
        """
        Enable strata.
        """
        self.strata = strata

    def enable_u(self, u):
        """
        Enable mode u. Map all reads.
        """
        self.u = u
    
    def enable_v(self, v=2):
        """
        Enable mode v
        """
        self.v = 2
        self.n = False
        
    def set_m(self, m):
        """
        Set allowed number of multiple hits per read for the hits to be reported.
        """
        self.m = m
        
    def set_sam(self, sam):
        """
        Enable sam output.
        """
        self.sam = sam
        
    def set_processors(self, processors):
        """
        Set number of processors to use.
        """
        self.processors = processors
        
    def set_mode_fasta(self):
        """
        Enable fasta mode.
        """
        self.mode_fasta = True
        self.enable_v()

    def set_mode_fastq(self):
        """
        Enable fastq mode.
        """
        self.mode_fasta = False
        self.enable_n()
        
    def enable_trim3(self, trim3_iter=None, trim3_step=None):
        """
        Enable iterative trimming of unmapped reads from 3' end.
        
        :param trim3_iter: number of iterations
        :param trim3_step: step of trimming (in nucleotides)
        """
        self.trim3 = True
        self.trim3_iter = trim3_iter if trim3_iter!=None else 5
        self.trim3_step = trim3_step if trim3_step!=None else 2
        
    def disable_trim3(self):
        """
        Disable iterative trimming of unmapped reads from 3' end.
        """
        self.trim3 = False
        
    def enable_un(self, un):
        self.un_enabled = un
        
    def enable_max(self, mx):
        self.max_enabled = mx
        
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
        
    def map(self, index, input, output, index_path=None, create_stats=True, bam=True, delete_temp=True, bam_include_unmapped = False, only_show_cmd = False):
        """
        Map reads to the reference.
        
        :param index: the name of the reference sequence index (dd/dp)
        :param input: full path to FASTA or FASTQ file with reads
        :param output: output folder
        :param index_path: specify full path to genome index, instead of using the indexes in biox.config.bowtie_index_folder folder
        
        """
        
        input_decompressed = biox.utils.decompress(input)
        
        sam_par = "--sam" if self.sam else ""
        bam_unmapped_par = "-F 4" if bam_include_unmapped else ""
        u_par = "-u %s" % self.u if self.u!=None else ""
        n_par = "-n %s" % self.n if self.n!=False else ""
        v_par = "-v %s" % self.v if self.v!=False else ""
        m_par = "-m %s" % self.m if self.m != None else ""
        strata_par = "--strata" if self.strata !=None else ""
        l_par = "-l %s" % self.l if self.l != None else ""
        fasta_par = "-f" if self.mode_fasta==True else ""
        index_par = pjoin(biox.bowtie_index_folder, index)
        if index_path!=None: # specify direct path to bowtie index
            index_par = index_path
        sam_files = []
        stat_files = []
        un_files = []
        bam_files = []
        if self.trim3==False:
            stats_par = "%s.stats" % (output) if create_stats else "/dev/null" 
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
            command = "{bowtie_exec} {strata_par} {l_par} -a {u_par} {processors} {un_par} {max_par} {n_par} {v_par} {sam_par} {m_par} {index_par} {fasta_par} {input_par} 1>{output_par} 2>{output_stats}".format \
            (bowtie_exec = self.bowtie_exec, fasta_par = fasta_par, index_par = index_par, input_par = input_par, output_par = output_par, \
            sam_par = sam_par, n_par = n_par, v_par = v_par, output_stats = stats_par, m_par = m_par, \
            un_par = un_par, max_par = max_par, processors = "-p %s" % self.processors, u_par = u_par, strata_par=strata_par, l_par = l_par)
            if not only_show_cmd:
                str, err = biox.utils.cmd(command)
            else:
                print command            
            sam_files.append(output_par)
            if self.un_enabled:
                un_files.append(output+".unmapped")
            if bam: # create bam file
                bam_output = "%s.bam" % output
                command = "{samtools_exec} view -F 4 {bam_unmapped_par} -bS {sam} > {bam}".format(samtools_exec = self.samtools_exec, sam = output_par, bam = bam_output, bam_unmapped_par = bam_unmapped_par)
                if not only_show_cmd:
                    str, err = biox.utils.cmd(command)
                else:
                    print command
            if create_stats and not only_show_cmd:
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
                stats_par = "%s.trim%s.stats" % (output, trim3) if create_stats else "/dev/null" 
                
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
                command = "{bowtie_exec} {strata_par} {l_par} -a {u_par} {processors} --un {un_par} {max_par} {trim3_par} {n_par} {v_par} {sam_par} {m_par} {index_par} {fasta_par} {input_par} 1>{output_par} 2>{output_stats}".format \
                (bowtie_exec = self.bowtie_exec, fasta_par = fasta_par, index_par = index_par, input_par = input_par, output_par = output_par, \
                sam_par = sam_par, n_par = n_par, v_par = v_par, output_stats = stats_par, m_par = m_par, trim3_par = trim3_par, \
                un_par = un_par, max_par = max_par, processors = "-p %s" % self.processors, u_par = u_par, strata_par=strata_par, l_par = l_par)
                if not only_show_cmd:
                    str, err = biox.utils.cmd(command)
                else:
                    print command                
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
                    command = "{samtools_exec} view -F 4 {bam_unmapped_par} -bS {sam} > {bam}".format(samtools_exec = self.samtools_exec, sam = output_par, bam = bam_file, bam_unmapped_par = bam_unmapped_par)
                    bam_files.append(output_par[:-3]+"bam")
                    if not only_show_cmd:
                        str, err = biox.utils.cmd(command)
                    else:
                        print command                    
                if create_stats and not only_show_cmd:
                    stat_files.append(stats_par)
                    reads[trim3], mapped[trim3] = self.read_statfile(stats_par)
            # merge bam files
            if bam:
                bam_output = "%s.bam" % output
                command = "{samtools_exec} merge -h {header_sam} -f {bam_output} {bam_files}".format(samtools_exec = self.samtools_exec, \
                header_sam = output+".trim0.sam", bam_output = bam_output, bam_files = " ".join(bam_files))
                if not only_show_cmd:
                    str, err = biox.utils.cmd(command)
                else:
                    print command                
            # create trim statistics
            if create_stats and not only_show_cmd:
                all_reads = 0
                mapped_reads = 0
                trims = reads.keys()
                trims.sort()
                stats_par_tab = "%s.stats.tab" % (output)
                f = open(stats_par_tab, "wt")
                f.write("trim3 size\treads processed\treads mapped\treads mapped (perc.)\n")
                for trim in trims:
                    perc = (float(mapped.get(trim, 0)) / max(float(reads.get(trim, 1)), 1) ) * 100
                    f.write("%s\t%s\t%s\t%.2f%%\n" % (trim, locale.format("%d", reads.get(trim, 0), True), locale.format("%d", mapped.get(trim, 0), True), perc))
                    all_reads = max(all_reads, float(reads.get(trim, 0)))
                    mapped_reads += float(mapped.get(trim, 0))
                f.write("all\t100%%\t%.2f%%\t\n" % (mapped_reads/max(all_reads, 1) * 100))
                f.close()
                
        if bam: # sort and index bam file
            bam_sorted = bam_output[:-4]+"_sorted"
            command = "{samtools_exec} sort {bam} {bam_sorted}".format(samtools_exec = self.samtools_exec, bam = bam_output, bam_sorted = bam_sorted)
            if not only_show_cmd:
                str, err = biox.utils.cmd(command)
            else:
                print command            
            command = "{samtools_exec} index {bam_sorted}".format(samtools_exec = self.samtools_exec, bam_sorted = bam_sorted+".bam")
            if not only_show_cmd:
                str, err = biox.utils.cmd(command)
            else:
                print command
            if not only_show_cmd:
                os.rename(bam_sorted+".bam", output+".bam")
                os.rename(bam_sorted+".bam.bai", output+".bam.bai")

        if delete_temp and not only_show_cmd:
            for stat_file in stat_files:
                if os.path.exists(stat_file):
                    os.remove(stat_file)
        if delete_temp and bam: # remove sam
            for sam_file in sam_files:
                if os.path.exists(sam_file):
                    os.remove(sam_file)
            for bam_file in bam_files:
                if os.path.exists(bam_file):
                    os.remove(bam_file)
        if delete_temp and self.un_enabled in [False, None]:
            for un_file in un_files:
                if os.path.exists(un_file):
                    os.remove(un_file)
        
        if input_decompressed!=input and not only_show_cmd:
            try:
                os.remove(input_decompressed)
            except:
                pass
        
        return True
        
    def make_index(self, fasta, index_name):
        output = pjoin(biox.bowtie_index_folder, index_name)
        command = "{bowtie_build_exec} {fasta} {output}".format \
        (bowtie_build_exec = self.bowtie_build_exec, fasta = fasta, output = output)
        print command
        std, err = biox.utils.cmd(command)
