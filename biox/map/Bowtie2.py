import biox
from os.path import join as pjoin
import os
import sys
import locale

#locale.setlocale(locale.LC_ALL, '')

class Bowtie2():

    """
    Mapping of reads to a reference genome
    """

    def __init__(self):
        self.bowtie_exec = pjoin(biox.bowtie2_folder, "bowtie2")
        self.bowtie_build_exec = pjoin(biox.bowtie2_folder, "bowtie2-build")
        self.samtools_exec = pjoin(biox.samtools_folder, "samtools")
        self.trim3 = False
        self.processors = 2
        self.mode_fasta = False
        self.un_enabled = False
        self.quality = "phred64-quals"
        self.mode_par = "--very-fast-local"
        self.u = None
        self.X = 1000

    def enable_u(self, u):
        """
        Enable mode u. Map all reads.
        """
        self.u = u

    def set_processors(self, processors):
        """
        Set number of processors to use.
        """
        self.processors = processors

    def set_X(self, X):
        self.X = X

    def set_mode_fasta(self):
        """
        Enable fasta mode.
        """
        self.mode_fasta = True

    def set_mode_fastq(self):
        """
        Enable fastq mode.
        """
        self.mode_fasta = False

    def set_mode(self, mode_par):
        """
        Set alignment mode
        """
        self.mode_par = mode_par

    def enable_un(self, un):
        self.un_enabled = un

    def set_quality(self, quality):
        """
        Set qualities of data from FASTQ. Possible options as in bowtie manual.
        :param quality: "solexa-quals", "phred33-quals", "phred64-quals", "integer-quals"
        """
        self.quality = quality

    def map(self, index, input, output, index_path=None, create_stats=True, bam=True, delete_temp=True, bam_include_unmapped = False, simulate = False, verbose=False, keep_unmapped = True):
        """
        Map reads to the reference.

        :param index: the name of the reference sequence index (dd/dp)
        :param input: full path to FASTA or FASTQ file with reads
        :param output: output folder
        :param index_path: specify full path to genome index, instead of using the indexes in biox.config.bowtie_index_folder folder        :param verbose: Print each command that is executed (verbose mode)
        :param simulate: Only simulate mappings and doesn't execute any commands (together with verbose, prints out all the mapping commands)
        """

        u_par = "-u %s" % self.u if self.u != None else ""

        executed_commands = []
        output_log = open("%s.log" % output, "wt")
        paired_input = False

        # input decompressed is the same as input (bowtie2 reads .gz)
        if type(input)==list:
            if len(input)==2:
                paired_input = True
                input_decompressed = [input[0], input[1]]
            else:
                input_decompressed = input[0]
        else:
            input_decompressed = input

        bam_unmapped_par = "-F 4" if bam_include_unmapped else ""
        fasta_par = "-f" if self.mode_fasta==True else ""
        if index_path!=None: # specify direct path to bowtie index
            index_par = index_path
        else:
            index_par = pjoin(biox.bowtie2_index_folder, index)
        sam_files = []
        stat_files = []
        un_files = []
        bam_files = []
        stats_par = "%s.stats" % (output) if create_stats else "/dev/null"
        if type(input_decompressed)==list:
            if len(input_decompressed)==2:
                input_par = "-1 %s -2 %s -X %s" % (input_decompressed[0], input_decompressed[1], self.X)
            else:
                input_par = "-U %s" % input_decompressed[0]
        else:
            input_par = "-U %s" % input_decompressed
            if input_decompressed.endswith(".fasta"):
                self.set_mode_fasta()
        output_par = "%s.sam" % output
        un_par = "--un-gz "+output+".unmapped.gz" if self.un_enabled else ""
        command = "{bowtie_exec} --{quality} {processors} {un_par} {u_par} {mode_par} {index_par} {fasta_par} {input_par} 1>{output_par} 2>{output_stats}".format \
        (bowtie_exec = self.bowtie_exec, fasta_par = fasta_par, index_par = index_par, input_par = input_par, output_par = output_par, u_par = u_par, \
        output_stats = stats_par, un_par = un_par, processors = "-p %s" % self.processors, quality = self.quality, mode_par = self.mode_par)
        output_log.write(command+"\n")
        executed_commands.append(command)
        if verbose:
            print command
        if not simulate:
            str, err = biox.utils.cmd(command)
        sam_files.append(output_par)
        if self.un_enabled:
            un_files.append(output+".unmapped.gz")
        if bam: # create bam file
            bam_output = "%s.bam" % output
            command = "{samtools_exec} view -F 4 {bam_unmapped_par} -bS {sam} > {bam}".format(samtools_exec = self.samtools_exec, sam = output_par, bam = bam_output, bam_unmapped_par = bam_unmapped_par)
            if verbose:
                print command
            if not simulate:
                str, err = biox.utils.cmd(command)
            output_log.write(command+"\n")
            executed_commands.append(command)
        if create_stats and not simulate:
            print "stats"

        if bam: # sort and index bam file
            bam_sorted = bam_output[:-4]+"_sorted"
            command = "{samtools_exec} sort -o {bam_sorted} {bam}".format(samtools_exec = self.samtools_exec, bam = bam_output, bam_sorted = bam_sorted)
            if verbose:
                print command
            if not simulate:
                str, err = biox.utils.cmd(command)
            output_log.write(command+"\n")
            executed_commands.append(command)
            command = "{samtools_exec} index {bam_sorted}".format(samtools_exec = self.samtools_exec, bam_sorted = bam_sorted+".bam")
            if verbose:
                print command
            if not simulate:
                str, err = biox.utils.cmd(command)
            output_log.write(command+"\n")
            executed_commands.append(command)
            if not simulate:
                if verbose:
                    print "mv %s %s\n" % (bam_sorted+".bam", output+".bam")
                    print "mv %s %s\n" % (bam_sorted+".bam.bai", output+".bam.bai")
                os.rename(bam_sorted+".bam", output+".bam")
                os.rename(bam_sorted+".bam.bai", output+".bam.bai")
                output_log.write("mv %s %s\n" % (bam_sorted+".bam", output+".bam"))
                executed_commands.append("mv %s %s" % (bam_sorted+".bam", output+".bam"))
                output_log.write("mv %s %s\n" % (bam_sorted+".bam.bai", output+".bam.bai"))
                executed_commands.append("mv %s %s\n" % (bam_sorted+".bam.bai", output+".bam.bai"))
                if verbose:
                    print "mv %s %s\n" % (bam_sorted+".bam", output+".bam")
                    print "mv %s %s\n" % (bam_sorted+".bam.bai", output+".bam.bai")

        # delete temp files
        if bam and os.path.exists(output+".sam"):
            os.remove(output+".sam")

        output_log.close()
        return executed_commands

    def make_index(self, fasta, index_name):
        output = pjoin(biox.bowtie_index_folder, index_name)
        command = "{bowtie_build_exec} {fasta} {output}".format \
        (bowtie_build_exec = self.bowtie_build_exec, fasta = fasta, output = output)
        print command
        std, err = biox.utils.cmd(command)
        # also make color index
        command = "{bowtie_build_exec} --color {fasta} {output}".format \
        (bowtie_build_exec = self.bowtie_build_exec, fasta = fasta, output = output+"_color")
        print command
        std, err = biox.utils.cmd(command)
