import os
import sys
import subprocess
import gzip

endings = [".gzip", ".gz", ".bz2"]

class Cmd():

    def __init__(self, command):
        self.command = command
        self.process = subprocess.Popen(['/bin/bash', '-cl', command], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.pid = self.process.pid
        
    def start(self):
        output, error = self.process.communicate()
        self.return_code = self.process.returncode
        return output, error

def cmd(command, shell=True):
    #process = subprocess.Popen(command, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.PIPE, executable='/bin/bash')
    # -cl : this loads the path (loads .bashrc) for the user running the PIPAx wsgi daemon
    process = subprocess.Popen(['/bin/bash', '-cl', command], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
    return output, error

def cmd_pipe(command, shell=True):
    #process = subprocess.Popen(command, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process = subprocess.Popen(['/bin/bash', '-cl', command], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return process.stdout
    
def decompress(source, dest=None):
    """
    Decompress gzip or bzip2 files. The original compressed file is left untouched.
    
    :param dest: the uncompressed output is written to this file. If omitted, the uncompressed content is written to source name but without ".gz" or ".bz2". The original compressed file is left untouched.
    """
    
    if dest==None:
        for end in endings:
            if source.endswith(end):
                dest = source[:-len(end)]
                break
    
    if source.endswith(".gzip") or source.endswith(".gz"):
        command = "gunzip -c %s > %s" % (source, dest)
        out, err = cmd(command)
        return dest
        
    if source.endswith(".bz2"):
        command = "bunzip2 -c %s > %s" % (source, dest)
        out, err = cmd(command)
        return dest

    return source # no decompression

def gzip(source):
    """
    Compress (gzip) input.
    """
    
    command = "gzip -f %s" % (source)
    out, err = cmd(command)
    return source+".gz"
