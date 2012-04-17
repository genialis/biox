import os
import sys
import subprocess
import gzip

endings = [".gzip", ".gz", ".bz2"]

def cmd(command, shell=True):
    process = subprocess.Popen(command, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
    return output, error
    
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
