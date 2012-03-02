import os
import sys
import subprocess
import gzip

endings = [".gzip", ".gz", ".bz2"]

def cmd(command):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
    return output, error
    
def decompress(source, dest=None):
    if dest==None:
        for end in endings:
            if source.endswith(end):
                dest = source[:-len(end)]
                break
    
    if source.endswith(".gzip") or source.endswith(".gz"):
        command = "gunzip -c %s > %s" % (source, dest)
        print command
        out, err = cmd(command)
        return dest
    
    return source
