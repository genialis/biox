import os
import sys
import subprocess
import gzip
import biox

def interval_overlap(start1, stop1, start2, stop2):
    if stop1 < start2 or stop2 < start1:
        return 0
    else:        
        return max(0, min(stop1, stop2) - max(start1, start2)) + 1

def merge_ints(ints):

    def take_next(ints):
        to_return = ints[0]
        del ints[0]
        return to_return      
        
    def add_new(ints, (start, stop)):
        if len(ints)==0:
            ints.append((start, stop))
            return
        (start_0, stop_0) = ints[-1]
        if start_0 <= start <= stop_0:
            del ints[-1]
            ints.append((start_0, max(stop_0, stop)))
        else:
            ints.append((start, stop))
            
    ints.sort()
    ints_merged = []            
    while len(ints)>0:
        (start, stop) = take_next(ints)
        add_new(ints_merged, (start, stop))
    ints_merged.sort()
    return ints_merged

    # ints = [(1,5), (6,10)]
    # assert(join_overlapping(ints)==[(1,5), (6,10)])

    # ints = [(6,10), (1,5)]
    # assert(join_overlapping(ints)==[(1,5), (6,10)])

    # ints = [(1,5), (2,6)]
    # assert(join_overlapping(ints)==[(1,6)])

    # ints = [(2,6), (1,5)]
    # assert(join_overlapping(ints)==[(1,6)])

    # ints = [(1,5), (1,5)]
    # assert(join_overlapping(ints)==[(1,5)])

    # ints = [(1,5), (1,5)]
    # assert(join_overlapping(ints)==[(1,5)])


#process = subprocess.Popen(command, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.PIPE, executable='/bin/bash')
# -cl : this loads the path (loads .bashrc) for the user running the PIPAx wsgi daemon
endings = [".gzip", ".gz", ".bz2"]

class Cmd():

    def __init__(self, command):
        self.command = command
        if biox.os_shell!="":
            self.process = subprocess.Popen(['/bin/bash', '-cl', command], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            self.process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.pid = self.process.pid
        
    def start(self):
        output, error = self.process.communicate()
        self.return_code = self.process.returncode
        return output, error

def cmd(command, shell=True):
    if biox.os_shell!="":
        process = subprocess.Popen(['/bin/bash', '-cl', command], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
    return output, error

def cmd_pipe(command, shell=True):
    if biox.os_shell!="":
        process = subprocess.Popen(['/bin/bash', '-cl', command], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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
    
def process_exists(pid, os="linux"):
    if os=="linux":
        output, error = biox.utils.cmd("ps -p %s" % pid)
        output = output.split("\n")
        return len(output)>=3
    if os=="windows":
        from win32com.client import GetObject
        WMI = GetObject('winmgmts:')
        processes = WMI.InstancesOf('Win32_Process')
        plist = [process.Properties_('ProcessID').Value for process in processes]
        return pid in plist

def web_2_txt(string):
    return string.replace("\"", "").replace("'", "").replace("%2C", ",").replace("%2F", "/").replace("%3C", "<").replace("%3D", "=").replace("%3E", ">").replace("%22", "\"").replace("%2B", "+").replace("%3B", ";").replace("%0A", "").replace("%25", "%").replace("%27", "'").replace("%28", "(").replace("%29", ")").replace("%26", "&").replace("%09", "").replace("%23", "#")
