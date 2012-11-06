import biox
import math
import os

pysam_enabled = False

try:
    import pysam
    pysam_enabled = True
except:
    pysam_enabled = False

def testBit(int_type, offset):
   mask = 1 << offset
   return(int_type & mask)

def region_expression_pileup(sam_filename, chr, start, stop):
    import pysam
    f = pysam.Samfile(sam_filename)
    sum = 0
    for rec in f.pileup(chr, start, stop):
        if not start<=rec.pos<=stop:
            continue
        sum += rec.n    
    return round(sum/float(stop-start+1))
   
def gene_expression_overlap(gtf_file, bam_file, quality = 30, genes = None):
    if pysam_enabled:
        pysam_bam = pysam.Samfile(bam_file)
    gtf = biox.data.Gtf(gtf_file)
    genes_exp = {}
    if genes==None:
        for gene_id in gtf.genes:
            genes_exp[gene_id] = 0
    else:
        for gene_id in genes:
            genes_exp[gene_id] = 0
            
    current = 0
    for gene_id, gene in gtf.genes.items():
        current += 1
        if current%300==0:
            print "%.2f" % (float(current)/len(gtf.genes)), bam_file
        if genes!=None and gene.id not in genes:
            continue
        for feature in gene.features:
            if feature.type!="exon":
                continue
            assert(feature.start<=feature.stop)
            if pysam_enabled:
                genes_exp[gene_id] = genes_exp.get(gene_id, 0) + pysam_bam.count(gene.chr, feature.start-1, feature.stop-1)
            else:
                command = "samtools view -F 4 -q {quality} -c {bam_file} {chr}:{start}-{stop}".format(bam_file = bam_file, quality = 30, chr=gene.chr, start=feature.start, stop=feature.stop)
                output, error = biox.utils.cmd(command)
                if output!="":
                    genes_exp[gene_id] = genes_exp.get(gene_id, 0) + int(output)
    return genes_exp
   
def gene_expression(gtf_file, bam_file, quality = 30, genes = None):
    gtf = biox.data.Gtf(gtf_file)
    genes_exp = {}
    if genes==None:
        for gene_id in gtf.genes:
            genes_exp[gene_id] = 0
    else:
        for gene_id in genes:
            genes_exp[gene_id] = 0
            
    command = "samtools view -F 4 -q {quality} -c {bam_file}".format(bam_file = bam_file, quality = quality)
    output, error = biox.utils.cmd(command)
    reads = int(output)            
    
    command = "samtools view -F 4 -q {quality} {bam_file}".format(bam_file = bam_file, quality = quality)
    current = 0
    for line in biox.utils.cmd_pipe(command):
        current += 1
        if current%200000==0:
            print "%.2f" % (current/float(reads)), bam_file
        line = line.split("\t")
        if len(line)>3:
            chr = line[2]
            flag = int(line[1])
            seq_len = len(line[9])
            strand = flag & 16 # is bit 5 set?
            strand = "+" if strand==0 else "-"
            pos = int(line[3]) if strand=="+" else int(line[3])+seq_len-1
            position_genes = gtf.get_genes(chr, pos)
            for gene_id in position_genes:
                if genes!=None and gene_id not in genes:
                    continue
                genes_exp[gene_id] = genes_exp.get(gene_id, 0) + 1
    return genes_exp
    
def gene_expression_promoters(gtf_file, bam_file, quality = 30, genes = None):
    gtf = biox.data.Gtf(gtf_file)
    genes_exp = {}
    if genes==None:
        for gene_id in gtf.genes:
            genes_exp[gene_id] = 0
    else:
        for gene_id in genes:
            genes_exp[gene_id] = 0
            
    command = "samtools view -F 4 -q {quality} -c {bam_file}".format(bam_file = bam_file, quality = quality)
    output, error = biox.utils.cmd(command)
    reads = int(output)            
    
    command = "samtools view -F 4 -q {quality} {bam_file}".format(bam_file = bam_file, quality = quality)
    current = 0
    for line in biox.utils.cmd_pipe(command):
        current += 1
        if current%200000==0:
            print "%.2f" % (current/float(reads)), bam_file
        line = line.split("\t")
        if len(line)>3:
            chr = line[2]
            flag = int(line[1])
            seq_len = len(line[9])
            strand = flag & 16 # is bit 5 set?
            strand = "+" if strand==0 else "-"
            pos = int(line[3]) if strand=="+" else int(line[3])+seq_len-1
            position_genes = gtf.get_genes(chr, pos)
            for gene_id in position_genes:
                if genes!=None and gene_id not in genes:
                    continue
                genes_exp[gene_id] = genes_exp.get(gene_id, 0) + 1
    return genes_exp    

def bam_chromosomes(bam_file):
    chrs = {}
    command = "samtools view -H {bam_file}".format(bam_file = bam_file)
    output, error = biox.utils.cmd(command)
    output = output.split("\n")
    for line in output:
        line = line.split("\t")
        if line[0]=="@SQ":
            chrs[line[1].split("SN:")[1]] = int(line[2].split("LN:")[1])
    return chrs
    
def write_bam_chr(bam_file, chr_file):
    chrs = bam_chromosomes(bam_file)
    f = open(chr_file, "wt")
    for chr, chr_len in chrs.items():
        f.write("%s\t%s\n" % (chr, chr_len))
    f.close()
    return chr_file
    
def write_fasta_chr(fasta_file, chr_file):
    fin = biox.data.Fasta(fasta_file)
    f = open(chr_file, "wt")
    while fin.read():
        f.write("%s\t%s\n" % (fin.id, len(fin.sequence)))
    f.close()
    
def bam_coverage(bam_file, chr="chr1", strand=None, start=1, stop=None, position = '5prime'):
    if strand=="+":
        strand_par = "-F 0x0010"
    elif strand=="-":
        strand_par = "-f 0x0010"
    else:
        strand_par = ""
    chrs = bam_chromosomes(bam_file)
    if start==None:
        start = 1
    if stop==None:
        stop = chrs[chr]
    start = int(start)
    stop = int(stop)
    command = "samtools view {strand_par} {bam_file} {chr}:{start}-{stop}".format(bam_file = bam_file, strand_par=strand_par, chr=chr, start=start, stop=stop)
    result = {}
    for line in biox.utils.cmd_pipe(command):
        line = line.split("\t")
        if len(line)>3:
            strand = "+" if int(line[1])==0 else "-"
            if position=='5prime':
                pos = int(line[3]) if strand=="+" else int(line[3])+len(line[9])-1
                if pos<start or pos>stop:
                    continue
                result[pos] = result.setdefault(pos, 0) + 1
            if position=='span':
                for pos in range(int(line[3]), int(line[3])+len(line[9])):
                    if pos<start or pos>stop:
                        continue                
                    result[pos] = result.setdefault(pos, 0) + 1
    result_list = []
    for pos in range(start, stop+1):
        result_list.append(result.get(pos, 0))
    return result_list
    
def bam2wig(bam_file, wig_file, strand=None, position = 'span', bigWig = True):
    f = open(wig_file, "wt")
    chrs = bam_chromosomes(bam_file)
    for chr, chr_len in chrs.items():
        f.write("variableStep chrom=%s span=1\n" % chr)
        z = bam_coverage(bam_file, chr, strand = strand, position = position)
        for pos, value in enumerate(z):
            if value==0:
                continue
            f.write("%s\t%s\n" % (pos+1, value))
    f.close()
    bam_chrs = write_bam_chr(bam_file, wig_file+".chrs")
    if bigWig:
        command = "wigToBigWig %s %s %s" % (wig_file, bam_chrs, wig_file.replace(".wig", ".bw"))
        output, error = biox.utils.cmd(command)
        os.remove(wig_file)
    os.remove(bam_chrs)
