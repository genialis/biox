import biox
import math
import os

def gene_expression(gtf_file, bam_file, quality = 30, genes = None):
    gtf = biox.data.Gtf(gtf_file)
    all = len(gtf.genes.keys())
    current = 0

    position_2_gene = {}
    for gene_id, gene in gtf.genes.items():
        if genes!=None and gene_id not in genes:
            continue
        for feature in gene.features:
            if feature.type!="exon":
                continue
            for pos in range(feature.start, feature.stop+1):
                # present = position_2_gene.get("%s_%s" % (gene.chr, pos), None)
                # if present!=None:
                    # print present, gene_id
                position_2_gene["%s_%s" % (gene.chr, pos)] = gene_id
    
    genes_data = {}
    for gene_id in gtf.genes.keys():
        genes_data[gene_id] = 0
 
    command = "samtools view -F 4 -q {quality} -F 0x0010 {bam_file}".format(bam_file = bam_file, quality = quality)
    for line in biox.utils.cmd_pipe(command):
        line = line.split("\t")
        if len(line)>3:
            pos = int(line[3])
            chr = line[2]
            gene_id = position_2_gene.get("%s_%s" % (chr, pos), None)
            if gene_id!=None:
                genes_data[gene_id] = genes_data[gene_id] + 1

    command = "samtools view -F 4 -q {quality} -f 0x0010 {bam_file}".format(bam_file = bam_file, quality = quality)
    for line in biox.utils.cmd_pipe(command):
        line = line.split("\t")
        if len(line)>3:
            seq_len = len(line[9])
            pos = int(line[3])+seq_len-1
            chr = line[2]
            gene_id = position_2_gene.get("%s_%s" % (chr, pos), None)
            if gene_id!=None:
                genes_data[gene_id] = genes_data[gene_id] + 1
    return genes_data

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
