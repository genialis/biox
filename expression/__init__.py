import biox
import math

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
    
def bam_coverage(bam_file, chr="chr1", strand=None, start=1, stop=None):
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
            pos = int(line[3]) if strand=="+" else int(line[3])+len(line[9])-1
            if pos<start or pos>stop:
                continue
            result[pos] = result.setdefault(pos, 0) + 1
    result_list = []
    for pos in range(start, stop+1):
        result_list.append(result.get(pos, 0))
    return result_list
