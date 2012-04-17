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
    print "done"
    
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
    print "done"

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
    print "done"
    
    return genes_data
