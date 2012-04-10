import biox

def gene_expression(gtf_file, bam_file, genes = None):
    result = {}
    gtf = biox.data.Gtf(gtf_file)
    for gene_id, gene in gtf.genes.items():
        if genes!=None and gene_id not in genes:
            continue
        raw_count = 0
        for feature in gene.features:
            if feature.type!="exon":
                continue
            command = "samtools view -F 4 -q 30 %s %s:%s-%s" % (bam_file, gene.chr, feature.start, feature.stop)
            output, error = biox.utils.cmd(command)
            output = output.split("\n")
            for line in output:
                line = line.split("\t")
                if len(line)>3:
                    if feature.start<=int(line[3])<=feature.stop:
                        raw_count += 1
        result[gene_id] = raw_count
    return result
