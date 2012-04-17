import biox
import math

def gene_expression(gtf_file, bam_file, quality = 30, genes = None):
    result = {}
    gtf = biox.data.Gtf(gtf_file)
    all = len(gtf.genes.keys())
    current = 0
    for gene_id, gene in gtf.genes.items():
        if genes!=None and gene_id not in genes:
            continue
        raw_count = 0
        coding_positions = []
        for feature in gene.features:
            if feature.type!="exon":
                continue
            coding_positions += range(feature.start, feature.stop + 1)
        coding_positions = set(coding_positions)
        
        # positive strand
        command = "samtools view -F 4 -q {quality} -F 0x0010 {bam_file} {chr}:{start}-{stop}".format(bam_file = bam_file, chr = gene.chr, start = gene.start, stop = gene.stop, quality = quality)
        output, error = biox.utils.cmd(command)
        output = output.split("\n")
        for line in output:
            line = line.split("\t")
            if len(line)>3:
                position = int(line[3])
                if position in coding_positions:
                    raw_count += 1

        # negative strand
        command = "samtools view -F 4 -q {quality} -f 0x0010 {bam_file} {chr}:{start}-{stop}".format(bam_file = bam_file, chr = gene.chr, start = gene.start, stop = gene.stop, quality = quality)
        output, error = biox.utils.cmd(command)
        output = output.split("\n")
        for line in output:
            line = line.split("\t")
            if len(line)>3:
                seq_len = len(line[9])
                position = int(line[3])+seq_len-1
                if position in coding_positions:
                    raw_count += 1
        result[gene_id] = raw_count
        current += 1
        if current%10==0:
            print "%.2f" % (float(current)*100/all)
    return result
