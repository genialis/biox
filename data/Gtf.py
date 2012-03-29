import biox
import os
import sys

class Gtf():

    def __init__(self, filename):
        self.genes = {}
        f = biox.data.TabReader(filename)
        while f.readline():
            chr = f.r[0]
            type = f.r[2]
            start = int(f.r[3])
            stop = int(f.r[4])
            strand = f.r[6]
            attrs = {}
            temp = f.r[-1].split("; ")
            for t in temp:
                t = t.replace(";", "")
                t = t.split(" ")
                attrs[t[0]] = " ".join(t[1:])
            if attrs.get("gene_id", None)==None:
                continue
            gene = self.genes.get(attrs["gene_id"], biox.data.Gene(attrs["gene_id"], chr, strand))
            feature = biox.data.GeneFeature(start, stop, type, gene)
            gene.add_feature(feature)
            self.genes[gene.id] = gene
    
    def return_genes(self):
        return self.genes
    
    def compute_overlap(self, start1, stop1, start2, stop2):
        if stop1 < start2 or stop2 < start1:
            return 0
        else:        
            return max(0, min(stop1, stop2) - max(start1, start2)) + 1
    
    def find_overlap(self, gene_source):
        """
        Return overlapping feature to given feature
        """
        overlapping_genes = []
        for gene_id, gene in self.genes.iteritems():
            if gene.strand != gene_source.strand or gene.chr != gene_source.chr:
                continue
            overlap = self.compute_overlap(gene.start, gene.stop, gene_source.start, gene_source.stop)
            if overlap>0:
                overlapping_genes.append((overlap, gene))
        return overlapping_genes

    def write_gff3(self, filename):
        f = open(filename, "wt")
        for gene_id, gene in self.genes.iteritems():
            row = [gene.chr, "biox", "gene", gene.start, gene.stop, "", gene.strand, ".", "ID=%s" % gene_id] # gene
            f.write("\t".join(str(x) for x in row) + "\n")
            row = [gene.chr, "biox", "mRNA", gene.start, gene.stop, "", gene.strand, ".", "ID=%s.t1;Parent=%s" % (gene_id, gene_id)] # mRNA
            f.write("\t".join(str(x) for x in row) + "\n")
            for feature in gene.features:
                row = [gene.chr, "biox", "CDS", feature.start, feature.stop, "", gene.strand, ".", "ID=%s.t1.cds;Parent=%s.t1" % (gene_id, gene_id)] # mRNA
                f.write("\t".join(str(x) for x in row) + "\n")
        f.close()