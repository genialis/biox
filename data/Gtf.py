import biox
import os
import sys
import cPickle

cache_data = {}
def cache_string(string):
    if cache_data.get(string, None)==None:
        cache_data[string] = string
        return string
    else:
        return cache_data[string]

class Gtf():

    def __init__(self, filename):
        self.genes = {}
        self.filename = filename
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
            gene = self.genes.get(attrs["gene_id"], biox.data.Gene(attrs["gene_id"], chr, strand, attrs=attrs))
            feature = biox.data.GeneFeature(start, stop, type, gene)
            gene.add_feature(feature)
            self.genes[gene.id] = gene
            self.bin_size = 1000
        self.load_index()
    
    def return_genes(self):
        return self.genes
        
    def load_index(self):
        if not os.path.exists(self.filename+".pindex"):
            pindex = {}
            # create index of exons
            c = 0
            for gene_id, gene in self.genes.items():
                c+=1
                print c
                pindex_chr = pindex.get(gene.chr, {})
                for feature in gene.features:
                    if feature.type!="exon":
                        continue
                    for pos in range(feature.start, feature.stop+1):
                        bin = pos/self.bin_size
                        GL = pindex_chr.get(bin, set())
                        GL.add(cache_string(gene_id))
                        pindex_chr[bin] = GL
                pindex[gene.chr] = pindex_chr
            cPickle.dump(pindex, open(self.filename+".pindex", "wb"), -1)
            self.pindex = pindex
        else:
            self.pindex = cPickle.load(open(self.filename+".pindex"))
    
    def get_genes(self, chr, pos):
        bin = pos/self.bin_size
        candidate_genes = self.pindex.get(chr, {}).get(bin, [])
        position_genes = set()
        for gene_id in candidate_genes:
            for feature in self.genes[gene_id].features:
                if feature.type!="exon":
                    continue
                if feature.start<=pos<=feature.stop:
                    position_genes.add(gene_id)
        return position_genes
    
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