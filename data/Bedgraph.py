class Bedgraph():
    def __init__(self, filename1=None, min_value=None, only_positions=None):
        self.trackname = {"type":"bedGraph", "color": "120,101,172", "altColor" : "200,120,59", "maxHeightPixels":"100:50:0", "visibility":"full", "priority":"20"}
        self.track_header = ""
        self.data = {}
        self.sites = 0
        self.data_sum = 0
        self.norm_factor = 1000000
        self.strand_file = {"+":"", "-": "-"}
        if filename1 != None:
            self.load_from_file(filename1, min_value=min_value, only_positions = only_positions)
            
    def load_from_file(self, filename1, min_value=0, only_positions=None):
        if filename1.endswith(".gz"):
            f = gzip.open(filename1, "rb")
        else:
            f = open(filename1, "rb")
        self.track_header = f.readline()
        r = f.readline()
        while r:
            r = r.rstrip("\r").rstrip("\n").split("\t")
            chr = r[0]
            start = int(r[1])
            stop = int(r[2])
            cDNA = float(r[3])
            if abs(cDNA) < min_value:
                r = f.readline()
                continue
            strand = "+" if cDNA>=0 else "-"
            self.data.setdefault(chr, {}).setdefault(strand, {})
            for p in xrange(start, stop):
                if only_positions!=None:
                    if p not in only_positions.get("%s%s" % (strand, chr), None):
                        continue
                self.data[chr][strand][p] = self.data[chr][strand].setdefault(p, 0) + abs(cDNA)
            self.data_sum += abs(cDNA) * (stop-start)
            r = f.readline()
        f.close()
        
    def chromosomes(self):
        return self.data.keys()
        
    def get_value(self, chr, strand, pos):
        return self.data.get(chr, {}).get(strand, {}).get(pos, 0)

    def peaks(self, flanking):
        data_peaks = {}
        self.data_sum = 0
        for chr, chr_data in self.data.iteritems():
            data_peaks[chr] = {}
            for strand, strand_data in chr_data.iteritems():
                data_peaks[chr][strand] = {}
                pos_sorted = []
                for (pos, val) in strand_data.iteritems():
                    pos_sorted.append((val, pos))
                pos_sorted.sort(reverse=True)
                for (value, pos) in pos_sorted:
                    if self.data[chr][strand].get(pos, 0)==0:
                        continue
                    sum_peak = 0
                    pos_from = max(0, pos-flanking)
                    pos_to = max(0, pos+flanking+1)
                    for i in range(pos_from, pos_to):
                        at_pos = self.get_value(chr, strand, i)
                        sum_peak += self.get_value(chr, strand, i)
                        if at_pos>0:
                            self.data[chr][strand][i] = 0
                    data_peaks[chr][strand][pos] = sum_peak
                    self.data_sum += sum_peak
        # replace original data with peak data
        self.data = data_peaks

    def region(self, chr, strand, pos_from, pos_to):
        return sum([self.get_value(chr, strand, i) for i in xrange(pos_from, pos_to+1)])
        
    def get_region_max(self, chr, strand, pos_from, pos_to):
        return max([self.get_value(chr, strand, i) for i in xrange(pos_from, pos_to+1)])
        
    def window_sum(self, chr, strand, pos, window_size=50):
        sum_window = 0
        pos_from = max(0, pos-int(window_size/2))
        pos_to = max(0, pos+int(window_size/2)+1)
        return sum([self.get_value(chr, strand, i) for i in range(pos_from, pos_to)])
		
    def window(self, window_size, smoothing=False):
        data_window = {}
        self.data_sum = 0
        for chr, chr_data in self.data.iteritems():
            print "window on ", chr
            data_window[chr] = {}
            for strand, strand_data in chr_data.iteritems():
                data_window[chr][strand] = {}
                for pos, value in strand_data.iteritems():
                    position_range = [pos] if not smoothing else range(max(0, pos-window_size/2), max(0, pos+window_size/2+1))
                    for pos_set in position_range:
                        data_window[chr][strand][pos_set] = self.window_sum(chr, strand, pos_set, window_size=window_size)
        # replace original data with peak data
        self.data = data_window		
                    
    def set_name(self, name):
        self.trackname["name"] = name
        
    def normalize(self):
        data_sum_new = 0
        for chr, strand_data in self.data.iteritems():
            for strand, position_data in strand_data.iteritems():
                for position, value in position_data.iteritems():
                    new_value = (value/float(self.data_sum)) * self.norm_factor
                    self.data[chr][strand][position] = new_value
                    data_sum_new += new_value
                    
    def addup(self, Bg):
        for chr, strand_data in Bg.data.iteritems():
            for strand, position_data in strand_data.iteritems():
                for position, value_new in position_data.iteritems():
                    value_current = self.data.setdefault(chr, {}).setdefault(strand, {}).setdefault(position, 0)
                    value_new = value_current + value_new
                    delta = value_new - value_current
                    self.data[chr][strand][position] = value_new
                    self.data_sum += delta                    
