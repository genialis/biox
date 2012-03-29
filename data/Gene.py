class Gene():

    def __init__(self, id = None, chr = None, strand = None, name = None, alias = None):
        self.id = id
        self.name = name
        self.chr = chr
        self.strand = strand
        self.alias = alias
        self.features = []
        self.start = ()
        self.stop = 0
    
    def add_feature(self, feature_new):
        if len(self.features)==0:
            self.features.append(feature_new)
            return True
        insert_at = 0
        for feature in self.features:
            if feature.start>=feature_new.start:
                break
            insert_at += 1
        self.start = min(feature_new.start, self.start)
        self.stop = max(feature_new.stop, self.stop)
        self.features.insert(insert_at, feature_new)