class TabReader():
    def __init__(self, filename="", header_at = 0, fast = False):
        self.f = None
        self.header = None
        self.fast = fast
        if filename!="":
            self.initialize(filename, header_at = header_at)
            
    def initialize(self, filename, header_at = 0):
        if filename.endswith(".gz"):
            self.f = gzip.open(filename)
        else:
            self.f = open(filename)
        while header_at>0:
            self.f.readline()
            header_at -= 1
        self.header = self.f.readline()
        self.header = self.header.replace("\"", "")
        self.header = self.header.rstrip("\r").rstrip("\n").split("\t")
    
    def readline(self):
        line = self.f.readline()
        if not line:
            return False
        line = line.rstrip("\r").rstrip("\n").split("\t")
        if not self.fast:
            for index, value in enumerate(line):
                try:
                    line[index] = line[index].replace("\"", "")
                    line[index] = float(value)
                except:
                    pass
        self.r = line
        if not self.fast:
            for index, name in enumerate(self.header):
                self.set_value(name, line[index])
                self.set_value("col%s" % index, line[index])
        return True
        
    def set_value(self, name, value):
        name = name.replace(" ", "_").rstrip("\r").rstrip("\n")
        try:
            value = float(value)
        except:
            pass
        setattr(self, name, value)
    