import biox

f = biox.data.TabReader("data.tab")
while f.read():
    print f.Gene_name # string:Gene_name, note that "Gene name" from the tab file was renamed to "Gene_name" (space to _)
    print f.Alias # string:Alias
    print f.Expression # float:Expression
    print f.r # this is a list of all columns at the current line
    print f.r[0] # same as f.Gene_name
