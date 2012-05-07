bc = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}

def rev_comp(seq):
    seq = [bc[s] for s in seq]
    return "".join(seq[::-1])