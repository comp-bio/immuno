# Table2Matrix
# Diagonal: Min:   4.0  Max:  11.0  Mean:  5.80  Sigma: 1.90
# Rest:     Min:  -4.0  Max:   3.0  Mean: -1.43  Sigma: 1.51
# Matrix:   Min:  -4.0  Max:  11.0  Mean: -0.74  Sigma: 2.63

class Blosum:
    def __init__(self, file):
        def splt(x): return [t.strip() for t in x.replace('\n', '').split(',')]
        with open(file, 'r') as bs:
            self.header = splt(bs.readline())
            matrix = [splt(ln) for ln in bs.readlines()]
        for l, ln in enumerate(matrix):
            for c, val in enumerate(ln):
                matrix[l][c] = matrix[c][l] if val == "" else int(val)
        self.matrix = matrix

    def get(self, a, b):
        if a in self.header and b in self.header:
            return self.matrix[self.header.index(a)][self.header.index(b)]
        return -10

    def max(self, seq):
        value = 0
        for a in seq:
            if a not in self.header:
                value -= 4
            else:
                idx = self.header.index(a)
                value += self.matrix[idx][idx]
        return value
