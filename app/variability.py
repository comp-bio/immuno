#!/usr/bin/env python3
import sys
from collections import Counter
from core.biogl import fasta_parse, translate, rev_comp
from config import *


def offset(seq: str, pos: int) -> int:
    """ Returns the index in the aligned sequence 
    corresponding to the pos-th letter, excluding the leading '-' """
    count = 0
    for i, ch in enumerate(seq):
        if ch == "-":
            continue
        if count == pos:
            return i
        count += 1
    return -1


def frame2pos(seq, frame, pos, pep):
    """ Finds a peptide in the original sequence, 
    taking into account the reading frame and direction """
    p1 = int(pos) * 3 + int(frame[1])
    p2 = p1 + len(pep) * 3
    if frame[3] != 'R':
        return (offset(seq, p1), offset(seq, p2))

    seq = seq[::-1]
    o1, o2 = (offset(seq, p1), offset(seq, p2))
    ro1, ro2 = (len(seq) - o2, len(seq) - o1)
    return ro1, ro2

    if o1 == -1 or o2 == -1:
        print("!!!!", seq, frame, pos, pep, p1, p2, o1, o2)


def main():
    prm = params.ParamsSim()

    sequences = {}
    for name, s in fasta_parse(prm.src, separator=""):
        sequences[name] = s.strip().upper()

    seqs = list(sequences.values())
    column_matches = []
    for i in range(len(seqs[0])):
        col = [s[i] for s in seqs]
        counts = Counter(col)
        column_matches.append(max(counts.values()))

    def math_score(o1, o2):
        total = len(sequences.keys())
        current = sum([column_matches[i] for i in range(o1, o2)])
        return current / ((o2 - o1) * total)

    scores = {}
    with open(prm.results, 'r') as res:
        head = res.readline().replace('\n', '').split('\t')
        i = head.index("Locations")
        for line in res.readlines():
            item = line.replace('\n', '').split('\t')
            for loc in item[i].split(', '):
                name, frame, pos = loc.split(':')
                seq = str(sequences[name])
                o1, o2 = frame2pos(seq, frame, pos, item[1])

                key = f"{o1}-{o2}-{frame[3]}"
                if key not in scores:
                    s = math_score(o1, o2)
                    scores[key] = f"{s:.5f}\t{item[1]}\t{frame[3]}\t{o1}-{o2}\n"

                # Debug
                """
                part = seq[o1:o2].replace('-', '')
                trans = translate(part).replace('X', '*')
                print(name, frame, pos,
                      f"[{part}] > [{trans}]", item[1], o1, o2, sep='\t')
                """

    with open(prm.out, 'w+') as out:
        out.write("Variability\tPepide\tFrame\tPosition\n")
        for key in scores:
            out.write(scores[key])
        echo("Done" + (" " * 20) + "\n", 36)
        log(f"Results:        {prm.out}")


if __name__ == "__main__":
    main()
