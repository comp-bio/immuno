#!/usr/bin/env python3
from config import *
from core import *


def main():
    prm = params.Params()
    stt = state.State(prm)

    # ------------------------------ #
    # 1. Read & Split
    if prm.src:
        stt.seq = reader.src_parse(prm.src)
        log(f"Total sequences: {len(stt.seq)}")

        chop.netchop(stt, prm)

        stt.all_peps = reader.splitter(stt.seq, prm.c)
        log(f"Total peptides:  {len(stt.all_peps)}")

        if len(stt.all_peps) == 0:
            help("Peptides not found!")

    # ------------------------------ #
    # 2. Predict Immunogenicity
    ordered = immunogenicity.predict(stt)

    # ------------------------------ #
    # 3. Predict Affinity
    page_size = 500
    min_score = 0
    iteration = 0

    for start in range(0, len(ordered), page_size):
        iteration += 1
        page = ordered[start:(start + page_size)]
        max_possibe_score = affinity.score(1, stt.imm[page[0]])
        if max_possibe_score <= min_score:
            stt.done()
            break
        affinity.predict(stt, prm, page, iteration)
        min_score, max_score = align.clustering(stt, prm)
        align.get_clusters(stt)

    # ------------------------------ #
    # 4. External DB search
    log(f"External DB search:")

    ieatlas_results = db.search(db.IEAtlas(prm), stt.all_peps)
    db.export_file(prm, "IEAtlas", ieatlas_results)

    iedb_results = db.search(db.IEDB(prm), stt.all_peps)
    db.export_file(prm, "IEDB", iedb_results)


if __name__ == "__main__":
    main()
