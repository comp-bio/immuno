from .affinity import get_affinity
import re


def smith_waterman(seq1, seq2, blosum, gap=-2):
    m, n = len(seq1), len(seq2)
    score_matrix = [[0] * (n + 1) for _ in range(m + 1)]

    max_score = 0
    max_pos = (0, 0)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            score = blosum.get(seq1[i - 1], seq2[j - 1])

            match_score = score_matrix[i - 1][j - 1] + score
            delete_score = score_matrix[i - 1][j] + gap
            insert_score = score_matrix[i][j - 1] + gap
            score_matrix[i][j] = max(
                0, match_score, delete_score, insert_score)

            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)

    aligned_seq1 = []
    aligned_seq2 = []
    i, j = max_pos

    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        current_score = score_matrix[i][j]
        diagonal_score = score_matrix[i - 1][j - 1]
        up_score = score_matrix[i - 1][j]
        left_score = score_matrix[i][j - 1]

        score = blosum.get(seq1[i - 1], seq2[j - 1])

        if current_score == diagonal_score + score:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif current_score == up_score + gap:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        elif current_score == left_score + gap:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    aligned_seq1 = ''.join(reversed(aligned_seq1))
    aligned_seq2 = ''.join(reversed(aligned_seq2))
    dist = 1 - max_score/max(blosum.max(seq1), blosum.max(seq2), 1)
    return dist, max_score, aligned_seq1, aligned_seq2


def clc(state, params):
    # complete linkage clustering
    peptides = [p for p in state.pep]

    n = len(peptides)
    dist = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            SW = smith_waterman(peptides[i], peptides[j], state.blosum62)
            dist[i][j] = dist[j][i] = SW[0]

    n = len(peptides)
    clusters = [{i} for i in range(n)]

    def cluster_distance(c1, c2):
        return max(dist[i][j] for i in c1 for j in c2)

    while True:
        best_pair = None
        best_dist = float('inf')

        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                d = cluster_distance(clusters[i], clusters[j])
                if d <= params.dist and d < best_dist:
                    best_dist = d
                    best_pair = (i, j)

        if not best_pair:
            break

        i, j = best_pair
        merged = clusters[i] | clusters[j]
        clusters.pop(j)
        clusters.pop(i)
        clusters.append(merged)

    return [[peptides[i] for i in cluster] for cluster in clusters]


def clustering(state, params):
    cls = clc(state, params)
    scores = [max([state.pep[p]['score'] for p in cl]) for cl in cls]
    ids = [i for i, _ in sorted(
        enumerate(scores), key=lambda x: x[1], reverse=True)]

    lim = min(params.n, len(cls)-1)
    state.cls = [cls[i] for i in ids[0:lim]]
    filtred = [j for sub in state.cls for j in sub]
    # -> Save peps from first N clusters only
    state.pep = {k: v for k, v in state.pep.items() if k in filtred}
    # Score of min and max clusters
    return (scores[ids[min(len(ids), lim) - 1]], scores[ids[0]])


def finder(state):
    for p in state.pep:
        state.pep[p]['places'] = []
        state.pep[p]['chops'] = []
        for one in state.seq:
            for m in re.finditer(p, one['seq']):
                name = one['name'].replace('___', ':')
                state.pep[p]['places'].append(name + ':' + str(m.start()))
                state.pep[p]['chops'].append(f"{one['chop'][m.end() - 1]:.5f}")


def get_clusters(state):
    finder(state)

    header = ['Cluster', 'Peptide', 'Score', 'Immunogenicity', 'TAP_score']
    header += state.allele_freq.keys()
    header += ['Locations', 'NetChop_score']
    result = "\t".join(header) + "\n"
    # state.places
    for k, cls in enumerate(state.cls):
        for p in cls:
            # Cluster, Pep, Score, Immuno
            result += f"{k}\t{p}\t{state.pep[p]['score']:.5f}\t{state.pep[p]['imm']}\t"
            # TAP
            if len(p) > 8:
                tap = state.tap[p[-9:]]
            else:
                if p not in state.tap:
                    print("!", p)
                tap = state.tap[p]
            result += f"{tap:.5f}\t"
            # Affinity
            aff = [str(get_affinity(state, p, a)) for a in state.allele_freq]
            result += "\t".join(aff) + "\t"
            # Places
            result += ", ".join(state.pep[p]['places']) + "\t"
            # Chops
            result += ", ".join(state.pep[p]['chops']) + "\n"

    with open(f"{state.wd}/results.tsv", 'w') as res:
        res.write(result)
