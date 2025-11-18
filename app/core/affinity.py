import subprocess
import sys
from collections import defaultdict
from core import exec
from config import error, echo


def score(x, y):
    return (x + y)/2 - abs(x - y)/(2 * 1.41421)


def get_affinity(state, p, allele):
    if allele in state.pep[p]['aff']:
        return state.pep[p]['aff'][allele][0]/100
    return 0


def calc_score(state):
    for p in state.pep:
        state.pep[p]['score'] = 0
        for allele in state.allele_freq:
            sc = score(get_affinity(state, p, allele), state.pep[p]['imm'])
            state.pep[p]['score'] += state.allele_freq[allele] * sc


def parse_mhci(state, results):
    for cmd, code, stdout, stderr in results:
        if code != 0:
            continue
        for line in stdout.split('\n')[1:]:
            cols = line.split('\t')
            if len(cols) <= 1:
                continue
            l = state.model(cols)
            all, pep, data = (l[0], l[1], l[2:])
            if pep not in state.pep:
                state.pep[pep] = {'imm': state.imm[pep], 'aff': {}}
            if all not in state.pep[pep]['aff']:
                state.pep[pep]['aff'][all] = data


def predict(state, params, page, iteration):
    MHC = f"{sys.executable} {state.pwd}/mhc_i/src/predict_binding.py"

    runs = []
    sizes = defaultdict(list)
    for p in page:
        sizes[len(p)].append(p)

    sizes = dict(sizes)
    depends = {}
    created = set()
    for s in sizes:
        src = f"{state.wd}/tmp/s{s}-{page[0]}.fa"
        with open(src, 'w+') as fa:
            created.add(src)
            fasta = [f">P{s}_{k}\n{p}\n" for k, p in enumerate(sizes[s])]
            fa.write("".join(fasta))
        for a in state.allele_freq:
            if s not in state.mhc[a]:
                error(f"Allele not found: {a}, size: {s}")
                continue
            cmd = f"{MHC} {params.model} \"{a}\" {s} {src}"
            depends[cmd] = src
            runs.append(cmd)

    results = exec.pool(runs, params.t, iteration)
    for cmd, code, stdout, stderr in results:
        if code != 0:
            created.remove(depends[cmd])
    subprocess.run(['rm', '-f'] + list(created))

    echo(f"➜ Parse MHCI results...\r", 36)
    parse_mhci(state, results)

    echo(f"➜ Score caclulation... \r", 36)
    calc_score(state)
