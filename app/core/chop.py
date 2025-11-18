import subprocess
import sys
from collections import defaultdict
from core import exec
from config import error, echo


def available_points(stdout):
    lines = stdout.split('\n')
    if len(lines) < 3:
        return False, []
    pep_index = int(lines[0][1:].split(':')[0])
    prob = [line.split('\t')[2] for line in lines[2:-1]]
    return (pep_index, prob)


def tap_values(stdout):
    lines = stdout.split('\n')
    if len(lines) < 3:
        return False
    for line in lines[2:-1]:
        d = line.split('\t')
        pep, score = (d[1], float(d[5]))
        yield (pep, score)
        yield (pep[-8:], score)


def netchop(state, params):
    EXEC = f"{sys.executable} {state.pwd}/netchop/predict.py"

    files, runs = ([], [])
    for i, obj in enumerate(state.seq):
        tmp = f"/tmp/prt_{i}.fa"
        with open(tmp, 'w+') as f:
            seq = obj['seq'].replace('*', 'X')
            f.write(f">P{i}\n{seq}\n")
        files.append(tmp)
        runs.append(f"{EXEC} -n -m netchop {tmp}")
        runs.append(f"{EXEC} -n -m netctl {tmp}")

    state.tap = {}
    results = exec.pool(runs, params.t)
    for cmd, code, stdout, stderr in results:
        if code != 0:
            error("NetChop ERROR:")
            error(stderr)
        else:
            if '-m netchop' in cmd:
                i, ptx = available_points(stdout)
                state.seq[i]['chop'] = list(map(float, ptx))
            if '-m netctl' in cmd:
                for pep, score in tap_values(stdout):
                    state.tap[pep] = score

    subprocess.run(['rm', '-f'] + list(files))
