import os
import sys
import time
from . import error, help, log, echo, blosum


class State:
    def __init__(self, params):
        self.params = params
        self.model = None
        self.pep = {}
        self.mhc = {}
        self.allele_freq = {}
        self.time = time.time()
        self.wd = ''
        self.dt = time.time()
        self.pwd = os.path.dirname(os.path.realpath(sys.argv[0]))
        self.blosum62 = blosum.Blosum(f'{self.pwd}/meta/blosum62.txt')

        self._load_model()
        self._load_allele_frequencies()
        self._prepare_output_dirs()

    def _load_model(self):
        def float_it(val):
            return '-' if val == '-' else float(val)

        models = {
            'netmhcpan_el': lambda x: [x[0], x[5], float(x[9]), {
                'score': float(x[8])
            }],
            'netmhcpan_ba': lambda x: [x[0], x[5], float(x[9]), {
                'ic50': float(x[8])
            }],
            'consensus': lambda x: [x[0], x[5], float(x[6]), {
                'ann_ic50': float_it(x[7]),
                'smm_ic50': float_it(x[9])
            }],
        }

        if self.params.model not in models:
            self.params.model = 'consensus'

        self.model = models[self.params.model]
        mhc = os.path.join(self.pwd, 'meta', f'{self.params.model}.tsv')
        with open(mhc) as f:
            for line in f.readlines()[2:]:
                org, matrix, length = line.strip().split('\t')
                matrix = matrix.strip()
                self.mhc.setdefault(matrix, []).append(int(length))

    def _add_allele(self, name, weight=1):
        try:
            if name not in self.mhc:
                error(f"Allele {name} not found in model {self.params.model}")
                return False
            self.allele_freq[name] = float(weight)
            return True
        except ValueError:
            error("Error reading .tsv file")
            return False

    def _load_allele_frequencies(self):
        freq_path = self.params.freq
        if freq_path and os.path.isfile(freq_path):
            with open(freq_path) as f:
                lines, added = (0, 0)
                for line in f:
                    parts = line.replace(" ", "\t").split("\t")
                    lines += 1
                    if len(parts) >= 2:
                        if self._add_allele(parts[0], parts[1]):
                            added += 1
            if lines > added:
                pass
                # error(f"Allele {name} not found in model {self.params.model}")

        if self.params.allele:
            self._add_allele(self.params.allele)

        if not self.allele_freq:
            # raise ValueError()
            err = [
                "Frequencies not found!",
                "  Use -freq <file> to load .tsv frequency file",
                "  or use -allele [allele name] to specify a single allele"]
            help("\n".join(err))

        total = sum(self.allele_freq.values())
        self.allele_freq = {k: v / total for k, v in self.allele_freq.items()}

    def _prepare_output_dirs(self):
        self.wd = os.path.join(self.params.dir, self.params.model)
        os.makedirs(os.path.join(self.wd, "tmp"), exist_ok=True)

    def done(self):
        tx = time.strftime("%H:%M:%S", time.gmtime(time.time() - self.dt))
        echo("Done" + (" " * 20) + "\n", 36)
        log(f"Execution time: {tx} (H:M:S)")
        log(f"Results:        {self.wd}/results.tsv")
