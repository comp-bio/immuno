import sys
import os
import time
import multiprocessing
from . import help, help_sim


class Params:
    def __init__(self):
        self.src = None
        self.freq = None
        self.allele = None
        self.organism = 'human'
        self.dist = 0.15
        self.model = 'consensus'
        self.c = 0.15
        self.n = 20
        self.dir = time.strftime("__OUT_%Y%m%d_%H%M%S")
        self.t = multiprocessing.cpu_count()
        self._parse_args()
        self._validate()

    def _parse_args(self):
        args = sys.argv[1:]
        if len(args) == 0:
            return help()
        for i in range(len(args)):
            if args[i].startswith("-") and i + 1 < len(args):
                key = args[i][1:]
                val = args[i + 1]
                if hasattr(self, key):
                    setattr(self, key, val)

    def _validate(self):
        self.c = min(max(0, float(self.c)), 1)
        self.n = min(max(1, int(self.n)), 1000)
        self.t = min(int(self.t), multiprocessing.cpu_count())
        self.dist = min(max(0, float(self.dist)), 1)

        if self.organism not in ['human', 'mouse']:
            self.organism = 'human'


class ParamsSim:
    def __init__(self):
        self.src = None
        self.results = None
        self._parse_args()
        self._validate()

    def _parse_args(self):
        args = sys.argv[1:]
        if len(args) == 0:
            return help_sim()
        for i in range(len(args)):
            if args[i].startswith("-") and i + 1 < len(args):
                key = args[i][1:]
                val = args[i + 1]
                if hasattr(self, key):
                    setattr(self, key, val)

    def _validate(self):
        if not self.src or not os.path.isfile(self.src):
            return help_sim("Input file not found!")
        if not self.results or not os.path.isfile(self.results):
            return help_sim("Results file not found!")
        self.out = os.path.dirname(
            os.path.realpath(self.results)) + "/variability.tsv"
