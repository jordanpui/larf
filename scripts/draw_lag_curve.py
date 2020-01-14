#!/usr/bin/env python3

from numpy import *
import math
import matplotlib.pyplot as plt
import sys

class Benchmark:
    def __init__(self, full_name, abbr_name):
        self.full_name = full_name
        self.abbr_name = abbr_name
    def __repr__(self):
        return self.full_name

class Benchmarks:
    def __init__(self):
        self.__bms = []

    def add(self, full_name_pat, abbr_name_pat, ids):
        for id in ids:
            self.__bms.append(Benchmark(full_name_pat.format(id), abbr_name_pat.format(id)))

    def get_choices(self):
        choices = []
        for bm in self.__bms:
            choices.append(bm.full_name)
            choices.append(bm.abbr_name)
        choices.append('all')
        return choices

    def get_selected(self, names):
        if 'all' in names:
            return self.__bms
        else:
            selected = []
            for name in names:
                for bm in self.__bms:
                    if name == bm.abbr_name or name == bm.full_name:
                        selected.append(bm)
                        break
            return selected

all_benchmarks = Benchmarks()
all_benchmarks.add('FPGA{}', 'f{}', ['01', '02', '03', '04','05','06','07','08','09','10','11','12'])

bms = all_benchmarks.get_selected([sys.argv[1]])

for bm in bms:
    print(bm)
    curve_file=open(bm.abbr_name+'/'+bm.abbr_name+'.curve', 'r')
    primal=[]
    dual=[]

    curve_file.readline()

    for line in curve_file:
        primal.append(float(line.split()[0]))
        dual.append(float(line.split()[1]))

    iter = len(primal)

    plt.xlabel('#Iteration')
    plt.ylabel('Primal/Dual Value')
    plt.plot(range(0,iter), primal, 'r')
    plt.plot(range(0,iter), dual, 'b')
    plt.savefig(bm.abbr_name + '/' + bm.abbr_name + '_' + str(len(primal)) + '.pdf')
    plt.clf() 
    # plt.show()
