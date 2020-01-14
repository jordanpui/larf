#!/usr/bin/env python3

from numpy import *
import math
import matplotlib.pyplot as plt
import sys

curve_file=open('usageError.txt', 'r')
act=[]
org=[]
ratio=[]

numVal = 1600

for line in curve_file:
    v1 = float(line.split()[0])
    v2 = float(line.split()[1])
    act.append(float(v1))
    org.append(float(v2))
    ratio.append((v1-v2)/v2)

plt.xlabel('TDM ratio $v$')
plt.ylabel('Wire usage error rate $u_{err}$')
plt.plot(range(1,numVal+1), ratio, 'r')
plt.show()