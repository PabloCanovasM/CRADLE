# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 15:34:48 2016

@author: leendert
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

def DalitzPlot(name1, name2, name3):
    range1 = [i for i in range(len(names)) if names[i] == name1]
    range2 = [i for i in range(len(names)) if names[i] == name2]
    range3 = [i for i in range(len(names)) if names[i] == name3]
    
    r = set(range1) ^ set(range2)
    r = set(r) ^ set(range3)

    data1 = data[list(range1),3:]
    data2 = data[list(range2),3:]
    data3 = data[list(range3),3:]
    
    s12 = data1 + data2
    s13 = data1 + data3
    
    s12 = [i[0]**2.0-sum(i[1:]**2.0) for i in s12]
    s13 = [i[0]**2.0-sum(i[1:]**2.0) for i in s13]
    
    plt.figure(figsize=(9, 8))
    plt.hist2d(s12, s13, bins=200)
    plt.title("Dalitz Plot", fontsize=18)
    plt.xlabel(r"$m_{\nu recoil}$", fontsize=18)
    plt.ylabel(r"$m_{\nu e}$", fontsize=18)
    plt.savefig(name2 + '_DalizPlot.pdf')


fileName = "/Users/victor/Documents/git/CRADLE (copie)/build/output.txt"


data = np.genfromtxt(fileName, dtype=None)
names = [str(i[2], 'utf-8') for i in data]
data = np.array([list(i)[:2] + list(i)[3:] for i in data])
data[np.isnan(data)]=0.


NAMES = list(set(names))
for i in NAMES:    
    if (str(i) == 'gamma') or ( str(i) == 'gammaRC'):
        NAMES.remove(str(i))
NAMES.sort()
print(NAMES)

DalitzPlot(NAMES[2], NAMES[0], NAMES[1])

plt.show(block=True)
