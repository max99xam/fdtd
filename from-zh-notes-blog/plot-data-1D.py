#!/usr/bin/env python
##
## file fdtd1d-GaussianBeamPlot.py
## brief plot fdtd1d-Gaussian results
## version 2013/08/09
## Author Zong-Han, Xie <icbm0926@gmail.com>
## All rights reserved

import matplotlib.pyplot as plt
import sys
import numpy as np
import time

def countLine(fileName):
    f = open(fileName,'r')
    line = f.readline()
    numberOfLines = 0
    if len(line) < 1:
        f.close()
        return 0
    while len(line) > 0:
        numberOfLines = numberOfLines + 1
        line = f.readline(100)
    f.close()
    return numberOfLines

def readData(fileName):
    numberOfLines = countLine(fileName)
    data=np.zeros((2,numberOfLines))
    f = open(fileName, 'r')
    line = f.readline(100)
    index = 0
    while len(line) > 0:
        g=line.split(',')
        data[0][index] = float(g[0])
        data[1][index] = float(g[1])
        line = f.readline(100)
        index = index +1
    f.close()
    return data

def main():
    if len(sys.argv) < 2:
        return
    mmain(sys.argv[1])

def mmain(file):
    print(file)
    if len(file) < 1:
        return
    data=readData(file)
    plt.plot(data[0], data[1])
    plt.show()

if __name__ == "__main__":
    main()