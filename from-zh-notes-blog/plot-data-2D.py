#!/usr/bin/env python
##
## file fdtd2d-GaussianBeamPlot.py
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
    ''' read Nx, Ny first'''
    f = open(fileName, 'r')
    Nx = int(f.readline(100))
    Ny = int(f.readline(100))
    data=np.zeros((Nx,Ny))

    line = f.readline(100)
    index = 0
    for j in range(Ny):
        for i in range(Nx):
            data[i][j] = float(line)
            line = f.readline(100)
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
    adjustPlotRatio = 0.5 ##to adjust plot size
    plot = plt.imshow(data, aspect='equal', extent=(0, 1000, 0, int(1000/adjustPlotRatio)))
    plot.set_cmap("jet")
    plt.colorbar()
    plt.show()

if __name__ == "__main__":
    main()