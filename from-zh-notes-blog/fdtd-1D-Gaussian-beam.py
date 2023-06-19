"""
    Purpose: This program is to demostrate one dimensional
            finite difference time domain electromagnetic simulation.
             This program os simply demostrate one dimensional electromagnetic
             wave propagates along z direction.
             This code is an implementation of examples in Sullivan's FDTD textbook.
"""
sd
import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import sys
import time


#output field data to a text file
class FieldOutput:
    def __init__(self, fileName, xpt, v):
        self.prefixName = fileName
        self.xPoints = xpt
        self.field = v
        self.dumpSteps = set()

    def addDumpStep(self, step):
        self.dumpSteps.add(step)

    #perform output actions
    def dump(self, nowStep):
        if not (nowStep in self.dumpSteps):
            return
        nm = self.prefixName + "_" + str(nowStep) + ".txt"
        f = open(nm, 'w')
        for i in range(len(self.xPoints)):
            f.write(str(self.xPoints[i]) + "," + str(self.field[i]) + "\n")
        f.close()


#update electric field and magnetic field according to FDTD algorithm
class FDTDUpdater:
    def __init__(self, exField, hyField, dtt, dtCoeff):
        self.ex = exField
        self.hy = hyField
        self.nowTime = 0.0
        self.dtCoeffecient = dtCoeff
        self.dt = dtt
        self.numberOfGrids = len(exField)
        self.src = None

    def update(self):
        if (self.src is None):
            print("source of the FDTD simulation is not given!n")
            sys.exit(-1)

        self.nowTime = self.nowTime + self.dt
        ##update electric field
        self.ex[1:len(self.ex) - 1] = self.ex[1:len(self.ex) - 1] + self.dtCoeffecient * (
                self.hy[0:len(self.hy) - 1] - self.hy[1:len(self.hy)])

        ##update source
        self.ex[0] = self.src(self.nowTime)

        ##update magnetic field
        self.hy[0:len(self.hy)] = self.hy[0:len(self.hy)] + self.dtCoeffecient * (
                self.ex[0:len(self.ex) - 1] - self.ex[1:len(self.ex)])

    def setSource(self, inp):
        self.src = inp

    def getNowTime(self):
        return self.nowTime

    def printNowTime(self):
        print("nowTime = " + str(self.nowTime) + "n")


##generate a time sequential 1D Gaussian beam
def GaussianBeam(e0, t0In, tau_in, wc_in):
    E0 = e0  ##amplitude
    t0 = t0In  ##time when center appeared
    tau = tau_in  ##Gaussian beam half-width
    wc = wc_in  ##radial frequency

    def inner(T):
        ##update source
        if (abs(t0 - T) < 2.0 * tau):
            return E0 * np.exp(-1.0 * np.power((t0 - T) / tau, 2.0)) * np.cos(wc * (t0 - T))
        return 0

    return inner


def main():
    startTime = time.time()
    #Define used parameters
    ##nitialize simulation and laser parameters
    wavelength = 810e-9  ## laser wavelength
    tau = 29.72e-15  ##half temporal width of simulated electromagnetic gaussian beam in second.

    ##radial frequency of the simulated electromagnetic wave, corresponding wavelength is 810 nm.
    wc = const.c / wavelength * const.pi * 2.0
    ##define simulation resolution parameters
    dtCoeffecient = 0.99  ##for Courant condition
    dz = wavelength / 32.0  ##spatial grid interval in meter
    dt = dz * dtCoeffecient / const.c  ##time step interval used in the sumulation
    numberOfGrids = int(
        np.floor(4.0 * tau * const.c / dz)) + 200  ##grid number is 4 times half temporal width and extra 200 cell.
    t0 = 2.0 * tau + 100.0 * dt  ##The time when pulse peak appeared in the simulation domain
    NSTEPS = int(np.floor(4.0 * tau / dt))  ##Maximum number of steps
    E0 = 1.0  ##maximum amplitude of the electric field in V/m

    ##declare and initialize field data
    ex = np.zeros(numberOfGrids)  ##electric field data
    hy = np.zeros(numberOfGrids - 1)  ##magnetic field data
    updater = FDTDUpdater(ex, hy, dt, dtCoeffecient)
    src = GaussianBeam(E0, t0, tau, wc)
    updater.setSource(src)

    ##setup dump related setting
    xGridPt = np.zeros(numberOfGrids)  ##grid coordinates
    xCellPt = np.zeros(numberOfGrids - 1)  ##grid coordinates
    ##Initialize x grid points data
    for k in range(numberOfGrids):
        xGridPt[k] = k * dz
    ##Initialize x cell points data
    for k in range(numberOfGrids - 1):
        xCellPt[k] = (k + 0.5) * dz

    exOut = FieldOutput("gaussian_Ex", xGridPt, ex)  ##electric field is on grid point.
    hyOut = FieldOutput("gaussian_Hy", xCellPt, hy)  ##magnetic field is at half cell.
    ##setup dump time
    dumpNumber = 5.0  ##total number of data dump
    interval = int(np.floor(NSTEPS / dumpNumber))  ##get interval between data dump according to dumpNUmber
    for i in range(int(np.floor(dumpNumber))):
        exOut.addDumpStep(interval * (i + 1))
        hyOut.addDumpStep(interval * (i + 1))

    ##main loop
    for n in range(NSTEPS):
        updater.update()
        ##updater.printNowTime()
        exOut.dump(n)
        hyOut.dump(n)

    endTime = time.time()
    print("Execution time = " + str(endTime - startTime) + " (s)")
    return 0


if __name__ == "__main__":
    main()
