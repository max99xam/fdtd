############################################################################################################
## Program: fdtd2d-GaussianBeam
## Author: Zong-Han, Xie
## Date: 2013/08/13
## Description: Vacuum 2d FDTD, modified 2dTM program from Sullivan's textbook.It solves Max Eq. in vacuum.
##          It generates a Gaussian beam from left BC, with all three other conductor BCs.
##          This version add w0 test.
## BCs: All decided by boundary d field.
## sqrt(mu_0/eps_0) is the intrinsic impedance of free space
## Units:E=sqrt(mu_0/eps_0)##e
## Units:D=sqrt(1/eps_0##mu_0)##d
## Units:d=e
## units:E=sqrt(mu_0/eps_0)##d
## units:h=H=B/mu_0
## All rights reserved!
############################################################################################################

import numpy as np
import scipy.constants as const
import sys
import time

InImp = np.sqrt(4.0 * np.pi * 1.0e-7 / const.epsilon_0)  ##the intrinsic impedance of free space


##field data of all electromagnetic field
class EMField:
    def __init__(self, nx, ny):
        self.Nx = nx
        self.Ny = ny
        self.Nx1 = self.Nx + 1
        self.Ny1 = self.Ny + 1
        self.dx = np.zeros((self.Nx, self.Ny + 1))  ## dx field size is (nx) X (ny+1)
        self.dy = np.zeros((self.Nx + 1, self.Ny))  ## dy field size is (nx+1) X (ny)
        self.dz = np.zeros((self.Nx + 1, self.Ny + 1))  ## dz field size is (nx+1) X (ny+1)
        self.hx = np.zeros((self.Nx, self.Ny))  ## hx field size is (nx) X (ny)
        self.hy = np.zeros((self.Nx, self.Ny))  ## hy field size is (nx) X (ny)
        self.hz = np.zeros((self.Nx, self.Ny))  ## hz field size is (nx) X (ny)

    def getNxNy(self):
        return (self.Nx, self.Ny)


##http://en.wikipedia.org/wiki/Gaussian_beam
def GaussianBeam(yLen_in, wc_in, w0_in, tau_in, t0_in, x0_in, amp):
    yLen = yLen_in
    wc = wc_in
    w0 = w0_in
    tau = tau_in
    t0 = t0_in
    x0 = x0_in
    E0 = amp
    y0 = yLen / 2.0
    kay = wc / const.c
    Z0 = w0 * w0 * kay / 2.0

    def inner(nowTime, x, yVec):
        phi = (x - x0) / Z0
        inv_phasefront = (x - x0) / ((x - x0) * (x - x0) + Z0 * Z0)
        ans = np.zeros(len(yVec))
        index = 0
        for y in yVec:
            if (abs((y - y0)) < 2.4 * w0):
                if (abs(t0 - nowTime) < 2.0 * tau):
                    ans[index] = E0 / InImp * np.exp(-1.0 * pow((t0 - nowTime) / tau, 2.0)) * np.cos(
                        kay * ((x - x0) * (x - x0) + (y - y0) * (y - y0)) * inv_phasefront * 0.5 - np.tan(phi) + kay * (
                                    0 - x0) - wc * (t0 - nowTime)) * np.exp(-1.0 * pow((y - y0) / w0, 2.0))
                else:
                    ans[index] = 0.0
            else:
                ans[index] = 0.0
            index = index + 1
        return ans

    return inner


class FDTDUpdater:
    def __init__(self, emField, ddx, ddy, CRatio):
        self.dyFieldSrc = 0  ##constants for specifying source polarization
        self.dzFieldSrc = 1  ##constants for specifying source polarization
        self.em = emField
        self.nowTime = 0.0
        self.deltaX = ddx
        self.deltaY = ddy
        self.src = None
        self.CourantRatio = CRatio
        self.dt = self.CourantRatio / const.c / np.sqrt((1.0 / ddx) * (1.0 / ddx) + (1.0 / ddy) * (1.0 / ddy))
        self.XdtCoefficient = const.c * self.dt / ddx  ## dt coefficient multiplied before x direction difference
        self.YdtCoefficient = const.c * self.dt / ddy  ## dt coefficient multiplied before y direction difference
        temp = self.em.getNxNy()
        self.Nx = temp[0]
        self.Ny = temp[1]
        self.srcDirection = None
        ## This program only launchs wave from lower x-direction boundary
        self.srcYCoord = np.zeros(self.Ny - 1)
        for i in range(int(self.Ny - 1)):
            self.srcYCoord[i] = (i + 1) * self.deltaY

    def update(self):
        em = self.em  ## for cleaner code
        nx = self.Nx  ## cell number along x direction
        ny = self.Ny  ## cell number along y direction

        ## Calculate the dz field, dz field size is (nx+1) X (ny+1)
        em.dz[1:nx, 1:ny] += (self.XdtCoefficient * (em.hy[1:nx, 1:ny] - em.hy[0:nx - 1, 1:ny])
                              - self.YdtCoefficient * (em.hx[1:nx, 1:ny] - em.hx[1:nx, 0:ny - 1]))

        ## Calculate the Hx
        em.hx[0:nx, 0:ny] += self.YdtCoefficient * (em.dz[0:nx, 0:ny] - em.dz[0:nx, 1: ny + 1])  ##norm dz = norm ez

        ## Calculate the Hy
        em.hy[0:nx, 0:ny] += self.XdtCoefficient * (em.dz[1:nx + 1, 0:ny] - em.dz[0:nx, 0:ny])  ##norm dz = norm ez

        ## Calculate the dx ##/
        em.dx[0:nx, 1:ny] += self.YdtCoefficient * (em.hz[0:nx, 1:ny] - em.hz[0:nx, 0:ny - 1])

        ## Calculate the dy ##/
        em.dy[1:nx, 0:ny] += self.XdtCoefficient * (em.hz[0:nx - 1, 0:ny] - em.hz[1:nx, 0:ny])

        ## Calculate the hz ##/
        em.hz[0:nx, 0:ny] += self.YdtCoefficient * (em.dx[0:nx, 1:ny + 1] - em.dx[0:nx, 0:ny]) - self.XdtCoefficient * (
                    em.dy[1:nx + 1, 0:ny] - em.dy[0:nx, 0:ny])  ##norm dz = norm ez

        self.nowTime += self.dt

        if (self.srcDirection == self.dzFieldSrc):
            self.em.dz[0, 1:ny] = self.src(self.nowTime, 0.0, self.srcYCoord)

        if (self.srcDirection == self.dyFieldSrc):
            self.em.dy[0, 1:ny] = self.src(self.nowTime, 0.0, self.srcYCoord)

    def setSourceDir(self, xx):
        self.srcDirection = xx

    def setSrc(self, xx):
        self.src = xx

    def getDt(self):
        return self.dt

    def getNowTime(self):
        return self.nowTime

    def printNowTime(self):
        print("nowTime = " + str(self.nowTime) + " (s)")


##This object is used to write electro magnecitc field data into txt files
class EMFieldOutput:
    def __init__(self, xx):
        self.em = xx
        self.dumpSteps = set()

    def writeFields(self, nowStep, suffix=""):
        if (not (nowStep in self.dumpSteps)):
            return 0
        nx = self.em.getNxNy()[0]
        ny = self.em.getNxNy()[1]
        em = self.em

        ##write simulation result to files.
        f = open("Ey_Step_" + str(nowStep) + suffix + ".txt", 'w')
        f.write(str(nx + 1) + "\n")
        f.write(str(ny) + "\n")
        for j in range(ny):
            for i in range(nx + 1):
                f.write(str(InImp * (em.dy[i, j])) + "\n")
        f.close()

        f = open("Ex_Step_" + str(nowStep) + suffix + ".txt", 'w')
        f.write(str(nx) + "\n")
        f.write(str(ny + 1) + "\n")
        for j in range(ny + 1):
            for i in range(nx):
                f.write(str(InImp * (em.dx[i, j])) + "\n")
        f.close()

        f = open("Ez_Step_" + str(nowStep) + suffix + ".txt", 'w')
        f.write(str(nx + 1) + "\n")
        f.write(str(ny + 1) + "\n")
        for j in range(ny + 1):
            for i in range(nx + 1):
                f.write(str(InImp * (em.dz[i, j])) + "\n")
        f.close()

        f = open("Bx_Step_" + str(nowStep) + suffix + ".txt", 'w')
        f.write(str(nx) + "\n")
        f.write(str(ny) + "\n")
        for j in range(ny):
            for i in range(nx):
                f.write(str(InImp * (em.hx[i, j])) + "\n")
        f.close()

        f = open("By_Step_" + str(nowStep) + suffix + ".txt", 'w')
        f.write(str(nx) + "\n")
        f.write(str(ny) + "\n")
        for j in range(ny):
            for i in range(nx):
                f.write(str(InImp * (em.hy[i, j])) + "\n")
        f.close()

        f = open("Bz_Step_" + str(nowStep) + suffix + ".txt", 'w')
        f.write(str(nx) + "\n")
        f.write(str(ny) + "\n")
        for j in range(ny):
            for i in range(nx):
                f.write(str(InImp * (em.hz[i, j])) + "\n")
        f.close()

    def addDumpStep(self, xx):
        self.dumpSteps.add(xx)


def main():
    startTime = time.time()
    w0 = 4.0e-6  ## transverse width of the gaussian beam at beam waist in meter
    tau = 29.72e-15  ##temporal width of the gaussina beam
    CourantRatio = 0.99  ##for Courant condition
    wc = 3.75e+14 * 3.1415926 * 2.0  ##carier frequency of the gaussian beam, it's wavelength is 810 nm.
    ddx = 810e-9 / 32.0  ##grid resolution, using 32 grids to resolve the carrier wavelength of gaussian beam.
    ddy = 4.0 * ddx  ##transverse grid resolution
    E0 = 1.0  ##amplitude of incomming Gaussian beam
    Nx = int(np.floor((tau * const.c * 6.0) / ddx) + 100)  ##longitudial lenght is slightly larger than 140 micrometer
    Ny = int(np.floor(6.0 * 1.414 * w0 / ddy) + 20)  ##transverse width is slightly larger than 5##sqrt(2)##w0
    x0 = Nx * ddx / 2.0  ##x-direction beamn waist at center of the simulation space.
    t0 = float(2.5 * tau)
    y0 = Ny * ddy / 2.0  ##y-direction beamn waist at center of the simulation space.

    em = EMField(Nx, Ny)  ##construct electromagnetic field

    updater = FDTDUpdater(em, ddx, ddy, CourantRatio)

    src = GaussianBeam(Ny * ddy, wc, w0, tau, t0, x0, E0)
    updater.setSrc(src)
    updater.setSourceDir(updater.dzFieldSrc)

    out = EMFieldOutput(em)  ##field data writer

    numberOfSteps = int(np.floor(Nx * ddx / const.c / updater.getDt()))
    print("numberOfSteps= " + str(numberOfSteps) + "n")
    numberOfSteps1 = numberOfSteps + 1

    dumpNumber = 5  ##intended number of dumps
    outputInterval = int(np.floor(numberOfSteps / dumpNumber))
    for i in range(dumpNumber + 1):  ## try to include the last step
        out.addDumpStep(i * outputInterval)

    out.addDumpStep(numberOfSteps - 1)  ## try to include the last step

    ##main loop of the simulations
    ##main loop start!!!!!!
    for i in range(numberOfSteps):
        updater.update()
        updater.printNowTime()
        out.writeFields(i)
    ##main loop end!!!!!!!
    endTime = time.time()
    print("Execution time = " + str(endTime - startTime) + " (s)")

if __name__ == "__main__":
    main()