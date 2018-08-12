import engine
import general
import gpb
from math import *
import time as tim
import matplotlib.pyplot as plt


def main():

    N = 50

    exhaust_port = gpb.Port(25, 0.75, width=25, height=25, fromtdc=34, R=12.5)
    transfer_a = gpb.Port(10, 0.65, width=10, height=15, fromtdc=44)
    intake = gpb.Port(20, 0.75, width=30, height=25, fromtdc=50)
    cyl = gpb.Plenum(601.43, 293, 101325, 580, 287, 1.4, 293, 1.5*101325, ports={'exhaust':exhaust_port, 'transfer':transfer_a, 'intake':intake}, deck=0)
    crankcase = gpb.Plenum(rho0=1.43, T0=293, p0=101325, a0=340, R=287, G=1.4, T=293, p=101325, ports={'transfer':transfer_a})
    atmos = gpb.Plenum(1.43, 293, 101325, 343, 287, 1.4, 293, 101324)
    pipe = gpb.Duct(d1=25, d2=25, rho0=1.43, T0=293, p0=101325, a0=343.1, R=287, G=1.4, L=5901, n=N, inHandle=cyl, outHandle=atmos)
    #intake = gpb.Duct(d1=30, d2=28, rho0=1.43, T0=293, p0=101325, a0=340, R=287, G=1.4, L=120, n=12, inHandle=atmos, outHandle=cyl)

    e = engine.Engine(bore=60, stroke=60, conrod=110, trapped_stroke=60, compression_ratio=8.1, n=1, RPM=8000,
                      intake=None, exhaust=pipe, cylinder=cyl, crankcase=crankcase, atmos=atmos, debug = True, piston_intake_h=50)


    p1, p2, p3, p4, p5 = [], [], [], [], []
    times = []
    time = 0
    angles = []

    idxx = int(3691/(5901/N))

    idxx = 20

    for i in range(200):
        e.step()
        #print("Crankshaft angle: {:.2f}, i: {}".format(e.angle, i))
        #print(pipe.meshes[0].XL1)
        #pressures.append(pipe.meshes[0].XR1**general.G7(1.4))
        p1.append(pipe.meshes[idxx].Xj ** general.G7(1.4))
        p2.append(pipe.meshes[idxx].XR ** general.G7(1.4))
        p3.append(pipe.meshes[idxx].XR1 ** general.G7(1.4))
        p4.append(pipe.meshes[idxx].XL ** general.G7(1.4))
        p5.append(pipe.meshes[idxx].XL1 ** general.G7(1.4))
        time += pipe.dt
        times.append(time)
        angles.append(e.angle)
        #print(time, e.angle)
        idx = -2
        #print("Xj: {:.2f}, XR: {:.2f}, XL: {:.2f}, XR1: {:.2f}, XL1: {:.2f}".format(pipe.meshes[idx].Xj, pipe.meshes[idx].XR, pipe.meshes[idx].XL, pipe.meshes[idx].XR1, pipe.meshes[idx].XL1))

    plt.plot(times, p1, 'r')
    plt.plot(times, p2, 'b', alpha=0.7)
    plt.plot(times, p3, 'g', alpha=0.7)
    plt.plot(times, p4, 'c', alpha=0.7)
    plt.plot(times, p5, 'm', alpha=0.7)
    plt.show()


main()