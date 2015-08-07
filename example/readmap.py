#!/usr/bin/python

import numpy as nm
import matplotlib.pyplot as plt
from time import sleep
from drift import drift

dt=1e-5
dl=10

drf=drift(mapfile="1D.map")

nx=drf.dimension()
dim=nx[0]*nx[1]*nx[2]

drf.timestep(dt)
drf.spacestep(dl)
drf.algorithm("basic")
source=drf.source()
con=drf.con()
cap=drf.cap()

for i in range(0,dim):
	source[i]=nm.cos(i*dl*nm.pi/(dim*dl)/2.0)
	con[i]=0.2497
	cap[i]=1.8088e-7

drf.setcell(50.0,cells="s")

drf.step()

for i in range(0,dim):
	source[i]=0.0

plt.ion()
for i in range(0,100):
    drf.step(250)
    A=drf.collect(drf.cell())
    plt.clf()
    plt.plot(A)
    plt.draw()
#    sleep(1)

plt.show()
