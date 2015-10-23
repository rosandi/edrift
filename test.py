#!/usr/bin/python

import numpy as nm
import matplotlib.pyplot as plt
from time import sleep
from drift import drift

dim=1000

symlist=("x ++---- "
         "s +----- ")

mapstr=""

for i in range(0,dim-1):
    mapstr+="x"
mapstr+="s"

print len(mapstr)
print mapstr

dim=1000
dt=1e-5
dl=10

drf=drift(dim,1,1,symlist,mapstr)
drf.timestep(dt)
drf.spacestep(dl)
drf.algorithm("basic")
source=drf.source()
con=drf.con()
cap=drf.cap()
cell=drf.cell()

for i in range(0,dim):
	source[i]=nm.cos(i*dl*nm.pi/(dim*dl)/2.0)
	con[i]=0.2497
	cap[i]=1.8088e-7

#@MAC numbers having more than 4 digit is treated as variable!
drf.link_equhead(cell, symbol='s', equ='RCELL(c,0,0)-(RCELL(0,0,0)-30)*DL*0.001;c=-1')

drf.step()

for i in range(0,dim):
	source[i]=0.0

drf.dump("dump.out")
plt.ion()
for i in range(0,10):
    drf.step(2500)
    A=drf.collect(drf.cell())
    plt.clf()
    plt.plot(A)
    plt.draw()

#    sleep(1)
