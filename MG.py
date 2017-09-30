#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 30 11:22:51 2017

@author: AndersonLab
"""
import pandas as pd
import pypower as pp
import pyomo.core
import numpy as np
from pypower.api import case30,makeYbus
from numpy.linalg import inv
from pyomo.environ import *

#%% System Parameters
nt = 24;
lmp = 5.3*np.ones(nt)
cg = np.array([[3.3, 3.3, 3.3]]) 
cd = 1.01*np.array([[3.3, 3.3, 3.3, 3.3, 3.3, 3.3]]) # slightly larger than cg to encourage cost
cb = 0.1 # needs to be small so that MG stores energy rather than export energy
ci = lmp # change
ce = 0.8*lmp #change

ndl = 6 # number of dispatchable load
ng = 3 # number of generators
nb = 1 # number of storage
model = ConcreteModel()
model.nb = RangeSet(0,nb-1)
model.ndl = RangeSet(0,ndl-1)
model.nt = RangeSet(0,nt-1)
model.ng = RangeSet(0,ng-1)
#%% Define Optimization Variables and Their Bounds
# 1. Dispatchable loads

one = np.ones(nt)
Pdmin = 0.1*np.array([0.5*one,4*one,2*one,5.5*one,1*one,7*one])
Pdmax = 0.03*np.array([10*one,16*one,15*one,20*one,27*one,32*one])
def b1(model,i,j):
    return (Pdmin[i,j],Pdmax[i,j])
model.pd=Var(model.ndl, model.nt, bounds = b1)

# 2. Non-controllable load    
NCL = 0.03 * np.array([[120.6, 115.8, 114.8, 112.6, 114.0, 113.4, 
                  117.1, 126.3, 130.7, 132.5, 135.6, 134.8, 
                  136.5, 137.7, 137.1, 138.0, 136.3, 133.3, 
                  131.7, 129.3, 128.2, 127.4, 125.6, 124.2]]) # sample load profile
NCL = NCL[:,0:nt];

# 3. Generators
Gmax = 0.3 * np.array([5*one,4.5*one,7*one])
Gmin = np.array([1*one, 0.8*one, 1.5*one])
def b2(model,i,j):
    return (Gmin[i,j],Gmax[i,j])
model.pg = Var(model.ng, model.nt, bounds = b2)
#Rgsmin = 0.3*Gmax
#Rup = 0.3*Gmax
#def b8(model,i,j):
#    return (Rdn[i,j],Rup[i,j])
model.rgs = Var(model.ng, model.nt)
model.z = Var(model.ng, model.nt)

# 4. Storage
# power to and from storage
Pbmin = -3*np.ones((nb,nt))
Pbmax = 3*np.ones((nb,nt))
def b3(model,i,j):
    return (Pbmin[i,j],Pbmax[i,j])
model.pb = Var(model.nb, model.nt, bounds=b3)

# state of charge
Bmin = np.zeros((nb,nt))
Bmax = 10*np.ones((nb,nt))
def b4(model,i,j):
    return (Bmin[i,j],Bmax[i,j])
model.b = Var(model.nb, model.nt, bounds=b4)
model.b[0,0] = 5 # initialize battery energy state

# 4. Import & Export
def b5(model,i):
    return (0,1000)
model.ex = Var(model.nt, bounds=b5)
def b6(model,i):
    return (0,1000)
model.im = Var(model.nt, bounds=b6)
Nmin = -20*one
Nmax = 20*one
def b7(model,i):
    return (Nmin[i],Nmax[i])
model.net = Var(model.nt, bounds=b7)
 
#%% Constraints
# DG Constraints
Rdn = 0.3*Gmax
Rup = 0.3*Gmax
model.cons = ConstraintList()
for g in xrange(ng):
    for h in xrange(nt-1):
        model.cons.add(-Rdn[g,h+1] <= model.pg[g,h+1]-model.pg[g,h] <= Rup[g,h+1]) # ramping constraints
        model.cons.add((model.pg[g,h+1]+model.rgs[g,h+1])-(model.pg[g,h]-model.rgs[g,h]) <= float(Rup[g,h+1])) # Up ramping Constraints with Rgs
        model.cons.add(float(-Rdn[g,h+1]) <= (model.pg[g,h+1]-model.rgs[g,h+1])-(model.pg[g,h]+model.rgs[g,h])) # Dn Ramping Constraints with Rgs

for g in xrange(ng):
    for h in xrange(nt):
        model.cons.add(Gmin[g,h] <= model.pg[g,h] + model.rgs[g,h])#Generator Bounds with Rgs
        model.cons.add(model.pg[g,h] + model.rgs[g,h] <= Gmax[g,h])
        model.cons.add(model.rgs[g,h] <= model.z[g,h])
        model.cons.add(-model.rgs[g,h] <= model.z[g,h])

for h in xrange(nt):
    model.cons.add(sum(model.rgs[r,h] for r in xrange(ng)) == sum(0.00*NCL[r,h] for r in xrange(1)))
    model.cons.add(sum(model.pg[r,h] for r in xrange(ng)) == sum(NCL[r,h] + model.pb[r,h] for r in xrange(1))+ model.net[h] + sum(model.pd[r,h] for r in xrange(ndl))) # Power balance equation forecast
#    model.cons.add(model.ex[h] == max(model.net[h],0))
#    model.cons.add(model.im[h] == max(-model.net[h],0))

for i in xrange(nb):
    for j in xrange(nt-1):
        model.cons.add(model.b[i,j+1] == model.b[i,j]+model.pb[i,j])
for i in xrange(nb):
    for j in xrange(nt):
        model.cons.add(-model.pb[i,j] <= model.b[i,j])

#%% Objective
#%% Objective
#def exps(model):
a = sum(cb * model.b[i,j] for i in xrange(nb) for j in xrange(nt))
b = sum(model.pd[i,j]*cd[0,i] for i in xrange(ndl) for j in xrange(nt))
c = sum(model.pg[i,j]*cg[0,i] for i in xrange(ng) for j in xrange(nt))
ex = sum(model.ex[i]*ce[i] for i in xrange(nt))
im = sum(model.im[i]*ci[i] for i in xrange(nt))
#    return a+b+c+d
#model.obj = Objective(rule = exps)
model.obj = Objective(expr = a+b+c+ex+im)

#%  Objective=sum(pg')*cg-sum(pd')*cd+sum(b')*cb+sum(im)*ci-sum(ex)*ce+sum(abs(rls))*crls+sum(abs(rgs'))*crgs+sum(abs(rgd'))*crgd;
#Objective=sum(pg')*cg-sum(pd')*cd   +sum(b')*cb  +im*ci'-ex*ce';
#      
opt = SolverFactory("gurobi")
results = opt.solve(model)
model.pg.display()
model.rgs.display()
model.ex.display()
model.pb.display()
model.b.display()
model.pd.display()
model.im.display()
model.obj.display()        