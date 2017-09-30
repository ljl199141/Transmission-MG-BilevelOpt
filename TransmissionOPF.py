# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import pypower as pp
import pyomo.core
import numpy as np
from pypower.api import case30

from numpy.linalg import inv

ppc = case30()

nt = 2 # time period
gennum = ppc["gencost"].shape[0]
busnum = ppc["bus"].shape[0]
linenum = ppc["branch"].shape[0]
 
# define cost
cg1 = ppc["gencost"][:, [5]].T # % 1st order generator cost
cg2 = ppc["gencost"][:, [4]].T # 2nd order generator cost
crgs = array([[6, 6.75, 7, 5.25, 5, 5]])
conoff = 2*ones((1,gennum)) # generator on cost

sampleL = array([[120.6, 115.8, 114.8, 112.6, 114.0, 113.4, 
                  117.1, 126.3, 130.7, 132.5, 135.6, 134.8, 
                  136.5, 137.7, 137.1, 138.0, 136.3, 133.3, 
                  131.7, 129.3, 128.2, 127.4, 125.6, 124.2]]) # sample load profile
sf = sampleL[0,:nt]/np.mean(sampleL[0,:nt])
sf = np.matlib.repmat(sf,busnum,1)
loads = np.matlib.repmat(ppc["bus"][:, [2]],1,nt)
loads = loads*sf

#     loads(8,:)=loads(8,:)+5*ones(1,nt);%add 5mw of load at bus 8 to create congestion
#     load ('newnetCavg.mat');
#     mg=net; % microgrid contribution, +:export, -:import
#     mbus=8; % bus for microgrid
#     wbus=25; 
#     wl=0; % wind penetration scaling factor


#%% Create B Matrix
#    make Y/B bus matrix and A matrix in Pline=A*Pinj
Ybus, Yf, Yt = makeYbus(100,ppc["bus"],ppc["branch"])
from pyomo.environ import *
B = (1j*Ybus).real.todense()
NB = -B
Bred = B[0:busnum-1,0:busnum-1] # Reduced B Matrix
#    
frm = ppc["branch"][:,0]
to = ppc["branch"][:,1]
M = np.zeros((linenum,busnum)) # from/to matrix
lineC = np.zeros((linenum,1))

for i in xrange(linenum):
    M[i,int(frm[i]-1)] = 1
    M[i,int(to[i]-1)]= -1
    lineC[i]=1*ppc["branch"][i,6]
    
m = M[:,0:busnum-1]    
x = ppc["branch"][:,[3]] # resistance
r = ppc["branch"][:,[2]] # reactance
b = x/(x*x+r*r) # susceptance
D = np.diag(b[:,0]) # diagnol matrix of the susceptance


#%% Define Optimization Variables and Their Bounds

# 3. Generators
one = ones(nt)
Gmax = array([80*one,80*one,40*one,50*one,30*one,55*one])
Gmin = array(np.zeros((6,nt)))

model = ConcreteModel()
model.bn = RangeSet(0,busnum-1)
model.bn1 = RangeSet(0,busnum-2)
model.nt = RangeSet(0,nt-1)
model.gn = RangeSet(0,gennum-1)
model.ln = RangeSet(0,linenum-1)

def b1(model,i,j):
    return (Gmin[i,j],Gmax[i,j])
model.pg = Var(model.gn, model.nt, bounds = b1)
model.onoff = Var(model.gn, model.nt, within = Binary)

Rgsmin = -0.5*Gmax;
Rgsmax = 0.5*Gmax;
def b2(model,i,j):
    return (Rgsmin[i,j],Rgsmax[i,j])
model.rgs = Var(model.gn, model.nt, bounds = b2)
model.z = Var(model.gn, model.nt) # slack variable for rgs

#%% Power Flow Constraints
model.Pinj = Var(model.bn, model.nt) #% bus nodal matrix with forecast wind
model.cons = ConstraintList()
for i in xrange(busnum):
    for j in xrange(nt):
        if i==0:
            model.cons.add(model.Pinj[i,j]==-loads[i,j]+model.pg[0,j])
        elif i==1:
            model.cons.add(model.Pinj[i,j]==-loads[i,j]+model.pg[1,j])
        elif i==12:
            model.cons.add(model.Pinj[i,j]==-loads[i,j]+model.pg[2,j])
        elif i==21:
            model.cons.add(model.Pinj[i,j]==-loads[i,j]+model.pg[3,j])
        elif i==22:
            model.cons.add(model.Pinj[i,j]==-loads[i,j]+model.pg[4,j])
        elif i==26:
            model.cons.add(model.Pinj[i,j]==-loads[i,j]+model.pg[5,j])
    #      elseif i==wbus
    #      constraints=[constraints,Pinj(i,:)==-loads(i,:)+wst(1,:)];
    #      elseif i==mbus
    #      constraints=[constraints,Pinj(i,:)==-loads(i,:)+mg];
        else:
            model.cons.add(model.Pinj[i,j]==-loads[i,j])

model.theta =  Var(model.bn1, model.nt) #bus angle
Bredinv = inv(Bred)
for r in xrange(busnum-1):
    for c in xrange(nt):
        model.cons.add(model.theta[r,c] == sum(Bredinv[r,i]*model.Pinj[i,c] for i in xrange(busnum-1)))

def b3(model,i,j):
    return (-lineC[i,0], lineC[i,0])
model.pflow =  Var(model.ln, model.nt, bounds=b3) #line flow
prod = np.dot(D,m)
for r in xrange(linenum):
    for c in xrange(nt):
        model.cons.add(model.pflow[r,c] == sum(prod[r,i]*model.theta[i,c] for i in xrange(busnum-1)))
#%% Constraints
  
#  %constraints=[constraints,sum(rlsup)+sum(rgsdn)==pvup,sum(rlsdn)+sum(rgsup)==pvdn]; 
#  constraints=[constraints,sum(rgs)==0.03*sum(loads)+0.05*wst(1,:)]; % secondary reserve constraint
for h in xrange(nt):
    model.cons.add(sum(model.rgs[r,h] for r in xrange(gennum)) == sum(0.03*loads[r,h] for r in xrange(busnum)))
    model.cons.add(sum(model.pg[r,h] for r in xrange(gennum)) == sum(loads[r,h] for r in xrange(busnum))) # Power balance equation forecast
#constraints=[constraints,sum(rgs)==0.03*sum(loads)+0.05*wst(1,:)]; % secondary reserve constraint
  
# Generator Constraints
Rdn = 0.3*Gmax;
Rup = 0.3*Gmax;
for g in xrange(gennum):
    for h in xrange(nt-1):
        model.cons.add(-Rdn[g,h+1] <= model.pg[g,h+1]-model.pg[g,h] <= Rup[g,h+1]) # ramping constraints
        model.cons.add((model.pg[g,h+1]+model.rgs[g,h+1])-(model.pg[g,h]-model.rgs[g,h]) <= float(Rup[g,h+1])) # Up ramping Constraints with Rgs
        model.cons.add(float(-Rdn[g,h+1]) <= (model.pg[g,h+1]-model.rgs[g,h+1])-(model.pg[g,h]+model.rgs[g,h])) # Dn Ramping Constraints with Rgs

for g in xrange(gennum):
    for h in xrange(nt):
        model.cons.add(Gmin[g,h] * model.onoff[g,h] <= model.pg[g,h] + model.rgs[g,h])#Generator Bounds with Rgs
        model.cons.add(model.pg[g,h] + model.rgs[g,h] <= Gmax[g,h] * model.onoff[g,h])
        model.cons.add(model.rgs[g,h] <= model.z[g,h])
        model.cons.add(-model.rgs[g,h] <= model.z[g,h])
        
#%% Objective
#def exps(model):
a = sum(crgs[0,i] * model.z[i,j] for i in xrange(gennum) for j in xrange(nt))
b = sum(model.onoff[i,j]*conoff[0,i] for i in xrange(gennum) for j in xrange(nt))
c = sum(model.pg[i,j]*cg1[0,i] for i in xrange(gennum) for j in xrange(nt))
d = sum(model.pg[i,j]*model.pg[i,j]*cg2[0,i] for i in xrange(gennum) for j in xrange(nt))
#    return a+b+c+d
#model.obj = Objective(rule = exps)
model.obj = Objective(expr = a+b+c+d)

opt = SolverFactory("gurobi")
results = opt.solve(model)
model.pg.display()
model.rgs.display()
model.onoff.display()
#model.Pinj.display()
model.pflow.display()
model.obj.display()

