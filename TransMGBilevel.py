# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import pypower as pp
import pyomo.core
import numpy as np
from pypower.api import case30,makeYbus
from numpy.linalg import inv
from pyomo.environ import *
from pyomo.bilevel import *
from pyomo.opt import SolverFactory




#%% Upper Level problem: Transmission network 
ppc = case30()

nt = 10 # time period
gennum = ppc["gencost"].shape[0]
busnum = ppc["bus"].shape[0]
linenum = ppc["branch"].shape[0]
 
# define cost
cg1 = ppc["gencost"][:, [5]].T # % 1st order generator cost
cg2 = ppc["gencost"][:, [4]].T # 2nd order generator cost
crgs = np.array([[6, 6.75, 7, 5.25, 5, 5]])
conoff = 2*np.ones((1,gennum)) # generator on cost

sampleL = np.array([[120.6, 115.8, 114.8, 112.6, 114.0, 113.4, 
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
one = np.ones(nt)
Gmax = np.array([80*one,80*one,40*one,50*one,30*one,55*one])
Gmin = np.array(np.zeros((6,nt)))

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
model.obj = Objective(expr = a+b+c+d, sense = minimize)








#%% Lower Level Problem：Microgrid Model
#%% System Parameters
lmp = 5.3*np.ones(nt)
cg = np.array([[3.3, 3.3, 3.3]]) 
cd = 1.01*np.array([[3.3, 3.3, 3.3, 3.3, 3.3, 3.3]]) # slightly larger than cg to encourage cost
cb = 0.1 # needs to be small so that MG stores energy rather than export energy
ci = lmp # change
ce = 0.8*lmp #change
#
ndl = 6 # number of dispatchable load
ng = 3 # number of generators
nb = 1 # number of storage
#
model.sub = SubModel()
model.sub.nb = RangeSet(0,nb-1)
model.sub.ndl = RangeSet(0,ndl-1)
model.sub.ntt = RangeSet(0,nt-1)
model.sub.ng = RangeSet(0,ng-1)
#%% Define Optimization Variables and Their Bounds
# 1. Dispatchable loads

one = np.ones(nt)
Pdmin = 0.1*np.array([0.5*one,4*one,2*one,5.5*one,1*one,7*one])
Pdmax = 0.03*np.array([10*one,16*one,15*one,20*one,27*one,32*one])

def b1(model,i,j):
    return (Pdmin[i,j],Pdmax[i,j])
model.sub.pd=Var(model.sub.ndl, model.sub.ntt, bounds = b1)

# 2. Non-controllable load    
nload = 0.03 * np.array([[120.6, 115.8, 114.8, 112.6, 114.0, 113.4, 
                  117.1, 126.3, 130.7, 132.5, 135.6, 134.8, 
                  136.5, 137.7, 137.1, 138.0, 136.3, 133.3, 
                  131.7, 129.3, 128.2, 127.4, 125.6, 124.2]]) # sample load profile
nload = nload[:,0:nt]

# 3. Generators
Gmax = 0.3 * np.array([5*one,4.5*one,7*one])
Gmin = np.array([1*one, 0.8*one, 1.5*one])
def b2(model,i,j):
    return (Gmin[i,j],Gmax[i,j])
model.sub.pg = Var(model.sub.ng, model.sub.ntt, bounds=b2)

Rgsmin = -0.3*Gmax
Rgsup = 0.3*Gmax
def b8(model,i,j):
    return (0,0)
model.sub.rgs = Var(model.sub.ng, model.sub.ntt, bounds=b8)
model.sub.z = Var(model.sub.ng, model.sub.ntt)

# 4. Storage
# power to and from storage
Pbmin = -3*np.ones((nb,nt))
Pbmax = 3*np.ones((nb,nt))
def b3(model,i,j):
    return (Pbmin[i,j],Pbmax[i,j])
model.sub.pb = Var(model.sub.nb, model.sub.ntt, bounds=b3)

# state of charge
Bmin = np.zeros((nb,nt))
Bmax = 10*np.ones((nb,nt))
def b4(model,i,j):
    return (Bmin[i,j],Bmax[i,j])
model.sub.b = Var(model.sub.nb, model.sub.ntt, bounds=b4)



# 4. Import & Export
def b5(model,i):
    return (0, 1000)
model.sub.ex = Var(model.sub.ntt, bounds=b5)

def b6(model,i):
    return (0,1000)
model.sub.im = Var(model.sub.ntt, bounds=b6)

Nmin = -20*one
Nmax = 20*one
def b7(model,i):
    return (Nmin[i],Nmax[i])
model.sub.net = Var(model.sub.ntt, bounds=b7)
 
#%% Constraints
# DG Constraints
Rdn = 0.3*Gmax
Rup = 0.3*Gmax
model.sub.cons = ConstraintList()
model.sub.cons.add(model.sub.b[0,0] <= 5) # initialize battery energy state
model.sub.cons.add(model.sub.b[0,0] >= 5) # initialize battery energy state

for g in xrange(ng):
    for h in xrange(nt-1):
        model.sub.cons.add(-Rdn[g,h+1] <= model.sub.pg[g,h+1]-model.sub.pg[g,h] <= Rup[g,h+1]) # ramping constraints
        model.sub.cons.add((model.sub.pg[g,h+1]+model.sub.rgs[g,h+1])-(model.sub.pg[g,h]-model.sub.rgs[g,h]) <= float(Rup[g,h+1])) # Up ramping Constraints with Rgs
        model.sub.cons.add(float(-Rdn[g,h+1]) <= (model.sub.pg[g,h+1]-model.sub.rgs[g,h+1])-(model.sub.pg[g,h]+model.sub.rgs[g,h])) # Dn Ramping Constraints with Rgs

for g in xrange(ng):
    for h in xrange(nt):
#        model.sub.cons.add(Gmin[g,h] <= model.sub.pg[g,h] + model.sub.rgs[g,h])#Generator Bounds with Rgs
        model.sub.cons.add(model.sub.pg[g,h] + model.sub.rgs[g,h] <= Gmax[g,h])
        model.sub.cons.add(model.sub.rgs[g,h] <= model.sub.z[g,h])
        model.sub.cons.add(-model.sub.rgs[g,h] <= model.sub.z[g,h])

for h in xrange(nt):
    model.sub.cons.add(sum(model.sub.pg[r,h] for r in xrange(ng)) >= sum(nload[r1,h] + model.sub.pb[r1,h] for r1 in xrange(1))+ model.sub.net[h] + sum(model.sub.pd[r2,h] for r2 in xrange(ndl))) # Power balance equation forecast 
    model.sub.cons.add(sum(model.sub.pg[r,h] for r in xrange(ng)) <= sum(nload[r1,h] + model.sub.pb[r1,h] for r1 in xrange(1))+ model.sub.net[h] + sum(model.sub.pd[r2,h] for r2 in xrange(ndl))) # Power balance equation forecast 
    model.sub.cons.add(model.sub.ex[h] - model.sub.im[h] <= model.sub.net[h])
    model.sub.cons.add(model.sub.ex[h] - model.sub.im[h] >= model.sub.net[h])
    
for i in xrange(nb):
    for j in xrange(nt-1):
        model.sub.cons.add(model.sub.b[i,j+1] >= model.sub.b[i,j]+model.sub.pb[i,j])
        model.sub.cons.add(model.sub.b[i,j+1] <= model.sub.b[i,j]+model.sub.pb[i,j])

for i in xrange(nb):
    for j in xrange(nt):
        model.sub.cons.add(-model.sub.pb[i,j] <= model.sub.b[i,j])

#%% Objective
aa = sum(cb * model.sub.b[i,j] for i in xrange(nb) for j in xrange(nt))
bb = -sum(model.sub.pd[i,j]*cd[0,i] for i in xrange(ndl) for j in xrange(nt))
cc = sum(model.sub.pg[i,j]*cg[0,i] for i in xrange(ng) for j in xrange(nt))
dd = -sum(model.sub.ex[i]*ce[i] for i in xrange(nt))
ee = sum(model.sub.im[i]*ci[i] for i in xrange(nt))

model.sub.obj = Objective(expr = aa+bb+cc+dd+ee, sense = minimize)
#model.sub.obj = Objective(expr = aa, sense = minimize)




#model.sub.y = Var(bounds=(0,None))
#model.sub.o = Objective(expr=model.sub.y, sense=minimize)
#model.sub.con = ConstraintList()
#model.sub.con.add ( model.sub.y >= 3)
#model.sub.c3 = Constraint(expr=model.sub.y <= 12)

opt = SolverFactory("bilevel_blp_global")
opt.options["solver"] = 'gurobi'
results = opt.solve(model,tee=True)


model.sub.o.display()
model.sub.y.display()
model.sub.pg.display()
model.sub.rgs.display()
model.sub.pb.display()
model.sub.b.display()
model.sub.pd.display()
model.sub.ex.display()
model.sub.im.display()
model.sub.net.display()
model.sub.obj.display()  

#model.pg.display()
#model.rgs.display()
#model.onoff.display()
##model.Pinj.display()
#model.pflow.display()
#model.obj.display()

