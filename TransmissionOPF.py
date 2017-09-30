# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import pypower as pp
import pyomo.core
import numpy as np

from numpy.linalg import inv

ppc = case30()

nt = 24 # time period
gennum = ppc["gencost"].shape[0]
busnum = ppc["bus"].shape[0]
linenum = ppc["branch"].shape[0]
 
# define cost
cg1 = ppc["gencost"][:, [5]].T # % 1st order generator cost
cg2 = ppc["gencost"][:, [4]].T # 2nd order generator cost
crgs = array([[6, 6.75, 7, 5.25, 5, 5]])
conoff = ones((1,gennum)) # generator on cost

sampleL = array([[120.6, 115.8, 114.8, 112.6, 114.0, 113.4, 
                  117.1, 126.3, 130.7, 132.5, 135.6, 134.8, 
                  136.5, 137.7, 137.1, 138.0, 136.3, 133.3, 
                  131.7, 129.3, 128.2, 127.4, 125.6, 124.2]]) # sample load profile
sf = sampleL/np.mean(sampleL)
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
Bred = B[0:29,0:29] # Reduced B Matrix
#    
frm = ppc["branch"][:,0]
to = ppc["branch"][:,1]
M = np.zeros((linenum,busnum)) # from/to matrix
lineC = np.zeros((linenum,1))

for i in xrange(41):
    M[i,int(frm[i]-1)] = 1
    M[i,int(to[i]-1)]= -1
    lineC[i]=1*ppc["branch"][i,6]
    
m = M[:,0:29]    
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

model.pg = Var(model.gn, model.nt)
model.onoff = Var(model.gn, model.nt, within = Binary)

Rgsmin = -0.5*Gmax;
Rgsmax = 0.5*Gmax;
model.rgs = Var(model.gn, model.nt)

#    % LB & UB 
#    constraints = [Gmin<=pg<=Gmax,Rgsmin<=rgs<=Rgsmax];   

#%% Power Flow Constraints

# power flow with forecast wind

model.Pinj = Var(model.bn, model.nt) #% bus nodal matrix with forecast wind
#def injection(model,i,j):
#    
#model.inj = Constraint(model.bn, model.nt, rule = injection)
model.cons = ConstraintList()
for i in xrange(30):
    for j in xrange(24):
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
for r in xrange(29):
    for c in xrange(24):
        model.cons.add(model.theta[r,c] == sum(Bredinv[r,i]*model.Pinj[i,c] for i in xrange(29)))

model.pflow =  Var(model.ln, model.nt) #line flow
prod = np.dot(D,m)
for r in xrange(41):
    for c in xrange(24):
        model.cons.add(model.pflow[r,c] == sum(prod[r,i]*model.theta[i,c] for i in xrange(29)))
        model.cons.add(-lineC[r,0] <= model.pflow[r,c] <= lineC[r,0])  
#%% Constraints
  
#  %constraints=[constraints,sum(rlsup)+sum(rgsdn)==pvup,sum(rlsdn)+sum(rgsup)==pvdn]; 
#  constraints=[constraints,sum(rgs)==0.03*sum(loads)+0.05*wst(1,:)]; % secondary reserve constraint
for h in xrange(24):
    model.cons.add(sum(model.rgs[r,h] for r in xrange(6)) == sum(0.03*loads[r,h] for r in xrange(30)))
    model.cons.add(sum(model.pg[r,h] for r in xrange(6)) == sum(loads[r,h] for r in xrange(30))) # Power balance equation forecast
#constraints=[constraints,sum(rgs)==0.03*sum(loads)+0.05*wst(1,:)]; % secondary reserve constraint
  
# Generator Constraints
Rdn = 0.3*Gmax;
Rup = 0.3*Gmax;
for g in xrange(6):
    for h in xrange(23):
        model.cons.add(-Rdn[g,h+1] <= model.pg[g,h+1]-model.pg[g,h] <= Rup[g,h+1]) # ramping constraints
        model.cons.add((model.pg[g,h+1]+model.rgs[g,h+1])-(model.pg[g,h]-model.rgs[g,h]) <= float(Rup[g,h+1])) # Up ramping Constraints with Rgs
        model.cons.add(float(-Rdn[g,h+1]) <= (model.pg[g,h+1]-model.rgs[g,h+1])-(model.pg[g,h]+model.rgs[g,h])) # Dn Ramping Constraints with Rgs
#  constraints=[constraints,-Rdn(:,2:nt)<=pg(:,2:nt)-pg(:,1:nt-1)<=Rup(:,2:nt)]; % ramping constraints
for g in xrange(6):
    for h in xrange(24):
        model.cons.add(Gmin[g,h] * model.onoff[g,h] <= model.pg[g,h] + model.rgs[g,h])#Generator Bounds with Rgs
        model.cons.add(model.pg[g,h] + model.rgs[g,h] <= Gmax[g,h] * model.onoff[g,h])
        
#  constraints=[constraints,Gmin.*onoff<=pg+rgs<=Gmax.*onoff]; % Generator Bounds with Rgs
#  constraints=[constraints,(pg(:,2:nt)+rgs(:,2:nt))-(pg(:,1:nt-1)-rgs(:,1:nt-1))<=Rup(:,2:nt)]; % Up ramping Constraints with Rgs
#  constraints=[constraints,-Rdn(:,2:nt)<=(pg(:,2:nt)-rgs(:,2:nt))-(pg(:,1:nt-1)+rgs(:,1:nt-1))]; % Dn Ramping Constraints with Rgs
#
#  constraints=[constraints,sum(pg)-sum(loads)==0];  % Power balance equation forecast
#  
#%% Objective
def exps(model):
    a = sum(crgs[0,i]*abs(model.rgs[i,j]) for i in xrange(6) for j in xrange(24))
    b = sum(model.onoff[i,j]*conoff[0,i] for i in xrange(6) for j in xrange(24))
    c = sum(model.pg[i,j]*cg1[0,i] for i in xrange(6) for j in xrange(24))
    d = sum(model.pg[i,j]*model.pg[i,j]*cg2[0,i] for i in xrange(6) for j in xrange(24))
    return a+c
model.obj = Objective(rule = exps)

opt = SolverFactory("gurobi")
results = opt.solve(model)
model.pg.display()
#     Objective=sum(pg')*cg1+sum(abs(rgs'))*crgs+sum(onoff')*Conoff';
#     for i=1:nt
#         Objective=Objective+pg(:,i)'*diag(cg2)*pg(:,i);
#     end

#%% Solve Problem
#
#    %sdpsettings('solver','gurobi');
#%     opts = sdpsettings('verbose',0,'saveduals',0);
#    solopt=optimize(constraints,Objective)
#    z=dual(constraints('lmp'));
#    sol = solvesdp(constraints,Objective,sdpsettings('savesolveroutput',1))
#sol.solveroutput
#%% Optimization Variables
#
#    pg=value(pg);
#    rgs=value(rgs);
#    onoff=value(onoff);
#    Pinj=value(Pinj);
#    pflow=value(pflow);
#    Objective=value(Objective);
