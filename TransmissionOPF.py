# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import pypower as pp
import pyomo.core
import numpy as np
from pyomo.environ import *

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
[Ybus, Yf, Yt] = makeYbus(100,ppc["bus"],ppc["branch"])
#    B=full(real(i.*Ybus));
#    NB=-B;
#    Bred=B(1:29,1:29); % Reduced B Matrix
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
ng = gennum
one = ones(nt)
Gmax = array([80*one,80*one,40*one,50*one,30*one,55*one])
Gmin = array(np.zeros((6,nt)))

model = ConcreteModel()

model.pg = Var([gennum,nt])
model.onoff = Var([gennum,nt], within = Binary)

Rgsmin = -0.5*Gmax;
Rgsmax = 0.5*Gmax;
model.rgs = Var([gennum,nt])

#    % LB & UB 
#    constraints = [Gmin<=pg<=Gmax,Rgsmin<=rgs<=Rgsmax];   

#%% Power Flow Constraints

# power flow with forecast wind
model.Pinj = Var([busnum,nt]) #% bus nodal matrix with forecast wind
#    for i=1:30
#      if i==1
#      constraints=[constraints,Pinj(i,:)==-loads(i,:)+pg(1,:)];
#      elseif i==2
#      constraints=[constraints,Pinj(i,:)==-loads(i,:)+pg(2,:)];
#      elseif i==13
#      constraints=[constraints,Pinj(i,:)==-loads(i,:)+pg(3,:)];
#      elseif i==22
#      constraints=[constraints,Pinj(i,:)==-loads(i,:)+pg(4,:)];
#      elseif i==23
#      constraints=[constraints,Pinj(i,:)==-loads(i,:)+pg(5,:)];
#      elseif i==27
#      constraints=[constraints,Pinj(i,:)==-loads(i,:)+pg(6,:)];
#      elseif i==wbus
#      constraints=[constraints,Pinj(i,:)==-loads(i,:)+wst(1,:)];
#      elseif i==mbus
#      constraints=[constraints,Pinj(i,:)==-loads(i,:)+mg];
#      else
#      constraints=[constraints,Pinj(i,:)==-loads(i,:)];
#      end
#    end
#   
#
theta =  Var([busnum-1,nt]) #bus angle
 #constraints = [constraints,(theta==inv(Bred)*Pinj(1:29,:)):'lmp'];
pflow =  Var([linenum,nt]) #line flow
#   constraints = [constraints,pflow==D*m*theta];
#   constraints = [constraints,-repmat(lineC,1,nt)<=pflow<=repmat(lineC,1,nt)];%  line constraints
#   
#%% Constraints
#  
#  %constraints=[constraints,sum(rlsup)+sum(rgsdn)==pvup,sum(rlsdn)+sum(rgsup)==pvdn]; 
#  constraints=[constraints,sum(rgs)==0.03*sum(loads)+0.05*wst(1,:)]; % secondary reserve constraint
#  
# Generator Constraints
Rdn = 0.3*Gmax;
Rup = 0.3*Gmax;
#  constraints=[constraints,-Rdn(:,2:nt)<=pg(:,2:nt)-pg(:,1:nt-1)<=Rup(:,2:nt)]; % ramping constraints
#  constraints=[constraints,Gmin.*onoff<=pg+rgs<=Gmax.*onoff]; % Generator Bounds with Rgs
#  constraints=[constraints,(pg(:,2:nt)+rgs(:,2:nt))-(pg(:,1:nt-1)-rgs(:,1:nt-1))<=Rup(:,2:nt)]; % Up ramping Constraints with Rgs
#  constraints=[constraints,-Rdn(:,2:nt)<=(pg(:,2:nt)-rgs(:,2:nt))-(pg(:,1:nt-1)+rgs(:,1:nt-1))]; % Dn Ramping Constraints with Rgs
#
#  constraints=[constraints,sum(pg)-sum(loads)+mg==0];  % Power balance equation forecast
#  
#%% Objective
#
#     Objective=sum(pg')*cg1+sum(abs(rgs'))*crgs+sum(onoff')*Conoff';
#     for i=1:nt
#         Objective=Objective+pg(:,i)'*diag(cg2)*pg(:,i);
#     end
#
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
