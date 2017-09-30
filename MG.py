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
lmp = 5.3*np.ones((1,nt))
cg = np.array([[3.3, 3.3, 3.3]]) 
cg = 1.01*np.array([[3.3, 3.3, 3.3, 3.3, 3.3, 3.3]]) # slightly larger than cg to encourage cost
cb = 0.1 # needs to be small so that MG stores energy rather than export energy
ci = lmp # change
ce = 0.8*lmp #change

ndl = 6 # number of dispatchable load
ng = 3 # number of generators
model = ConcreteModel()
model.ndl = RangeSet(0,ndl-1)
model.nt = RangeSet(0,nt-1)
model.ng = RangeSet(0,ng-1)
model.ln = RangeSet(0,linenum-1)
#%% Define Optimization Variables and Their Bounds
# 1. Dispatchable loads

one = np.ones(nt)
Pdmin = 0.1*np.array([0.5*one,4*one,2*one,5.5*one,1*one,7*one])
Pdmax = 0.03*np.array([10*one,16*one,15*one,20*one,27*one,32*one])
def b1(model,i,j):
    return (Pdmin[i,j],Pdmax[i,j])
pd=sdpvar(model.nd, model.nt, bounds = b1)

# 2. Non-controllable load    
NCL = 0.03 * np.array([[120.6, 115.8, 114.8, 112.6, 114.0, 113.4, 
                  117.1, 126.3, 130.7, 132.5, 135.6, 134.8, 
                  136.5, 137.7, 137.1, 138.0, 136.3, 133.3, 
                  131.7, 129.3, 128.2, 127.4, 125.6, 124.2]]) # sample load profile
NCL = NCL[:,0:nt];

# 3. Generators
Gmax = 0.3*[5*ones(1,nt);4.5*ones(1,nt);7*ones(1,nt)];
Gmin=[1*ones(1,nt);0.8*ones(1,nt);1.5*ones(1,nt)];
pg=sdpvar(ng,nt,'full');
rgs=sdpvar(ng,nt,'full');

% 4. Storage
nb=1;
% power to and from storage
Pbmin=-3*ones(nb,nt);
Pbmax=3*ones(nb,nt);
pb=sdpvar(nb,nt,'full');

% state of charge
Bmin=zeros(nb,nt);
Bmax=10*ones(nb,nt);
b=sdpvar(nb,nt,'full');
b(1)=5;% initialize battery energy state

% 4. Import & Export
ex=sdpvar(1,nt,'full');
im=sdpvar(1,nt,'full');
net=sdpvar(1,nt,'full');
Nmin=-20*ones(1,nt);
Nmax=20*ones(1,nt);

% LB & UB 
constraints = [Pdmin<=pd<=Pdmax,Gmin<=pg<=Gmax,Pbmin<=pb<=Pbmax,Bmin<=b<=Bmax,0<=ex,0<=im,Nmin<=net<=Nmax]; 
              