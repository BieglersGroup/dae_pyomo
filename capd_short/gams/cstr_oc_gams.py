# -*- coding: utf-8 -*-
#

from __future__ import division
from __future__ import print_function
#: pyomo imports
from pyomo.environ import *

__author__ = "David Thierry"  #: May 2018


#: Number of finite elements
nfe = 5
#: Number of collocation points
ncp = 3
m = ConcreteModel()

m.i = Set(initialize=range(0, nfe))
m.j = Set(initialize=range(0, ncp))

m.C = Var(m.i, m.j, initialize=1)
m.T = Var(m.i, m.j, initialize=1)
m.u = Var(m.i, m.j, initialize=1)

a_init = {}
a_init[0, 0] = 0.0
a_init[0, 1] = 0.0
a_init[0, 2] = 0.0
a_init[0, 3] = 0.0
a_init[1, 0] = 0.0
a_init[2, 0] = 0.0
a_init[3, 0] = 0.0

a_init[1, 1] = 0.19681547722366 
a_init[1, 2] = 0.39442431473909 
a_init[1, 3] = 0.37640306270047
a_init[2, 1] = -0.06553542585020 
a_init[2, 2] = 0.29207341166523 
a_init[2, 3] = 0.51248582618842
a_init[3, 1] = 0.02377097434822 
a_init[(3, 2)] = -0.04154875212600 
a_init[3, 3] = 0.11111111111111

def ainit_rule(m, i, j):
    return a_init[i, j]
m.a = Param(m.j, m.j, initialize=ainit_rule)


m.pprint()
