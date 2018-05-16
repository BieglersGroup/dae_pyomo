# -*- coding: utf-8 -*-
#

from __future__ import division
from __future__ import print_function
#: pyomo imports
from pyomo.environ import *
from pyomo.core.kernel.expr import exp
from pyomo.core.kernel.numvalue import value
import sys
__author__ = "David Thierry"  #: May 2018


#: Number of finite elements
nfe = 5
#: Number of collocation points
ncp = 3
m = ConcreteModel()

m.nfe = nfe
m.ncp = ncp


m.i = Set(initialize=range(1, nfe + 1))
m.j = Set(initialize=range(1, ncp + 1))

m.C = Var(m.i, m.j, initialize=0.0)
m.T = Var(m.i, m.j, initialize=0.0)
m.u = Var(m.i, m.j, initialize=0.0)
m.tt = Var(m.i, m.j, initialize=0.0)

m.cdot = Var(m.i, m.j, initialize=0.0)
m.tdot = Var(m.i, m.j, initialize=0.0)
m.c0 = Var(m.i, initialize=0.0)
m.t0 = Var(m.i, initialize=0.0)
m.tt0 = Var(m.i, initialize=0.0)
m.phi = Var()

a_init = {}
# a_init[0, 0] = 0.0
# a_init[0, 1] = 0.0
# a_init[0, 2] = 0.0
# a_init[0, 3] = 0.0
# a_init[1, 0] = 0.0
# a_init[2, 0] = 0.0
# a_init[3, 0] = 0.0

a_init[1, 1] = 0.19681547722366 
a_init[1, 2] = 0.39442431473909 
a_init[1, 3] = 0.37640306270047
a_init[2, 1] = -0.06553542585020 
a_init[2, 2] = 0.29207341166523 
a_init[2, 3] = 0.51248582618842
a_init[3, 1] = 0.02377097434822 
a_init[3, 2] = -0.04154875212600
a_init[3, 3] = 0.11111111111111

alpha_init = {0: 1.95e-04, 1: 1e+06, 2: 2e+03, 3: 1e-03}


#: Scalars
m.cinit = Param(initialize=0.1367)
m.tinit = Param(initialize=0.7293)
m.uinit = Param(initialize=390.0)
m.cdes = Param(initialize=0.0944)
m.tdes = Param(initialize=0.7766)
m.udes = Param(initialize=340)
m.k10 = Param(initialize=300)
m.n = Param(initialize=5)


m.cf = Param(initialize=7.6)
m.tf = Param(initialize=300)
m.tc = Param(initialize=290)



m.theta = Param(initialize=20)
m.yf = Param(initialize=0.3947)
m.yc = Param(initialize=0.3816)
m.time = Param(initialize=10)
m.point = Param(initialize=0)

# m.nfe = Param(initialize=100)
# m.ncp = Param(initialize=3)
m.slopec = Param()
m.slopet = Param()
m.slopeu = Param()
m.ii = Param()
m.jj = Param()





m.a = Param(m.j, m.j, initialize=a_init)
m.alpha = Param([0, 1, 2, 3], initialize=alpha_init)
m.h = Param(m.i, initialize=1/value(m.nfe))

m.cguess = Param(m.i, m.j)
m.tguess = Param(m.i, m.j)
m.ttguess = Param(m.i, m.j)
m.uguess = Param(m.i, m.j)


def _fecolc_rule(m, i, j):
    #: Note the ambiguity between m.j and j inside the function.
    if i <= m.nfe:
        return m.C[i, j] == m.c0[i] + m.time * m.h[i] * sum(m.a[k, j] * m.cdot[i, k] for k in m.j)
    else:
        return Constraint.Skip


def _fecolt_rule(m, i, j):
    if i <= m.nfe:
        return m.T[i, j] == m.t0[i] + m.time * m.h[i] * sum(m.a[k, j] * m.tdot[i, k] for k in m.j)
    else:
        return Constraint.Skip


def _fecoltt_rule(m, i, j):
    if i <= m.nfe:
        return m.tt[i, j] == m.tt0[i] + m.time * m.h[i] * sum(m.a[k, j] for k in m.j)
    else:
        return Constraint.Skip


def _conc_rule(m, i):
    if 1 < i <= m.nfe:
        return m.c0[i] == m.c0[i-1] + m.time * m.h[i-1] * sum(m.cdot[i-1, j] * m.a[j, 3] for j in m.j)
    else:
        return Constraint.Skip


def _cont_rule(m, i):
    if 1 < i <= m.nfe:
        return m.t0[i] == m.t0[i-1] + m.time * m.h[i-1] * sum(m.tdot[i-1, j] * m.a[j, 3] for j in m.j)
    else:
        return Constraint.Skip


def _contt_rule(m, i):
    if 1 < i <= m.nfe:
        return m.tt0[i] == m.tt0[i-1] + m.time * m.h[i-1] * sum(m.a[j, 3] for j in m.j)
    else:
        return Constraint.Skip


def _odec_rule(m, i, j):
    if i <= m.nfe:
        return m.cdot[i, j] == (1 - m.C[i, j])/m.theta - m.k10 * exp(-m.n/m.T[i, j]) * m.C[i, j]


def _odet_rule(m, i, j):
    if i <= m.nfe:
        return m.tdot[i, j] == \
               (m.yf - m.T[i, j])/m.theta + m.k10 * exp(-m.n/m.T[i, j]) * m.C[i, j] - \
               m.alpha[0] * m.u[i, j] * (m.T[i, j] - m.yc)
    else:
        return Constraint.Skip


def _ic_rule(m):
    return m.c0[0] == m.cinit


def _it_rule(m):
    return m.t0[0] - m.init


m.FeColC = Constraint(m.i, m.j, rule=_fecolc_rule)
m.FeColT = Constraint(m.i, m.j, rule=_fecolt_rule)
m.FeColTt = Constraint(m.i, m.j, rule=_fecoltt_rule)

m.ConC = Constraint(m.i, rule=_conc_rule)
m.ConT = Constraint(m.i, rule=_cont_rule)
m.ConTt = Constraint(m.i, rule=_contt_rule)

m.OdeTt = Constraint(m.i, m.j, rule=_odec_rule)
m.OdeC = Constraint(m.i, m.j, rule=_odet_rule)


