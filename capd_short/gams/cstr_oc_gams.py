#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

from __future__ import division
from __future__ import print_function
#: pyomo imports
from pyomo.environ import *
from pyomo.core.kernel.expr import exp
from pyomo.core.kernel.numvalue import value
from pyomo.opt import SolverFactory, SolverStatus
import matplotlib.pyplot as plt
import sys
__author__ = "David Thierry"  #: May 2018

def main():
    #: Number of finite elements
    nfe = 100
    #: Number of collocation points
    ncp = 3
    m = ConcreteModel()

    m.nfe = nfe
    m.ncp = ncp


    m.i = Set(initialize=range(1, nfe + 1))
    m.j = Set(initialize=range(1, ncp + 1))

    a_init = {}
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

    # m.cguess = Param(m.i, m.j)
    # m.tguess = Param(m.i, m.j)
    # m.ttguess = Param(m.i, m.j)
    # m.uguess = Param(m.i, m.j)

    point = 0
    slopec = (value(m.cdes) - value(m.cinit)) / (m.nfe * m.ncp)
    slopet = (value(m.tdes) - value(m.tinit)) / (m.nfe * m.ncp)
    slopeu = (value(m.udes) - value(m.uinit)) / (m.nfe * m.ncp)

    cguess = {}
    tguess = {}
    ttguess = {}
    uguess = {}

    for i in range(1, m.nfe + 1):
        for j in range(1, m.ncp + 1):
            point += 1
            cguess[i, j] = slopec * point + value(m.cinit)
            tguess[i, j] = slopet * point + value(m.tinit)
            ttguess[i, j] = value(m.time) * point
            uguess[i, j] = slopeu * point + value(m.uinit)

    #: One can pass the dictionary, as long as the keys are defined within the index set.
    m.C = Var(m.i, m.j, initialize=cguess)
    m.T = Var(m.i, m.j, initialize=tguess)
    m.u = Var(m.i, m.j, initialize=uguess)
    m.tt = Var(m.i, m.j, initialize=ttguess)

    m.cdot = Var(m.i, m.j, initialize=0.0)
    m.tdot = Var(m.i, m.j, initialize=0.0)
    m.c0 = Var(m.i, initialize=0.0)
    m.t0 = Var(m.i, initialize=0.0)
    m.tt0 = Var(m.i, initialize=0.0)
    m.phi = Var()

    def _fecolc_rule(mod, i, j):
        #: Note the ambiguity between m.j and j inside the function.
        if i <= mod.nfe:
            return mod.C[i, j] == mod.c0[i] + mod.time * mod.h[i] * sum(mod.a[k, j] * mod.cdot[i, k] for k in mod.j)
        else:
            return Constraint.Skip


    def _fecolt_rule(mod, i, j):
        if i <= mod.nfe:
            return mod.T[i, j] == mod.t0[i] + mod.time * mod.h[i] * sum(mod.a[k, j] * mod.tdot[i, k] for k in mod.j)
        else:
            return Constraint.Skip


    def _fecoltt_rule(mod, i, j):
        if i <= mod.nfe:
            return mod.tt[i, j] == mod.tt0[i] + mod.time * mod.h[i] * sum(mod.a[k, j] for k in mod.j)
        else:
            return Constraint.Skip


    def _conc_rule(mod, i):
        if 1 < i <= mod.nfe:
            return mod.c0[i] == mod.c0[i-1] + mod.time * mod.h[i-1] * sum(mod.cdot[i-1, j] * mod.a[j, 3] for j in mod.j)
        else:
            return Constraint.Skip


    def _cont_rule(mod, i):
        if 1 < i <= mod.nfe:
            return mod.t0[i] == mod.t0[i-1] + mod.time * mod.h[i-1] * sum(mod.tdot[i-1, j] * mod.a[j, 3] for j in mod.j)
        else:
            return Constraint.Skip


    def _contt_rule(mod, i):
        if 1 < i <= mod.nfe:
            return mod.tt0[i] == mod.tt0[i-1] + mod.time * mod.h[i-1] * sum(mod.a[j, 3] for j in mod.j)
        else:
            return Constraint.Skip


    def _odec_rule(mod, i, j):
        if i <= mod.nfe:
            return mod.cdot[i, j] == (1 - mod.C[i, j])/mod.theta - mod.k10 * exp(-mod.n/mod.T[i, j]) * mod.C[i, j]
        else:
            return Constraint.Skip


    def _odet_rule(mod, i, j):
        if i <= mod.nfe:
            return mod.tdot[i, j] == \
                   (mod.yf - mod.T[i, j])/mod.theta + mod.k10 * exp(-mod.n/mod.T[i, j]) * mod.C[i, j] - \
                   mod.alpha[0] * mod.u[i, j] * (mod.T[i, j] - mod.yc)
        else:
            return Constraint.Skip


    def _ic_rule(mod):
        return mod.c0[0] == mod.cinit


    def _it_rule(mod):
        return mod.t0[0] - mod.init


    m.FeColC = Constraint(m.i, m.j, rule=_fecolc_rule)
    m.FeColT = Constraint(m.i, m.j, rule=_fecolt_rule)
    m.FeColTt = Constraint(m.i, m.j, rule=_fecoltt_rule)

    m.ConC = Constraint(m.i, rule=_conc_rule)
    m.ConT = Constraint(m.i, rule=_cont_rule)
    m.ConTt = Constraint(m.i, rule=_contt_rule)

    m.OdeT = Constraint(m.i, m.j, rule=_odet_rule)
    m.OdeC = Constraint(m.i, m.j, rule=_odec_rule)

    for var in m.C.itervalues():
        var.setlb(0)
        var.setub(1)

    for var in m.T.itervalues():
        var.setlb(0.1)
        var.setub(1)

    for var in m.c0.itervalues():
        var.setlb(0)
        var.setub(value(m.cinit))

    for var in m.t0.itervalues():
        var.setlb(0.1)
        var.setub(1)


    def objective_rule(mod):
        return sum((mod.alpha[1] * (mod.cdes - mod.C[i, j]) ** 2 +
                                             mod.alpha[2] * (mod.tdes - mod.T[i, j]) ** 2 +
                                             mod.alpha[3] * (mod.udes - mod.u[i, j]) ** 2)
                   for i in m.i for j in m.j)

    # def objective_rule(mod):
    #     return sum(mod.h[i] * mod.a[j, 3] * (mod.alpha[1] * (mod.cdes - mod.C[i, j]) ** 2 +
    #                                          mod.alpha[2] * (mod.tdes - mod.T[i, j]) ** 2 +
    #                                          mod.alpha[3] * (mod.udes - mod.u[i, j]) ** 2)
    #                for i in m.i for j in m.j)

    m.fobj = Objective(sense=minimize, rule=objective_rule)
    #: call ipopt to solve the model
    ipopt = SolverFactory('ipopt')
    results = ipopt.solve(m, tee=True)

    if results.solver.status == SolverStatus.ok:
        print("Okay")
        templ = []
        cl = []
        tl = []
        ul = []
        for key in m.u.keys():
            var = m.u[key]
            if var.stale:
                continue
            templ.append(value(m.T[key]))
            cl.append(value(m.C[key]))
            tl.append(value(m.tt[key]))
            ul.append(value(m.u[key]))
        plt.subplot(3, 1, 1)
        plt.plot(tl, templ)
        plt.subplot(3, 1, 2)
        plt.plot(tl, cl)
        plt.subplot(3, 1, 3)
        plt.plot(tl, ul)
        plt.show()


if __name__ == '__main__':
    main()
