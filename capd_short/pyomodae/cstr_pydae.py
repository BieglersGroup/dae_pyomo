#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

from __future__ import division
from __future__ import print_function

import matplotlib.pyplot as plt
from pyomo.core.kernel.expr import exp
from pyomo.core.kernel.numvalue import value
from pyomo.dae import *
#: pyomo imports
from pyomo.environ import *
from pyomo.opt import SolverFactory, SolverStatus

__author__ = "David Thierry"  #: May 2018

def main():
    with_plots = False
    #: Number of finite elements
    nfe = 100
    #: Number of collocation points
    ncp = 3
    m = ConcreteModel()

    m.nfe = nfe
    m.ncp = ncp

    m.t = ContinuousSet(bounds=(0, 1))

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

    m.alpha = Param([0, 1, 2, 3], initialize=alpha_init)

    # m.cguess = Param(m.fe, m.cp)
    # m.tguess = Param(m.fe, m.cp)
    # m.ttguess = Param(m.fe, m.cp)
    # m.uguess = Param(m.fe, m.cp)

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
    m.C = Var(m.t, initialize=cguess)
    m.T = Var(m.t, initialize=tguess)
    m.u = Var(m.t, initialize=uguess)
    # m.tt = Var(m.t, initialize=ttguess)

    m.dC_dt = DerivativeVar(m.C)
    m.dT_dt = DerivativeVar(m.T)

    m.phi = Var()


    def _odec_rule(mod, t):
        if t > 0:
            return mod.dC_dt[t] == \
                   ((1 - mod.C[t]) / mod.theta - mod.k10 * exp(-mod.n / mod.T[t]) * mod.C[t])
        else:
            return Constraint.Skip


    def _odet_rule(mod, t):
        if t > 0:
            return mod.dT_dt[t] == \
                   ((mod.yf - mod.T[t]) / mod.theta + mod.k10 * exp(-mod.n / mod.T[t]) * mod.C[t] -
                    mod.alpha[0] * mod.u[t] * (mod.T[t] - mod.yc))
        else:
            return Constraint.Skip


    def _ic_rule(mod):
        return mod.C[0] == mod.cinit


    def _it_rule(mod):
        return mod.T[0] == mod.tinit


    m.OdeT = Constraint(m.t, rule=_odet_rule)
    m.OdeC = Constraint(m.t, rule=_odec_rule)

    m.IC = Constraint(rule=_ic_rule)
    m.IT = Constraint(rule=_it_rule)


    def objective_rule(mod):
        return sum((mod.alpha[1] * (mod.cdes - mod.C[t]) ** 2 +
                    mod.alpha[2] * (mod.tdes - mod.T[t]) ** 2 +
                    mod.alpha[3] * (mod.udes - mod.u[t]) ** 2)
                   for t in m.t if t > 0)


    m.fobj = Objective(sense=minimize, rule=objective_rule)

    dae = TransformationFactory('dae.collocation')
    dae.apply_to(m, nfe=m.nfe, ncp=m.ncp, scheme='LAGRANGE-RADAU')


    for var in m.C.itervalues():
        var.setlb(0)
        var.setub(1)

    for var in m.T.itervalues():
        var.setlb(0.1)
        var.setub(1)

    for var in m.u.itervalues():
        var.setlb(0)
        var.setub(500)

    for t in m.t:
        m.C[t].set_value(slopec * t + value(m.cinit))
        m.T[t].set_value(slopet * t + value(m.tinit))
        m.u[t].set_value(slopeu * t + value(m.uinit))

    ipopt = SolverFactory('ipopt')
    results = ipopt.solve(m, tee=True)
    #
    # if results.solver.status == SolverStatus.ok:
    #     print("Okay")
    #     ul = []
    #     tl = []
    #     for key in m.t:
    #         var = m.u[key]
    #         if var.stale:
    #             continue
    #         ul.append(value(var))
    #         tl.append(key)
    #     plt.plot(tl, ul)
    #     plt.show()

    if results.solver.status == SolverStatus.ok and with_plots:
        print("Okay")
        templ = []
        cl = []
        tl = []
        ul = []
        for key in m.t:
            var = m.u[key]
            if var.stale:
                continue
            templ.append(value(m.T[key]))
            cl.append(value(m.C[key]))
            tl.append(key)
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