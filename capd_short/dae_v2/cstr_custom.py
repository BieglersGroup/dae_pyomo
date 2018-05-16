#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

from __future__ import division
from __future__ import print_function

import matplotlib.pyplot as plt
from pyomo.core.kernel.expr import exp
from pyomo.core.kernel.numvalue import value
#: pyomo imports
from pyomo.environ import *
from pyomo.opt import SolverFactory, SolverStatus

from capd_short.dae_v2.collocation_funcs.cpoinsc import collptsgen
from capd_short.dae_v2.collocation_funcs.lagrange_f import lgr, lgrdot

__author__ = "David Thierry"  #: May 2018

def main():
    #: Number of finite elements
    nfe = 100
    #: Number of collocation points
    ncp = 3
    m = ConcreteModel()

    m.nfe = nfe
    m.ncp = ncp


    m.fe = Set(initialize=range(0, nfe))
    m.cp = Set(initialize=range(0, ncp + 1))


    m.tau_t = collptsgen(m.ncp, 1, 0)

    # start at zero
    m.tau_i_t = {0: 0.}
    # create a list

    for ii in range(1, m.ncp + 1):
        m.tau_i_t[ii] = m.tau_t[ii - 1]

    m.taucp_t = Param(m.cp, initialize=m.tau_i_t)

    m.ldot_t = Param(m.cp, m.cp,
                     initialize=(lambda m, j, k: lgrdot(k, m.taucp_t[j], m.ncp, 1, 0)))
    m.l1_t = Param(m.cp,
                   initialize=(lambda m, j: lgr(j, 1, m.ncp, 1, 0)))

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
    m.h = Param(m.fe, initialize=1/value(m.nfe))

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
    m.C = Var(m.fe, m.cp, initialize=cguess)
    m.T = Var(m.fe, m.cp, initialize=tguess)
    m.u = Var(m.fe, m.cp, initialize=uguess)
    # m.tt = Var(m.fe, m.cp, initialize=ttguess)

    m.dC_dt = Var(m.fe, m.cp, initialize=0.0)
    m.dT_dt = Var(m.fe, m.cp, initialize=0.0)

    m.c0 = Var(m.fe, initialize=0.0)
    m.t0 = Var(m.fe, initialize=0.0)
    m.tt0 = Var(m.fe, initialize=0.0)
    m.phi = Var()

    def _c_coll(mod, i, j):
        if j > 0:
            return mod.dC_dt[i, j] == sum(mod.ldot_t[j, k] * mod.C[i, k] for k in mod.cp)
        else:
            return Constraint.Skip

    def _c_cont(mod, i):
        if i < mod.nfe - 1:
            return mod.C[i + 1, 0] == sum(mod.l1_t[j] * mod.C[i, j] for j in mod.cp)
        else:
            return Constraint.Skip

    def _t_coll(mod, i, j):
        if j > 0:
            return mod.dT_dt[i, j] == sum(mod.ldot_t[j, k] * mod.T[i, k] for k in mod.cp)
        else:
            return Constraint.Skip

    def _t_cont(mod, i):
        if i < mod.nfe - 1:
            return mod.T[i + 1, 0] == sum(mod.l1_t[j] * mod.T[i, j] for j in mod.cp)
        else:
            return Constraint.Skip

    def _odec_rule(mod, i, j):
        if j > 0:
            return mod.dC_dt[i, j] == \
                   ((1 - mod.C[i, j]) / mod.theta - mod.k10 * exp(-mod.n / mod.T[i, j]) * mod.C[i, j]) * mod.h[i]
        else:
            return Constraint.Skip


    def _odet_rule(mod, i, j):
        if j > 0:
            return mod.dT_dt[i, j] == \
                   ((mod.yf - mod.T[i, j]) / mod.theta + mod.k10 * exp(-mod.n / mod.T[i, j]) * mod.C[i, j] -
                    mod.alpha[0] * mod.u[i, j] * (mod.T[i, j] - mod.yc)) * mod.h[i]
        else:
            return Constraint.Skip


    def _ic_rule(mod):
        return mod.C[0, 0] == mod.cinit


    def _it_rule(mod):
        return mod.T[0, 0] - mod.init


    m.Col_C = Constraint(m.fe, m.cp, rule=_c_coll)
    m.Col_T = Constraint(m.fe, m.cp, rule=_t_coll)

    m.Con_C = Constraint(m.fe, rule=_c_cont)
    m.Con_T = Constraint(m.fe, rule=_t_cont)

    m.OdeT = Constraint(m.fe, m.cp, rule=_odet_rule)
    m.OdeC = Constraint(m.fe, m.cp, rule=_odec_rule)

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
                   for i in m.fe for j in m.cp if j > 0)

    m.fobj = Objective(sense=minimize, rule=objective_rule)


    ipopt = SolverFactory('ipopt')
    results = ipopt.solve(m, tee=True)
    tij = {}
    tij[(0, 0)] = 0
    for j in m.cp:
        tij[(0, j)] = m.h[0] * m.tau_i_t[j]
    for i in m.fe:
        if i > 0:
            tij[(i, 0)] = tij[i - 1, 0] + m.h[i]
            for j in m.cp:
                if j > 0:
                    tij[(i, j)] = tij[i, 0] + m.h[i] * m.tau_i_t[j]

    print(tij)
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
            tl.append(tij[key])
            ul.append(value(m.u[key]))
            print(tij[key])
        plt.subplot(3, 1, 1)
        plt.plot(tl, templ)
        plt.subplot(3, 1, 2)
        plt.plot(tl, cl)
        plt.subplot(3, 1, 3)
        plt.plot(tl, ul)
        plt.show()


if __name__ == '__main__':
    main()