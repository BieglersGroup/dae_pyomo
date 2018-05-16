#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
from __future__ import print_function

from pyomo.environ import *
from pyomo.dae import *
from pyomo.opt import *
import matplotlib.pyplot as plt

#: URL: https://github.com/BieglersGroup/dae_pyomo

__author__ = 'David Thierry'

#: Example 10.2 Larry's book
#: dzdt = z**2  - 2*z +1, z(0) = -3

def solvemodel(nfe, ncp):
    m = ConcreteModel()
    m.t = ContinuousSet(bounds=(0,1))
    m.z = Var(m.t, initialize=1)
    m.dzdt = DerivativeVar(m.z)

    m.c = Constraint(m.t,
                     rule=lambda d, t: d.dzdt[t] == (d.z[t]**2) - 2*d.z[t] + 1 if t>0 else Constraint.Skip)
    m.z[0].fix(-3)
    m2 = m.clone()

    #: Discretize
    dae = TransformationFactory('dae.collocation')
    dae.apply_to(m2, nfe=nfe, ncp=ncp)

    ipopt = SolverFactory('ipopt')
    ipopt.solve(m2, tee=True)

    #:Some plots
    tl = []
    zl = []
    for key in m2.t:
        tl.append(key)
        zl.append(value(m2.z[key]))
    return tl, zl

def plot_data(tl, zl):
    fig = plt.plot(tl, zl, 'o-')
    plt.show()


def main(i, j):
    tl, zl = solvemodel(i, j)
    plot_data(tl, zl)


if __name__ == '__main__':
    main(1, 3)
