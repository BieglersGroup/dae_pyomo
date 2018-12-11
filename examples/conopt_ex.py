#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function
from pyomo.environ import *
from pyomo.opt import SolverFactory, ProblemFormat

"""Example taken from the sipopt manual
please check
https://github.com/coin-or/Ipopt/blob/master/Ipopt/contrib/sIPOPT/examples/redhess_ampl/red_hess.run"""

__author__ = 'David Thierry'  #: @2018
def main():
    #: Declare Model
    m = ConcreteModel()

    m.i = Set(initialize=[1, 2, 3])

    init_vals = {1:25E+00, 2:0.0, 3:0.0}
    #: Variables
    m.x = Var(m.i, initialize=init_vals)
    #: Objective
    m.oF = Objective(expr=(m.x[1] - 1.0)**2 + exp(m.x[2] - 2.0)**2 + (m.x[3] - 3.0)**2, sense=minimize)
    #: Constraints
    m.c1 = Constraint(expr=m.x[1] + 2 * m.x[2] + 3 * m.x[3] == 0.0)

    #: ipopt suffixes  REQUIRED FOR K_AUG!
    m.dual = Suffix(direction=Suffix.IMPORT_EXPORT)
    m.ipopt_zL_out = Suffix(direction=Suffix.IMPORT)
    m.ipopt_zU_out = Suffix(direction=Suffix.IMPORT)
    m.ipopt_zL_in = Suffix(direction=Suffix.EXPORT)
    m.ipopt_zU_in = Suffix(direction=Suffix.EXPORT)

    #: sipopt suffix

    ipopt = SolverFactory('ipopt')
    conopt = SolverFactory('conopt')
    m.write(filename='mynl.nl', format=ProblemFormat.nl)
    conopt.options['outlev'] = 3
    conopt.solve(m, tee=True)
    #: works with Pyomo 5.5.1 (CPython 2.7.15 on Linux 4.15.0-42-generic)
if __name__ == '__main__':
    main()
