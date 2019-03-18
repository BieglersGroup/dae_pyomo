#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

from pyomo.environ import *
from pyomo.opt import SolverFactory, ProblemFormat

"""Example taken from the sipopt manual
please check
https://github.com/coin-or/Ipopt/blob/master/Ipopt/contrib/sIPOPT/examples/"""
#: Please check:
#: https://projects.coin-or.org/Ipopt/wiki/sIpopt

__author__ = 'David Thierry'  #: @2019

#: Declare Model
m = ConcreteModel()

m.i = Set(initialize=[1, 2, 3])

init_vals = {1:0.15, 2:0.15, 3:0.0}
#: Variables
m.x = Var(m.i, initialize=init_vals)
m.p1 = Var()   #: dummy 1
m.p2 = Var()   #: dummy 2

#: Objective
m.oF = Objective(expr=m.x[1]**2 + m.x[2]**2 + m.x[3]**2, sense=minimize)
#: Constraints
m.c1 = Constraint(expr=6.0 * m.x[1] + 3.0 * m.x[2] + 2 * m.x[3] - m.p1 == 0.0)
m.c2 = Constraint(expr=m.p2 * m.x[1] + m.x[2] - m.x[3] - 1.0 == 0.0)
#: Dummy Constraints
m.r1 = Constraint(expr=m.p1 == 5.0)
m.r2 = Constraint(expr=m.p2 == 1.0)

#: ipopt suffixes  REQUIRED FOR K_AUG!
m.dual = Suffix(direction=Suffix.IMPORT_EXPORT)
m.ipopt_zL_out = Suffix(direction=Suffix.IMPORT)
m.ipopt_zU_out = Suffix(direction=Suffix.IMPORT)
m.ipopt_zL_in = Suffix(direction=Suffix.EXPORT)
m.ipopt_zU_in = Suffix(direction=Suffix.EXPORT)

#: sipopt suffix
m.sens_state_0 = Suffix(direction=Suffix.EXPORT)
m.sens_state_1 = Suffix(direction=Suffix.EXPORT)
m.sens_state_value_1 = Suffix(direction=Suffix.EXPORT)
m.sens_sol_state_1 = Suffix(direction=Suffix.IMPORT) #: this one contains the results
m.sens_init_constr = Suffix(direction=Suffix.EXPORT)
m.sens_sol_state_1_z_U = Suffix(direction=Suffix.IMPORT)
m.sens_sol_state_1_z_L = Suffix(direction=Suffix.IMPORT)

#: set suffixes
m.p1.set_suffix_value(m.sens_state_0, 1)
m.p1.set_suffix_value(m.sens_state_1, 1)
m.p1.set_suffix_value(m.sens_state_value_1, 4.5)

m.p2.set_suffix_value(m.sens_state_0, 2)
m.p2.set_suffix_value(m.sens_state_1, 2)
m.p2.set_suffix_value(m.sens_state_value_1, 1.0)


m.r1.set_suffix_value(m.sens_init_constr, 1)
m.r2.set_suffix_value(m.sens_init_constr, 1)

sipopt = SolverFactory('ipopt_sens')

#: write some options for ipopt sens
with open('ipopt.opt', 'w') as f:
    f.write('run_sens yes\n')  #: computes the reduced hessian (sens_ipopt)
    f.close()
#: Solve
sipopt.solve(m, tee=True)

#ipopt.solve(m, tee=True)

# print results.
for i in m.component_objects(Suffix):
    i.pprint()
