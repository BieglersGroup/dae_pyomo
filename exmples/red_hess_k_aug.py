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

#: Declare Model
m = ConcreteModel()

m.i = Set(initialize=[1, 2, 3])

init_vals = {1:25E+07, 2:0.0, 3:0.0}
#: Variables
m.x = Var(m.i, initialize=init_vals)
#: Objective
m.oF = Objective(rule=(m.x[1] - 1.0)**2 + (m.x[2] - 2.0)**2 + (m.x[3] - 3.0)**2, sense=minimize)
#: Constraints
m.c1 = Constraint(expr=m.x[1] + 2 * m.x[2] + 3 * m.x[3] == 0.0)

#: ipopt suffixes  REQUIRED FOR K_AUG!
m.dual = Suffix(direction=Suffix.IMPORT_EXPORT)
m.ipopt_zL_out = Suffix(direction=Suffix.IMPORT)
m.ipopt_zU_out = Suffix(direction=Suffix.IMPORT)
m.ipopt_zL_in = Suffix(direction=Suffix.EXPORT)
m.ipopt_zU_in = Suffix(direction=Suffix.EXPORT)

#: sipopt suffix
m.red_hessian = Suffix(direction=Suffix.EXPORT)

#: K_AUG SUFFIXES
m.dof_v = Suffix(direction=Suffix.EXPORT)  #: SUFFIX FOR K_AUG
m.rh_name = Suffix(direction=Suffix.IMPORT)  #: SUFFIX FOR K_AUG AS WELL
#: rh_name will tell us which position the corresponding variable has on the reduced hessian text file.


#: be sure to declare the suffix value (order)
m.x[2].set_suffix_value(m.red_hessian, 1)
m.x[3].set_suffix_value(m.red_hessian, 2)
#: be sure to have ipopt_sens in the path variable
ipopt = SolverFactory('ipopt')
sipopt = SolverFactory('ipopt_sens')
kaug = SolverFactory('k_aug')

#: be sure to declare the suffix value (order)
m.x[2].set_suffix_value(m.dof_v, 1)
m.x[3].set_suffix_value(m.dof_v, 1)

kaug.options["compute_inv"] = ""  #: if the reduced hessian is desired.

#: write some options for ipopt sens
with open('ipopt.opt', 'w') as f:
    f.write('compute_red_hessian yes\n')  #: computes the reduced hessian
    f.write('output_file my_ouput.txt\n')  #: you probably want this file to parse the RH values
    f.close()
#: Solve
sipopt.solve(m, tee=True)
with open('ipopt.opt', 'w') as f:
    f.close()

ipopt.solve(m, tee=True)

m.ipopt_zL_in.update(m.ipopt_zL_out)
m.ipopt_zU_in.update(m.ipopt_zU_out)
#: k_aug
#m.write('problem.nl', format=ProblemFormat.nl)
kaug.solve(m, tee=True)

