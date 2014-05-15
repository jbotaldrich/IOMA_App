from gurobipy import *

m = Model("mip1")

x = m.addVar(vtype=GRB.SEMICONT, name = "x")
y = m.addVar(vtype=GRB.SEMICONT, name = "y")
z = m.addVar(vtype=GRB.SEMICONT, name = "z")

m.update()

m.setObjective(x + y + 2 *z, GRB.MAXIMIZE)

m.addConstr(x + 2 * y - 3 * z <= 4, "c0")
m.addConstr(x + y >= 1, "c1")

m.optimize()
m.printAttr('X')
