import general
import flowsolver

test = flowsolver.P2C(0.03, 0.02, 0.75, 0.65, 1.0, 273+70, 273+70, 101325, 1.0, 1.0, 287, 1.40, True)
test.solve()